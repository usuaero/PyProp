import math
import sys
import warnings
import json
import os
import copy
import time

from random import randint,seed

import multiprocessing as mp
import matplotlib.pyplot as plt
import sqlite3 as sql
import numpy as np

from .electronics import Battery, Motor, ESC
from .propulsion_unit import PropulsionUnit
from .propellers import DatabaseFitProp, DatabaseDataProp, BladeElementProp
from .exceptions import MaxCurrentExceededError, TorquesNotMatchedError, ThrottleNotFoundError, InvalidRuntimeError
from .special_functions import create_component_from_database

class Optimizer:
    """A class for optimizing propulsion units. Will act as a wrapper for lower-level classes
    to provide easy analysis and optimization functionality to the user.
    """

    def __init__(self):

        # Get database location
        self._db_file = os.path.join(os.path.dirname(__file__), "components.db")


    def random_search(self, **kwargs):
        """Randomly searches the design space of components to maximize flight time.
        Plots these flight times against component parameters. User can select specific points in
        the design space to view component parameters and thrust curves.

        Parameters
        ----------
        goal_type : str
            Type of target required. Can be "thrust", "thrust_to_weight_ratio", or "power_to_weight_ratio".
            Note that thrust to weight ratio is dimensionless and power to weight ratio is given in 
            Watts per pound.

        goal_val : float
            Target value.

        airspeed : float
            The cruise velocity at which to evaluate each propulsion system in feet per second.
        
        altitude : float, optional
            The altitude at which to evaluate each propulsion system in feet. Defaults to 0.0.

        airframe_weight : float, optional
            Weight of the aircraft minus the propulsion system in lbf. Defaults to 0.0.

        n_units : int, optional
            The number of random combinations to try. Defaults to 1000.

        n_processes : int, optional
            The maximum number of computer processes that may be spawned at a time. Defaults to 8.

        filename : str, optional
            A file to export the results of the search to, sorted by flight time. Defaults to no file.

        plot : bool, optional
            Whether to plot the results of the search. Defaults to True. A figure will be displayed
            containing 6 plots which describe (in part) the design space. Each point on a plot represents a possible
            design; each design is reflected in each plot. The y axis of each plot is the flight time given by a
            design, and the x axis of each plot is a defining parameter of the design (currently: prop diameter, prop
            pitch, motor Kv constant, battery voltage, battery capacity, and total unit weight). The user may select
            any design in any one of the plots to see a plot of its thrust at various airspeeds and throttle settings.
            A corresponding plot of propeller speeds is also shown and all details of the design are printed to the
            terminal. Selecting a design will also highlight that design in each of the 6 plots, so that general
            patterns in the design space can be opserved.

        prop_model_type : str, optional
            Determines which type of prop model is to be used. Can be "data", "fit", or "BET". "data" props are defined
            by tabulated experimental data for COTS parts. "fit" props are defined by polynomial fits of the 
            experimental data. "BET" props are defined by blade element theory, a numerical model for determining
            prop behavior. **"fit" props are not guaranteed to be realistic. It is recommended these not be used.**
            Defaults to "data".

        prop_constraints : dict, optional
            This can be used to constrain the propellers to a single prop or manufacturer. Formatted as

                {
                    "name" : "<PROP_NAME>"
                }

            or

                {
                    "manufacturer" : "<MANUFACTURER_NAME>"
                }

            These names should be given here exactly as they appear in the database. Defaults to allowing
            any propeller in the database.

        motor_constraints : dict, optional
            Same as "prop_constraints".

        battery_constraints : dict, optional
            Same as "prop_constraints".

        esc_constraints : dict, optional
            Same as "prop_constraints".
        """

        # Get computation parameters
        N_proc_max = kwargs.get("n_processes", 8)
        N_units = kwargs.get("n_units", 1000)
        V_req = kwargs.get("airspeed")
        h = kwargs.get("altitude", 0.0)
        W_frame = kwargs.get("airframe_weight", 0.0)
        prop_model_type = kwargs.get("prop_model_type", "data")

        seed(time.time())

        # Determine goal
        goal_type = kwargs.get("goal_type")
        try :
            goal_code = ["thrust", "thrust_to_weight_ratio", "power_to_weight_ratio"].index(goal_type)
        except ValueError:
            raise IOError("{0} is not a valid goal type.".format(goal_type))

        goal_val = kwargs.get("goal_val")

        # Print out optimization settings
        print("Flight conditions:\nAirspeed {0} ft/s, Altitude {1} ft, Airframe weight: {2} lbs".format(V_req, h, W_frame))
        if goal_code == 1:
            print("Optimizing for a thrust to weight ratio of {0}.".format(goal_val))
        elif goal_code == 0:
            print("Optimizing for a required thrust of {0} lbf.".format(goal_val))
        else:
            print("Optimizing for a power to weight ratio of {0} W/lbf.".format(goal_val))

        # Determine constraints
        names = []
        manufacturers = []
        print("Optimization constrained as follows:")
        for component in ["prop", "motor", "esc", "battery"]:
            self._append_constraints(names, manufacturers, kwargs.get("{0}_constraints".format(component), {}), component)

        # Distribute work
        with mp.Pool(processes=N_proc_max) as pool:
            args = [(V_req, goal_val, h, goal_code, W_frame, names, manufacturers, prop_model_type) for i in range(N_units)]
            data = pool.map(self._get_random_unit, args)

        # Package results
        t_flight, throttles, units, P_e, P_p, eff, cruise_thrust = map(list,zip(*data))
        t_flight = np.asarray(t_flight)
        throttles = np.asarray(throttles)
        units = np.asarray(units)
        P_e = np.asarray(P_e)
        P_p = np.asarray(P_p)
        eff = np.asarray(eff)
        cruise_thrust = np.asarray(cruise_thrust)

        # Determine optimum
        max_ind = np.argmax(t_flight)
        t_max = t_flight[max_ind]
        best_unit = units[max_ind]
        throttle_at_max = throttles[max_ind]
        if goal_code == 1:
            T_req = goal_val*(best_unit.get_weight()+W_frame)
        elif goal_code == 0:
            T_req = goal_val
        else:
            T_req = self._get_thrust_from_power_ratio(goal_val, best_unit.get_weight()+W_frame, V_req)

        # Print information to the terminal
        print(best_unit)
        print("Flight characteristics:")
        print("    Maximum flight time: {0} min".format(t_max))
        print("    Throttle setting: {0}".format(throttles[max_ind]))
        print("    Current draw: {0} A".format(P_e[max_ind]/best_unit.batt.V0))
        print("    Power drawn: {0} W".format(P_e[max_ind]))
        print("    Power developed: {0} lbf*ft/s ({1} W)".format(P_p[max_ind], P_p[max_ind]/0.7375621494575464))
        print("    Propulsive efficiency: {0}%".format(eff[max_ind]*100.0))
        print("    Cruise thrust: {0} lbf".format(cruise_thrust[max_ind]))

        # Plot results
        if kwargs.get("plot", True):

            # Define picker
            def on_pick(event):
                # Handles when the user picks a plotted point in the design space. Highlights that point and plots that unit's thrust curves.

                # Get plot info
                artist = event.artist
                fig = plt.figure(plt.get_fignums()[0])
                ax = fig.axes

                # Determine index
                ind = int(event.ind[0])
                selected_unit = units[ind]

                # Alter title
                fig.suptitle("SELECTED Prop: {0} Motor: {1} Battery: {2} ESC: {3}".format(selected_unit.prop.name, selected_unit.motor.name, selected_unit.batt.name, selected_unit.esc.name))

                # Highlight points in all plots
                ax[0].plot(selected_unit.prop.diameter,t_flight[ind],'*')
                ax[1].plot(selected_unit.prop.pitch,t_flight[ind],'*')
                ax[2].plot(selected_unit.motor.Kv,t_flight[ind],'*')
                ax[3].plot(selected_unit.batt.V0,t_flight[ind],'*')
                ax[4].plot(selected_unit.batt.capacity,t_flight[ind],'*')
                ax[5].plot(selected_unit.get_weight()+W_frame,t_flight[ind],'*')
                ax[6].plot(throttles[ind],t_flight[ind],'*')
                ax[7].plot(eff[ind],t_flight[ind],'*')

                # Print out info
                print(selected_unit)
                print("Flight characteristics:")
                print("    Flight Time: {0} min".format(t_flight[ind]))
                print("    Throttle setting: {0}".format(throttles[ind]*100.0))
                print("    Current draw: {0} A".format(P_e[ind]/selected_unit.batt.V0))
                print("    Power drawn: {0} W".format(P_e[ind]))
                print("    Power developed: {0} lbf*ft/s ({1} W)".format(P_p[ind], P_p[ind]/0.7375621494575464))
                print("    Propulsive efficiency: {0}%".format(eff[ind]*100.0))
                print("    Cruise thrust: {0} lbf".format(cruise_thrust[ind]))

                # Plot performance curves
                selected_unit.plot_thrust_curves([0.0, V_req*2.0+10.0], 11, 51)
                selected_unit.prop.plot_coefs()

            # Plot design space
            plt.ion()
            fig,((ax1,ax2,ax3,ax4),(ax5,ax6,ax7,ax8)) = plt.subplots(nrows=2,ncols=4)
            fig.suptitle("OPTIMUM Prop: {0} Motor: {1} Battery: {2} ESC: {3}".format(best_unit.prop.name, best_unit.motor.name, best_unit.batt.name, best_unit.esc.name))

            ax1.scatter([units[i].prop.diameter for i in range(N_units)], t_flight, picker=3, c=eff)
            ax1.plot(best_unit.prop.diameter,t_max,'r*')
            ax1.set_xlabel("Prop Diameter [in]")
            ax1.set_ylabel("Flight Time [min]")

            ax2.scatter([units[i].prop.pitch for i in range(N_units)],t_flight, c=eff,picker=3)
            ax2.plot(best_unit.prop.pitch,t_max,'r*')
            ax2.set_xlabel("Prop Pitch [in]")
            ax2.set_ylabel("Flight Time [min]")

            ax3.scatter([units[i].motor.Kv for i in range(N_units)],t_flight, c=eff,picker=3)
            ax3.plot(best_unit.motor.Kv,t_max,'r*')
            ax3.set_xlabel("Motor Kv [rpm/V]")
            ax3.set_ylabel("Flight Time [min]")

            ax4.scatter([units[i].batt.V0 for i in range(N_units)],t_flight, c=eff,picker=3)
            ax4.plot(best_unit.batt.V0,t_max,'r*')
            ax4.set_xlabel("Battery Voltage [V]")
            ax4.set_ylabel("Flight Time [min]")

            ax5.scatter([units[i].batt.capacity for i in range(N_units)],t_flight, c=eff,picker=3)
            ax5.plot(best_unit.batt.capacity,t_max,'r*')
            ax5.set_xlabel("Battery Capacity [mAh]")
            ax5.set_ylabel("Flight Time [min]")

            ax6.scatter([units[i].get_weight()+W_frame for i in range(N_units)],t_flight, c=eff,picker=3)
            ax6.plot(best_unit.get_weight()+W_frame,t_max,'r*')
            ax6.set_xlabel("Total Unit Weight [lb]")
            ax6.set_ylabel("Flight Time [min]")

            ax7.scatter(throttles,t_flight, c=eff,picker=3)
            ax7.plot(throttle_at_max,t_max,'r*')
            ax7.set_xlabel("Throttle Setting at Max Flight Time")
            ax7.set_ylabel("Flight Time [min]")

            ax8.scatter(eff, t_flight,  c=eff, picker=3)
            ax8.plot(eff[max_ind], t_max, 'r*')
            ax8.set_xlabel("Efficiency at Max Flight Time")
            ax8.set_ylabel("Flight Time [min]")

            fig.canvas.mpl_connect('pick_event',on_pick)
            plt.show(block=True)
            plt.ioff()

        # Export results
        filename = kwargs.get("filename", None)
        if filename is not None:

            # Open file
            with open(filename, 'w') as export_file:

                # Write header
                header = "{:<25}{:<25}{:<25}{:<25}{:<25}{:<25}{:<25}{:<25}{:<25}{:<25}{:<25}{:<45}{:<25}{:<25}{:<35}{:<25}{:<35}{:<25}{:<25}{:<25}{:<25}"
                header = header.format("Flight Time [min]", "Throttle Setting [%]", "Current Draw [A]", "Power Drawn [W]", "Power Produced [lbf*ft/s]", "Efficiency [%]",
                                       "Cruise Thrust [lbf]", "Prop Name", "Prop Manufacturer", "Prop Diameter [in]", "Prop Pitch [in]", "Motor Name",
                                       "Motor Manufacturer", "Motor Kv [rpm/V]", "ESC Name", "ESC Manufacturer", "Battery Name", "Battery Manufacturer",
                                       "Battery Voltage [V]", "Battery Capacity [mAh]", "Total Weight [lb]")
                print(header, file=export_file)

                # Loop through flight times from greater to least
                sorted_indices = np.argsort(t_flight)
                for i in sorted_indices[::-1]:
                    curr_unit = units[i]
                    row = "{:<25.10}{:<25.10}{:<25.10}{:<25.10}{:<25.10}{:<25.10}{:<25.10}{:<25}{:<25}{:<25.10}{:<25.10}{:<45}{:<25}{:<25.10}{:<35}{:<25}{:<35}{:<25}{:<25.10}{:<25.10}{:<25.10}"
                    row = row.format(t_flight[i], throttles[i]*100.0, P_e[i]/curr_unit.batt.V0, P_e[i], P_p[i], eff[i]*100.0, cruise_thrust[i], curr_unit.prop.name,
                                     curr_unit.prop.manufacturer, curr_unit.prop.diameter, curr_unit.prop.pitch, curr_unit.motor.name, curr_unit.motor.manufacturer,
                                     curr_unit.motor.Kv, curr_unit.esc.name, curr_unit.esc.manufacturer, curr_unit.batt.name, curr_unit.batt.manufacturer, curr_unit.batt.V0,
                                     curr_unit.batt.capacity, curr_unit.get_weight())
                    print(row, file=export_file)


    def _append_constraints(self, names, manufacturers, cnstr_dict, component):
        # Adds the constraints given by the user to the lists

        # Check for name constraint
        name = cnstr_dict.get("name", None)
        names.append(name)

        # Check for manufacturer constraint
        manufacturer = cnstr_dict.get("manufacturer", None)
        manufacturers.append(manufacturer)

        # Check we only have one
        if name is not None and manufacturer is not None:
            raise RuntimeError("A manufacturer and name constraint cannot both be applied to {0}.".format(component))
        
        # Print constraint information
        print("    {0}".format(component.title()))
        if name is not None:
            print("        Name:", name)
        elif manufacturer is not None:
            print("        Manufacturer:", manufacturer)
        else:
            print("        Not constrained.")


    def _get_random_unit(self, args):
        #Selects a propulsion unit and calculates its flight time.

        # Parse arguments
        V_req = args[0]
        goal_val = args[1]
        h = args[2]
        goal_code = args[3]
        W_frame = args[4]
        names = args[5]
        manufacturers = args[6]
        prop_model_type = args[7]

        # Loop through random components until we get a valid one
        t_flight_curr = None
        while t_flight_curr is None or math.isnan(t_flight_curr):

            #Fetch prop data
            if prop_model_type == "data":
                props = os.listdir(os.path.join(os.path.dirname(__file__), "props"))
                prop_ind = randint(0, len(props)-1)
                prop_name = props[prop_ind].replace(".ppdat", "").replace(".ppinf", "")
                prop = DatabaseDataProp(prop_name)
            elif prop_model_type == "fit":
                prop = create_component_from_database(component="prop", name=names[0], manufacturer=manufacturers[0])
            elif prop_model_type == "BET":
                raise RuntimeError("BET prop model is not currently available.")
            else:
                raise IOError("{0} is not a valid prop model type.".format(prop_model_type))

            #Fetch motor data
            motor = create_component_from_database(component="motor", name=names[1], manufacturer=manufacturers[1])

            #Fetch ESC data
            esc = create_component_from_database(component="ESC", name=names[2], manufacturer=manufacturers[2])

            #Fetch battery data
            batt = create_component_from_database(component="battery", name=names[3], manufacturer=manufacturers[3])

            # Determine required thrust
            curr_unit = PropulsionUnit(prop, motor, batt, esc, h)
            if goal_code == 1:
                T_req = (curr_unit.get_weight()+W_frame)*goal_val
            elif goal_code == 0:
                T_req = goal_val
            else:
                T_req = self._get_thrust_from_power_ratio(goal_val, curr_unit.get_weight()+W_frame, V_req)

            # Determine performance
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    t_flight_curr = curr_unit.calc_batt_life(V_req, T_req)
                    thr_curr = curr_unit.calc_cruise_throttle(V_req, T_req)
                    P_e_curr = curr_unit.get_electric_power()
                    P_p_curr = V_req*T_req
                    eff = P_p_curr/(P_e_curr*0.7375621494575464)

            # Anything that goes wrong just means this particular combination isn't up to snuff
            except ZeroDivisionError:
                continue
            except ThrottleNotFoundError:
                continue
            except TorquesNotMatchedError:
                continue
            except MaxCurrentExceededError:
                continue
            except InvalidRuntimeError:
                continue

        # Return params
        return t_flight_curr, thr_curr, curr_unit, P_e_curr, P_p_curr, eff, T_req


    def _get_thrust_from_power_ratio(self, power_ratio, weight, V):
        # Returns the thrust needed based on W/lbf
        return 0.7375621494575464*power_ratio*weight/V