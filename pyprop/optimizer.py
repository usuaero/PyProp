import math
import sys
import warnings
import datetime
import json
import os
import copy

from random import randint,seed

import multiprocessing as mp
import matplotlib.pyplot as plt
import sqlite3 as sql
import numpy as np

from .electronics import Battery, Motor, ESC
from .propulsion_unit import PropulsionUnit
from .propellers import DatabaseFitProp, DatabaseDataProp, BladeElementProp

class Optimizer:
    """A class for optimizing propulsion units. Will act as a wrapper for lower-level classes
    to provide easy analysis and optimization functionality to the user.
    """

    def __init__(self):

        # Get database location
        self._db_file = os.path.join(os.path.dirname(__file__), "components.db")


    def random_search(self, input_dict):
        """Randomly searches the design space of components to maximize flight time.
        Plots these flight times against component parameters. User can select specific points in
        the design space to view component parameters and thrust curves.

        The only argument is a dictionary defining the search parameters, an example of which is
        given below (explanatory comments given within // //):

        {
            "computation":{
                "units":1000, //Number of propulsion units to find in the design space.//
                "processes":8, //Maximum number of processes to be used in parallel computation.//
                "outlierStdDevs":5 //Number of standard deviations of the half-normal distribution within which designs are considered feasible.//
            },
            "condition":{
                "altitude":0, //Flight altitude.//
                "airspeed":10 //Flight cruise speed.//
            },
            "goal":{ //One and only one of these parameters must be specified. Set for cruise condition.//
                "thrust":0, //Thrust required from the propulsion unit.//
                "thrustToWeightRatio":0.3 //Thrust to weight ratio required (requires emptyWeight to be defined.//
            },
            "aircraft":{
                "emptyWeight":1, //Weight of the aircraft minus the propulsion system.//
                "components":{ //These parameters are optional, but only one for each component may be specified.//
                    "propeller":{
                        "name":"",
                        "manufacturer":""
                    },
                    "motor":{
                        "name":"",
                        "manufacturer":""
                    },
                    "esc":{
                        "name":"",
                        "manufacturer":""
                    },
                    "battery":{
                        "name":"",
                        "manufacturer":""
                    }
                }
            }
        }

        The component parameters in the dict file are all optional. Specifying a component name limits
        the search to propulsion units including that specific component. Specifying a component manufacturer
        limits the search to a single manufacturer for that component. Please note that some component
        manufacturers have very few components in our current database, and specifying this may limit the
        search more than desirable. Only one of these parameters may be specified for each component at most.

        Once the search is complete (i.e. the specified number of designs has been considered), the propulsion
        unit which has the longest flight time will be output to the terminal. A figure will also be displayed
        containing 6 plots which describe (in part) the design space. Each point on a plot represents a possible
        design; each design is reflected in each plot. The y axis of each plot is the flight time given by a
        design, and the x axis of each plot is a defining parameter of the design (currently: prop diameter, prop
        pitch, motor Kv constant, battery voltage, battery capacity, and total unit weight). The user may select
        any design in any one of the plots to see a plot of its thrust at various airspeeds and throttle settings.
        A corresponding plot of propeller speeds is also shown and all details of the design are printed to the
        terminal. Selecting a design will also highlight that design in each of the 6 plots, so that general
        patterns in the design space can be opserved.

        The developer is of the opinion that outliers should be ignored. Testing has shown that these arise from
        error in the component models and should not be trusted as feasible, high-performance designs. Realistic
        designs will be found closer to the main cluster of designs.
        """

        # Load input
        if isinstance(input_dict, str):
            with open(input_dict, 'r') as input_handle:
                input_dict = json.load(input_handle)

        # Get parameters
        N_proc_max = input_dict["computation"].get("processes", 8)
        N_units = input_dict["computation"]["units"]
        V_req = input_dict["condition"]["airspeed"]
        h = input_dict["condition"]["airspeed"]
        W_frame = input_dict["aircraft"]["emptyWeight"]

        # Determine goal
        if input_dict["goal"]["thrust"] is 0:
            if input_dict["goal"]["thrustToWeightRatio"] is 0:
                raise RuntimeError("No goal specified!")
            optimize_thrust_to_weight = True
            thrust_param = input_dict["goal"]["thrustToWeightRatio"]
            R_tw_req = thrust_param
        else:
            optimize_thrust_to_weight = False
            thrust_param = input_dict["goal"]["thrust"]
            T_req = thrust_param

        # Print out optimization settings
        print("Flight conditions:\nAirspeed {0} ft/s, Altitude {1} ft, Airframe weight: {2} lbs".format(V_req, h, W_frame))
        if optimize_thrust_to_weight:
            print("Optimizing for a thrust to weight ratio of",thrust_param)
        else:
            print("Optimizing for a required thrust of",thrust_param)

        # Determine constraints
        names = []
        manufacturers = []
        print("Optimization constrained as follows:")
        for component in input_dict["aircraft"]["components"]:
            name = input_dict["aircraft"]["components"][component]["name"]
            if len(name) is 0:
                name = None
            names.append(name)

            manufacturer = input_dict["aircraft"]["components"][component]["manufacturer"]
            if len(manufacturer) is 0:
                manufacturer = None
            manufacturers.append(manufacturer)

            if name is not None and manufacturer is not None:
                raise RuntimeError("Component: "+component+" is overconstrained!")
            
            print(component.title())
            if name is not None:
                print("Name:",name)
            elif manufacturer is not None:
                print("Manufacturer:",manufacturer)
            else:
                print("Not constrained.")

        # Distribute work
        with mp.Pool(processes=N_proc_max) as pool:
            args = [(V_req, thrust_param, h, optimize_thrust_to_weight, W_frame, names, manufacturers) for i in range(N_units)]
            data = pool.map(self._get_random_unit, args)

        # Package results
        t_flight, throttles, units = map(list,zip(*data))
        t_flight = np.asarray(t_flight)
        throttles = np.asarray(throttles)
        units = np.asarray(units)

        # Determine optimum
        max_ind = np.argmax(t_flight)
        t_max = t_flight[max_ind]
        best_unit = units[max_ind]
        throttle_at_max = throttles[max_ind]
        if optimize_thrust_to_weight:
            T_req = thrust_param*(best_unit.get_weight()+W_frame)
        else:
            T_req = thrust_param

        # Print information to the terminal
        print("Maximum flight time found: {0} min".format(t_max))
        print(best_unit)
        print("Throttle setting for max flight: {0}".format(best_unit.calc_cruise_throttle(V_req, T_req)))
        print("Current draw: {0} A".format(best_unit.I_motor))

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
            fig.suptitle("SELECTED Prop: "+str(selected_unit.prop.name)+"  Motor: "+str(selected_unit.motor.name)+"  Battery: "+str(selected_unit.batt.name)+"  ESC: "+str(selected_unit.esc.name))

            # Highlight points in all plots
            ax[0].plot(selected_unit.prop.diameter,t_flight[ind],'o')
            ax[1].plot(selected_unit.prop.pitch,t_flight[ind],'o')
            ax[2].plot(selected_unit.motor.Kv,t_flight[ind],'o')
            ax[3].plot(selected_unit.batt.V0,t_flight[ind],'o')
            ax[4].plot(selected_unit.batt.cell_cap,t_flight[ind],'o')
            ax[5].plot(selected_unit.get_weight()+W_frame,t_flight[ind],'o')
            ax[6].plot(throttles[ind],t_flight[ind],'o')

            # Print out info
            print(selected_unit)
            print("Flight Time:",t_flight[ind],"min")
            if optimize_thrust_to_weight:
                print("    at {:4.2f}% throttle".format(selected_unit.calc_cruise_throttle(V_req,(selected_unit.get_weight()+W_frame)*R_tw_req)*100))
            else:
                print("    at {:4.2f}% throttle".format(selected_unit.calc_cruise_throttle(V_req,T_req)*100))

            # Plot performance curves
            selected_unit.plot_thrust_curves([0.0, V_req*2.0+10.0], 11, 51)
            selected_unit.prop.plot_coefs()

        # Plot design space
        plt.ion()
        fig,((ax1,ax2,ax3,ax4),(ax5,ax6,ax7,ax8)) = plt.subplots(nrows=2,ncols=4)
        fig.suptitle("OPTIMUM Prop: "+str(best_unit.prop.name)+"  Motor: "+str(best_unit.motor.name)+"  Battery: "+str(best_unit.batt.name)+"  ESC: "+str(best_unit.esc.name))

        ax1.plot([units[i].prop.diameter for i in range(N_units)],t_flight,'b*',picker=3)
        ax1.plot(best_unit.prop.diameter,t_max,'r*')
        ax1.set_xlabel("Prop Diameter [in]")
        ax1.set_ylabel("Flight Time [min]")

        ax2.plot([units[i].prop.pitch for i in range(N_units)],t_flight,'b*',picker=3)
        ax2.plot(best_unit.prop.pitch,t_max,'r*')
        ax2.set_xlabel("Prop Pitch [in]")
        ax2.set_ylabel("Flight Time [min]")

        ax3.plot([units[i].motor.Kv for i in range(N_units)],t_flight,'b*',picker=3)
        ax3.plot(best_unit.motor.Kv,t_max,'r*')
        ax3.set_xlabel("Motor Kv [rpm/V]")
        ax3.set_ylabel("Flight Time [min]")

        ax4.plot([units[i].batt.V0 for i in range(N_units)],t_flight,'b*',picker=3)
        ax4.plot(best_unit.batt.V0,t_max,'r*')
        ax4.set_xlabel("Battery Voltage [V]")
        ax4.set_ylabel("Flight Time [min]")

        ax5.plot([units[i].batt.cell_cap for i in range(N_units)],t_flight,'b*',picker=3)
        ax5.plot(best_unit.batt.cell_cap,t_max,'r*')
        ax5.set_xlabel("Cell Capacity [mAh]")
        ax5.set_ylabel("Flight Time [min]")

        ax6.plot([units[i].get_weight()+W_frame for i in range(N_units)],t_flight,'b*',picker=3)
        ax6.plot(best_unit.get_weight()+W_frame,t_max,'r*')
        ax6.set_xlabel("Total Unit Weight [lb]")
        ax6.set_ylabel("Flight Time [min]")

        ax7.plot(throttles,t_flight,'b*',picker=3)
        ax7.plot(throttle_at_max,t_max,'r*')
        ax7.set_xlabel("Throttle Setting at Max Flight Time")
        ax7.set_ylabel("Flight Time [min]")

        fig.canvas.mpl_connect('pick_event',on_pick)
        plt.show(block=True)
        plt.ioff()

    def _get_random_unit(self, args):
        #Selects a propulsion unit and calculates its flight time.

        # Parse arguments
        V_req = args[0]
        T = args[1]
        h = args[2]
        optimize_thrust_to_weight = args[3]
        W_frame = args[4]
        names = args[5]
        manufacturers = args[6]

        # Determine thrust requirement
        if optimize_thrust_to_weight:
            R_tw_req = T
        else:
            T_req = T

        # Loop through random components until we get a valid one
        t_flight_curr = None
        while t_flight_curr is None or math.isnan(t_flight_curr):

            #Fetch prop data
            prop = DatabaseFitProp(name=names[0], manufacturer=manufacturers[0])

            #Fetch motor data
            motor = Motor(name=names[1], manufacturer=manufacturers[1])

            #Fetch ESC data
            esc = ESC(name=names[2], manufacturer=manufacturers[2])

            #Fetch battery data
            batt = Battery(name=names[3], manufacturer=manufacturers[3])

            if batt.R == 0 and esc.R == 0 and motor.R == 0:
                continue

            # Determine performance
            curr_unit = PropulsionUnit(prop,motor,batt,esc,h)
            if optimize_thrust_to_weight:
                T_req = (curr_unit.get_weight()+W_frame)*R_tw_req
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                t_flight_curr = curr_unit.calc_batt_life(V_req,T_req)
                thr_curr = curr_unit.calc_cruise_throttle(V_req,T_req)

        # Return params
        return t_flight_curr, thr_curr, curr_unit