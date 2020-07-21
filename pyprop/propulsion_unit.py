"""A class that defines an entire electric propulsion unit."""

import numpy as np
import matplotlib.pyplot as plt

from .propellers import DatabaseFitProp, BladeElementProp
from .electronics import Battery, ESC, Motor
from .std_atmos import statee, statsi
from .helpers import to_rpm

class PropulsionUnit:
    """Defines a full electric propulsion unit.

    Parameters
    ----------
    prop : DatabaseFitProp, DatabaseDataProp, or BladeElementProp
        Propeller object.

    motor : Motor
        Motor object.

    battery : Battery
        Battery object.

    esc : ESC
        ESC object.

    altitude : float, optional
        Altitude in feet. Defaults to 0.0.
    """

    def __init__(self, prop, motor, battery, esc, altitude=0.0):
        
        self.prop = prop
        self.motor = motor
        self.batt = battery
        self.esc = esc

        _,_,_,self._rho = statee(altitude) # Converts kg/m^3 to slug/ft^3
        

    def calc_motor_torque(self, throttle, revs):
        """Computes motor torque (ft*lbf) given throttle setting and revolutions (rpm).

        Parameters
        ----------
        throttle : float
            Throttle setting.

        revs : float
            Angular velocity of the motor in rpm.

        Returns
        -------
        float
            Torque produced by the motor in Nm.
        """

        # Determine motor current
        etaS = 1 - 0.078*(1 - throttle)
        self.I_motor = (etaS*throttle*self.batt.V0 - (self.motor.Gr/self.motor.Kv)*revs)/(etaS*throttle*self.batt.R + self.esc.R + self.motor.R)

        # Determine torque
        # Note: the 7.0432 constant converts units [(Nm/ftlb)(min/s)(rad/rev)]^-1
        return 7.0432*self.motor.Gr/self.motor.Kv * (self.I_motor - self.motor.I0)
    
    def calc_cruise_thrust(self, v_cruise, throttle):
        """Computes thrust produced at a given cruise speed and throttle setting.

        Parameters
        ----------
        v_cruise : float
            Freestream velocity in feet per second.

        throttle : float
            Throttle setting.

        Returns
        -------
        float
            Thrust in lbf.
        """

        # Check for zero inputs
        if v_cruise == 0.0 and throttle == 0.0:
            return 0.0

        #Determine the shaft angular velocity at which the motor torque and propeller torque are matched
        #Uses a secant method
        err_max = 0.000001
        err_aprx = 1 + err_max #So that it executes at least once
        w_0 = 950 #An initial guess of the prop's angular velocity
        w_max = self.motor.Kv*self.batt.V0*throttle*(2*np.pi/60) # Theoretically the upper limit
        Cl_prop = self.prop.get_torque_coef(w_0, v_cruise)
        f_0 = self.calc_motor_torque(throttle, to_rpm(w_0))-Cl_prop*self._rho*(w_0/(2*np.pi))**2*(self.prop.diameter/12)**5
        w_1 = w_0 * 1.1
        iterations = 0
        
        while err_aprx >= err_max and iterations < 1000:
            iterations = iterations + 1
            Cl_prop = self.prop.get_torque_coef(w_1, v_cruise)
            T_motor = self.calc_motor_torque(throttle, to_rpm(w_1))
            T_prop = Cl_prop*self._rho*(w_1/(2*np.pi))**2*(self.prop.diameter/12)**5
            f_1 = T_motor - T_prop
            
            w_2 = w_1 - (f_1*(w_0 - w_1))/(f_0 - f_1)
            if w_2 < 0: # Prop angular velocity will never be negative even if windmilling
                w_2 = 0.00001

            err_aprx = abs((w_2 - w_1)/w_2)
            
            w_0 = w_1
            f_0 = f_1
            w_1 = w_2
    
        if False: #iterations >= 1000:
            w = np.linspace(0,30000,10000)
            T_motor = np.zeros(10000)
            T_prop = np.zeros(10000)
            for i,w_i in enumerate(w):
                self.prop.w = w_i
                Cl_prop = self.prop.get_torque_coef(w_i, v_cruise)
                T_motor[i] = self.calc_motor_torque(throttle, to_rpm(w_i))
                T_prop[i] = Cl_prop*self._rho*(w_i/(2*np.pi))**2*(self.prop.diameter/12)**5
            plt.plot(w,T_motor)
            plt.plot(w,T_prop)
            plt.title("Torques vs Angular Velocity")
            plt.legend(["Motor Torque","Prop Torque"])
            plt.show()
        
        Ct = self.prop.get_thrust_coef(w_2, v_cruise)
        self._w = w_2
        self.calc_motor_torque(throttle, to_rpm(w_2)) # To make sure member variables are fully updated

        return Ct*self._rho*(w_2/(2*np.pi))**2*(self.prop.diameter/12)**4
    
    #Computes required throttle setting for a given thrust and cruise speed
    def calc_cruise_throttle(self, v_cruise, T_req):
        #Uses a secant method
        err_max = 0.000001
        err_aprx = 1 + err_max
        t_0 = 0.5
        T_0 = self.calc_cruise_thrust(v_cruise, t_0)
        t_1 = t_0*1.1
        iterations = 0
        
        while err_aprx >= err_max and iterations < 1000:
            
            iterations = iterations + 1
            T_1 = self.calc_cruise_thrust(v_cruise, t_1) - T_req
            
            t_2 = t_1 - (T_1*(t_0 - t_1))/(T_0 - T_1)
            
            err_aprx = abs((t_2 - t_1)/t_2)
            
            if t_2 > 10:
                t_2 = 1.1
            elif t_2 < -10:
                t_2 = -0.1
            t_0 = t_1
            T_0 = T_1
            t_1 = t_2

        #if iterations == 1000:
        #    t = np.linspace(0,1.0,100)
        #    T = np.zeros(100)
        #    for i in range(100):
        #        T[i] = self.calc_cruise_thrust(v_cruise, t[i]) - T_req
        #    plt.plot(t,T) 
        #    plt.show()

        if t_2 > 1 or t_2 < 0:
            return None
        
        self.calc_cruise_thrust(v_cruise, t_2) # To make sure member variables are fully updated
        return t_2
        
    def plot_thrust_curves(self, v_lims, n_vel=10, n_thr=10):
        """Plots thrust curves for the propulsion unit through the range of specified velocities.

        Parameters
        ----------
        v_lims : list
            Limits for the range of airspeeds to be plotted.

        n_vel : int
            Number of velocity values to plot. Defaults to 10.

        n_thr : int
            Number of throttle values to plot. Throttle is ranged from 0 to 1. Defaults to 10.
        """
        
        # Get ranges
        vel = np.linspace(v_lims[0], v_lims[1], n_vel)
        thr = np.linspace(0.0, 1.0, n_thr)
        thrust = np.zeros((n_vel, n_thr))
        rpm = np.zeros((n_vel,n_thr))
        
        # Loop through to calculate thrust and rpm values
        for i, v in enumerate(vel):
            for j, t in enumerate(thr):
                
                thrust[i][j] = self.calc_cruise_thrust(v, t)
                rpm[i][j] = to_rpm(self._w)

        # Initialize figure
        fig = plt.figure()
        fig.suptitle("Components: " + str(self.prop.name) + ", " + str(self.motor.name) + ", and " + str(self.batt.name))

        # Plot thrust
        ax0 = fig.add_subplot(1,2,1)
        for i in range(n_vel):
            ax0.plot(thr, thrust[i])
        ax0.set_title("Thrust")
        ax0.set_ylabel("Thrust [lbf]")
        ax0.set_xlabel("Throttle Setting")
        ax0.legend(list(vel), title="Airspeed [ft/s]")

        # Plot rpm
        ax1 = fig.add_subplot(1,2,2)
        for i in range(n_vel):
            ax1.plot(thr, rpm[i])
        ax1.set_title("Prop Speed")
        ax1.set_ylabel("Speed [rpms]")
        ax1.set_xlabel("Throttle Setting")
        plt.show()


    def calc_batt_life(self, v_cruise, T_req):
        """Determines how long the battery w_ill last based on a required thrust and cruise speed.
        This assumes nominal cell capacity and constant battery discharge voltage.

        Parameters
        ----------
        v_cruise : float
            Cruise velocity in feet per second.

        T_req : float
            Required thrust in lbf.

        Returns
        -------
        float
            Flight time in minutes
        """

        # Calculate required throttle
        throttle = self.calc_cruise_throttle(v_cruise, T_req)

        # Check we have a result
        if(throttle==None or self.I_motor > self.esc.I_max or self.I_motor > self.batt.I_max):
            return None

        # Determine the run time
        run_time = (self.batt.cell_cap/1000)/self.I_motor*60 # Gives run time in minutes, assuming nominal cell capacity and constant battery votlage
        if run_time < 0:
            return None
        return run_time


    def get_weight(self):
        """Returns weight of electrical components in pounds.

        Returns
        -------
        float
            Weight of the propulsion system (minus the prop) in pounds.
        """
        return (self.batt.weight + self.motor.weight + self.esc.weight)/16.0

    def __str__(self):
        string = "----Propulsion Unit----"
        string += "\n{0}".format(str(self.prop))
        string += "\n{0}".format(str(self.motor))
        string += "\n{0}".format(str(self.esc))
        string += "\n{0}".format(str(self.batt))
        return string