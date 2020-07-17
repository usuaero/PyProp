"""A class that defines an entire electric propulsion unit."""

import numpy as np
import matplotlib.pyplot as plt

from .propellers import DatabaseProp, BladeElementProp
from .electronics import Battery, ESC, Motor
from .std_atmos import statee, statsi
from .helpers import to_rpm

class PropulsionUnit:

    #Initialize the class from subclasses which are previously initialized
    def __init__(self, prop, motor, battery, esc, altitude):
        
        self.prop = prop
        self.motor = motor
        self.batt = battery
        self.esc = esc

        _,_,_,self.airDensity = statee(altitude) # Converts kg/m^3 to slug/ft^3
        
        #Initialize exterior parameters to be set later
        self.prop.v_inf = 0
        self.prop.w = 0
        self.I_motor = 0 #Instantaneous current being drawn through the motor

    #Computes motor torque (ft*lbf) given throttle setting and revolutions (rpm)
    def calc_motor_torque(self, throttle, revs):
        etaS = 1 - 0.078*(1 - throttle)
        self.I_motor = (etaS*throttle*self.batt.V0 - (self.motor.Gr/self.motor.Kv)*revs)/(etaS*throttle*self.batt.R + self.esc.R + self.motor.R)
        # Note: the 7.0432 constant converts units [(Nm/ftlb)(min/s)(rad/rev)]^-1
        return 7.0432*self.motor.Gr/self.motor.Kv * (self.I_motor - self.motor.I0)
    
    #Computes thrust produced at a given cruise speed and throttle setting
    def calc_cruise_thrust(self, v_cruise, throttle):
        if v_cruise == 0 and throttle == 0:
            self.prop.w = 0
            return 0 #Don't even bother

        self.prop.v_inf = v_cruise

        #Determine the shaft angular velocity at which the motor torque and propeller torque are matched
        #Uses a secant method
        err_max = 0.000001
        err_aprx = 1 + err_max #So that it executes at least once
        w_0 = 950 #An initial guess of the prop's angular velocity
        w_max = self.motor.Kv*self.batt.V0*throttle*(2*np.pi/60) # Theoretically the upper limit
        self.prop.w = w_0
        self.prop.CalcTorqueCoef()
        f_0 = self.calc_motor_torque(throttle, to_rpm(w_0)) - self.prop.Cl*self.airDensity*(w_0/(2*np.pi))**2*(self.prop.diameter/12)**5
        w_1 = w_0 * 1.1
        iterations = 0
        
        while err_aprx >= err_max and iterations < 1000:
            iterations = iterations + 1
            self.prop.w = w_1
            self.prop.CalcTorqueCoef()
            T_motor = self.calc_motor_torque(throttle, to_rpm(w_1))
            T_prop = self.prop.Cl*self.airDensity*(w_1/(2*np.pi))**2*(self.prop.diameter/12)**5
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
                self.prop.CalcTorqueCoef()
                T_motor[i] = self.calc_motor_torque(throttle, to_rpm(w_i))
                T_prop[i] = self.prop.Cl*self.airDensity*(w_i/(2*np.pi))**2*(self.prop.diameter/12)**5
            plt.plot(w,T_motor)
            plt.plot(w,T_prop)
            plt.title("Torques vs Angular Velocity")
            plt.legend(["Motor Torque","Prop Torque"])
            plt.show()
        
        self.prop.w = w_2
        self.prop.CalcThrustCoef()
        _ = self.calc_motor_torque(throttle, to_rpm(w_2)) # To make sure member variables are fully updated

        return self.prop.Ct*self.airDensity*(w_2/(2*np.pi))**2*(self.prop.diameter/12)**4
    
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
        
        self.calc_cruise_thrust(v_cruise,t_2) # To make sure member variables are fully updated
        return t_2
        
    #Plots thrust curves for propulsion unit up to a specified airspeed
    def plot_thrust_curves(self, v_min, v_max, numVels, numThrSets):
        
        vel = np.linspace(v_min, v_max, numVels)
        thr = np.linspace(0, 1, numThrSets)
        thrust = np.zeros((numVels, numThrSets))
        rpm = np.zeros((numVels,numThrSets))
        
        for i in range(numVels):
            for j in range(numThrSets):
                
                #print("Freestream Velocity: ", vel[i])
                #print("Throttle Setting: ", thr[j])
                thrust[i][j] = self.calc_cruise_thrust(vel[i], thr[j])
                rpm[i][j] = to_rpm(self.prop.w)

        fig = plt.figure()
        fig.suptitle("Components: " + str(self.prop.name) + ", " + str(self.motor.name) + ", and " + str(self.batt.name))

        ax0 = fig.add_subplot(1,2,1)
        for i in range(numVels):
            ax0.plot(thr, thrust[i])
        ax0.set_title("Thrust")
        ax0.set_ylabel("Thrust [lbf]")
        ax0.set_xlabel("Throttle Setting")
        ax0.legend(list(vel), title="Airspeed [ft/s]")

        ax1 = fig.add_subplot(1,2,2)
        for i in range(numVels):
            ax1.plot(thr, rpm[i])
        ax1.set_title("Prop Speed")
        ax1.set_ylabel("Speed [rpms]")
        ax1.set_xlabel("Throttle Setting")
        plt.show()

    #Determines how long the battery w_ill last based on a required thrust and cruise speed
    def calc_batt_life(self, v_cruise, T_req):
        throttle = self.calc_cruise_throttle(v_cruise, T_req)
        if(throttle==None or self.I_motor > self.esc.I_max or self.I_motor > self.batt.I_max):
            return None
        #print("Throttle Setting:",throttle)
        #print("Current Draw:",self.I_motor)
        run_time = (self.batt.cell_cap/1000)/self.I_motor*60 # Gives run time in minutes, assuming nominal cell capacity and constant battery votlage
        if run_time < 0:
            return None
        return run_time

    def get_weight(self):#Returns weight of electrical components in pounds
        return (self.batt.weight + self.motor.weight + self.esc.weight)/16

    def __str__(self):
        string = "----Propulsion Unit----"
        string += "\n{0}".format(str(self.prop))
        string += "\n{0}".format(str(self.motor))
        string += "\n{0}".format(str(self.esc))
        string += "\n{0}".format(str(self.batt))
        return string