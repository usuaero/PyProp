import numpy as np
import copy
from scipy import integrate
import scipy.interpolate as interp
import os
from os import path
import matplotlib.pyplot as plt
import polyFit as fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from std_atmos import *
import sqlite3 as sql
from random import randint

#Classes in this file are defined such that their information is retrieved from the database (a database cursor must be given).
#If the component's exact name or id are given, that component w_ill be selected. If the manufacturer is given,
#a random component from that manufacturer w_ill be selected. If nothing is specified, a random component is selected.
#The number of battery cells should be specified. If not, it w_ill be randomly selected.

#Converts rads per second to rpms
def toRPM(rads):
    return rads*30/np.pi
    
#A class that defines a battery
class Battery:

    #Initialize the class from database
    def __init__(self, dbcur, name=None, manufacturer=None, dbid=None, numCells=None, capacity=None):

        command = "select * from Batteries"
        if name is not None:
            if manufacturer is not None or dbid is not None:
                raise ValueError("Too many battery parameters specified.")
            command = command+" where Name = '"+name+"'"
        elif manufacturer is not None:
            if dbid is not None:
                raise ValueError("Too many battery parameters specified.")
            command = command+" where manufacturer = '"+manufacturer+"'"
        elif dbid is not None:
            command = command+" where id = "+str(dbid)
        if capacity is not None:
            command = command+" order by abs("+str(capacity)+"-Capacity)"
        command = command+" order by RANDOM() limit 1"

        dbcur.execute(command)
        record = np.asarray(dbcur.fetchall())[0]

        if numCells is None:
            numCells = randint(1,8)

        #Define members from inputs
        self.n = int(numCells)
        self.cellCap = float(record[4])
        self.cellR = float(record[6])
        self.name = record[1]
        self.manufacturer = record[2]
        self.cellWeight = float(record[5])
        self.iMax = float(record[3])
        self.cellV = float(record[7])

        #Members derived from inputs
        self.V0 = self.cellV * self.n
        self.R = self.cellR * self.n
        self.weight = self.cellWeight*self.n

    def printInfo(self):
        print("Battery:",self.name)
        print("\tManufacturer:",self.manufacturer)
        print("\tCapacity:",self.cellCap)
        print("\tNum Cells:",self.n)
        print("\tVoltage:",self.V0)
        print("\tWeight:",self.weight)

#A class that defines an ESC (Electronic Speed Controller)
class ESC:

    #Initialization of the class from database
    def __init__(self, dbcur, name=None, manufacturer=None, dbid=None, I_max=None):

        command = "select * from ESCs"
        if name is not None:
            if manufacturer is not None or dbid is not None:
                raise ValueError("Too many esc parameters specified.")
            command = command+" where Name = '"+name+"'"
        elif manufacturer is not None:
            if dbid is not None:
                raise ValueError("Too many ESC parameters specified.")
            command = command+" where manufacturer = '"+manufacturer+"'"
        elif dbid is not None:
            command = command+" where id = "+str(dbid)
        if I_max is not None:
            command = command+" order by abs("+str(I_max)+"-I_motorax)"
        command = command+" order by RANDOM() limit 1"

        dbcur.execute(command)
        record = np.asarray(dbcur.fetchall())[0]

        self.R = float(record[6])
        self.name = record[1]
        self.manufacturer = record[2]
        self.iMax = float(record[3])
        self.weight = float(record[5])

    def printInfo(self):
        print("ESC:",self.name)
        print("\tManufacturer:",self.manufacturer)
        print("\tMax Current:",self.iMax)
        print("\tWeight:",self.weight)
        
#A class that defines an electric motor.
class Motor:

    #Initialization of the class from the database
    def __init__(self, dbcur, name=None, manufacturer=None, dbid=None, Kv=None):

        command = "select * from Motors"
        if name is not None:
            if manufacturer is not None or dbid is not None:
                raise ValueError("Too many motor parameters specified.")
            command = command+" where Name = '"+name+"'"
        elif manufacturer is not None:
            if dbid is not None:
                raise ValueError("Too many motor parameters specified.")
            command = command+" where manufacturer = '"+manufacturer+"'"
        elif dbid is not None:
            command = command+" where id = "+str(dbid)
        if Kv is not None:
            command = command+" order by abs("+str(Kv)+"-kv)"
        command = command+" order by RANDOM() limit 1"

        dbcur.execute(command)
        record = np.asarray(dbcur.fetchall())[0]

        self.Kv = float(record[3])
        self.Gr = float(record[4])
        self.I0 = float(record[6])
        self.R = float(record[5])
        self.name = record[1]
        self.manufacturer = record[2]
        self.weight = float(record[7])

    def printInfo(self):
        print("Motor:",self.name)
        print("\tManufacturer:",self.manufacturer)
        print("\tKv:",self.Kv)
        print("\tWeight:",self.weight)

#A class of propellers defined by database test files
class Propeller:
    
    #Initializes the prop from the database
    def __init__(self, dbcur, name=None, manufacturer=None, dbid=None, diameter=None, pitch=None):

        command = "select * from Props"
        if name is not None:
            if manufacturer is not None or dbid is not None:
                raise ValueError("Too many prop parameters specified.")
            command = command+" where Name = '"+name+"'"
        elif manufacturer is not None:
            if dbid is not None:
                raise ValueError("Too many prop parameters specified.")
            command = command+" where manufacturer = '"+manufacturer+"'"
        elif dbid is not None:
            command = command+" where id = "+dbid
        if diameter is not None:
            command = command+" order by abs("+str(diameter)+"-Diameter)"
        if pitch is not None:
            command = command+" order by abs("+str(pitch)+"-Pitch)"
        command = command+" order by RANDOM() limit 1"

        dbcur.execute(command)
        record = np.asarray(dbcur.fetchall())[0]

        self.name = record[1]
        self.manufacturer = record[2]
        self.diameter = float(record[3])
        self.pitch = float(record[4])
        self.thrustFitOrder = int(record[5])
        self.fitOfThrustFitOrder = int(record[6])
        self.powerFitOrder = int(record[7])
        self.fitOfPowerFitOrder = int(record[8])

        numThrustCoefs = (self.thrustFitOrder+1)*(self.fitOfThrustFitOrder+1)
        self.thrustCoefs = record[9:numThrustCoefs+9].reshape((self.thrustFitOrder+1,self.fitOfThrustFitOrder+1)).astype(np.float)
        self.powerCoefs = record[numThrustCoefs+9:].reshape((self.powerFitOrder+1,self.fitOfPowerFitOrder+1)).astype(np.float)

        #These parameters w_ill be set by later functions
        self.v_inf = 0.0
        self.angVel = 0.0
        
    def printInfo(self):
        print("Propeller:",self.name)
        print("\tManufacturer:",self.manufacturer)
        print("\tDiameter:",self.diameter)
        print("\tPitch:",self.pitch)

    def CalcTorqueCoef(self):
        self.rpm = toRPM(self.angVel)
        self.rps = self.rpm/60
        if abs(self.rps)<1e-10:
            self.J = 10000 #To prevent errors. Since angular velocity is 0, actual value w_ill also be 0.
        else:
            self.J = self.v_inf/(self.rps*self.diameter/12)
        a = fit.poly_func(self.powerCoefs.T, self.rpm)
        if(a[-1]>0):#Quadratic coefficient should always be non-positive
            a[-1] = 0
        self.Cl = fit.poly_func(a, self.J)/2*np.pi
        

    def CalcThrustCoef(self):
        self.rpm = toRPM(self.angVel)
        self.rps = self.rpm/60
        if abs(self.rps)<1e-10:
            self.J = 10000 #To prevent errors. Since angular velocity is 0, actual value w_ill also be 0.
        else:
            self.J = self.v_inf/(self.rps*self.diameter/12)
        a = fit.poly_func(self.thrustCoefs.T, self.rpm)
        if(a[-1]>0):#Quadratic coefficient should always be non-positive
            a[-1] = 0
        self.Ct = fit.poly_func(a, self.J)

    def PlotCoefs(self):
        #Plot thrust and torque coefficients
        rpms = np.linspace(0,35000,10)
        Js = np.linspace(0,1.4,10)
        fig = plt.figure(figsize=plt.figaspect(1.))
        fig.suptitle(self.name)
        ax = fig.add_subplot(1,2,1, projection='3d')

        for rpm in rpms:
            a = fit.poly_func(self.thrustCoefs.T, rpm)
            if(a[-1]>0):#Quadratic coefficient should always be non-positive
                a[-1] = 0
            thrust = fit.poly_func(a, Js)
            rpmForPlot = np.full(len(thrust),rpm)
            ax.plot(Js,rpmForPlot,thrust, 'r-')

        ax.set_title("Predicted Thrust")
        ax.set_xlabel("Advance Ratio")
        ax.set_ylabel("RPM")
        ax.set_zlabel("Thrust Coefficient")

        ax = fig.add_subplot(1,2,2, projection='3d')

        for rpm in rpms:
            a = fit.poly_func(self.powerCoefs.T, rpm)
            if(a[-1]>0):#Quadratic coefficient should always be non-positive
                a[-1] = 0
            power = fit.poly_func(a, Js)
            rpmForPlot = np.full(len(power),rpm)
            ax.plot(Js,rpmForPlot,power, 'r-')

        ax.set_title("Predicted Power")
        ax.set_xlabel("Advance Ratio")
        ax.set_ylabel("RPM")
        ax.set_zlabel("Power Coefficient")
        plt.show()

#A class that defines an entire electric propulsion unit
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
        self.prop.angVel = 0
        self.I_motor = 0 #Instantaneous current being drawn through the motor

    #Computes motor torque (ft*lbf) given throttle setting and revolutions (rpm)
    def CalcMotorTorque(self, throttle, revs):
        etaS = 1 - 0.078*(1 - throttle)
        self.I_motor = (etaS*throttle*self.batt.V0 - (self.motor.Gr/self.motor.Kv)*revs)/(etaS*throttle*self.batt.R + self.esc.R + self.motor.R)
        # Note: the 7.0432 constant converts units [(Nm/ftlb)(min/s)(rad/rev)]^-1
        return 7.0432*self.motor.Gr/self.motor.Kv * (self.I_motor - self.motor.I0)
    
    #Computes thrust produced at a given cruise speed and throttle setting
    def CalcCruiseThrust(self, v_cruise, throttle):
        if v_cruise == 0 and throttle == 0:
            self.prop.angVel = 0
            return 0 #Don't even bother

        self.prop.v_inf = v_cruise

        #Determine the shaft angular velocity at which the motor torque and propeller torque are matched
        #Uses a secant method
        err_max = 0.000001
        err_aprx = 1 + err_max #So that it executes at least once
        w_0 = 950 #An initial guess of the prop's angular velocity
        w_max = self.motor.Kv*self.batt.V0*throttle*(2*np.pi/60) # Theoretically the upper limit
        self.prop.angVel = w_0
        self.prop.CalcTorqueCoef()
        f_0 = self.CalcMotorTorque(throttle, toRPM(w_0)) - self.prop.Cl*self.airDensity*(w_0/(2*np.pi))**2*(self.prop.diameter/12)**5
        w_1 = w_0 * 1.1
        iterations = 0
        
        while err_aprx >= err_max and iterations < 1000:
            iterations = iterations + 1
            self.prop.angVel = w_1
            self.prop.CalcTorqueCoef()
            T_motor = self.CalcMotorTorque(throttle, toRPM(w_1))
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
                self.prop.angVel = w_i
                self.prop.CalcTorqueCoef()
                T_motor[i] = self.CalcMotorTorque(throttle, toRPM(w_i))
                T_prop[i] = self.prop.Cl*self.airDensity*(w_i/(2*np.pi))**2*(self.prop.diameter/12)**5
            plt.plot(w,T_motor)
            plt.plot(w,T_prop)
            plt.title("Torques vs Angular Velocity")
            plt.legend(["Motor Torque","Prop Torque"])
            plt.show()
        
        self.prop.angVel = w_2
        self.prop.CalcThrustCoef()
        _ = self.CalcMotorTorque(throttle, toRPM(w_2)) # To make sure member variables are fully updated

        return self.prop.Ct*self.airDensity*(w_2/(2*np.pi))**2*(self.prop.diameter/12)**4
    
    #Computes required throttle setting for a given thrust and cruise speed
    def CalcCruiseThrottle(self, v_cruise, T_req):
        #Uses a secant method
        err_max = 0.000001
        err_aprx = 1 + err_max
        t_0 = 0.5
        T_0 = self.CalcCruiseThrust(v_cruise, t_0)
        t_1 = t_0*1.1
        iterations = 0
        
        while err_aprx >= err_max and iterations < 1000:
            
            iterations = iterations + 1
            T_1 = self.CalcCruiseThrust(v_cruise, t_1) - T_req
            
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
        #        T[i] = self.CalcCruiseThrust(v_cruise, t[i]) - T_req
        #    plt.plot(t,T) 
        #    plt.show()

        if t_2 > 1 or t_2 < 0:
            return None
        
        self.CalcCruiseThrust(v_cruise,t_2) # To make sure member variables are fully updated
        return t_2
        
    #Plots thrust curves for propulsion unit up to a specified airspeed
    def PlotThrustCurves(self, v_min, v_max, numVels, numThrSets):
        
        vel = np.linspace(v_min, v_max, numVels)
        thr = np.linspace(0, 1, numThrSets)
        thrust = np.zeros((numVels, numThrSets))
        rpm = np.zeros((numVels,numThrSets))
        
        for i in range(numVels):
            for j in range(numThrSets):
                
                #print("Freestream Velocity: ", vel[i])
                #print("Throttle Setting: ", thr[j])
                thrust[i][j] = self.CalcCruiseThrust(vel[i], thr[j])
                rpm[i][j] = toRPM(self.prop.angVel)

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
    def CalcBattLife(self, v_cruise, T_req):
        throttle = self.CalcCruiseThrottle(v_cruise, T_req)
        if(throttle==None or self.I_motor > self.esc.iMax or self.I_motor > self.batt.iMax):
            return None
        #print("Throttle Setting:",throttle)
        #print("Current Draw:",self.I_motor)
        runTime = (self.batt.cellCap/1000)/self.I_motor*60 # Gives run time in minutes, assuming nominal cell capacity and constant battery votlage
        if runTime < 0:
            return None
        return runTime

    def GetWeight(self):#Returns weight of electrical components in pounds
        return (self.batt.weight + self.motor.weight + self.esc.weight)/16

    def printInfo(self):
        print("----Propulsion Unit----")
        self.prop.printInfo()
        self.motor.printInfo()
        self.esc.printInfo()
        self.batt.printInfo()
