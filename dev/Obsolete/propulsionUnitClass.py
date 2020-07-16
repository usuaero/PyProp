import numpy as np
import matplotlib.pyplot as plt
import sqlite3 as sql
import supportClasses as s
from skaero.atmosphere import coesa

#Note that dimensionless coefficients assume that angular velocity is given in revs/sec.

#Converts rads per second to rpms
def toRPM(rads):
    return rads*30/np.pi

#A class that defines an entire electric propulsion unit
class PropulsionUnit:

    #Initialize the class from subclasses which are previously initialized
    def __init__(self, prop, motor, battery, esc, altitude):
        
        self.prop = prop
        self.motor = motor
        self.batt = battery
        self.esc = esc

        _,_,_,self.airDensity = coesa.table(altitude*0.3048)
        self.airDensity = self.airDensity*0.0019403203 # Converts kg/m^3 to slug/ft^3
        
        #Initialize exterior parameters to be set later
        self.prop.vInf = 0
        self.prop.angVel = 0
        self.Im = 0 #Instantaneous current being drawn through the motor

    #Computes motor torque (ft*lbf) given throttle setting and revolutions (rpm)
    def CalcTorque(self, throttle, revs):
        etaS = 1 - 0.078*(1 - throttle)
        self.Im = (etaS*throttle*self.batt.V0 - (self.motor.Gr/self.motor.Kv)*revs)/(etaS*throttle*self.batt.R + self.esc.R + self.motor.R)
        # Note: the 7.0432 constant converts units [(Nm/ftlb)(min/s)(rad/rev)]^-1
        return 7.0432*self.motor.Gr/self.motor.Kv * (self.Im - self.motor.I0)
    
    #Computes thrust produced at a given cruise speed and throttle setting
    def CalcCruiseThrust(self, cruiseSpeed, throttle):
        if cruiseSpeed == 0 and throttle == 0:
            return 0 #Don't even bother

        self.prop.vInf = cruiseSpeed

        #Determine the shaft angular velocity at which the motor torque and propeller torque are matched
        #Uses a secant method
        errorBound = 0.000001
        approxError = 1 + errorBound #So that it executes at least once
        w0 = 300 #An initial guess of the prop's angular velocity
        self.prop.angVel = w0
        self.prop.CalcTorqueCoef()
        f0 = self.CalcTorque(throttle, toRPM(w0)) - self.prop.Cl*self.airDensity*(w0/(2*np.pi))**2*(self.prop.diameter/12)**5
        w1 = w0 * 1.1
        iterations = 0
        
        while approxError >= errorBound and iterations < 1000:
            iterations = iterations + 1
            self.prop.angVel = w1
            self.prop.CalcTorqueCoef()
            motorTorque = self.CalcTorque(throttle, toRPM(w1))
            propTorque = self.prop.Cl*self.airDensity*(w1/(2*np.pi))**2*(self.prop.diameter/12)**5
            f1 = motorTorque - propTorque
            
            w2 = w1 - (f1*(w0 - w1))/(f0 - f1)
            if w2 < 0: # Prop angular velocity will never be negative even if windmilling
                w2 = 0.01
            elif w2 > self.motor.Kv*self.batt.V0: #Theoretically the upper limit
                w2 = w2/3;

            approxError = abs((w2 - w1)/w2)
            
            w0 = w1
            f0 = f1
            w1 = w2
    
        if iterations == 1000:
            print("Problem speed:",cruiseSpeed)
            print("Problem throttle:",throttle)
            w = np.linspace(0,3000000,10000)
            Tm = np.zeros(10000)
            Tp = np.zeros(10000)
            for i in range(1000):
                self.prop.angVel = w[i]
                self.prop.CalcTorqueCoef()
                Tm[i] = self.CalcTorque(throttle, toRPM(w[i]))
                Tp[i] = self.prop.Cl*self.airDensity*(w[i]/(2*np.pi))**2*(self.prop.diameter/12)**5
            plt.plot(w,Tm)
            plt.plot(w,Tp)
            plt.plot(w,Tm-Tp)
            plt.title("Torque Balance vs Angular Velocity")
            plt.show()
        
        self.prop.angVel = w2
        self.prop.CalcThrustCoef()
        _ = self.CalcTorque(throttle, toRPM(w2)) # To make sure member variables are fully updated

        return self.prop.Ct*self.airDensity*(w2/(2*np.pi))**2*(self.prop.diameter/12)**4
    
    #Computes required throttle setting for a given thrust and cruise speed
    def CalcCruiseThrottle(self, cruiseSpeed, reqThrust):
        #Uses a secant method
        errorBound = 0.000001
        approxError = 1 + errorBound
        t0 = 0.5
        T0 = self.CalcCruiseThrust(cruiseSpeed, t0)
        t1 = t0*1.1
        iterations = 0
        
        while approxError >= errorBound and iterations < 1000:
            
            iterations = iterations + 1
            T1 = self.CalcCruiseThrust(cruiseSpeed, t1) - reqThrust
            
            t2 = t1 - (T1*(t0 - t1))/(T0 - T1)
            
            approxError = abs((t2 - t1)/t2)
            
            if t2 > 10:
                t2 = 1.1
            elif t2 < -10:
                t2 = -0.1
            t0 = t1
            T0 = T1
            t1 = t2

        #if iterations == 1000:
        #    t = np.linspace(0,1.0,100)
        #    T = np.zeros(100)
        #    for i in range(100):
        #        T[i] = self.CalcCruiseThrust(cruiseSpeed, t[i]) - reqThrust
        #    plt.plot(t,T) 
        #    plt.show()

        if t2 > 1 or t2 < 0:
            return None
        
        self.CalcCruiseThrust(cruiseSpeed,t2) # To make sure member variables are fully updated
        return t2
        
    #Plots thrust curves for propulsion unit up to a specified airspeed
    def PlotThrustCurves(self, maxAirspeed, numVels, numThrSets):
        
        vel = np.linspace(0, maxAirspeed, numVels)
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

    #Determines how long the battery will last based on a required thrust and cruise speed
    def CalcBattLife(self, cruiseSpeed, reqThrust):
        throttle = self.CalcCruiseThrottle(cruiseSpeed, reqThrust)
        if(throttle==None or self.Im > self.esc.iMax or self.Im > self.batt.iMax):
            return None
        #print("Throttle Setting:",throttle)
        #print("Current Draw:",self.Im)
        runTime = (self.batt.cellCap/1000)/self.Im*60 # Gives run time in minutes, assuming nominal cell capacity and constant battery votlage
        if runTime < 0:
            return None
        return runTime

    def GetWeight(self):#Returns weight of electrical components in pounds
        return (self.batt.weight + self.motor.weight + self.esc.weight)/16
