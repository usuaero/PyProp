"""Classes for modelling propellers."""

import os

import numpy as np
import matplotlib.pyplot as plt
import sqlite3 as sql

from .poly_fit import poly_func
from .helpers import to_rpm
from .base_classes import DatabaseComponent

class DatabaseProp(DatabaseComponent):
    """Defines a propeller by database test files."""
    
    def __init__(self, **kwargs):
        super().__init__()

        # Get database record
        record = self.get_database_record("prop", **kwargs)

        # Store params
        self.name = record[1]
        self.manufacturer = record[2]
        self.diameter = float(record[3])
        self.pitch = float(record[4])
        self.thrustFitOrder = int(record[5])
        self.fitOfThrustFitOrder = int(record[6])
        self.powerFitOrder = int(record[7])
        self.fitOfPowerFitOrder = int(record[8])

        # Determine fit properties
        numThrustCoefs = (self.thrustFitOrder+1)*(self.fitOfThrustFitOrder+1)
        self.thrustCoefs = record[9:numThrustCoefs+9].reshape((self.thrustFitOrder+1,self.fitOfThrustFitOrder+1)).astype(np.float)
        self.powerCoefs = record[numThrustCoefs+9:].reshape((self.powerFitOrder+1,self.fitOfPowerFitOrder+1)).astype(np.float)

        #These parameters will be set by later functions
        self.v_inf = 0.0
        self.w = 0.0

        
    def __str__(self):
        string = "Propeller: {0}".format(self.name)
        string += "\n\tManufacturer: {0}".format(self.manufacturer)
        string += "\n\tDiameter: {0}".format(self.diameter)
        string += "\n\tPitch: {0}".format(self.pitch)
        return string


    def CalcTorqueCoef(self):
        self.rpm = to_rpm(self.w)
        self.rps = self.rpm/60
        if abs(self.rps)<1e-10:
            self.J = 10000 #To prevent errors. Since angular velocity is 0, actual value will also be 0.
        else:
            self.J = self.v_inf/(self.rps*self.diameter/12)
        a = poly_func(self.powerCoefs.T, self.rpm)
        if(a[-1]>0):#Quadratic coefficient should always be non-positive
            a[-1] = 0
        self.Cl = poly_func(a, self.J)/2*np.pi
        

    def CalcThrustCoef(self):
        self.rpm = to_rpm(self.w)
        self.rps = self.rpm/60
        if abs(self.rps)<1e-10:
            self.J = 10000 #To prevent errors. Since angular velocity is 0, actual value will also be 0.
        else:
            self.J = self.v_inf/(self.rps*self.diameter/12)
        a = poly_func(self.thrustCoefs.T, self.rpm)
        if(a[-1]>0):#Quadratic coefficient should always be non-positive
            a[-1] = 0
        self.Ct = poly_func(a, self.J)

    def PlotCoefs(self):
        #Plot thrust and torque coefficients
        rpms = np.linspace(0,35000,10)
        Js = np.linspace(0,1.4,10)
        fig = plt.figure(figsize=plt.figaspect(1.))
        fig.suptitle(self.name)
        ax = fig.add_subplot(1,2,1, projection='3d')

        for rpm in rpms:
            a = poly_func(self.thrustCoefs.T, rpm)
            if(a[-1]>0):#Quadratic coefficient should always be non-positive
                a[-1] = 0
            thrust = poly_func(a, Js)
            rpmForPlot = np.full(len(thrust),rpm)
            ax.plot(Js,rpmForPlot,thrust, 'r-')

        ax.set_title("Predicted Thrust")
        ax.set_xlabel("Advance Ratio")
        ax.set_ylabel("RPM")
        ax.set_zlabel("Thrust Coefficient")

        ax = fig.add_subplot(1,2,2, projection='3d')

        for rpm in rpms:
            a = poly_func(self.powerCoefs.T, rpm)
            if(a[-1]>0):#Quadratic coefficient should always be non-positive
                a[-1] = 0
            power = poly_func(a, Js)
            rpmForPlot = np.full(len(power),rpm)
            ax.plot(Js,rpmForPlot,power, 'r-')

        ax.set_title("Predicted Power")
        ax.set_xlabel("Advance Ratio")
        ax.set_ylabel("RPM")
        ax.set_zlabel("Power Coefficient")
        plt.show()


class BladeElementProp:
    """Defines the performance of a propeller using blade element theory."""

    def __init__(self):
        pass