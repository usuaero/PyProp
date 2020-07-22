"""Classes for modelling propellers."""

import os

import numpy as np
import matplotlib.pyplot as plt
import sqlite3 as sql

from .poly_fit import poly_func
from .helpers import to_rpm
from .base_classes import DatabaseComponent

class DatabaseFitProp(DatabaseComponent):
    """Defines a propeller by database fits."""
    
    def __init__(self, **kwargs):
        super().__init__()

        # Get database record
        record = self.get_database_record("prop", **kwargs)

        # Store params
        self.name = record[1]
        self.manufacturer = record[2]
        self.diameter = float(record[3])
        self.pitch = float(record[4])

        # Determine fit properties
        thrust_fit_order = int(record[5])
        fit_of_thrust_fit_order = int(record[6])
        power_fit_order = int(record[7])
        fit_of_power_fit_order = int(record[8])
        num_coefs = (thrust_fit_order+1)*(fit_of_thrust_fit_order+1)
        self.thrust_coefs = record[9:num_coefs+9].reshape((thrust_fit_order+1,fit_of_thrust_fit_order+1)).astype(np.float)
        self.power_coefs = record[num_coefs+9:].reshape((power_fit_order+1,fit_of_power_fit_order+1)).astype(np.float)

        
    def __str__(self):
        string = "Propeller: {0}".format(self.name)
        string += "\n\tManufacturer: {0}".format(self.manufacturer)
        string += "\n\tDiameter: {0} inches".format(self.diameter)
        string += "\n\tPitch: {0} inches".format(self.pitch)
        return string


    def get_torque_coef(self, w, V):
        """Returns the torque coefficient for the prop at the given angular velocity and freestream velocity.

        Parameters
        ----------
        w : float
            Angular velocity of the prop in radians per second.

        V : float
            Freestream velocity in feet per second.

        Returns
        -------
        float
            Torque coefficient.
        """

        # Get revs per minute and second
        rpm = to_rpm(w)
        rps = rpm/60

        # Check for angular velocity equal to 0, actual value will also be 0
        if abs(rps)<1e-10:
            J = 10000
        else:
            J = V/(rps*self.diameter/12)

        # Evaluate function to get fit coefficients
        a = poly_func(self.power_coefs.T, rpm)
        if(a[-1]>0): # Quadratic coefficient should always be non-positive
            a[-1] = 0

        # Evaluate function to get torque coefficient
        return poly_func(a, J)/2*np.pi
        

    def get_thrust_coef(self, w, V):
        """Returns the thrust coefficient for the prop at the given angular velocity and freestream velocity.

        Parameters
        ----------
        w : float
            Angular velocity of the prop in radians per second.

        V : float
            Freestream velocity in feet per second.

        Returns
        -------
        float
            Thrust coefficient.
        """

        # Get revs per minute and second
        rpm = to_rpm(w)
        rps = rpm/60

        # Check for angular velocity equal to 0, actual value will also be 0
        if abs(rps)<1e-10:
            J = 10000
        else:
            J = V/(rps*self.diameter/12)
        
        # Evaluate function to get fit coefficients
        a = poly_func(self.thrust_coefs.T, rpm)

        # Quadratic coefficient should always be non-positive
        if(a[-1]>0):
            a[-1] = 0

        return poly_func(a, J)


    def plot_coefs(self, rpm_lims=[0.0, 35000.0], J_lims=[0.0, 1.4]):
        """Plot thrust and torque coefficients as functions of rpm and advance ratio.

        Parameters
        ----------
        rpm_lims : list, optional
            Limits for plotting in rpm. Defaults to [0.0, 35000.0].

        J_lims : list, optional
            Limits for plotting in advance ratio. Defaults to [0.0, 1.4].
        """

        # Get distributions
        rpms = np.linspace(rpm_lims[0], rpm_lims[1], 10)
        Js = np.linspace(J_lims[0], J_lims[1], 10)

        # Initialize figure
        fig = plt.figure(figsize=plt.figaspect(1.))
        fig.suptitle(self.name)
        ax = fig.add_subplot(1,2,1, projection='3d')

        # Loop through rpm values to plot thrust
        for rpm in rpms:
            a = poly_func(self.thrust_coefs.T, rpm)
            if(a[-1]>0):#Quadratic coefficient should always be non-positive
                a[-1] = 0
            thrust = poly_func(a, Js)
            rpms_to_plot = np.full(len(thrust),rpm)
            ax.plot(Js, rpms_to_plot, thrust, 'r-')

        # Format figure
        ax.set_title("Predicted Thrust")
        ax.set_xlabel("Advance Ratio")
        ax.set_ylabel("RPM")
        ax.set_zlabel("Thrust Coefficient")

        ax = fig.add_subplot(1,2,2, projection='3d')

        # Loop through rpm values to plot torque
        for rpm in rpms:
            a = poly_func(self.power_coefs.T, rpm)
            if(a[-1]>0):#Quadratic coefficient should always be non-positive
                a[-1] = 0
            power = poly_func(a, Js)
            rpms_to_plot = np.full(len(power),rpm)
            ax.plot(Js, rpms_to_plot, power, 'b-')

        # Format
        ax.set_title("Predicted Power")
        ax.set_xlabel("Advance Ratio")
        ax.set_ylabel("RPM")
        ax.set_zlabel("Power Coefficient")
        plt.show()


class DatabaseDataProp:
    """Defines a propeller by tabulated data."""

    def __init__(self):
        pass


class BladeElementProp:
    """Defines the performance of a propeller using blade element theory."""

    def __init__(self):
        pass