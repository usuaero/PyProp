"""Classes for modelling propellers."""

import os
import json

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp

from abc import abstractmethod

from .poly_fit import poly_func
from .helpers import to_rpm, to_rads

class BaseProp:
    """Defines a propeller."""

    def __init__(self):
        self.name = None
        self.manufacturer = None
        self.diameter = None
        self.pitch = None

        
    def __str__(self):
        string = "Propeller: {0}".format(self.name)
        string += "\n\tManufacturer: {0}".format(self.manufacturer)
        string += "\n\tDiameter: {0} inches".format(self.diameter)
        string += "\n\tPitch: {0} inches".format(self.pitch)
        return string


    def get_advance_ratio(self, w, V):
        """Determines the advance ratio from angular velocity in rad/s and velocity in ft/s.

        Parameters
        ----------
        w : float
            Angular velocity in rad/s.

        V : float
            Velocity in ft/s.

        Returns
        -------
        float
            Advance ratio.
        """
        rps = to_rpm(w)/60.0
        if abs(rps) < 1e-10:
            return np.inf
        else:
            return 12.0*V/(rps*self.diameter)


    def get_velocity(self, rpm, J):
        """Determines the velocity from the rpm and advance ratio.

        Parameters
        ----------
        rpm : float
            Angular velocity in rpm.

        J : float
            Advance ratio.

        Returns
        -------
        float
            Velocity in ft/s.
        """
        rps = rpm/60.0
        return J*rps*self.diameter/12.0


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
            thrust = np.zeros_like(Js)
            for i, J in enumerate(Js):
                w = to_rads(rpm)
                V = self.get_velocity(rpm, J)
                thrust[i] = self.get_thrust_coef(w, V)
            rpms_to_plot = np.full(len(thrust), rpm)
            ax.plot(Js, rpms_to_plot, thrust, 'r-')

        # Format figure
        ax.set_title("Predicted Thrust")
        ax.set_xlabel("Advance Ratio")
        ax.set_ylabel("RPM")
        ax.set_zlabel("Thrust Coefficient")

        ax = fig.add_subplot(1,2,2, projection='3d')

        # Loop through rpm values to plot torque
        for rpm in rpms:
            power = np.zeros_like(Js)
            for i, J in enumerate(Js):
                w = to_rads(rpm)
                V = self.get_velocity(rpm, J)
                power[i] = self.get_torque_coef(w, V)*2.0*np.pi
            rpms_to_plot = np.full(len(power), rpm)
            ax.plot(Js, rpms_to_plot, power, 'r-')

        # Format
        ax.set_title("Predicted Power")
        ax.set_xlabel("Advance Ratio")
        ax.set_ylabel("RPM")
        ax.set_zlabel("Power Coefficient")
        plt.show()


    @abstractmethod
    def get_torque_coef(self, w, V):
        pass


    @abstractmethod
    def get_thrust_coef(self, w, V):
        pass


class DatabaseFitProp(BaseProp):
    """Defines a propeller by database fits.

    Parameters
    ----------
    record : list
        SQL database record defining this prop. This class should
        be initialized using the pyprop.create_component_from_database()
        function.

    """
    
    def __init__(self, record):

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

        # Get revs per minute and J
        rpm = to_rpm(w)
        J = self.get_advance_ratio(w, V)

        # Evaluate function to get fit coefficients
        a = poly_func(self.power_coefs.T, rpm)
        if(a[-1]>0): # Quadratic coefficient should always be non-positive
            a[-1] = 0

        # Evaluate function to get torque coefficient
        return poly_func(a, J)/2.0*np.pi
        

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

        # Get revs per minute and J
        rpm = to_rpm(w)
        J = self.get_advance_ratio(w, V)
        
        # Evaluate function to get fit coefficients
        a = poly_func(self.thrust_coefs.T, rpm)

        # Quadratic coefficient should always be non-positive
        if(a[-1]>0):
            a[-1] = 0

        return poly_func(a, J)


class DatabaseDataProp(BaseProp):
    """Defines a propeller by tabulated data. The data for these
    props are contained in the pyprop/props/ directory. PyProp will
    import that data automatically to intialize the prop.

    Parameters
    ----------
    name : str
        The name of the prop as stored in the props/ directory. This
        is the same as the name of the prop in the components database.
        The name should be given without the ".ppdat" or ".ppinf" file
        extension.
    """

    def __init__(self, name):
        
        # Load data file
        prop_dir = os.path.join(os.path.dirname(__file__), "props")
        data_filename = os.path.join(prop_dir, name+".ppdat")
        with open(data_filename, 'r') as data_file:
            self._data = np.genfromtxt(data_file, skip_header=1)

        # Load info
        info_filename = os.path.join(prop_dir, name+".ppinf")
        with open(info_filename, 'r') as info_file:
            self._info_dict = json.load(info_file)

        # Store params
        self.name = self._info_dict["name"]
        self.manufacturer = self._info_dict["manufacturer"]
        self.diameter = self._info_dict["diameter"]
        self.pitch = self._info_dict["pitch"]


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

        # Get revs per minute and J
        rpm = to_rpm(w)
        J = self.get_advance_ratio(w, V)

        return interp.griddata(self._data[:,:2], self._data[:,3], np.array([[rpm, J]]), method='linear').item()/(2.0*np.pi)
        

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

        # Get revs per minute and J
        rpm = to_rpm(w)
        J = self.get_advance_ratio(w, V)
        
        return interp.griddata(self._data[:,:2], self._data[:,2], np.array([[rpm, J]]), method='linear').item()


class BladeElementProp(BaseProp):
    """Defines the performance of a propeller using blade element theory."""

    def __init__(self, **kwargs):
        pass