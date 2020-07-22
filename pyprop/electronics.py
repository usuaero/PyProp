"""Classes defining electronic components for the propulsion unit."""

import os

import numpy as np
import sqlite3 as sql

from random import randint

    
class Battery:
    """Defines a battery.

    Parameters
    ----------
    name : str, optional
        Name of the battery to be used.

    manufacturer : str, optional
        Manufacturer of the battery.

    num_cells : int
        Number of cells in the battery.

    capacity : float, optional
        Capacity of the battery in mAh.

    resistance : float, optional
        Internal resistance of the battery in Ohms.
    
    voltage : float
        Nominal voltage of the battery in Volts.

    weight : float
        Weight  of the battery in ounces.

    I_max : float
        Maximum current draw of the battery in Amps. Defaults to infinity.
    """

    def __init__(self, **kwargs):

        # Get parameters from kwargs
        self.capacity = kwargs.get("capacity")
        self.R = kwargs.get("resistance")
        self.V0 = kwargs.get("voltage")
        self.weight = kwargs.get("weight")
        self.name = kwargs.get("name", "Generic Battery")
        self.manufacturer = kwargs.get("manufacturer", "PyProp")
        self.I_max = kwargs.get("I_max", np.inf)
        self.n = kwargs.get("num_cells")


    def __str__(self):
        string = "Battery: {0}".format(self.name)
        string += "\n\tManufacturer: {0}".format(self.manufacturer)
        string += "\n\tCapacity: {0} mAh".format(self.capacity)
        string += "\n\tNum Cells: {0}".format(self.n)
        string += "\n\tVoltage: {0} V".format(self.V0)
        string += "\n\tWeight: {0} oz".format(self.weight)
        return string


class ESC:
    """Defines an electronic speed controller (ESC).

    Parameters
    ----------
    name : str, optional
        Name of the ESC to be used.

    manufacturer : str, optional
        Manufacturer of the ESC.

    resistance : float, optional
        Equivalent "on" resistance of the ESC in Ohms.
    
    I_max : float
        Maximum current that can be sourced by the ESC in amps. Defaults to infinity.

    weight : float
        Weight of the ESC in ounces.
    """

    def __init__(self, **kwargs):

        # Set params
        self.R = kwargs.get("resistance")
        self.name = kwargs.get("name", "Generic ESC")
        self.manufacturer = kwargs.get("manufacturer", "PyProp")
        self.I_max = kwargs.get("I_max", np.inf)
        self.weight = kwargs.get("weight")


    def __str__(self):
        string = "ESC: {0}".format(self.name)
        string += "\n\tManufacturer: {0}".format(self.manufacturer)
        string += "\n\tMax Current: {0} A".format(self.I_max)
        string += "\n\tWeight: {0} oz".format(self.weight)
        return string

        
class Motor:
    """Defines an electric motor.

    Parameters
    ----------
    name : str, optional
        Name of the motor to be used.

    manufacturer : str, optional
        Manufacturer of the motor.

    Kv : float, optional
        Kv rating (also called the speed constant) of the motor.

    resistance : float, optional
        DC resistance of the motor in Ohms.
    
    I_no_load : float
        No-load current of the motor in amps. Defaults to 0.0
    
    I_max : float
        Maximum current that can be handled by the motor in amps. Defaults
        to infinity.

    gear_ratio : float, optional
        Gear ratio of the motor. A value greater than unity will correspond
        to the output shaft of the motor turning slower than the motor itself.
        Defaults to 1.0.

    weight : float
        Weight of the motor in ounces.
    """

    def __init__(self, **kwargs):

        # Set params
        self.Kv = kwargs.get("Kv")
        self.Gr = kwargs.get("gear_ratio", 1.0)
        self.I0 = kwargs.get("I_no_load", 0.0)
        self.I_max = kwargs.get("I_max", np.inf)
        self.R = kwargs.get("resistance")
        self.name = kwargs.get("name", "Generic Motor")
        self.manufacturer = kwargs.get("manufacturer", "PyProp")
        self.weight = kwargs.get("weight")

        # Determine torque constant
        self.Kt = 7.0432/self.Kv # Converts to ft*lb/A


    def __str__(self):
        string = "Motor: {0}".format(self.name)
        string += "\n\tManufacturer: {0}".format(self.manufacturer)
        string += "\n\tKv: {0}".format(self.Kv)
        string += "\n\tWeight: {0} oz".format(self.weight)
        return string