"""Classes defining electronic components for the propulsion unit."""

import os

import numpy as np
import sqlite3 as sql

from random import randint

from .base_classes import DatabaseComponent
    
class Battery(DatabaseComponent):
    """Defines a battery.

    Parameters
    ----------
    type : str, optional
        Can be "database" or "user_defined". Defaults to "database".

        A "database" component is one selected from the database of COTS
        components available within PyProp. Any parameters given will limit
        the component to match that parameter. E.g. if "name" is given for
        a "database" component, the component with that name in the database
        will be selected (must be exactly as appears in the database).
        Otherwise, the component will be selected randomly. If "manufacturer" 
        is given, a random component from that manufacturer will be selected.
        A numerical value provided (such as battery capacity or motor Kv) will
        cause the component in the database with the value of that parameter
        closest to the given value to be selected.

        A "user_defined" component simply has its parameters given by the 
        user. This allows the user the freedom to specify components not in
        the database.

    name : str, optional
        Name of the battery to be used.

    manufacturer : str, optional
        Manufacturer of the battery.

    dbid : int, optional
        ID of the component in the database. Valid only for "database" type.

    capacity : float, optional
        Capacity of a *single cell* of the battery in mAh. Not optional for
        "user_defined" type.

    resistance : float, optional
        Internal resistance of a *single cell* of the battery in Ohms. Valid
        only for "user_defined" type. Defaults to 0.0.
    
    voltage : float
        Nominal voltage of a *single cell* of the battery in Volts. Valid only
        for "user_defined" type.

    weight : float
        Weight of a *single cell* of the battery in ounces. Valid only for
        "user_defined" type.

    I_max : float
        Maximum current draw of a *single cell* of the battery in Amps. Valid
        only for "user_defined" type.
    """

    def __init__(self, **kwargs):

        # Get type
        self._type = kwargs.get("type", "database")

        if self._type == "database":
            
            # Get database info
            super().__init__()
            record = self.get_database_record("battery", **kwargs)

            # If the number of cells has not been specified, set a random number
            if kwargs.get("num_cells", None) is None:
                num_cells = randint(1,8)

            #Define members from inputs
            self.n = int(num_cells)
            self.cell_cap = float(record[4])
            self.cell_R = float(record[6])
            self.cell_V = float(record[7])
            self.cell_weight = float(record[5])
            self.name = record[1]
            self.manufacturer = record[2]
            self.I_max = float(record[3])

        elif self._type == "user_defined":

            # Get parameters from kwargs
            self.n = kwargs.get("num_cells")
            self.cell_cap = kwargs.get("capacity")
            self.cell_R = kwargs.get("resistance", 0.0)
            self.cell_V = kwargs.get("voltage")
            self.cell_weight = kwargs.get("weight")
            self.name = kwargs.get("name", "user_battery")
            self.manufacturer = kwargs.get("manufacturer", "user")
            self.I_max = kwargs.get("I_max")

        #Members derived from inputs
        self.V0 = self.cell_V * self.n
        self.R = self.cell_R * self.n
        self.weight = self.cell_weight*self.n


    def __str__(self):
        string = "Battery: {0}".format(self.name)
        string += "\n\tManufacturer: {0}".format(self.manufacturer)
        string += "\n\tCapacity: {0}".format(self.cell_cap)
        string += "\n\tNum Cells: {0}".format(self.n)
        string += "\n\tVoltage: {0}".format(self.V0)
        string += "\n\tWeight: {0}".format(self.weight)
        return string


class ESC(DatabaseComponent):
    """Defines an electronic speed controller (ESC).

        Parameters
        ----------
        type : str, optional
            Can be "database" or "user_defined". For more details on this, see
            the docstrings for the Battery object. Defaults to "database".

        name : str, optional
            Name of the ESC to be used.

        manufacturer : str, optional
            Manufacturer of the ESC.

        dbid : int, optional
            ID of the component in the database. Valid only for "database" type.

        resistance : float, optional
            Equivalent resistance of the ESC in Ohms. Valid only for "user_defined"
            type. Defaults to 0.0.
        
        I_max : float
            Maximum current that can be sourced by the ESC in amps.

        weight : float
            Weight of the ESC in ounces. Valid only for "user_defined" type.
    """

    def __init__(self, **kwargs):

        # Get type
        self._type = kwargs.get("type", "database")

        if self._type == "database":
            
            # Get database info
            super().__init__()
            record = self.get_database_record("ESC", **kwargs)

            # Set params
            self.R = float(record[6])
            self.name = record[1]
            self.manufacturer = record[2]
            self.I_max = float(record[3])
            self.weight = float(record[5])

        elif self._type == "user_defined":

            # Set params
            self.R = kwargs.get("resistance", 0.0)
            self.name = kwargs.get("name", "user_ESC")
            self.manufacturer = kwargs.get("manufacturer", "user")
            self.I_max = kwargs.get("I_max")
            self.weight = kwargs.get("weight")


    def __str__(self):
        string = "ESC: {0}".format(self.name)
        string += "\n\tManufacturer: {0}".format(self.manufacturer)
        string += "\n\tMax Current: {0}".format(self.I_max)
        string += "\n\tWeight: {0}".format(self.weight)
        return string

        
class Motor(DatabaseComponent):
    """Defines an electric motor.

        Parameters
        ----------
        type : str, optional
            Can be "database" or "user_defined". For more details on this, see
            the docstrings for the Battery object. Defaults to "database".

        name : str, optional
            Name of the motor to be used.

        manufacturer : str, optional
            Manufacturer of the motor.

        dbid : int, optional
            ID of the component in the database. Valid only for "database" type.

        Kv : float, optional
            Kv rating (sometimes called the speed constant) of the motor.

        resistance : float, optional
            DC resistance of the motor in Ohms. Valid only for "user_defined"
            type. Defaults to 0.0.
        
        I_no_load : float
            No-load current of the motor in amps. Valid only for "user_defined".

        gear_ratio : float, optional
            Gear ratio of the motor. A valud greater than unity will correspond
            to the output shaft of the motor turning slower than the motor itself.
            Defaults to 1.0.

        weight : float
            Weight of the motor in ounces. Valid only for "user_defined" type.
    """

    def __init__(self, **kwargs):

        # Get type
        self._type = kwargs.get("type", "database")

        if self._type == "database":
            
            # Get database info
            super().__init__()
            record = self.get_database_record("motor", **kwargs)

            # Set params
            self.Kv = float(record[3])
            self.Gr = float(record[4])
            self.I0 = float(record[6])
            self.R = float(record[5])
            self.name = record[1]
            self.manufacturer = record[2]
            self.weight = float(record[7])

        elif self._type == "user_defined":

            # Set params
            self.Kv = kwargs.get("Kv")
            self.Gr = kwargs.get("gear_ratio", 1.0)
            self.I0 = kwargs.get("I_no_load")
            self.R = kwargs.get("resistance", 0.0)
            self.name = kwargs.get("name", "user_motor")
            self.manufacturer = kwargs.get("manufacturer", "user")
            self.weight = kwargs.get("weight")


    def __str__(self):
        string = "Motor: {0}".format(self.name)
        string += "\n\tManufacturer: {0}".format(self.manufacturer)
        string += "\n\tKv: {0}".format(self.Kv)
        string += "\n\tWeight: {0}".format(self.weight)
        return string