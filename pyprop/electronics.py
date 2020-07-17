"""Classes defining electronic components for the propulsion unit."""

import os

import numpy as np
import sqlite3 as sql

from random import randint

from .base_classes import DatabaseComponent
    
#A class that defines a battery
class Battery(DatabaseComponent):

    #Initialize the class from database
    def __init__(self, **kwargs):
        """A battery.

        Parameters
        ----------
        type : str, optional
            Can be "database" or "user_defined". Defaults to "database".

        """
        super().__init__()

        # Get type
        self._type = kwargs.get("type", "database")

        if self._type == "database":
            
            # Get database info
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
            pass

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

#A class that defines an ESC (Electronic Speed Controller)
class ESC(DatabaseComponent):

    #Initialization of the class from database
    def __init__(self, **kwargs):
        super().__init__()

        # Get database record
        record = self.get_database_record("ESC", **kwargs)

        self.R = float(record[6])
        self.name = record[1]
        self.manufacturer = record[2]
        self.I_max = float(record[3])
        self.weight = float(record[5])

    def __str__(self):
        string = "ESC: {0}".format(self.name)
        string += "\n\tManufacturer: {0}".format(self.manufacturer)
        string += "\n\tMax Current: {0}".format(self.I_max)
        string += "\n\tWeight: {0}".format(self.weight)
        return string
        
#A class that defines an electric motor.
class Motor(DatabaseComponent):

    #Initialization of the class from the database
    def __init__(self, **kwargs):
        super().__init__()

        # Get database record
        record = self.get_database_record("motor", **kwargs)

        self.Kv = float(record[3])
        self.Gr = float(record[4])
        self.I0 = float(record[6])
        self.R = float(record[5])
        self.name = record[1]
        self.manufacturer = record[2]
        self.weight = float(record[7])

    def __str__(self):
        string = "Motor: {0}".format(self.name)
        string += "\n\tManufacturer: {0}".format(self.manufacturer)
        string += "\n\tKv: {0}".format(self.Kv)
        string += "\n\tWeight: {0}".format(self.weight)
        return string