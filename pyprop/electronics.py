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

    chemistry : str
        Chemistry type of the battery.

    cell_arrangement : tuple
        The number of series and parallel cells in the battery. The first
        element should be the number of series cells and the second element
        should be the number of parallel cells. Defaults to (0, 0).
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
        self.chemistry = kwargs.get("chemistry", "NULL")
        self.S, self.P = kwargs.get("cell_arrangement", (0, 0))


    def __str__(self):
        string = "Battery: {0}".format(self.name)
        string += "\n\tManufacturer: {0}".format(self.manufacturer)
        string += "\n\tCapacity: {0} mAh".format(self.capacity)
        string += "\n\tCells: {0}S{1}P".format(self.S, self.P)
        string += "\n\tVoltage: {0} V".format(self.V0)
        string += "\n\tWeight: {0} oz".format(self.weight)
        return string


    def write_to_database(self):
        """Saves this component's information to the user database.
        Doing this allows the component to be used in optimization
        schemes that query the database.
        
        This function will replace duplicate components (i.e.
        same name). The database is stored in your Python site-packages
        folder. The name is ```user_components.db```.
        """

        # Locate database file
        db_file = os.path.join(os.path.dirname(__file__), "user_components.db")

        # Connect to database
        connection = sql.connect(db_file)
        cursor = connection.cursor()

        # Check if the table has been created
        try:
            cursor.execute("select * from Batteries")
        except:
            cursor.execute("""create table Batteries (id INTEGER PRIMARY KEY, 
                                                  Name VARCHAR(40), 
                                                  manufacturer VARCHAR,
                                                  Imax FLOAT, 
                                                  Capacity FLOAT, 
                                                  Weight FLOAT,
                                                  Ri FLOAT,
                                                  Volt FLOAT,
                                                  Chem VARCHAR,
                                                  S INTEGER default 0,
                                                  P INTEGER default 0);""")

        # Check for duplicate
        cursor.execute("delete from Batteries where Name = '{0}';".format(self.name))

        # Store component
        command = """insert into Batteries (Name, manufacturer, Imax, Capacity, Weight, Ri, Volt, Chem, S, P)
                  values ("{0}", "{1}", {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9});""".format(self.name,
                  self.manufacturer.replace(" ", "_"), self.I_max/self.P, self.capacity/self.P, self.weight/(self.S*self.P),
                  self.R*self.P/self.S, self.V0/self.S, self.chemistry, self.S, self.P)

        cursor.execute(command)
        cursor.close()
        connection.commit()
        connection.close()


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


    def write_to_database(self):
        """Saves this component's information to the user database.
        Doing this allows the component to be used in optimization
        schemes that query the database.
        
        This function will replace duplicate components 
        (i.e. same name). The database is stored in your Python site-packages
        folder. The name is ```user_components.db```.
        """

        # Locate database file
        db_file = os.path.join(os.path.dirname(__file__), "user_components.db")

        # Connect to database
        connection = sql.connect(db_file)
        cursor = connection.cursor()

        # Check if the table has been created
        try:
            cursor.execute("select * from ESCs")
        except:
            cursor.execute("""create table ESCs (id INTEGER PRIMARY KEY, 
                                                  Name VARCHAR(40), 
                                                  manufacturer VARCHAR,
                                                  Imax FLOAT, 
                                                  Ipeak FLOAT, 
                                                  Weight FLOAT,
                                                  Ri FLOAT);""")

        # Check for duplicate
        cursor.execute("delete from ESCs where Name = '{0}';".format(self.name))

        # Store component
        command = """insert into ESCs (Name, manufacturer, Imax, Ipeak, Weight, Ri)
                  values ("{0}", "{1}", {2}, {3}, {4}, {5});""".format(self.name,
                  self.manufacturer.replace(" ", "_"), self.I_max, self.I_max, self.weight, self.R)

        cursor.execute(command)
        cursor.close()
        connection.commit()
        connection.close()

        
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


    def write_to_database(self):
        """Saves this component's information to the user database.
        Doing this allows the component to be used in optimization
        schemes that query the database.
        
        This function will replace duplicate components
        (i.e. same name). The database is stored in your Python site-packages
        folder. The name is ```user_components.db```.
        """

        # Locate database file
        db_file = os.path.join(os.path.dirname(__file__), "user_components.db")

        # Connect to database
        connection = sql.connect(db_file)
        cursor = connection.cursor()

        # Check if the table has been created
        try:
            cursor.execute("select * from Motors")
        except:
            cursor.execute("""create table Motors (id INTEGER PRIMARY KEY, 
                                                  Name VARCHAR(40), 
                                                  manufacturer VARCHAR,
                                                  kv FLOAT, 
                                                  gear_ratio FLOAT default 1.0, 
                                                  no_load_current FLOAT,
                                                  weight FLOAT,
                                                  resistance FLOAT);""")

        # Check for duplicate
        cursor.execute("delete from Motors where Name = '{0}';".format(self.name))

        # Store component
        command = """insert into Motors (Name, manufacturer, kv, gear_ratio, no_load_current, weight,
                  resistance) values ("{0}", "{1}", {2}, {3}, {4}, {5}, {6});""".format(self.name,
                  self.manufacturer.replace(" ", "_"), self.Kv, self.Gr, self.I0, self.weight, self.R)

        cursor.execute(command)
        cursor.close()
        connection.commit()
        connection.close()