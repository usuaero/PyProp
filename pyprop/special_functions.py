"""Special functions used in the PyProp package."""

import os

import sqlite3 as sql
import numpy as np

from .propellers import DatabaseFitProp, DatabaseDataProp
from .electronics import Battery, Motor, ESC
from .exceptions import DatabaseRecordNotFoundError

from random import randint


def create_component_from_database(**kwargs):
    """Fetches the necessary information from the database to assemble a component.

    This function will randomly select a component that satisfies all the parameters
    given. So if none are given except what type of component is needed, that
    component will be selected at random from the database. If the manufacturer is 
    given, then a component made by that manufacturer will be selected at random. If
    the name is given, then that exact component will be selected.

    For numerical values, such as Kv for motors or capacity for batteries, the 
    component with the value of this closest to the value given will be
    selected. In case of a tie, the selection is made randomly.

    Alternatively, a two-element list can be given, as in

        "Kv" : [300.0, 500.0]

    This specifies that a motor with a Kv between these two values should be chosen.

    The following components are affected by the following parameters:

        prop: pitch, diameter, prop_type
        motor: Kv
        battery: capacity, num_cells
        ESC: I_max

    Only one exact numeric parameter may be given (i.e. you are not allowed to specify
    { "pitch" : 12, "diameter" : 3}), as doing otherwise leaves the program with an
    unsolvable problem, namely, which parameter am I to match more closely? Yes, there
    are cases where both can be matched exactly, but this is not guaranteed. Multiple
    parameter ranges may be specified (i.e. { "pitch" : [3, 5], "diameter" : [12, 14] }).

    Parameters
    ----------
    component : str
        Can be "prop", "battery", "ESC", or "motor".

    name : str, optional
        Name of the component. This should be exactly as the name appears in 
        pyprop/components.db.

    manufacturer : str, optional
        Component manufacturer.

    dbid : int, optional
        ID of the component in the database.

    prop_type : str, optional
        Type of the propeller model to use. Can either be "fit" or "data". "data" type
        props will calculate the torque and thrust coefficients based off of linearly
        interpolating experimental data. "fit" type props calculate the coefficients
        based off of fits of the same data. Defaults to "data".

    capacity : float or list, optional
        Battery capacity in mAh.

    num_cells : int or list, optional
        Number of series cells in the battery. Note that batteries are stored
        as single cells in the database, so all batteries are available with
        an arbitrary number of cells (this may not reflect commercial
        availability). If not specified, this defaults to a random integer
        between 1 and 8.

    I_max : float or list, optional
        Maximum ESC current draw.

    Kv : float or list, optional
        Motor speed constant.

    diameter : float or list, optional
        Propeller diameter in inches.

    pitch : float or list, optional
        Propeller pitch in inches.

    Returns
    -------
    Motor, Battery, ESC, DatabaseDataProp, or DatabaseFitProp
        The component specified, pulled from the database.
    """

    # Locate database file
    db_file = os.path.join(os.path.dirname(__file__), "components.db")

    # Set up params
    table_names = {
        "battery" : "Batteries",
        "ESC" : "ESCs",
        "motor" : "Motors"
    }

    # Determine where to get the info from
    component_type = kwargs.get("component")

    # Get kwargs
    name = kwargs.get("name", None)
    manufacturer = kwargs.get("manufacturer", None)
    dbid = kwargs.get("dbid", None)
    prop_type = kwargs.get("prop_type", "data")
    capacity = kwargs.get("capacity", None)
    I_max = kwargs.get("I_max", None)
    Kv = kwargs.get("Kv", None)
    diameter = kwargs.get("diameter", None)
    pitch = kwargs.get("pitch", None)
    num_cells = kwargs.get("num_cells", [1, 8])
    if isinstance(num_cells, list):
        num_cells = randint(num_cells[0], num_cells[1])

    # Get database location and connection
    with sql.connect(db_file) as conn:
        db_cur = conn.cursor()

        # Determine database table
        if component_type == "prop":
            if prop_type == "data":
                table_name = "DataProps"
            elif prop_type == "fit":
                table_name = "Props"
            else:
                raise IOError("{0} is not a valid prop type.".format(prop_type))
        else:
            table_name = table_names[component_type]

        # Format command generically
        command = "select * from {0}".format(table_name)
        if name is not None:
            if manufacturer is not None or dbid is not None:
                raise ValueError("Too many {0} parameters specified.".format(component_type))
            command = command+" where Name = '"+name+"'"
        elif manufacturer is not None:
            if dbid is not None:
                raise ValueError("Too many {0} parameters specified.".format(component_type))
            command = command+" where manufacturer = '"+manufacturer+"'"
        elif dbid is not None:
            command = command+" where id = "+str(dbid)

        # Add component-specific commands
        # Battery
        if component_type == "battery":
            if capacity is not None:
                if isinstance(capacity, list):
                    if "where" in command:
                        command = command+" and Capacity between {0} and {1} order by RANDOM() limit 1".format(min(capacity), max(capacity))
                    else:
                        command = command+" where Capacity between {0} and {1} order by RANDOM() limit 1".format(min(capacity), max(capacity))
                else:
                    command = command+" order by abs({0}-Capacity)".format(capacity)

        # ESC
        elif component_type == "ESC":
            if I_max is not None:
                if isinstance(I_max, list):
                    if "where" in command:
                        command = command+" and I_motorax between {0} and {1} order by RANDOM() limit 1".format(min(I_max), max(I_max))
                    else:
                        command = command+" where I_motorax between {0} and {1} order by RANDOM() limit 1".format(min(I_max), max(I_max))
                else:
                    command = command+" order by abs({0}-I_motorax)".format(I_max)

        # Motor
        elif component_type == "motor":
            if Kv is not None:
                if isinstance(Kv, list):
                    if "where" in command:
                        command = command+" and kv between {0} and {1} order by RANDOM() limit 1".format(min(Kv), max(Kv))
                    else:
                        command = command+" where kv between {0} and {1} order by RANDOM() limit 1".format(min(Kv), max(Kv))
                else:
                    command = command+" order by abs({0}-kv)".format(Kv)

        # Prop
        elif component_type == "prop":
            specific = False

            # Range of diameters
            if diameter is not None and isinstance(diameter, list):
                if "where" in command:
                    command = command+" and Diameter between {0} and {1}".format(min(diameter), max(diameter))
                else:
                    command = command+" where Diameter between {0} and {1}".format(min(diameter), max(diameter))

            # Range of pitches
            if pitch is not None and isinstance(pitch, list):
                if "where" in command:
                    command = command+" and Pitch between {0} and {1}".format(min(pitch), max(pitch))
                else:
                    command = command+" where Pitch between {0} and {1}".format(min(pitch), max(pitch))

            # Specific diameter
            if diameter is not None and not isinstance(diameter, list):
                specific = True
                command = command+" order by abs({0}-Diameter)".format(diameter)

            # Specific pitch
            if pitch is not None and not isinstance(pitch, list):
                if specific:
                    raise IOError("Exact values for both pitch and diameter may not be given when selecting a propeller from the database.")
                else:
                    specific = True
                    command = command+" order by abs({0}-Pitch)".format(pitch)

            if not specific:
                command = command+" order by RANDOM() limit 1"

        if "order" not in command:
            # Get in random order
            command = command+" order by RANDOM() limit 1"

        # Get record
        db_cur.execute(command)
        try:
            record = np.asarray(db_cur.fetchall())[0]
        except IndexError:
            raise DatabaseRecordNotFoundError(command)

    # Initialize object
    if component_type == "motor":
        obj = Motor(name=record[1],
                    manufacturer=record[2],
                    Kv=float(record[3]),
                    Gr=float(record[4]),
                    I_no_load=float(record[6]),
                    resistance=float(record[5]),
                    weight=float(record[7]))

    elif component_type == "battery":
        obj = Battery(name=record[1],
                      manufacturer=record[2],
                      num_cells=num_cells,
                      capacity=float(record[4]),
                      voltage=float(record[7])*num_cells,
                      resistance=float(record[6])*num_cells,
                      I_max=float(record[3]),
                      weight=float(record[5])*num_cells)

    elif component_type == "ESC":
        obj = ESC(name=record[1],
                  manufacturer=record[2],
                  resistance=float(record[6]),
                  I_max=float(record[3]),
                  weight=float(record[5]))

    elif component_type == "prop":
        if prop_type == "data":
            obj = DatabaseDataProp(record)
        else:
            obj = DatabaseFitProp(record)
        
    return obj