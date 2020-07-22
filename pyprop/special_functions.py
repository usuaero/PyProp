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

    Parameters
    ----------
    component : str
        Can be "fit_prop", "data_prop", "battery", "ESC", or "motor".

    Returns
    -------
    Motor, Battery, ESC, DatabaseFitProp, or DatabaseDataProp
        The component specified, pulled from the database.
    """

    # Locate database file
    db_file = os.path.join(os.path.dirname(__file__), "components.db")

    # Set up params
    table_names = {
        "battery" : "Batteries",
        "ESC" : "ESCs",
        "motor" : "Motors",
        "fit_prop" : "Props"
    }

    # Determine where to get the info from
    component_type = kwargs.get("component")

    # Raw data prop
    if component_type == "data_prop":
        pass

    # Other database components
    else:

        # Get kwargs
        name = kwargs.get("name", None)
        manufacturer = kwargs.get("manufacturer", None)
        dbid = kwargs.get("dbid", None)
        capacity = kwargs.get("capacity", None)
        I_max = kwargs.get("I_max", None)
        Kv = kwargs.get("Kv", None)
        diameter = kwargs.get("diameter", None)
        pitch = kwargs.get("pitch", None)
        num_cells = kwargs.get("num_cells", randint(1, 8))

        # Get database location and connection
        with sql.connect(db_file) as conn:
            db_cur = conn.cursor()

            # Format command generically
            table_name = table_names[component_type]
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
            if component_type == "battery" and capacity is not None:
                command = command+" order by abs("+str(capacity)+"-Capacity)"
            if component_type == "ESC" and I_max is not None:
                command = command+" order by abs("+str(I_max)+"-I_motorax)"
            if component_type == "motor" and Kv is not None:
                command = command+" order by abs("+str(Kv)+"-kv)"
            if component_type == "fit_prop" and diameter is not None:
                command = command+" order by abs("+str(diameter)+"-Diameter)"
            if component_type == "fit_prop" and pitch is not None:
                command = command+" order by abs("+str(pitch)+"-Pitch)"

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

        elif component_type == "fit_prop":
            obj = DatabaseFitProp(record)
        
    return obj