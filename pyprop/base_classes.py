"""Defines base classes for the module."""

import os

import sqlite3 as sql
import numpy as np

from .exceptions import DatabaseRecordNotFoundError

class DatabaseComponent:
    """A component defined in the database."""

    def __init__(self):

        # Get database file
        self._db_file = os.path.join(os.path.dirname(__file__), "components.db")

        # Set up params
        self._table_names = {
            "battery" : "Batteries",
            "ESC" : "ESCs",
            "motor" : "Motors",
            "prop" : "Props"
        }


    def get_database_record(self, component_type, **kwargs):
        """Extracts a record from the database."""

        # Get kwargs
        name = kwargs.get("name", None)
        manufacturer = kwargs.get("manufacturer", None)
        dbid = kwargs.get("dbid", None)
        capacity = kwargs.get("capacity", None)
        I_max = kwargs.get("I_max", None)
        Kv = kwargs.get("Kv", None)
        diameter = kwargs.get("diameter", None)
        pitch = kwargs.get("pitch", None)

        # Get database location and connection
        with sql.connect(self._db_file) as conn:
            db_cur = conn.cursor()

            # Format command generically
            table_name = self._table_names[component_type]
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
            if component_type == "prop" and diameter is not None:
                command = command+" order by abs("+str(diameter)+"-Diameter)"
            if component_type == "prop" and pitch is not None:
                command = command+" order by abs("+str(pitch)+"-Pitch)"

            # Get in random order
            command = command+" order by RANDOM() limit 1"

            # Get record
            db_cur.execute(command)
            try:
                record = np.asarray(db_cur.fetchall())[0]
            except IndexError:
                raise DatabaseRecordNotFoundError(command)

        return record