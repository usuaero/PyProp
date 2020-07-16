#used to organize ESCs from various sources into the component database
import os
import sqlite3 as sql
from dbfread import DBF

def isNum(x):
    try:
        float(x)
        return True
    except ValueError:
        return False
    except TypeError:
        return False

databaseFile = os.getcwd() + "/components.db"
connection = sql.connect(databaseFile)
cursor = connection.cursor()

cursor.execute("drop table ESCs")
cursor.execute("""create table ESCs (id INTEGER PRIMARY KEY, 
                                      Name VARCHAR(40), 
                                      manufacturer VARCHAR,
                                      Imax FLOAT, 
                                      Ipeak FLOAT, 
                                      Weight FLOAT,
                                      Ri FLOAT);""")

print("Reading MotoCalc Database")
escFilePath = os.getcwd() + "/ESCs/ESC8.DBF"
escFile = DBF(escFilePath)

for record in escFile:
    print(record)
    if record["MAXCURRENT"] == 0 or record["MAXCURRENT"] == None:
        continue

    formatStr = """INSERT INTO ESCs (Name, manufacturer, Weight, Imax, Ri) VALUES ("{name}", "{manu}", {weight},  {iMax}, {Ri});"""
    command = formatStr.format(name = record["ESCNAME"].strip(), manu = record["ESCNAME"].split(" " )[0].upper(), weight = record["WEIGHT"], iMax = record["MAXCURRENT"], Ri = record["RESISTANCE"])
    cursor.execute(command)

print("Reading Database after MotoCalc")
cursor.execute("SELECT * FROM ESCs")
result = cursor.fetchall()
for r in result:
    print(r)

print("Reading DriveCalc database")

inDatabaseFile = os.getcwd() + "/ESCs/DCbase.dcd"
inConnection = sql.connect(inDatabaseFile)
inCursor = inConnection.cursor()

inCursor.execute("SELECT * FROM ESC")
escs = inCursor.fetchall()

for esc in escs:

    if esc[4] == 0 or esc[4] == None:
        continue
    formatStr = """INSERT INTO ESCs (Name, manufacturer, Imax, Ipeak, Weight, Ri) VALUES ("{name}", "{manu}", {iMax}, {iPeak}, {weight}, {res});"""
    command = formatStr.format(name = esc[2].strip(), manu = esc[2].split(" ")[0].upper(), iMax = esc[4], iPeak = esc[5], weight = esc[7]*0.035274, res  = esc[6])
    cursor.execute(command)
    
print("Reading Database after DriveCalc")
cursor.execute("SELECT * FROM ESCs")
result = cursor.fetchall()
for r in result:
    print(r)

inCursor.close()
connection.commit()
connection.close()
