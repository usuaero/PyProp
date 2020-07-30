#Used to organize motors into the database
#There are three sources of motor parameters: DriveCalc, MotoCalc, csv from Josh
import sqlite3 as sql
import numpy as np
import os
from os import path
import csv
from dbfread import DBF

def isNum(x):
    try:
        float(x)
        return True
    except ValueError:
        return False

#Kv should be in rpms/V
#current should be in amps
#resistance should be in ohms
#weight should be in ounces

connection = sql.connect("components.db")
motoCalcFile = DBF("Motors/MOTOR8.DBF") 
cursor = connection.cursor()
cursor.execute("drop table Motors")# The database is refreshed every time

cursor.execute("""CREATE TABLE Motors (id INTEGER PRIMARY KEY, 
                                      name VARCHAR(40), 
                                      manufacturer VARCHAR(40),
                                      kv FLOAT, 
                                      gear_ratio FLOAT default 1.0, 
                                      resistance FLOAT, 
                                      no_load_current FLOAT default 0.0,
                                      weight FLOAT default -1.0);""")

print("Reading MotoCalc Database")
for record in motoCalcFile:
    
    print(record)
    if record["MOTORSTANT"] == 0:
        continue
    if "Car Motor" in record["MOTORNAME"]:
        manufacturer = "UNKNOWN"
    else:
        manufacturer = str(record["MOTORNAME"].split(" ")[0]).upper()
    
    formatStr = """INSERT INTO Motors (name, manufacturer, kv, resistance, no_load_current, weight) VALUES ("{motorName}", "{motorManu}", {motorKv}, {motorResistance}, {motorNoLoadCurrent}, {motorWeight});"""
    command = formatStr.format(motorName = record['MOTORNAME'], motorManu = manufacturer, motorKv = record['MOTORSTANT'], motorResistance = record['ARMATTANCE'], motorNoLoadCurrent = record['IDLECRRENT'], motorWeight = record['WEIGHT'])
    cursor.execute(command)

print("Reading Database after MotoCalc")
cursor.execute("SELECT * FROM Motors")
result = cursor.fetchall()
for r in result:
    print(r)

print("Reading DriveCalc database")

inDatabaseFile = os.getcwd() + "/Motors/DCbase.dcd"
inConnection = sql.connect(inDatabaseFile)
inCursor = inConnection.cursor()

inCursor.execute("SELECT * FROM motors")
motors = inCursor.fetchall()

for motor in motors:
    if motor[16] == 0:
        continue
    name = motor[2].replace("\""," ")
    manufacturer = "UNKNOWN"
    for word in name.split(" "):
        print(word)
        if len(word) > 0 and word[0].isalpha():
            manufacturer = word.upper()
            break

    formatStr = """INSERT INTO Motors (name, manufacturer, kv, resistance, weight, gear_ratio) VALUES ("{motorName}", "{manu}", {motorKv}, {motorResistance}, {motorWeight}, {motorGearRatio});"""
    command = formatStr.format(motorName = name, manu = manufacturer, motorKv = motor[16], motorResistance = motor[17], motorWeight = motor[13]*0.035274, motorGearRatio = motor[14])
    cursor.execute(command)
    
print("Reading Database after DriveCalc")
cursor.execute("SELECT * FROM Motors")
result = cursor.fetchall()
for r in result:
    print(r)
    
print(cursor.description)
inCursor.close()

print("Reading motors from csv")
motorCsv = open(os.getcwd()+"/Motors/motor.csv")
motorCsvReader = csv.reader(motorCsv,delimiter=',')
firstLine = True
for row in motorCsvReader:
    if firstLine:
        firstLine = False
        continue
    print(row)
    name = " ".join(row[0:1])
    kv = row[2]
    resistance = row[3]
    I0 = row[4]

    formatStr = """INSERT INTO Motors (name, manufacturer, kv, resistance, no_load_current) VALUES ("{motorName}", "{manufacturer}", {motorKv}, {motorResistance}, {motorI0});"""
    command = formatStr.format(motorName = name, manufacturer = row[0].upper(), motorKv = kv, motorResistance = resistance, motorI0 = I0)
    cursor.execute(command)

print("Reading Database after CSV")
cursor.execute("SELECT * FROM Motors")
result = cursor.fetchall()
for r in result:
    print(r)

#Code for cleaning newly created database:
cursor.execute("delete from Motors where weight = -1.0") # This deletes motors whose weight was not updated.

print("Reading Database after cleaning")
cursor.execute("SELECT * FROM Motors")
result = cursor.fetchall()
for r in result:
    print(r)

connection.commit()
cursor.close()
