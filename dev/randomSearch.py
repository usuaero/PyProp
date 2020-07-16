import sqlite3 as sql
import supportClasses as s
import numpy as np
from random import randint
import multiprocessing as mp
import math

dbFile = "Database/components.db"

def setGlobalCursor():
    global dbcur
    dbcur = sql.connect(dbFile).cursor()

def getCombination(args):

    reqCruiseSpeed = args[0]
    reqThrust = args[1]
    altitude = args[2]

    currFlightTime = None
    while currFlightTime is None or math.isnan(currFlightTime):

        #Fetch prop data
        prop = s.Propeller(dbcur)

        #Fetch motor data
        motor = s.Motor(dbcur)

        #Fetch ESC data
        esc = s.ESC(dbcur)

        #Fetch battery data
        batt = s.Battery(dbcur)

        if batt.R == 0 and esc.R == 0 and motor.R == 0:
            continue

        currUnit = s.PropulsionUnit(prop,motor,batt,esc,altitude)
        currFlightTime = currUnit.CalcBattLife(reqCruiseSpeed,reqThrust)
    return currFlightTime, currUnit


combinations = 100
maxProcesses = 8

# Will attempt to maximize flight time based on these parameters
# Randomly selects combinations of components
reqCruiseSpeed = 10
reqThrust = 0.5
altitude = 2000

with mp.Pool(processes=maxProcesses,initializer=setGlobalCursor,initargs=()) as pool:
    args = [(reqCruiseSpeed,reqThrust,altitude,i) for i in range(combinations)]
    data = pool.map(getCombination,args)
sql.connect(dbFile).close()

flightTimes,units = map(list,zip(*data))

bestFlightTime = max(flightTimes)
bestUnit = units[flightTimes.index(bestFlightTime)]


print("Maximum flight time found:",bestFlightTime,"min")
print("Prop:",bestUnit.prop.name)
print("Motor:",bestUnit.motor.name,"(Kv =",bestUnit.motor.Kv,")")
print("Battery:",bestUnit.batt.name,"(Capacity =",bestUnit.batt.cellCap,", Voltage =",bestUnit.batt.V0,")")
print("ESC:",bestUnit.esc.name)
print("Throttle setting for max flight:",bestUnit.CalcCruiseThrottle(reqCruiseSpeed,reqThrust))
print("Current draw:",bestUnit.Im,"A")
bestUnit.PlotThrustCurves(30,11,51)

