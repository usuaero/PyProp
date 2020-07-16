import sqlite3 as sql
import numpy as np
import supportClasses as s

props = ["gwsdd_2.5x0.8"]
motors = ["Mabuchi RC-280RA-2865"]
batteries = ["E-Tec 1200 Li-Poly"]
numCells = 6
escs = ["Kontronik Beat 60-8-24"]
altitude = 0
thrusts = [0.3]
speeds = [10]
weights = [1]

#Open database and read records from database
db = sql.connect("Database/components.db")
dbcur = db.cursor()

for prop,motor,battery,esc,thrust,cruiseSpeed,weight in zip(props,motors,batteries,escs,thrusts,speeds,weights):

    #Fetch prop data
    propC = s.Propeller(dbcur,name=prop)
    propC.printInfo()
    propC.PlotCoefs()

    #Fetch motor data
    motorC = s.Motor(dbcur,name=motor)
    motorC.printInfo()

    #Fetch battery data
    battC = s.Battery(dbcur)
    battC.printInfo()

    #Fetch ESC data
    escC = s.ESC(dbcur,name=esc)
    escC.printInfo()

    test = s.PropulsionUnit(propC, motorC, battC, escC, altitude)
    print("Initialization complete. Plotting thrust curves.")
    minAirspeed = 0
    maxAirspeed = 30
    numVelocities = 11
    numThrottles = 1001
    thrust = test.GetWeight()+weight*thrust
    test.PlotThrustCurves(minAirspeed,maxAirspeed, numVelocities, numThrottles)
    print("Flight time (thrust=",thrust,", cruiseSpeed=",cruiseSpeed,")",test.CalcBattLife(thrust,cruiseSpeed),"min")

db.close()
