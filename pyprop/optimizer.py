####################################################################
# plotDesignSpace.py
#
# Randomly searches possible combinations of components and determines viable flight times.
# Plots these flight times against component parameters. User can select specific points in
# the design space to view component parameters and thrust curves.
#
# The only argument is a .json file defining the search parameters, an example of which is
# given below (explanatory comments given within // //):
#
# ----sampleSearch.json----
#
# {
#     "computation":{
#         "units":1000, //Number of propulsion units to find in the design space.//
#         "processes":8, //Maximum number of processes to be used in parallel computation.//
#         "outlierStdDevs":5 //Number of standard deviations of the half-normal distribution within which designs are considered feasible.//
#     },
#     "condition":{
#         "altitude":0, //Flight altitude.//
#         "airspeed":10 //Flight cruise speed.//
#     },
#     "goal":{ //One and only one of these parameters must be specified. Set for cruise condition.//
#         "thrust":0, //Thrust required from the propulsion unit.//
#         "thrustToWeightRatio":0.3 //Thrust to weight ratio required (requires emptyWeight to be defined.//
#     },
#     "aircraft":{
#         "emptyWeight":1, //Weight of the aircraft minus the propulsion system.//
#         "components":{ //These parameters are optional, but only one for each component may be specified.//
#             "propeller":{
#                 "name":"",
#                 "manufacturer":""
#             },
#             "motor":{
#                 "name":"",
#                 "manufacturer":""
#             },
#             "esc":{
#                 "name":"",
#                 "manufacturer":""
#             },
#             "battery":{
#                 "name":"",
#                 "manufacturer":""
#             }
#         }
#     }
# }
#
# The component parameters in the .json file are all optional. Specifying a component name limits
# the search to propulsion units including that specific component. Specifying a component manufacturer
# limits the search to a single manufacturer for that component. Please note that some component
# manufacturers have very few components in our current database, and specifying this may limit the
# search more than desirable. Only one of these parameters may be specified for each component at most.
#
# Once the search is complete (i.e. the specified number of designs has been considered), the propulsion
# unit which has the longest flight time will be output to the terminal. A figure will also be displayed
# containing 6 plots which describe (in part) the design space. Each point on a plot represents a possible
# design; each design is reflected in each plot. The y axis of each plot is the flight time given by a
# design, and the x axis of each plot is a defining parameter of the design (currently: prop diameter, prop
# pitch, motor Kv constant, battery voltage, battery capacity, and total unit weight). The user may select
# any design in any one of the plots to see a plot of its thrust at various airspeeds and throttle settings.
# A corresponding plot of propeller speeds is also shown and all details of the design are printed to the
# terminal. Selecting a design will also highlight that design in each of the 6 plots, so that general
# patterns in the design space can be opserved.
#
# The developer is of the opinion that outliers should be ignored. Testing has shown that these arise from
# error in the component models and should not be trusted as feasible, high-performance designs. Realistic
# designs will be found closer to the main cluster of designs.
#
####################################################################

import matplotlib.pyplot as plt
import sqlite3 as sql
import supportClasses as s
import numpy as np
from random import randint,seed
import multiprocessing as mp
import math
import sys
import warnings
import datetime
import json
import os

dbFile = "Database/components.db"

#Defines what happens when the user picks a plotted point in the design space. Highlights that point and 
#plots that unit's thrust curves.
def on_pick(event):
    artist = event.artist
    fig = plt.figure(plt.get_fignums()[0])
    ax = fig.axes
    ind = int(event.ind[0])
    selUnit = units[ind]
    fig.suptitle("SELECTED Prop: "+str(selUnit.prop.name)+"  Motor: "+str(selUnit.motor.name)+"  Battery: "+str(selUnit.batt.name)+"  ESC: "+str(selUnit.esc.name))
    ax[0].plot(selUnit.prop.diameter,t_flight[ind],'o')
    ax[1].plot(selUnit.prop.pitch,t_flight[ind],'o')
    ax[2].plot(selUnit.motor.Kv,t_flight[ind],'o')
    ax[3].plot(selUnit.batt.V0,t_flight[ind],'o')
    ax[4].plot(selUnit.batt.cellCap,t_flight[ind],'o')
    ax[5].plot(selUnit.GetWeight()+W_frame,t_flight[ind],'o')
    ax[6].plot(throttles[ind],t_flight[ind],'o')
    selUnit.printInfo()
    print("Flight Time:",t_flight[ind],"min")
    if optimizeForRatio:
        print("    at {:4.2f}% throttle".format(selUnit.CalcCruiseThrottle(v_req,(selUnit.GetWeight()+W_frame)*R_tw_req)*100))
    else:
        print("    at {:4.2f}% throttle".format(selUnit.CalcCruiseThrottle(v_req,T_req)*100))
    selUnit.PlotThrustCurves(0,v_req*2+10,11,51)
    selUnit.prop.PlotCoefs()

#Defines a global database cursor giving all processes a connection to the database.
def setGlobalCursor():
    global dbcur
    dbcur = sql.connect(dbFile).cursor()
    seed(datetime.time.microsecond)

#Selects a propultion unit and calculates its flight time.
def getCombination(args):

    v_req = args[0]
    T = args[1]
    h = args[2]
    optimizeForRatio = args[3]
    W_frame = args[4]
    names = args[5]
    manufacturers = args[6]

    if optimizeForRatio:
        R_tw_req = T
    else:
        T_req = T

    t_flight_curr = None
    while t_flight_curr is None or math.isnan(t_flight_curr):

        #Fetch prop data
        prop = s.Propeller(dbcur,name=names[0],manufacturer=manufacturers[0])

        #Fetch motor data
        motor = s.Motor(dbcur,name=names[1],manufacturer=manufacturers[1])

        #Fetch ESC data
        esc = s.ESC(dbcur,name=names[2],manufacturer=manufacturers[2])

        #Fetch battery data
        batt = s.Battery(dbcur,name=names[3],manufacturer=manufacturers[3])

        if batt.R == 0 and esc.R == 0 and motor.R == 0:
            continue

        currUnit = s.PropulsionUnit(prop,motor,batt,esc,h)
        if optimizeForRatio:
            T_req = (currUnit.GetWeight()+W_frame)*R_tw_req
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            t_flight_curr = currUnit.CalcBattLife(v_req,T_req)
            thr_curr = currUnit.CalcCruiseThrottle(v_req,T_req)
    return t_flight_curr, thr_curr, currUnit

#----------------------BEGINNING OF COMPUTATION------------------------------------

if len(sys.argv) is not 2:
    raise RuntimeError("plotDesignSpace takes only one argument (the .json configuration filename)!")
configFile = sys.argv[1]
with open(configFile) as filename:
    settings = json.load(filename)

N_proc_max = settings["computation"]["processes"]
N_units = settings["computation"]["units"]
v_req = settings["condition"]["airspeed"]
h = settings["condition"]["airspeed"]
W_frame = settings["aircraft"]["emptyWeight"]

if settings["goal"]["thrust"] is 0:
    if settings["goal"]["thrustToWeightRatio"] is 0:
        raise RuntimeError("No goal specified!")
    optimizeForRatio = True
    thrustParam = settings["goal"]["thrustToWeightRatio"]
    R_tw_req = thrustParam
else:
    optimizeForRatio = False
    thrustParam = settings["goal"]["thrust"]
    T_req = thrustParam

print("Flight conditions: airspeed",v_req,"ft/s, altitude",h,"ft, airframe weight",W_frame,"lbs")
if optimizeForRatio:
    print("Optimizing for a thrust to weight ratio of",thrustParam)
else:
    print("Optimizing for a required thrust of",thrustParam)

names = []
manufacturers = []
print("Optimization constrained as follows:")
for component in settings["aircraft"]["components"]:
    name = settings["aircraft"]["components"][component]["name"]
    if len(name) is 0:
        name = None
    names.append(name)

    manufacturer = settings["aircraft"]["components"][component]["manufacturer"]
    if len(manufacturer) is 0:
        manufacturer = None
    manufacturers.append(manufacturer)

    if name is not None and manufacturer is not None:
        raise RuntimeError("Component: "+component+" is overconstrained!")

    print(component.title())
    if name is not None:
        print("Name:",name)
    elif manufacturer is not None:
        print("Manufacturer:",manufacturer)
    else:
        print("Not constrained.")

# Distribute work
with mp.Pool(processes=N_proc_max,initializer=setGlobalCursor,initargs=()) as pool:
    args = [(v_req,thrustParam,h,optimizeForRatio,W_frame,names,manufacturers) for i in range(N_units)]
    data = pool.map(getCombination,args)
sql.connect(dbFile).close()

t_flight,throttles,units = map(list,zip(*data))
t_flight = np.asarray(t_flight)
throttles = np.asarray(throttles)
units = np.asarray(units)

# Determine optimum
t_max = max(t_flight)
bestUnit = np.asscalar(units[np.where(t_flight==t_max)])
throttle_at_max = np.asscalar(throttles[np.where(t_flight==t_max)])
if optimizeForRatio:
    T_req = thrustParam*(bestUnit.GetWeight()+W_frame)
else:
    T_req = thrustParam

print("Maximum flight time found:",t_max,"min")
bestUnit.printInfo()
print("Throttle setting for max flight:",bestUnit.CalcCruiseThrottle(v_req,T_req))
print("Current draw:",bestUnit.I_motor,"A")

# Plot design space
plt.ion()
fig,((ax1,ax2,ax3,ax4),(ax5,ax6,ax7,ax8)) = plt.subplots(nrows=2,ncols=4)
fig.suptitle("OPTIMUM Prop: "+str(bestUnit.prop.name)+"  Motor: "+str(bestUnit.motor.name)+"  Battery: "+str(bestUnit.batt.name)+"  ESC: "+str(bestUnit.esc.name))

ax1.plot([units[i].prop.diameter for i in range(N_units)],t_flight,'b*',picker=3)
ax1.plot(bestUnit.prop.diameter,t_max,'r*')
ax1.set_xlabel("Prop Diameter [in]")
ax1.set_ylabel("Flight Time [min]")

ax2.plot([units[i].prop.pitch for i in range(N_units)],t_flight,'b*',picker=3)
ax2.plot(bestUnit.prop.pitch,t_max,'r*')
ax2.set_xlabel("Prop Pitch [in]")
ax2.set_ylabel("Flight Time [min]")

ax3.plot([units[i].motor.Kv for i in range(N_units)],t_flight,'b*',picker=3)
ax3.plot(bestUnit.motor.Kv,t_max,'r*')
ax3.set_xlabel("Motor Kv [rpm/V]")
ax3.set_ylabel("Flight Time [min]")

ax4.plot([units[i].batt.V0 for i in range(N_units)],t_flight,'b*',picker=3)
ax4.plot(bestUnit.batt.V0,t_max,'r*')
ax4.set_xlabel("Battery Voltage [V]")
ax4.set_ylabel("Flight Time [min]")

ax5.plot([units[i].batt.cellCap for i in range(N_units)],t_flight,'b*',picker=3)
ax5.plot(bestUnit.batt.cellCap,t_max,'r*')
ax5.set_xlabel("Cell Capacity [mAh]")
ax5.set_ylabel("Flight Time [min]")

ax6.plot([units[i].GetWeight()+W_frame for i in range(N_units)],t_flight,'b*',picker=3)
ax6.plot(bestUnit.GetWeight()+W_frame,t_max,'r*')
ax6.set_xlabel("Total Unit Weight [lb]")
ax6.set_ylabel("Flight Time [min]")

ax7.plot(throttles,t_flight,'b*',picker=3)
ax7.plot(throttle_at_max,t_max,'r*')
ax7.set_xlabel("Throttle Setting at Max Flight Time")
ax7.set_ylabel("Flight Time [min]")

fig.canvas.mpl_connect('pick_event',on_pick)
plt.show(block=True)
plt.ioff()
