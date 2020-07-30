#Used to create motor files for the optimization database
import sqlite3 as sql
import csv

def isNum(x):
    try:
        float(x)
        return True
    except ValueError:
        return False
    
motorFile = open("C:/Program Files (x86)/MotoCalc 8/Initial/MOTOR8.DBF")
databaseFile = "C:/Users/Cory/Google Drive/AeroLab/PropulsionUnitOptimization/components.db"

connection = sql.connect(databaseFile)
cursor = connection.cursor()

#cursor.execute("""CREATE TABLE motor (id INTEGER PRIMARY KEY, 
#                                      name VARCHAR(40), 
#                                      kv FLOAT, 
#                                      gear_ratio FLOAT, 
#                                      resistance FLOAT, 
#                                      no_load_current FLOAT,
#                                      weight FLOAT);""")

print("Reading MotoCalc Database")
firstLine = True
for line in motorFile:
    
    if firstLine:
        firstLine = False
        continue
    
    entries = line
    
print(entries)
totalLength = len(entries)
entryLength = len(" Neu F3A-1 1513/2Y                          1300.00   1.50000   0.01500  18.00000TF")
numMotors = len(entries)//entryLength
print(numMotors)

motors = [entries[i:i+entryLength] for i in range(0, totalLength, entryLength)]
motors.pop()

for motor in motors:
    print(motor)
    name = motor[0:40].replace("  ", "")
    if isNum(motor[41:51]):
        Kv = float(motor[41:51])
    else:
        Kv = None
    if isNum(motor[51:62]):
        I0 = float(motor[51:62])
    else:
        I0 = None
    if isNum(motor[63:73]):
        R = float(motor[63:73])
    else:
        R = None
    if isNum(motor[73:81]):
        weight = float(motor[73:81])
    else:
        weight = None
    
    print([name, Kv, I0, R, weight])

    formatStr = """INSERT INTO Motors (name, kv, resistance, no_load_current, weight) VALUES ("{motorName}", "{motorKv}", "{motorResistance}", "{motorNoLoadCurrent}", "{motorWeight}");"""
    command = formatStr.format(motorName = name, motorKv = Kv, motorResistance = R, motorNoLoadCurrent = I0, motorWeight = weight)
    cursor.execute(command)

print("Reading Database after CSV")
cursor.execute("SELECT * FROM Motors")
result = cursor.fetchall()
for r in result:
    print(r)

connection.commit()
connection.close()

#outDatabaseFile = "motors.db"
#outConnection = sql.connect(outDatabaseFile)
#outCursor = outConnection.cursor()
#
#inDatabaseFile = "DCbase.dcd"
#inConnection = sql.connect(inDatabaseFile)
#inCursor = inConnection.cursor()
#
#inCursor.execute("SELECT * FROM motors")
#motors = inCursor.fetchall()
#
#for motor in motors:
#
#    formatStr = """INSERT INTO motor (name, kv, resistance, weight, gear_ratio) VALUES ("{motorName}", "{motorKv}", "{motorResistance}", "{motorWeight}", "{motorGearRatio}");"""
#    command = formatStr.format(motorName = str(motor[2]).replace("\"", " "), motorKv = str(motor[16]), motorResistance = str(motor[17]), motorWeight = str(motor[13]), motorGearRatio = str(motor[14]))
#    outCursor.execute(command)
#    
#print("Reading Destination Database")
#outCursor.execute("SELECT * FROM motor")
#result = outCursor.fetchall()
#for r in result:
#    print(r)
#    
#print(outCursor.description)
#
#outConnection.commit()
#outConnection.close()
#inConnection.close()
