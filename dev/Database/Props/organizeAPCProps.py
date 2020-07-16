#Used to organize the APC prop database into subdirectories
import os
from os import path
import shutil

propDataFolder = "C:/Users/Cory/Google Drive/AeroLab/PropulsionUnitOptimization/Props/PERFILES_WEB/PERFILES2"
propDestFolder = "C:/Users/Cory/Google Drive/AeroLab/PropulsionUnitOptimization/Props"

for filename in os.listdir(path.join(propDataFolder)):
    
    dataFile = open(propDataFolder + "/" + filename)
    firstLine = dataFile.readline().split()
    dataFile.close()
    
    propName = firstLine[0]
    identifier = ""
    newPropName = ""
    prevDigit = ""
    
    xReached = False
    decimaled = False
    parameterEndReached = False
    for i in range(len(propName)):
        
        if (not (propName[i].isdigit() or propName[i] == ".")) or parameterEndReached:
            if not xReached:
                xReached = True
                newPropName = newPropName + propName[i]
            else:
                identifier = identifier + propName[i].lower()
                parameterEndReached = True
        else:
            newPropName = newPropName + propName[i]

    newFileName = "apc" + identifier + "_" + newPropName
    print(newFileName)
    
    foundFolder = False
    
    if not path.exists(path.join(propDestFolder + "/" + newFileName)):
        os.makedirs(path.join(propDestFolder + "/" + newFileName))
    
    src = propDataFolder + "/" + filename
    dest = propDestFolder + "/" + newFileName
        
    shutil.copy2(src, dest)