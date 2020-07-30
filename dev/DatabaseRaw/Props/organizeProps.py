#Used to organize Selig's prop database into subdirectories
import os
from os import path
import shutil

#Volume 1

propDataFolder = "C:/Users/Cory/Google Drive/AeroLab/UIUC-propDB/volume-1/data"
propDestFolder = "C:/Users/Cory/Google Drive/AeroLab/PropulsionUnitOptimization/Props"

for filename in os.listdir(path.join(propDataFolder)):
    
    filenameParts = filename.split("_")
    propName = filenameParts[0] + "_" + filenameParts[1]
    print(propName)
    foundFolder = False
    
    if not path.exists(path.join(propDestFolder + "/" + propName)):
        os.makedirs(path.join(propDestFolder + "/" + propName))
    
    src = propDataFolder + "/" + filename
    dest = propDestFolder + "/" + propName
        
    shutil.copy2(src, dest)

#Volume 2
    
propDataFolder = "C:/Users/Cory/Google Drive/AeroLab/UIUC-propDB/volume-2/data"
propDestFolder = "C:/Users/Cory/Google Drive/AeroLab/PropulsionUnitOptimization/Props"

for filename in os.listdir(path.join(propDataFolder)):
    
    filenameParts = filename.split("_")
    propName = filenameParts[0] + "_" + filenameParts[1]
    print(propName)
    foundFolder = False
    
    if not path.exists(path.join(propDestFolder + "/" + propName)):
        os.makedirs(path.join(propDestFolder + "/" + propName))
    
    src = propDataFolder + "/" + filename
    dest = propDestFolder + "/" + propName
        
    shutil.copy2(src, dest)