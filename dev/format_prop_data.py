"""Used to format propeller data into easily managed data files."""

import re
import os
import json

import sqlite3 as sql
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D


if __name__=="__main__":

    # Get location of props
    dev_dir = os.path.dirname(__file__)
    prop_source_dir = os.path.join(dev_dir, "Database/Props")
    prop_store_dir = os.path.join(dev_dir, "../pyprop/props")

    # Get files
    prop_data_dirs = os.listdir(prop_source_dir)

    # Initialize count
    prop_count = 0

    for prop_dir in prop_data_dirs:

        # Check for not a data folder
        if ".py" in str(prop_dir):
            continue
        if ".xls" in str(prop_dir):
            continue
        
        # Get path to current folder
        prop_dir_path = os.path.join(prop_source_dir, prop_dir)
        print(prop_dir)

        # Determine whether we have a Selig or APC prop (Selig is preferred)
        data_type = None
        for filename in os.listdir(prop_dir_path):
            if "PER" in filename:
                if data_type is None:
                    data_type = "apc"
                    data_file_name = filename
            if ("geom" in filename) or ("static" in filename):
                data_type = "selig"
                break

        # Increment count
        prop_count += 1

        # APC prop
        if data_type == "apc":

            # Open data file
            with open(os.path.join(prop_dir_path, data_file_name), 'r') as data_file:

                # Get pitch and diameter
                first_line = data_file.readline().split()
                dia_pitch = first_line[0].split("x")
                diameter = float(dia_pitch[0])
                pitch = float(re.search(r'\d+', dia_pitch[1]).group())
                manufacturer = "APC"

                # Read through the file to get how many sets of measurements were taken
                data_point_count = 0
                max_rpm = 30000
                if "apcc_7.4x7.5" in prop_dir:
                    max_rpm = 15000

                # Loop through lines
                for line in data_file:

                    line = line.replace("-", " ")
                    entries = line.split()
                    if len(entries) != 0:

                        # Line is the start of new rpm
                        if entries[0] == "PROP":

                            # Don't go past the max rpm
                            if float(entries[3]) > max_rpm:
                                break
                        
                        # Line has data
                        elif entries[0][0].isdigit() and "NaN" not in entries:
                            data_point_count += 1

                # Initialize storage array
                # We include an extra point because for any propeller, all coefficients are zero at zero J and rpm
                data = np.zeros((data_point_count+1, 4))

            # Go back to the beginning
            with open(os.path.join(prop_dir_path, data_file_name), 'r') as data_file:
                first_line = data_file.readline() # Get rid of the first line
                i = 1
                curr_rpm = None

                # Loop through lines
                for line in data_file:

                    # Split out line
                    line = line.replace("-", " ")
                    entries = line.split()

                    # Non-empty lines
                    if len(entries) != 0:

                        # New RPM value
                        if entries[0] == "PROP":

                            # Don't go past the max rpm
                            rpm = float(entries[3])
                            if rpm > max_rpm:
                                break

                            # Store current rpm
                            curr_rpm = rpm

                        # New advance ratio
                        elif entries[0][0].isdigit() and "NaN" not in entries:
                            data[i,0] = rpm
                            data[i,1] = float(entries[1])
                            data[i,2] = float(entries[3])
                            data[i,3] = float(entries[4])
                            i += 1

        #Data files from U of I U-C
        elif data_type == "selig":

            # Determine diameter
            dia_pitch = prop_dir.split("_")[1].split("x")
            diameter = float(dia_pitch[0])

            # Convert diameter from metric if needed
            met = False
            if diameter>50:
                diameter = diameter/25.4
                met = True

            # Determine pitch
            if "deg" in dia_pitch[1]:
                pitch = 2*np.pi*0.525*diameter*np.tan(np.radians(float(re.search(r'\d+', dia_pitch[1]).group())))
            else:
                pitch = float(dia_pitch[1])
                if met:
                    pitch = pitch/25.4

            # Find out manufacturer
            manufacturer = "UNKNOWN"
            manufacturers = [["an","AERONAUT"],
                             ["apc","APC"],
                             ["da","UIUC"],
                             ["ef","E-FLITE"],
                             ["gr","GRAUPNER"],
                             ["gws","GWS"],
                             ["kav","KAVON"],
                             ["kp","KP"],
                             ["kyosho","KYOSHO"],
                             ["ma","MASTERAIRSCREW"],
                             ["mi","MICROINVENT"],
                             ["nr","UIUC"],
                             ["pl","PLANTRACO"],
                             ["ru","REVUP"],
                             ["union","UNION"],
                             ["vp","VAPOR"],
                             ["zin","ZINGALI"]]
            for manu in manufacturers:
                if manu[0] in prop_dir:
                    manufacturer = manu[1]
                    break                 

            # Initialize count
            data_point_count = 0

            #Loop through files to count data points
            for data_file_name in os.listdir(prop_dir_path):

                # Get data files
                if not ("geom" in data_file_name or "PER" in data_file_name or "spec2" in data_file_name):

                    # Open file
                    with open(os.path.join(prop_dir_path, data_file_name)) as data_file:

                        # Throw away first line
                        first_line = data_file.readline()

                        # Increment count
                        for line in data_file:
                            entries = line.split()
                            if len(entries) != 0:
                                data_point_count += 1

            # Initialize storage array
            data = np.zeros((data_point_count+1, 4))

            # Loop thorugh files to read in measurements
            i = 1
            for data_file_name in os.listdir(prop_dir_path):

                # Get dynamic data
                if not ("static" in data_file_name or "geom" in data_file_name or "PER" in data_file_name or "spec2" in data_file_name):

                    # Get the rpm from the filename
                    data_file_name_parts = data_file_name.split("_")
                    curr_rpm = float(data_file_name_parts[-1].replace(".txt", ""))
                    with open(os.path.join(prop_dir_path, data_file_name), 'r') as data_file:

                        # Throw away first line
                        first_line = data_file.readline()

                        # Read in data points
                        for line in data_file:
                            entries = line.split()
                            if len(entries) != 0:
                                data[i,0] = curr_rpm
                                data[i,1] = float(entries[0])
                                data[i,2] = float(entries[1])
                                data[i,3] = float(entries[2])
                                i += 1
    
                # Get static data
                elif "static" in data_file_name and "spec2" not in data_file_name:
                    with open(os.path.join(prop_dir_path, data_file_name), 'r') as data_file:

                        # Throw away first line
                        first_line = data_file.readline()

                        # Read in data points
                        for line in data_file:
                            entries = line.split()
                            if len(entries) != 0:
                                data[i,0] = float(entries[0])
                                data[i,1] = 0.0
                                data[i,2] = float(entries[1])
                                data[i,3] = float(entries[2])
                                i += 1

        # Export data
        data_filename = os.path.join(prop_store_dir, prop_dir+".ppdat")
        header = "{:<25}{:<25}{:<25}{:<25}".format("RPM", "Adv Ratio", "Thrust Coef", "Power Coef")
        with open(data_filename, 'w') as data_file:
            np.savetxt(data_file, data, '%25.10E', header=header)

        # Parse info
        info_dict = {
            "name" : prop_dir,
            "manufacturer" : manufacturer,
            "diameter" : diameter,
            "pitch" : pitch,
            "data_file" : data_filename
        }
        info_filename = os.path.join(prop_store_dir, prop_dir+".ppinf")
        with open(info_filename, 'w') as info_file:
            json.dump(info_dict, info_file, indent=4)

    print("Sorted through {0} propellers.".format(prop_count))