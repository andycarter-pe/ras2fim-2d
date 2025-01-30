# Script 03 - Create Copies of base HEC-RAS model for each mainstem in the hydrofabric
# Also, create a precipitation file for each of the inflow points
#
# Created by: Andy Carter, PE
# Revised - 2025.01.07
# ************************************************************


# ************************************************************
import os
import re
import shutil

import h5py
import numpy as np
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon, LineString
import math
import configparser
from datetime import datetime, timedelta

import argparse
import time
import warnings
from tqdm import tqdm

import re
import xml.etree.ElementTree as ET

from pyproj import Transformer
from pyproj import CRS
# ************************************************************


# ---------------
def fn_format_flow_as_string(flt_flow):
    # format a flow value to no more than 5 characters
    # if bigger than 100000 ... returns 1.0e5
    
    # works up to 990 million cfs
    if flt_flow >= 100000:
        formatted_flow = "{:.1e}".format(flt_flow)
        parts = formatted_flow.split('e')
        formatted_flow = "{}e{}".format(parts[0], int(parts[1]))  # Reconstructing the notation
        if len(formatted_flow) > 8:  # Ensuring the total length is 8 characters
            formatted_flow = "{:.1e}".format(flt_flow).replace("e", "e+")
    else:
        formatted_flow = str(flt_flow)
    
    # returns a string
    return(formatted_flow)
# ---------------


# -------------------
# Function to find nearest boundary_int
def fn_find_nearest_boundary_int(value, boundary_values):
    nearest_boundary_int = None
    for boundary_int in boundary_values:
        if boundary_int <= value:
            nearest_boundary_int = boundary_int
        else:
            break
    return nearest_boundary_int
# -------------------


# --------------
def fn_build_plan_names(flt_delta_q, dict_mainstem):
    str_run_name, str_short_plan = None, None
    
    # Example:  str_run_name = '1884413_wb-2410249_wb-2410261_29-hr_14100-cfs_to_45000-cfs_step_2350-cfs'
    
    # systematic naming of plans (if only one flow, flt_upper_flow should be None)
    is_single_flow = False
    
    flt_lower_flow = dict_mainstem['q_lower_limit']
    flt_upper_flow = dict_mainstem['q_upper_limit']
    
    try:
        if flt_upper_flow == None:
            is_single_flow = True
            
        str_mainstem = str(int(dict_mainstem['mainstem']))
        str_start_node = dict_mainstem['id']
        str_end_node = dict_mainstem['id_end_node']
        int_firehose_time = int(dict_mainstem['time_sim_hr'])

        str_run_name = str_mainstem + "_" + str_start_node + "_"  + str_end_node + "_"
        str_run_name += str(int_firehose_time) + "-" + "hr"
        str_run_name += "_" + str(int(flt_lower_flow)) + "-cfs"

        if not is_single_flow:  
            # add the upper range flow to the description
            str_run_name += "_to_" + str(int(flt_upper_flow)) + "-cfs"
            
            # add the delta_q to the name of the run
            str_run_name += "_step_" + str(int(flt_delta_q)) + "-cfs"

        str_short_plan = str_start_node + '-' + fn_format_flow_as_string(flt_lower_flow)
    except:
        pass
    
    return(str_run_name, str_short_plan)
# --------------


# ----------------
def fn_list_of_ras_projects(str_ras_path):

    list_proj_title_files = []

    # locate the valid HEC-RAS project
    try:
        list_prj_files = [os.path.join(str_ras_path, f) for f in os.listdir(str_ras_path) if f.endswith('.prj')]

        if len(list_prj_files) > 0:
            list_proj_title_files = []

            for file_path in list_prj_files:
                with open(file_path, 'r') as file:
                    contents = file.read()
                    if 'Proj Title' in contents:
                        list_proj_title_files.append(file_path)

            for file_path in list_proj_title_files:
                #print(f"Valid HEC-RAS Project: {file_path}")
                pass

        else:
            print(f"No HEC-RAS PRJ found in {str_ras_path}")

    except Exception as e:
        print(f"An error occurred: {e}")
        
    return(list_proj_title_files)
# ----------------


# --------------------
def fn_list_filename_from_ras_prj(str_ras_prj_path, str_line_header):
    # Open the file
    with open(str_ras_prj_path, 'r') as file:
        # Read lines
        lines = file.readlines()

        # Initialize a list to store File values
        list_names = []

        # Iterate through each line
        for line in lines:
            # Check if the line contains str_line_header
            if str_line_header in line:
                # Extract the value after str_line_header
                value = line.split(str_line_header)[1].strip()
                # Add the value to the list
                list_names.append(value)
        return(list_names)
# --------------------


# ----------
def fn_extract_numbers_from_strings(lst):
    numbers = []
    for item in lst:
        number = re.search(r'\d+', item).group()
        numbers.append(int(number))
    return numbers
# ----------

# --------------
def fn_list_of_file_exists(list_filepaths):
    
    list_b_return = []
    for filepath in list_filepaths:
        if os.path.exists(filepath):
            list_b_return.append(True)
        else:
            list_b_return.append(False)
    return(list_b_return)
# --------------


# ---------------
def fn_generate_simulation_date_string(start_date_str, hours):
    
    # The first hour is the "zero-th" hour
    hours = hours - 1
    
    # Convert start_date_str to a datetime object
    start_date = datetime.strptime(start_date_str, '%d%b%Y')

    # Calculate end date by adding hours to the start date
    end_date = start_date + timedelta(hours=hours)

    # Format dates to the required string format
    start_date_formatted = start_date.strftime('%d%b%Y,%H%M')
    end_date_formatted = end_date.strftime('%d%b%Y,%H%M')

    # Concatenate start and end date strings
    date_string = f"{start_date_formatted},{end_date_formatted}"

    return date_string
# ---------------


# ***************
def fn_set_plan_simulation_date_range(str_copy_to_p_full_path,
                                      str_simulation_start_date,
                                      int_sim_time_hr):
    
    # In the plan file, set the simulation date range

    # Read the file content
    with open(str_copy_to_p_full_path, "r") as file:
        lines = file.readlines()

    # Find indices and lines that start with "Geom File="
    list_tup_match = [(index, line.strip()) for index, line in enumerate(lines) if re.match(r'^Simulation Date=', line)]

    str_simulation_range = fn_generate_simulation_date_string(str_simulation_start_date, int_sim_time_hr)

    # Replace lines in the list
    for index, _ in list_tup_match:
        lines[index] = "Simulation Date=" + str_simulation_range + "\n"  # Replace Simulation Date

    # Write the modified lines back to the file
    with open(str_copy_to_p_full_path, "w") as file:
        file.writelines(lines)
# ***************


# -----------------
def fn_generate_step_integer(flt_lower_flow, flt_upper_flow, flt_delta_q):
    integers = []
    current = flt_lower_flow
    if flt_upper_flow != None and flt_delta_q != None:
        if flt_upper_flow > flt_lower_flow:
            while current < flt_upper_flow:
                integers.append(int(current))
                current += flt_delta_q

            # Add flt_upper_flow back if it's not in the list
            if integers[-1] != flt_upper_flow:
                integers.append(int(flt_upper_flow))
        else:
            integers.append(int(flt_lower_flow))
    else:
        integers.append(int(flt_lower_flow))
    return (integers)
# -----------------


# ---------------------
def fn_copy_ras_file(tup_of_file_params,
                     str_project_path,
                     str_run_name,
                     str_new_folder):
    
    
    # tup_of_file_params: Example ... ('Plan File=', 'p', 1)
    # str_project_path --- path to project to copy from Ex: 'E:\\HECRAS_2D_12070205\\base_model_20240414_copy\\BLE_LBSG_501.prj'
    # str_run_name:  new models name:  Ex: 1884895_wb-2410559_wb-2410559_11-hr_500-cfs
    # str_new_folder:  folder to place the new spawned copies... will contain multiple models
    
    # ---------
    str_file_type_header = tup_of_file_params[0]
    str_header = tup_of_file_params[1]
    int_to_copy = tup_of_file_params[2]
    # ---------
    
    str_copy_to_full_path = None

    # Splitting the path into folder and filename
    str_prj_folder, str_filename = os.path.split(str_project_path)

    # Splitting the filename and its extension
    str_file_only, str_extension = os.path.splitext(str_filename)

    # determine the files of type 'str_file_type_header' that are in the project
    list_names = fn_list_filename_from_ras_prj(str_project_path, str_file_type_header)
    list_int = fn_extract_numbers_from_strings(list_names)

    list_fullpath = []

    for item in list_names:
        str_file = str_file_only + "." + item

        # Combine the folder path and filename
        str_full_path = os.path.join(str_prj_folder, str_file)

        list_fullpath.append(str_full_path)

    list_b_exists = fn_list_of_file_exists(list_fullpath)

    # Assuming all lists have the same length
    data = {
        'int': list_int,
        'name': list_names,
        'path': list_fullpath,
        'exists': list_b_exists,
    }

    df = pd.DataFrame(data)

    # Filter the DataFrame
    df_filtered = df[(df['int'] == int_to_copy) & (df['exists'])]

    # Check if any rows match the condition
    if not df_filtered.empty:
        # If there's at least one row matching the condition
        #print('---------------------------')
        #print('Valid "' + str_file_type_header + '" Found')

        # copy file
        str_copy_from = df_filtered.iloc[0]['path']
        str_copy_to = str_run_name + '.' + str_header + "01"
        str_copy_to_full_path = os.path.join(str_new_folder , str_copy_to)

        shutil.copy(str_copy_from, str_copy_to_full_path)
        
        # --- try to copy the hdf versions of the same file ---
        str_copy_from_hdf = str_copy_from + ".hdf"
        str_copy_to_full_path_hdf = str_copy_to_full_path + ".hdf"
        
        # Copy only the geom and unsteady hdf (gXX.hdf and uXX.hdf)
        if str_header == 'g' or str_header == 'u':
            shutil.copy(str_copy_from_hdf, str_copy_to_full_path_hdf)
        # ---
        
        # ---
        if str_header == 'g':
            str_c_geom = str_copy_from[:-3] + 'c' + str_copy_from[-2:]
            str_x_geom = str_copy_from[:-3] + 'x' + str_copy_from[-2:]
            
            str_copy_c_full_path = str_copy_to_full_path[:-3] + 'c' + str_copy_to_full_path[-2:]
            str_copy_x_full_path = str_copy_to_full_path[:-3] + 'x' + str_copy_to_full_path[-2:]
            
            shutil.copy(str_c_geom, str_copy_c_full_path)
            shutil.copy(str_x_geom, str_copy_x_full_path)
        
    else:
        #print('---------------------------')
        print('Invalid "' + str_file_type_header + '" NOT FOUND')
# ---------------------


# ----------------------
def fn_create_prj(str_run_name, str_project_path, str_new_folder):

    # ~~~~~~~~~~
    # Get the current date
    current_date = datetime.now()

    # Format the date as "23 July 2022"
    str_formatted_date = current_date.strftime("%d %B %Y")

    str_description = 'RAS2FIM-2D Run: ' + str_run_name + '\n'
    str_description += 'Created from source model: ' + str_project_path + '\n'
    str_description += 'Created for RAS2FIM-2D analysis \n'
    str_description += 'Andy Carter - University of Texas at Austin \n'
    str_description += 'Center for Water and the Environment \n'
    str_description += 'Created: ' + str_formatted_date + '\n'
    # ~~~~~~~~~~

    # Create a new PRJ file for the new model.

    str_file_body = 'Proj Title=' + str_run_name + '\n'
    str_file_body += 'Current Plan=p01\n' # assumed to be the only plan
    str_file_body += 'Default Exp/Contr=0.3,0.1\n'
    str_file_body += 'English Units\n'

    # assuming the g,u and p are 01
    str_file_body += 'Geom File=g01\n'
    str_file_body += 'Unsteady File=u01\n'
    str_file_body += 'Plan File=p01\n'

    str_file_body += 'Y Axis Title=Elevation\n'
    str_file_body += 'X Axis Title(PF)=Main Channel Distance\n'
    str_file_body += 'X Axis Title(XS)=Station\n'
    str_file_body += 'BEGIN DESCRIPTION:\n'

    str_file_body += str_description

    str_file_body += 'END DESCRIPTION:\n'
    str_file_body += 'DSS Start Date=\n'
    str_file_body += 'DSS Start Time=\n'
    str_file_body += 'DSS End Date=\n'
    str_file_body += 'DSS End Time=\n'
    str_file_body += 'DSS File=dss\n'
    str_file_body += 'DSS File= \n'
    str_file_body += 'DSS Export Filename= \n'
    str_file_body += 'DSS Export Rating Curves= 0 \n'
    str_file_body += 'DSS Export Rating Curve Sorted= 0 \n'
    str_file_body += 'DSS Export Volume Flow Curves= 0 \n'
    str_file_body += 'DXF Filename= \n'
    str_file_body += 'DXF OffsetX= 0 \n'
    str_file_body += 'DXF OffsetY= 0 \n'
    str_file_body += 'DXF ScaleX= 1 \n'
    str_file_body += 'DXF ScaleY= 10 \n'
    str_file_body += 'GIS Export Profiles= 0 \n'


    str_prj_file_name = str_run_name + ".prj" 

    str_prj_fullpath = os.path.join(str_new_folder, str_prj_file_name)
    with open(str_prj_fullpath, 'w') as file:
        file.write(str_file_body)
# ----------------------


# ~~~~~~~~~~~~~~~~~
def fn_adjust_plan(str_run_name,
                   str_new_folder,
                   dict_flows,
                   int_sim_time_hr,
                   str_simulation_start_date):
    
    # *************************************
    # Adjust the plan file
    # *************************************

    str_plan_file_name = str_run_name + ".p01" 
    str_plan_fullpath = os.path.join(str_new_folder, str_plan_file_name)

    # ---------------------
    # Within the plan file, adjust the selected geometry file

    # Read the file content
    with open(str_plan_fullpath, "r") as file:
        lines = file.readlines()

    # Find indices and lines that start with "Geom File="
    list_tup_match = [(index, line.strip()) for index, line in enumerate(lines) if re.match(r'^Geom File=', line)]

    # Replace lines in the list
    for index, _ in list_tup_match:
        lines[index] = "Geom File=g01\n"

    # Write the modified lines back to the file
    with open(str_plan_fullpath, "w") as file:
        file.writelines(lines)

    # ---------------------
    # Within the plan file, adjust the selected unsteady file

    # Read the file content
    with open(str_plan_fullpath, "r") as file:
        lines = file.readlines()

    # Find indices and lines that start with "Flow File="
    list_tup_match = [(index, line.strip()) for index, line in enumerate(lines) if re.match(r'^Flow File=', line)]

    # Replace lines in the list
    for index, _ in list_tup_match:
        lines[index] = "Flow File=u01\n"

    # Write the modified lines back to the file
    with open(str_plan_fullpath, "w") as file:
        file.writelines(lines)

    # ---------------------
    # Within the plan file, change the plan 

    # Read the file content
    with open(str_plan_fullpath, "r") as file:
        lines = file.readlines()

    # --- Append the plan title in the plan file
    # Find indices and lines that start with "Flow File="
    list_tup_match = [(index, line.strip()) for index, line in enumerate(lines) if re.match(r'^Plan Title=', line)]

    # Replace lines in the list
    for index, _ in list_tup_match:
        lines[index] = "Plan Title=" + dict_flows['run_name'] + "_plan" + "\n"

    # --- Append the short identifier in the plan file
    # Find indices and lines that start with "Flow File="
    list_tup_match = [(index, line.strip()) for index, line in enumerate(lines) if re.match(r'^Short Identifier=', line)]  

    # Replace lines in the list
    for index, _ in list_tup_match:
        lines[index] = "Short Identifier=" + dict_flows['short_name'] + "\n"

    # Write the modified lines back to the file
    with open(str_plan_fullpath, "w") as file:
        file.writelines(lines)

    fn_set_plan_simulation_date_range(str_plan_fullpath,
                                      str_simulation_start_date,
                                      int_sim_time_hr)
# ~~~~~~~~~~~~~~~~~


# .................
def fn_adjust_geom(str_run_name, str_new_folder, dict_flows):

    # *************************************
    # Adjust the geometry file
    # *************************************

    str_geom_file_name = str_run_name + ".g01" 
    str_geom_fullpath = os.path.join(str_new_folder, str_geom_file_name)

    # --------------------- 
    # Change the name of the Geom Title
    # Read the file content
    with open(str_geom_fullpath, "r") as file:
        lines = file.readlines()

    str_run_name = dict_flows['run_name']

    # Find indices and lines that start with "Geom Title="
    list_tup_match = [(index, line.strip()) for index, line in enumerate(lines) if re.match(r'^Geom Title=', line)]

    # Replace lines in the list
    for index, _ in list_tup_match:
        lines[index] = "Geom Title=" + str_run_name + '_geom' + "\n"

    # Write the modified lines back to the file
    with open(str_geom_fullpath, "w") as file:
        file.writelines(lines)
# .................


# ----------------------
def fn_adjust_unsteady_flow(str_run_name, str_new_folder, dict_flows):

    # *************************************
    # Adjust the unsteady flow file
    # *************************************

    str_flow_file_name = str_run_name + ".u01" 
    str_flow_fullpath = os.path.join(str_new_folder, str_flow_file_name)

    # Enforce a 1HOUR Interval ("Inverval=1HOUR") in the unsteady flow file

    # Read the file content
    with open(str_flow_fullpath, "r") as file:
        lines = file.readlines()

    # Find indices and lines that start with "Boundary Location="
    list_tup_match = [(index, line.strip()) for index, line in enumerate(lines) if re.match(r'^Interval=', line)]

    # Replace lines in the list
    for index, _ in list_tup_match:
        lines[index] = "Interval=1HOUR\n"  # Replace the line with "Interval=1HOUR"

    # Write the modified lines back to the file
    with open(str_flow_fullpath, "w") as file:
        file.writelines(lines)

    # --------------------- 
    # Change the name of the Unsteady Flow Title
    # Read the file content
    with open(str_flow_fullpath, "r") as file:
        lines = file.readlines()

    str_run_name = dict_flows['run_name']

    # Find indices and lines that start with "Flow Title="
    list_tup_match = [(index, line.strip()) for index, line in enumerate(lines) if re.match(r'^Flow Title=', line)]

    # Replace lines in the list
    for index, _ in list_tup_match:
        lines[index] = "Flow Title=" + str_run_name + '_unsteady_flow' + "\n"

    # Write the modified lines back to the file
    with open(str_flow_fullpath, "w") as file:
        file.writelines(lines)
# ----------------------


# ===========================
def fn_copy_source_terrain(str_folder_for_copies, str_new_terain_folder_name, str_source_geom_hdf):
    
    print('+++++++++++++++')
    print('str_folder_for_copies:'  + str_folder_for_copies)
    print('str_new_terain_folder_name:'  + str_new_terain_folder_name)
    print('str_source_geom_hdf:'  + str_source_geom_hdf)
    print('+++++++++++++++')
    
    # Define the path for the new subfolder
    source_terrain_path = os.path.join(str_folder_for_copies, str_new_terain_folder_name)

    # Remove the existing source_terrain folder if it exists
    if os.path.exists(source_terrain_path):
        try:
            shutil.rmtree(source_terrain_path)
            #print(f"Removed existing directory: {source_terrain_path}")
        except Exception as e:
            print(f"Error removing existing directory: {e}")
            # Exit the function if the directory removal fails
            exit(1)

    # Create the folder where the copies will be placed
    #os.makedirs(source_terrain_path)

    # Check if the source geometry HDF exists
    if os.path.exists(str_source_geom_hdf):
        # Read the HDF file
        with h5py.File(str_source_geom_hdf, 'r') as hdf_file:
            # Navigate to the "/Geometry" group
            geometry_group = hdf_file["/Geometry"]

            # Retrieve and print the attribute value of "Terrain Filename"
            str_terrain_filename = geometry_group.attrs.get("Terrain Filename", b"Attribute not found").decode('utf-8')
            
            print('str_terrain_filename: ' + str_terrain_filename)

            # Get the directory of the HDF file and the absolute terrain path
            hdf_dir = os.path.dirname(str_source_geom_hdf)
            absolute_terrain_path = os.path.abspath(os.path.join(hdf_dir, str_terrain_filename))

            # Get the absolute folder and the filename of the terrain
            absolute_terrain_dir = os.path.dirname(absolute_terrain_path)
            absolute_terrain_basename = os.path.basename(absolute_terrain_path)

            # Check if the terrain directory exists
            if os.path.exists(absolute_terrain_dir):
                # Copy everything from the terrain directory to the new source_terrain_path
                try:
                    print('absolute_terrain_dir: '+ absolute_terrain_dir)
                    print('source_terrain_path: '+ source_terrain_path)
                    print('+++++++++++++++')
                    
                    shutil.copytree(absolute_terrain_dir, source_terrain_path)
                    #print(f"Copied terrain files from {absolute_terrain_dir} to {source_terrain_path}")
                except Exception as e:
                    print(f"Error copying terrain files: {e}")
            else:
                print(f"Terrain directory not found: {absolute_terrain_dir}")
    else:
        print(f"Source geom HDF not found: {str_source_geom_hdf}")
        
    return(absolute_terrain_basename)
# ===========================


# ---------------------
def fn_adjust_geom_hdf(str_run_name, str_new_folder, str_terrain_relativepath):

    str_geom_hdf_path = os.path.join(str_new_folder,str_run_name + ".g01.hdf")
    str_length = len(str_terrain_relativepath)

    # Convert the string to a numpy array with the required shape and datatype
    terrain_array = np.array(str_terrain_relativepath.ljust(str_length, '\0'), dtype=h5py.string_dtype(encoding='ascii', length=str_length))

    with h5py.File(str_geom_hdf_path, 'r+') as f:
        # Access the '/Geometry' group and set the attribute
        geometry_group = f['/Geometry']
        geometry_group.attrs['Terrain Filename'] = terrain_array
        
        # *** hard coding the 2D area
        # Specify the HDF5 file path and group path
        str_hdf_geom_path = '/Geometry/2D Flow Areas/'
        
        # Get names of HDF5 Group objects in the specified group
        list_group_names = fn_get_group_names(str_geom_hdf_path, str_hdf_geom_path)
        
        geometry_group = f['/Geometry/2D Flow Areas/' + list_group_names[0]]
        geometry_group.attrs['Terrain Filename'] = terrain_array
# ---------------------


# ------------------------
def fn_get_group_names(hdf5_file_path, group_path):
    """
    Retrieve the names of groups within a specified HDF5 file under a given group path.

    Parameters:
    hdf5_file_path (str): The file path to the HDF5 file.
    group_path (str): The path to the group whose subgroups' names are to be retrieved.

    Returns:
    list or None: A list containing the names of groups found under the specified group path. 
                  Returns None if the group path does not exist in the HDF5 file.
    """
    try:
        with h5py.File(hdf5_file_path, 'r') as hdf_file:
            # Check if the specified group path exists
            if group_path in hdf_file:
                group = hdf_file[group_path]

                # Extract names of HDF5 Group objects
                group_names = [name for name in group if isinstance(group[name], h5py.Group)]

                return group_names
            else:
                print(f"Group '{group_path}' not found in the HDF5 file.")
                return None
    except Exception as e:
        print(f"An error occurred: {e}")
# ------------------------


# ------------------------
def fn_get_projection_wkt(hdf_file_path):

    # Specify the HDF5 file path and group path
    str_hdf_geom_path = '/Geometry/2D Flow Areas/'

    # Get names of HDF5 Group objects in the specified group
    list_group_names = fn_get_group_names(hdf_file_path, str_hdf_geom_path)

    b_has_2d_area = False

    if len(list_group_names) > 1:
        print('Multiple 2D areas found -- Using first area:', list_group_names[0])
        b_has_2d_area = True
    elif len(list_group_names) == 0:
        print('Error: No 2D areas found')
    else:
        # Only one 2D area found
        b_has_2d_area = True

    if b_has_2d_area:
        # Open the HDF file
        with h5py.File(hdf_file_path, 'r') as hdf_file:
            # Extract the projection
            projection_wkt = hdf_file.attrs['Projection'].decode('utf-8')

        return(projection_wkt)
    else:
        pass
        # return nothing as there is nothing to return
# ------------------------


# ----------
def fn_find_rasmap_projection_filenames(str_ras_path):

    # Determine if source folder (str_ras_path) contains any *.rasmap files
    # If is does, read the rasmap file (xml) and look for "RASProjectionFilename" element
    # and the "Filename" attribute

    list_proj_filenames = []
    list_rasmap_filenames = []

    # Loop through files in the folder
    for file_name in os.listdir(str_ras_path):
        # Check for files with .rasmap extension
        if file_name.endswith('.rasmap'):
            rasmap_path = os.path.join(str_ras_path, file_name)
            list_rasmap_filenames.append(rasmap_path)

            # Parse the XML file
            try:
                tree = ET.parse(rasmap_path)
                root = tree.getroot()

                # Find RASProjectionFilename element
                ras_proj = root.find('RASProjectionFilename')
                if ras_proj is not None:
                    # Get the Filename attribute if it exists
                    filename_attr = ras_proj.get('Filename')
                    if filename_attr:
                        # Check if filename_attr is a relative path
                        if not os.path.isabs(filename_attr):
                            # Convert to an absolute path relative to rasmap_path
                            filename_attr = os.path.abspath(os.path.join(os.path.dirname(rasmap_path), filename_attr))
                        list_proj_filenames.append(filename_attr)
            except ET.ParseError:
                print(f"Error parsing {rasmap_path}")

    return list_proj_filenames, list_rasmap_filenames
# ----------


# ---------
def fn_replace_rasmap_elements(str_rasmap_copyto, str_run_name, str_header_val):
    
    # replace all attributes where "Filename" attributes end with 
    #".uXX.hdf" (XX is wildcard) with ".\" + str_run_name + ".u01.hdf"
    # In this case, str_header_val is 'u'
    
    # Load the XML file
    tree = ET.parse(str_rasmap_copyto)
    root = tree.getroot()

    # Regex pattern to match filenames ending with ".uXX.hdf" where XX is a wildcard of two digits
    pattern = re.compile(r'\.' + str_header_val + '\d{2}\.hdf$')

    # Modify attributes where "Filename" matches the pattern
    for elem in root.iter():
        filename = elem.attrib.get("Filename")
        if filename and pattern.search(filename):
            # Update the filename as specified
            elem.set("Filename", f".\\{str_run_name}.{str_header_val}01.hdf")

    # Save the modified XML back to file or overwrite the original
    tree.write(str_rasmap_copyto)
# ---------

# ~~~~~~~~~~~~~~~~~
def fn_update_rasmap_xml(str_rasmap_copyto, str_proj_relativepath, str_terrain_relativepath):
    # Load the XML file
    tree = ET.parse(str_rasmap_copyto)
    root = tree.getroot()

    # Update the RASProjectionFilename attribute
    projection_element = root.find(".//RASProjectionFilename")
    if projection_element is not None:
        projection_element.set("Filename", str_proj_relativepath)

    # Update all Layer Filename attributes within the Terrains element
    for layer in root.findall(".//Terrains/Layer"):
        layer.set("Filename", str_terrain_relativepath)

    # Save the modified XML file
    tree.write(str_rasmap_copyto)
# ~~~~~~~~~~~~~~~~~


# -------
def fn_convert_point_to_output_crs(x,y,str_input_crs, str_output_crs):

    # Note: lon, lat is in a x,y format
    #transformer = Transformer.from_crs("EPSG:4326", "EPSG:2277", always_xy=True)
    transformer = Transformer.from_crs(str_input_crs, str_output_crs, always_xy=True)

    # Perform the transformation
    x_transformed, y_transformed = transformer.transform(x, y)

    # return a tuple
    return((x_transformed, y_transformed))
# -------


# .............
# Function to calculate the Euclidean distance between two points
def euclidean_distance(point1, point2):
    x1, y1 = point1
    x2, y2 = point2
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
# .............


# ~~~~~~~~~~~~~
# Function to determine the largest pixel dimension (one point per pixel)
def fn_largest_square_pixel_dimension(points):
    min_distance = float('inf')  # Start with a very large value
    
    # Compare each point with every other point
    for i in range(len(points)):
        for j in range(i+1, len(points)):
            dist = euclidean_distance(points[i], points[j])
            if dist < min_distance:
                min_distance = dist
    
    # The largest pixel size will now be exactly the minimum distance
    return min_distance
# ~~~~~~~~~~~~~


# ------------
# Function to determine the cell in which each point is located
def fn_get_cell_for_point(point, pixel_size, upper_left_x, upper_left_y):
    x, y = point
    
    # Calculate the cell column and row based on the pixel size
    # Adding a small offset to ensure separation
    col = int((x - upper_left_x) // pixel_size)
    row = int((upper_left_y - y) // pixel_size)

    # Ensure we account for potential negative indices (i.e., if the point is outside the grid)
    col = max(0, col)
    row = max(0, row)

    return (row, col)
# ------------


# ------------
# Use this function to get cell coordinates for each point
def fn_assign_points_to_cells(points, pixel_size):
    # Calculate the bounding box, similar to the plotting function
    min_x = min(point[0] for point in points)
    max_x = max(point[0] for point in points)
    min_y = min(point[1] for point in points)
    max_y = max(point[1] for point in points)

    # Adjust for half the pixel size to center the grid
    min_x = min_x - pixel_size / 2
    max_x = max_x + pixel_size / 2
    min_y = min_y  - pixel_size / 2
    max_y = max_y + pixel_size / 2
    
    # The upper-left corner of the grid
    upper_left_x = min_x
    upper_left_y = max_y
    
    # Assign each point to a cell
    cell_assignments = [fn_get_cell_for_point(point, pixel_size, upper_left_x, upper_left_y) for point in points]
    
    return cell_assignments
# ------------


# ~~~~~~~~~~~~~~~
# Function to plot cells along with points
def fn_raster_parameters(points, pixel_size):
    min_x = min(point[0] for point in points)
    max_x = max(point[0] for point in points)
    min_y = min(point[1] for point in points)
    max_y = max(point[1] for point in points)

    # Calculate grid boundaries
    min_x = min_x - pixel_size / 2
    max_x = max_x + pixel_size / 2
    min_y = min_y - pixel_size / 2
    max_y = max_y + pixel_size / 2
    
    upper_left_x = min_x
    upper_left_y = max_y
    
    # Calculate the number of rows and columns
    int_width = int(np.ceil((max_x - min_x) / pixel_size))
    int_height = int(np.ceil((max_y - min_y) / pixel_size))
    
    return(upper_left_x, upper_left_y, int_width, int_height)
# ~~~~~~~~~~~~~~~


# ------------
def fn_create_time_list(str_start_time, int_delta_hr, list_peak_flows):
    int_len = len(list_peak_flows) + 1

    # Convert the start time to a datetime object
    start_time = datetime.strptime(str_start_time, "%Y-%m-%d %H:%M:%S")

    # Create the list of time strings
    list_time_strings = [(start_time + timedelta(hours=i * int_delta_hr)).strftime("%Y-%m-%d %H:%M:%S") for i in range(int_len)]
    
    # Return a list of datetime strings
    return list_time_strings
# ------------


# ~~~~~~~~~~~
def fn_list_of_peak_flow(flt_low_q, flt_high_q, flt_q_step):
    list_peak_flows = []
    flt_q_val = flt_low_q

    # Increment by flt_q_step using a while loop
    while flt_q_val < flt_high_q:
        list_peak_flows.append(flt_q_val)
        flt_q_val += flt_q_step

    # Include the last value if it isn't covered by the loop
    if list_peak_flows[-1] < flt_high_q:
        list_peak_flows.append(flt_high_q)

    return list_peak_flows
# ~~~~~~~~~~~


# ................
def fn_create_cumulative_rainfall_array(dict_params):
    
    list_time_strings = dict_params['list_time_strings']
    list_peak_flows = dict_params['list_peak_flows']
    list_scaled_flows = dict_params['list_scaled_flows']
    int_delta_hr = dict_params['int_delta_hr']
    list_cell_assignments = dict_params['list_cell_assignments']
    flt_largest_pixel_size = dict_params['flt_largest_pixel_size']
    
    list_flow_per_time = []
    
    for idx_flow, item_time in enumerate(list_time_strings):
        list_scaled_flowrate_per_timestep = [scaled_flow * list_peak_flows[idx_flow] for scaled_flow in list_scaled_flows]
        list_flow_per_time.append(list_scaled_flowrate_per_timestep)

    # convert to numpy array
    array_flow_per_time = np.array(list_flow_per_time)

    # convert the list of flows to cummulative precip by timestep duration
    flt_pixel_area = flt_largest_pixel_size * flt_largest_pixel_size

    # convert from cfs to inches of depth per time step
    # Depth (inches) = (Q/A) * 12 * 3600 * duration of rain
    convert_to_depth = 12 * 3600 * int_delta_hr / flt_pixel_area
    array_flow_per_time_depth = array_flow_per_time * convert_to_depth

    # transpose the array
    array_flow_per_point = array_flow_per_time_depth.T

    # create the rainfall array (as zero values)
    array_flows = np.zeros((len(list_peak_flows), dict_params['int_height'], dict_params['int_width']))
    #print(array_flows.shape)  # This will print the dimensions of the array
    #print(array_flow_per_point)

    # Apply cumulative sum along each row (axis=1)
    accumulated_array = np.cumsum(array_flow_per_point, axis=1)
    #print(accumulated_array)
    #print(list_cell_assignments)

    # put the flow data into the array
    for i in range(len(list_cell_assignments)):
        for j in range(len(list_peak_flows)):
            array_flows[j, list_cell_assignments[i][0], list_cell_assignments[i][1]] = accumulated_array[i][j]

    # reshape to what HEC-RAS unsteady HDF uses a 2D array
    array_flows_reshaped = array_flows.reshape(len(list_peak_flows), -1)  # -1 lets NumPy calculate the size automatically
    #print(array_flows_reshaped.shape)
    
    # ------Test Added : 2024.11.22 -----
    # extend the accumulated rainfall one more time step at the same gradient as between the last two rows
    if array_flows_reshaped.shape[0] >= 2:
        # Extract the last row and second-to-last row
        last_row = array_flows_reshaped[-1, :]
        second_to_last_row = array_flows_reshaped[-2, :]
        
        # Compute the new row
        new_row = 2 * last_row - second_to_last_row
        
        # Add the new row to the array
        array_flows_reshaped = np.vstack([array_flows_reshaped, new_row])
    # ------Test Added : 2024.11.22 -----
    
    return(array_flows_reshaped)
# ................


# =====================
def fn_overwrite_rainfall_hdf(dict_params, array_flows, hdf_path, str_dataset_path):
    # function to overwrite the precipitation hdf file
    
    # Open the HDF5 file in read/write mode
    with h5py.File(hdf_path, 'r+') as hdf_file:
        # Navigate to the dataset where the attribute exists

        # Check if the dataset exists
        if str_dataset_path in hdf_file:
            dataset = hdf_file[str_dataset_path]

            # Backup attributes before replacing the dataset
            attrs_backup = {key: dataset.attrs[key] for key in dataset.attrs}

            # delete the dataset and create a new one
            del hdf_file[str_dataset_path]

            # Create new dataset with specified settings (float32, GZIP, unlimited dims)
            new_dataset = hdf_file.create_dataset(
                str_dataset_path, 
                data=array_flows,
                dtype=np.float32,                # Ensure 32-bit float precision
                compression="gzip",              # Enable GZIP compression
                maxshape=(None, None)      # Set max dimensions to unlimited
            )

            # Restore the attributes
            for key, value in attrs_backup.items():
                new_dataset.attrs[key] = value

            # --- Set new attributes values ---
            # --- revise the projection ---
            # Define the length of the projection string
            str_length = len(dict_params['str_projection'])  # Adjust this length based on your attribute requirements

            # Convert the projection string to a numpy array with the required shape and datatype
            projection_array = np.array(dict_params['str_projection'].ljust(str_length, '\0'),
                                        dtype=h5py.string_dtype(encoding='ascii', length=str_length))

            # Create or overwrite the projection attribute with the new string
            new_dataset.attrs.create('Projection', projection_array)

            # --- adjust the times ---
            list_time_strings = dict_params['list_time_strings']
            list_time_strings.append(dict_params['str_last_time'])
            
            #arr_binary_time = np.array(dict_params['list_time_strings'])
            arr_binary_time = np.array(list_time_strings)

            # Create a dtype for the times - ASCII strings of length 19, null-terminated
            time_dtype = h5py.string_dtype(encoding='ascii', length=19)

            # Convert the times array to this dtype
            arr_binary_time_h5 = np.array(arr_binary_time, dtype=time_dtype)

            new_dataset.attrs['Times'] = arr_binary_time_h5

            str_units = 'in' # hard coded units
            units_array = np.array(str_units.ljust(len(str_units), '\0'),
                                   dtype=h5py.string_dtype(encoding='ascii',length=len(str_units)))
            new_dataset.attrs.create('Units', units_array)

            # ----
            # Change or create the cell size attribute as a 64-bit float
            new_dataset.attrs.create('Raster Cellsize', np.array(dict_params['flt_largest_pixel_size'], dtype='float64'))
            new_dataset.attrs.create('Raster Cols', np.array(dict_params['int_width'], dtype='int32'))
            new_dataset.attrs.create('Raster Left', np.array(dict_params['upper_left_x'], dtype='float64'))
            new_dataset.attrs.create('Raster Rows', np.array(dict_params['int_height'], dtype='int32'))
            new_dataset.attrs.create('Raster Top', np.array(dict_params['upper_left_y'], dtype='float64'))

        else:
            print(f"Dataset '{str_dataset_path}' not found in the HDF file.")
# =====================


# ~~~~~~~~~~~~~~~~~~~~~~~
def fn_revise_rainfall_hdf(dict_for_precip, projection_wkt):

    # Convert the Point to list of tuples of (x,y)
    point = dict_for_precip['geometry']
    list_rainfall_points = [(point.x, point.y)]

    # ** Hard coded - no flow scaling of input precip **
    list_scaled_flows = [1.0 for _ in range(len(list_rainfall_points))]

    # Timestep variables
    str_start_time = dict_for_precip['str_date_for_hdf']
    int_delta_hr = dict_for_precip['time_sim_hr']

    str_hdf_to_edit = dict_for_precip['str_hdf_to_edit']
    # -----------
    
    flt_cell_edge_length = dict_for_precip['cell_edge_length']

    # pixel area in square feet
    ftl_cell_area = flt_cell_edge_length * flt_cell_edge_length # value in square feet

    # list of the stepped flows... 0,250,500 cfs... etc.
    list_peak_flows = fn_list_of_peak_flow(dict_for_precip['lower_flow'],
                                           dict_for_precip['upper_flow'],
                                           dict_for_precip['delta_q'])

    # List of converted points (projected to output CRS)
    list_converted_points_tup = [fn_convert_point_to_output_crs(x, y, dict_for_precip['input_crs'], dict_for_precip['output_crs']) for x, y in list_rainfall_points]

    # Calculate the largest square pixel side dimension where each square contains only one point
    flt_largest_pixel_size = fn_largest_square_pixel_dimension(list_converted_points_tup)
    
    # Set the pixel size to the either user specified size or just
    # a little smaller than the "cell per point" raster... whichever is smaller
    if flt_largest_pixel_size > flt_cell_edge_length:
        flt_largest_pixel_size = flt_cell_edge_length
    else:
        flt_largest_pixel_size = flt_largest_pixel_size * 0.90
        
    # Get the cell assignments for each point
    list_cell_assignments = fn_assign_points_to_cells(list_converted_points_tup, flt_largest_pixel_size)
    
    # Call the function to get raster parameters
    upper_left_x,upper_left_y,int_width,int_height = fn_raster_parameters(list_converted_points_tup, flt_largest_pixel_size)

    # Check for non-unique cell assignments
    if len(list_cell_assignments) != len(set(list_cell_assignments)):
        b_all_cells_unique = False
        #print("Warning: There are non-unique cell assignments.")
    else:
        #print("All cell assignments are unique.")
        b_all_cells_unique = True

    if b_all_cells_unique:
        list_time_strings = fn_create_time_list(str_start_time, int_delta_hr, list_peak_flows)

        # now create an array to stuff into the unsteady hdf file
        #print('/Event Conditions/Meteorology/Precipitation/Imported Raster Data/Values')

        dict_rainfall_params = {
            "flt_largest_pixel_size": flt_largest_pixel_size,
            "int_width": int_width,
            "upper_left_x": upper_left_x,
            "int_height": int_height,
            "upper_left_y": upper_left_y,
            "int_delta_hr": int_delta_hr,
            "list_rainfall_points": list_rainfall_points,
            "list_cell_assignments": list_cell_assignments,
            "list_time_strings": list_time_strings[:-1],
            "str_last_time": list_time_strings[-1],
            "list_peak_flows": list_peak_flows,
            "list_scaled_flows": list_scaled_flows,
            "str_projection": projection_wkt}

        #print(dict_rainfall_params)
        array_flows_reshaped = fn_create_cumulative_rainfall_array(dict_rainfall_params)
        
        # Extend the list_time_strings in dict_rainfall_params one additional timestep
        # Added - 2024.11.22
        
        # Overwrite the "Values" tables
        str_dataset_path = '/Event Conditions/Meteorology/Precipitation/Imported Raster Data/Values'
        fn_overwrite_rainfall_hdf(dict_rainfall_params, array_flows_reshaped, str_hdf_to_edit, str_dataset_path)

        # Overwrite the "Values (Vertical)" tables
        str_dataset_path = '/Event Conditions/Meteorology/Precipitation/Imported Raster Data/Values (Vertical)'
        fn_overwrite_rainfall_hdf(dict_rainfall_params, array_flows_reshaped, str_hdf_to_edit, str_dataset_path)

        #print(f'Precipitation HDF updated: {str_hdf_to_edit}')
    else:
        print('Error: Unique points not created... No update')
# ~~~~~~~~~~~~~~~~~~~~~~~


# ------------
def fn_determine_ras_hdf(str_ras_path, int_num, str_type):
    """
    Determines the .hdf file corresponding to a .prj file and checks if it exists.

    Args:
        str_ras_path (str): The path to search for the .prj file.
        int_geom (int): The geometry index to construct the expected .hdf filename.
        str_type (str): String of desired file type (g,p,u...)

    Returns:
        str: The full path to the .hdf file if it exists, otherwise None.
    """
    # Iterate through files in the root directory
    for file in os.listdir(str_ras_path):
        if file.endswith('.prj'):
            # Extract the base name of the .prj file (without the extension)
            str_prj_file = file[:-4]

            # Construct the expected .hdf filename
            str_file_suffix = f".{str_type}{int_num:02d}.hdf"
            str_expected_geom = os.path.join(str_ras_path, str_prj_file + str_file_suffix)

            # Check if the file exists
            if os.path.exists(str_expected_geom):
                return str_expected_geom
    
    # If no matching .prj file or .hdf file is found, return None
    return None
# ------------


# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist" % arg)
    else:
        # File exists so return the directory
        return arg
        return open(arg, 'r')  # return an open file handle
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


# -------------
def fn_parse_tuple(value):
    # Remove any spaces and parentheses, then split by commas
    return tuple(map(int, value.strip('()').split(',')))
# -------------


# ----------------
def fn_str_to_bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in {'true', 't', '1'}:
        return True
    elif value.lower() in {'false', 'f', '0'}:
        return False
    else:
        raise argparse.ArgumentTypeError(f"Boolean value expected. Got '{value}'.")
# ----------------

# ---------
def fn_normalize_mixed_path(mixed_path: str) -> str:
    # Check if path starts with a Linux-style root
    if mixed_path.startswith("/"):
        # Convert backslashes to forward slashes
        return mixed_path.replace("\\", "/")
    return mixed_path  # Return as is if it's not Linux-style at the start
# ---------


# .........................................................
def fn_spawn_hecras_copies(str_config_file_path,
                           str_output_dir,
                           tup_hecras_files_to_copy,
                           b_print_output):
    
    
    # supress all warnings
    warnings.filterwarnings("ignore", category=UserWarning )
    
    print(" ")
    if b_print_output:
        print("+=================================================================+")
        print("|                SPAWN COPIES OF HEC-RAS MODELS                   |")
        print("|                Created by Andy Carter, PE of                    |")
        print("|             Center for Water and the Environment                |")
        print("|                 University of Texas at Austin                   |")
        print("+-----------------------------------------------------------------+")
    
        
        print("  ---(c) INPUT GLOBAL CONFIGURATION FILE: " + str_config_file_path)
        print("  ---(o) OUTPUT PATH: " + str_output_dir)
        print("  ---(f) FILES TO COPY (geom, unsteady flow, plan): " + ', '.join(map(str, tup_hecras_files_to_copy)))
        print("  ---[r] PRINT OUTPUT: " + str(b_print_output))
        
        print("===================================================================")
    else:
        print('Step 2: Spawn copies of HEC-RAS models')
    
    str_ras_path = os.path.join(str_output_dir, '00_base_hecras')
    str_model_hydrofabric_gpkg = os.path.join(str_output_dir, '01_model_hydrofabric', 'model_hydrofabric.gpkg')
    
    # create a folder
    str_folder_for_copies = os.path.join(str_output_dir, '02_model_copies')
    # If the folder does not exist, create it
    os.makedirs(str_folder_for_copies, exist_ok=True)
    
    # get the path to the required gXX.hdf
    str_source_geom_hdf = fn_determine_ras_hdf(str_ras_path, tup_hecras_files_to_copy[0], 'g')
    
    # Create the folder were the copies will be placed
    if not os.path.exists(str_folder_for_copies):
        os.makedirs(str_folder_for_copies)
    
    # --- Read variables from config.ini ---
    # Create a ConfigParser object
    config = configparser.ConfigParser()
    
    # Read the config.ini file
    config.read(str_config_file_path)
    # ---
    
    # +++ Copy the projection file +++
    # Get a list of filepaths to projection files in rasmap (xml)
    list_proj_filenames, list_rasmap_filenames = fn_find_rasmap_projection_filenames(str_ras_path)
    
    if len(list_proj_filenames) > 0:
        # Copy the first found prj file
        str_proj = list_proj_filenames[0]
        
        str_new_proj_folder_name = "projection"  #the folder place projection
        str_proj_folder_path = os.path.join(str_folder_for_copies,str_new_proj_folder_name)
        
        # Create the folder were the copies will be placed
        if not os.path.exists(str_proj_folder_path):
            os.makedirs(str_proj_folder_path)
        
        str_proj = fn_normalize_mixed_path(str_proj)
        
        str_copy_to_filepath = os.path.join(str_proj_folder_path, os.path.basename(str_proj))
        
        # Copy from str_proj to str_copy_to_filepath
        shutil.copy(str_proj, str_copy_to_filepath)
        
        # get the relative path of the new projection file
        str_proj_relativepath = "..\\" + str_new_proj_folder_name + "\\" + os.path.basename(str_proj)
        #print(str_proj_relativepath)
        
    else:
        print('No projection found in rasmap')
        
    str_rasmap_copy_from = list_rasmap_filenames[0]
    # +++
    
    
    # ---
    # Copy the base model terrain to a 02_model_copies folder
    str_new_terain_folder_name = "source_terrain"  #the folder to create and place terrain
    str_out = fn_copy_source_terrain(str_folder_for_copies, str_new_terain_folder_name, str_source_geom_hdf)
    str_terrain_relativepath = "..\\" + str_new_terain_folder_name + "\\" + str_out
    # ---
    
    # ---
    # Get the variables in the [02_copy_geom] section
    if '02_copy_geom' in config:
        section = config['02_copy_geom']
    
        # Read variables in the section
        int_time_rolling_avg = section.getint('int_time_rolling_avg', 0)
        int_buffer_time = section.getint('int_buffer_time', 0)
        str_simulation_start_date = section.get('str_simulation_start_date', '')
        flt_delta_q = section.getfloat('flt_delta_q', 0.0)
        flt_cell_edge_length = section.getfloat('flt_cell_edge_length', 0.0)
    else:
        print("[02_copy_geom] section not found in the config file.")
    # ---
    
    # --- Read the geopackage ---
    gdf_area = gpd.read_file(str_model_hydrofabric_gpkg, layer='00_area_2d')
    gdf_mainstems = gpd.read_file(str_model_hydrofabric_gpkg, layer='02_flow_points')
    
    dict_area_2d = gdf_area.iloc[0].to_dict()
    str_2d_area_name = dict_area_2d['area_2d_name']
    
    # Get a list of the HEC-RAS projects in a given directory
    list_proj_title_files = fn_list_of_ras_projects(str_ras_path)
    
    if len(list_proj_title_files) > 0:
        # select the first HEC-RAS Project found
        str_project_path = list_proj_title_files[0]
    
    # Define the corresponding headers and identifiers
    file_headers = ['Geom File=', 'Unsteady File=', 'Plan File=']
    file_identifiers = ['g', 'u', 'p']
    
    # Create the list of tuples
    list_hecras_files = [
        (file_headers[i], file_identifiers[i], tup_hecras_files_to_copy[i])
        for i in range(len(tup_hecras_files_to_copy))
    ]
    
    # projection parameters
    str_input_crs = f"EPSG:{gdf_mainstems.crs.to_epsg()}"
    
    
    # TODO -- 2025.01.10 -- this need parallel and TQDM
    
    # Loop through all the mainstems
    # Set up the tqdm bar
    with tqdm(total=len(gdf_mainstems), 
              desc='   -- Spawn HEC-RAS',
              bar_format="{desc}:({n_fmt}/{total_fmt})|{bar}| {percentage:.1f}%",
              ncols=75) as pbar:
    
        for index, row in gdf_mainstems.iterrows():
            
            # ---- example input: gdf_mainstems.iloc[17] ----
            dict_mainstem = gdf_mainstems.iloc[index].to_dict()
        
            dict_area_2d = gdf_area.iloc[0].to_dict()
            str_2d_area_name = dict_area_2d['area_2d_name']
        
            # Set the run time for stepped flows (this is the time for each constant flow)
            int_hour_count = int(dict_mainstem['travel_time_hr']) + 1 + int_time_rolling_avg + int_buffer_time
        
            # Add the constant flow to the dict_mainstem dictionary
            # this is so a run name can be created
            dict_mainstem['time_sim_hr'] = int_hour_count
        
            # create the run names (note: set flt_upper_flow to None for single flow simulation)
            str_run_name,str_short_name = fn_build_plan_names(flt_delta_q, dict_mainstem)
        
            str_new_folder = os.path.join(str_folder_for_copies,str_run_name)
        
            # Create the folder
            if not os.path.exists(str_new_folder):
                os.makedirs(str_new_folder)
        
            # Convert the flow data to dictionary
            dict_flows = {
                "lower_flow": dict_mainstem['q_lower_limit'],
                "upper_flow": dict_mainstem['q_upper_limit'],
                "delta_q": flt_delta_q,
                "hour_count": int_hour_count,
                "run_name": str_run_name,
                "short_name": str_short_name
            }
        
            # From  'flt_lower_flow,flt_upper_flow,flt_delta_q' determine 
            # the number of simulations... then determine the int_sim_time_hr ... using int_hour_count
            
            int_num_of_flows = len(fn_generate_step_integer(dict_flows['lower_flow'], dict_flows['upper_flow'], flt_delta_q))
            
            # This is one time step too long as last step has no precip.  Fixing this (2024.11.22)
            int_sim_time_hr = int_hour_count * (int_num_of_flows - 1)
            
            # spawn a copy of the desired files
            for tup_of_file_params in list_hecras_files:
                fn_copy_ras_file(tup_of_file_params, str_project_path, str_run_name, str_new_folder)
            
            # create a project file for each run
            fn_create_prj(str_run_name, str_project_path, str_new_folder)
            
            # Adjusting the plan
            fn_adjust_plan(str_run_name,str_new_folder,dict_flows,int_sim_time_hr,str_simulation_start_date)
            
            # Adjusting the unsteady flow
            fn_adjust_unsteady_flow(str_run_name, str_new_folder, dict_flows)
            
            # Adjusting the geometry
            fn_adjust_geom(str_run_name, str_new_folder, dict_flows)
            
            # Adjusting the geometry.hdf
            fn_adjust_geom_hdf(str_run_name, str_new_folder, str_terrain_relativepath)
            
            # --- rasmap copy and adjustments ---
            # File and folder name of rasmap to create/copy
            str_rasmap_copyto = os.path.join(str_new_folder, str_run_name + ".rasmap")
            
            # Copy from str_proj to str_rasmap_copyto
            shutil.copy(str_rasmap_copy_from, str_rasmap_copyto)
            
            # Augment this run's .rasmap
            fn_replace_rasmap_elements(str_rasmap_copyto, str_run_name, 'g')
            fn_replace_rasmap_elements(str_rasmap_copyto, str_run_name, 'p')
            fn_replace_rasmap_elements(str_rasmap_copyto, str_run_name, 'u')
            
            # Replace the .rasmap terrain and projection relative paths
            fn_update_rasmap_xml(str_rasmap_copyto, str_proj_relativepath, str_terrain_relativepath)
            # ---
            
            str_geom_hdf = os.path.join(str_new_folder,str_run_name + '.g01.hdf')
            projection_wkt = fn_get_projection_wkt(str_geom_hdf)
            crs = CRS.from_wkt(projection_wkt)
            epsg_code = crs.to_epsg()
            str_output_crs = f"EPSG:{epsg_code}"
            str_hdf_to_edit = os.path.join(str_new_folder,str_run_name + '.u01.hdf')
            
            # --- Format the start time to HDF compliant format ---
            str_start_time_rev = str_simulation_start_date + ",0000"
            
            # Convert to desired format
            str_formatted_date_for_hdf = datetime.strptime(str_start_time_rev, "%d%b%Y,%H%M").strftime("%Y-%m-%d %H:%M:%S")
            # ---
            
            # Create combined flat dictionary for precip revisions
            dict_for_precip = {
            **dict_mainstem,
            **dict_flows,
            "input_crs": str_input_crs,
            "output_crs": str_output_crs,
            "cell_edge_length": flt_cell_edge_length,
            "str_hdf_to_edit": str_hdf_to_edit,
            "str_date_for_hdf": str_formatted_date_for_hdf}
            
            # revise the rainfall (precip in u01.hdf)
            fn_revise_rainfall_hdf(dict_for_precip, projection_wkt)
            
            # Update the tqdm bar after processing each row
            pbar.update(1)
    print("+-----------------------------------------------------------------+")
# .........................................................



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()
    
    parser = argparse.ArgumentParser(description='========= SPAWN COPIES OF HEC-RAS MODELS =========')
    
    parser.add_argument('-c',
                        dest = "str_config_file_path",
                        help=r'REQUIRED: Global configuration filepath Example: C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini',
                        required=True,
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-o',
                        dest = "str_output_dir",
                        help=r'REQUIRED: output directory of RAS2FIM-2D output Example: E:\ras2fim_test_20250107',
                        required=True,
                        metavar='DIR',
                        type=str)
    
    parser.add_argument('-f',
                        dest = "tup_hecras_files_to_copy",
                        help=r'REQUIRED: Files to copy from source HEC-RAS (geom, unsteady flow, plan) Example: "(1,1,1)"',
                        required=True,
                        metavar='TUPLE',
                        type=fn_parse_tuple)
    
    parser.add_argument('-r',
                        dest = "b_print_output",
                        help=r'OPTIONAL: Print output messages Default: True',
                        required=False,
                        default=True,
                        metavar='T/F',
                        type=fn_str_to_bool)
    
    args = vars(parser.parse_args())
    
    str_config_file_path = args['str_config_file_path']
    str_output_dir = args['str_output_dir']
    tup_hecras_files_to_copy = args['tup_hecras_files_to_copy']
    b_print_output = args['b_print_output']
    
    fn_spawn_hecras_copies(str_config_file_path,
                           str_output_dir,
                           tup_hecras_files_to_copy,
                           b_print_output)

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)
    
    print('Compute Time: ' + str(time_pass))
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~