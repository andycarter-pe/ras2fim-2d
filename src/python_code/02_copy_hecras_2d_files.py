# Script 02 - Copy an existing HEC-RAS Geometry and set the "Infow" boundary 
# condition to a desired cell location - Establish a "firehose" of constant
# flow from this cell.  Create a copy of a 2D HEC-RAS geometry, unsteady flow
# and plan file.  Flows can be a single flow or a range of flows.
#
#
# Created by: Andy Carter, PE
# Created - 2024.05.11
# Last Revised - 2024.08.01

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
# ************************************************************


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
# argparse custom types for tuple processing
# Custom type to parse a tuple of integers
def tuple_of_ints(arg):
    try:
        values = tuple(map(int, arg.split(',')))
        return values
    except Exception as e:
        raise argparse.ArgumentTypeError(f"Invalid input: {e}")
        
def tuple_of_ints_or_none(arg):
    try:
        values = tuple(int(x) if x.strip().lower() != 'none' else None for x in arg.split(','))
        return values
    except Exception as e:
        raise argparse.ArgumentTypeError(f"Invalid input: {e}")
# -------------


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


# --------------
def fn_build_plan_names(flt_lower_flow, flt_upper_flow, flt_delta_q, dict_mainstem):
    str_run_name, str_short_plan = None, None
    
    # Example:  str_run_name = '1884413_wb-2410249_wb-2410261_29-hr_14100-cfs_to_45000-cfs_step_2350-cfs'
    
    # systematic naming of plans (if only one flow, flt_upper_flow should be None)
    is_single_flow = False
    
    try:
        if flt_upper_flow == None:
            is_single_flow = True
            
        str_mainstem = str(int(dict_mainstem['mainstem']))
        str_start_node = dict_mainstem['id_start_node']
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
                print(f"Valid HEC-RAS Project: {file_path}")

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


# ....................
def fn_copy_geom_file(list_proj_title_files, int_geom_to_copy, dict_flows):
    
    str_copy_to_geom_full_path = None
    str_copy_to_geom_hdf_full_path = None

    # Splitting the path into folder and filename
    str_prj_folder, str_filename = os.path.split(list_proj_title_files[0])

    # Splitting the filename and its extension
    str_file_only, str_extension = os.path.splitext(str_filename)

    # determine the geometry files that are in the project
    list_geom_names = fn_list_filename_from_ras_prj(list_proj_title_files[0], "Geom File=")
    list_geom_int = fn_extract_numbers_from_strings(list_geom_names)

    list_geom_hdf_fullpath = []
    list_geom_fullpath = []
    list_binary_fullpath = []
    list_config_fullpath = []

    for geom_item in list_geom_names:
        str_hdf_geom_file = str_file_only + "." + geom_item + ".hdf"
        str_geom_file = str_file_only + "." + geom_item
        str_binary_preprocessor_file = str_file_only + "." + "c" + geom_item[1:]
        str_config_preprocessor_file = str_file_only + "." + "x" + geom_item[1:]

        # Combine the folder path and filename
        str_hdf_full_path = os.path.join(str_prj_folder, str_hdf_geom_file)
        str_geom_full_path = os.path.join(str_prj_folder, str_geom_file)
        str_binary_preprocessor_full_path = os.path.join(str_prj_folder, str_binary_preprocessor_file)
        str_config_preprocessor_full_path = os.path.join(str_prj_folder, str_config_preprocessor_file)

        list_geom_hdf_fullpath.append(str_hdf_full_path)
        list_geom_fullpath.append(str_geom_full_path)
        list_binary_fullpath.append(str_binary_preprocessor_full_path)
        list_config_fullpath.append(str_config_preprocessor_full_path)

    list_b_hdf_exists = []

    list_b_hdf_exists = fn_list_of_file_exists(list_geom_hdf_fullpath)
    list_b_geom_exists = fn_list_of_file_exists(list_geom_fullpath)
    list_b_binary_exists = fn_list_of_file_exists(list_binary_fullpath)
    list_b_config_exists = fn_list_of_file_exists(list_config_fullpath)
    

    # Assuming all lists have the same length
    data = {
        'geom_int': list_geom_int,
        'geom_name': list_geom_names,
        'hdf_path': list_geom_hdf_fullpath,
        'hdf_exists': list_b_hdf_exists,
        'geom_path': list_geom_fullpath,
        'geom_exists': list_b_geom_exists,
        'binary_path': list_binary_fullpath,
        'binary_exists': list_b_binary_exists,
        'config_path': list_config_fullpath,
        'config_exists': list_b_config_exists
    }

    df = pd.DataFrame(data)
    
    # Filter the DataFrame
    df_filtered = df[(df['geom_int'] == int_geom_to_copy) & (df['hdf_exists']) & (df['geom_exists']) & (df['binary_exists']) & (df['config_exists'])]
    
    # Check if any rows match the condition
    if not df_filtered.empty:
        # If there's at least one row matching the condition
        print('---------------------------')
        #print('Valid Geometry Match Found')

        # Determine the highest number in geom_int
        highest_geom_int = df['geom_int'].max()

        # Add one to the highest number
        next_geom_int = highest_geom_int + 1

        # Convert next_geom_int to string with leading zero padding if necessary
        next_geom_int_str = '{:02d}'.format(next_geom_int)

        # copy geom file
        str_copy_from = df_filtered.iloc[0]['geom_path']
        str_copy_to_geom = str_file_only + ".g" + next_geom_int_str
        str_copy_to_geom_full_path = os.path.join(str_prj_folder, str_copy_to_geom)

        shutil.copy(str_copy_from, str_copy_to_geom_full_path)
        print(f"Copied {str_copy_from} to {str_copy_to_geom_full_path}")

        # copy the hdf geom file
        str_copy_from_hdf = df_filtered.iloc[0]['hdf_path']
        str_copy_to_geom = str_file_only + ".g" + next_geom_int_str + ".hdf"
        str_copy_to_geom_hdf_full_path = os.path.join(str_prj_folder, str_copy_to_geom)

        shutil.copy(str_copy_from_hdf, str_copy_to_geom_hdf_full_path)
        print(f"Copied {str_copy_from_hdf} to {str_copy_to_geom_hdf_full_path}")
        
        # copy the binary precrocessed geom file
        str_copy_from_binary = df_filtered.iloc[0]['binary_path']
        str_copy_to_binary = str_file_only + ".c" + next_geom_int_str
        str_copy_to_binary_full_path = os.path.join(str_prj_folder, str_copy_to_binary)

        shutil.copy(str_copy_from_binary, str_copy_to_binary_full_path)
        print(f"Copied {str_copy_from_binary} to {str_copy_to_binary_full_path}")
        
        # copy the binary precrocessed geom file
        str_copy_from_config = df_filtered.iloc[0]['config_path']
        str_copy_to_config = str_file_only + ".x" + next_geom_int_str
        str_copy_to_config_full_path = os.path.join(str_prj_folder, str_copy_to_config)

        shutil.copy(str_copy_from_config, str_copy_to_config_full_path)
        print(f"Copied {str_copy_from_config} to {str_copy_to_config_full_path}")

        # ~~~~~~~~~~~~~~~~~~~~~
        # add the geom to the project file

        # Read the contents of the file
        with open(list_proj_title_files[0], 'r') as file:
            lines = file.readlines()

        # Find the index of the last occurrence of a line starting with "Geom File="
        last_geom_index = -1
        for i in range(len(lines)):
            if lines[i].strip().startswith("Geom File="):
                last_geom_index = i

        # Insert a new line after the last occurrence of "Geom File="
        if last_geom_index != -1:
            lines.insert(last_geom_index + 1, "Geom File=g" + next_geom_int_str + "\n")

        # Write the modified content back to the file
        with open(list_proj_title_files[0], 'w') as file:
            file.writelines(lines)
            
        # ---------------------
        # Within the geom file, change the plan 

        # Read the file content
        with open(str_copy_to_geom_full_path, "r") as file:
            lines = file.readlines()

        # --- Append the plan title in the plan file
        # Find indices and lines that start with "Flow File="
        list_tup_match = [(index, line.strip()) for index, line in enumerate(lines) if re.match(r'^Geom Title=', line)]

        # Replace lines in the list
        for index, _ in list_tup_match:
            lines[index] = "Geom Title=" + dict_flows['run_name'] + "_geom" + "\n"
            
        # Write the modified lines back to the file
        with open(str_copy_to_geom_full_path, "w") as file:
            file.writelines(lines)

    else:
        # Otherwise, print an error statement
        print("Error: Geometry HDF, gXX, cXX and/or xXX not available.")
        # TODO - 2024.08.01 - Could add additional detail here - are we missing preprocess files?
        
    return(str_copy_to_geom_full_path, str_copy_to_geom_hdf_full_path)
# ....................


# ++++++++++++++++++++++++++++
def fn_get_gdf_of_cells_from_list(hdf_file_path, list_unique_indices_sorted, str_2darea_name):

    # Specify the HDF5 file path and group path
    area_2D_path = '/Geometry/2D Flow Areas/'

    str_hdf_folder_2darea = area_2D_path + str_2darea_name + '/'

    # Location of Face Point Coordinates in HDF5
    str_facepoint_coords = str_hdf_folder_2darea + 'FacePoints Coordinate'

    # Open the HDF5 file
    with h5py.File(hdf_file_path, 'r') as hdf_file:
        # Extract X and Y coordinates
        x_coordinates = hdf_file[str_facepoint_coords][:, 0]
        y_coordinates = hdf_file[str_facepoint_coords][:, 1]

    # Create a pandas DataFrame
    df_facepoints = pd.DataFrame({'X': x_coordinates, 'Y': y_coordinates})

    # Location of Indices of face points making up the cells
    str_cells_facepoint_indexes = str_hdf_folder_2darea + 'Cells FacePoint Indexes'

    # Open the HDF5 file
    with h5py.File(hdf_file_path, 'r') as hdf_file:
        # Extract face points coordinate data
        facepoints_data = hdf_file[str_cells_facepoint_indexes][:]

        # Extract the projection
        projection_wkt = hdf_file.attrs['Projection'].decode('utf-8')

    # Create a GeoDataFrame to store the polygons
    geometry = []
    indices = []

    for row_idx, row in enumerate(facepoints_data):
        polygon_coords = []

        if row_idx in list_unique_indices_sorted:
            for idx in row:
                if idx != -1:
                    x = df_facepoints.loc[idx, 'X']
                    y = df_facepoints.loc[idx, 'Y']
                    polygon_coords.append((x, y))
            # Check if the polygon has at least 3 points (needed to create a polygon)
            if len(polygon_coords) >= 3:
                # Connect to the first point to close the polygon
                polygon_coords.append(polygon_coords[0])
                geometry.append(Polygon(polygon_coords))
                indices.append(row_idx)  # Append the row index as the cell index

    # Create a GeoDataFrame
    gdf_cells = gpd.GeoDataFrame(geometry=geometry, index=indices, columns=['geometry'], crs=projection_wkt)
    
    return(gdf_cells)
# ++++++++++++++++++++++++++++

# ---------------------------
def fn_create_line_inside_polygon(shp_polygon):
    
    # Define your original polygon
    #original_polygon = gdf_firehose_cell.iloc[0]['geometry']
    original_polygon = shp_polygon

    # Calculate the length of the shortest side
    shortest_side_length = min(original_polygon.length for side in original_polygon.exterior.coords[:-1])

    # Offset the polygon internally by 1% of the length of the shortest side
    offset_distance = 0.01 * shortest_side_length
    offset_polygon = original_polygon.buffer(-offset_distance)

    # Extract the exterior boundary of the offset polygon
    offset_exterior = offset_polygon.exterior

    # Find the shortest side of the offset polygon
    shortest_side_length = float('inf')
    shortest_side = None
    for i in range(len(offset_exterior.coords) - 1):
        p1 = offset_exterior.coords[i]
        p2 = offset_exterior.coords[i + 1]
        length = LineString([p1, p2]).length
        if length < shortest_side_length:
            shortest_side_length = length
            shortest_side = (p1, p2)

    # Create a LineString representing the shortest side
    shortest_side_line = LineString(shortest_side)
    
    return(shortest_side_line)
# ---------------------------


# --------------
def fn_number_round_digits(number, int_requested_digits):
    integer_digits = int(math.log10(abs(number))) + 1 if number != 0 else 1
    decimal_places = max(0, int_requested_digits - integer_digits - 1)  # Ensure decimal_places is non-negative
    formatted_number = '{:.{}f}'.format(number, decimal_places)
    return formatted_number
# --------------


# --------------
def fn_format_coords(coords, int_requested_digits):
    formatted_coords = []
    for pair in coords:
        formatted_pair = []
        for num in pair:
            formatted_num = fn_number_round_digits(num, int_requested_digits)
            formatted_pair.append(formatted_num)
        formatted_coords.append(formatted_pair)
    return formatted_coords
# --------------


# --------------
def fn_midpoint(coords):
    # this assumes only two points
    # Extracting coordinates
    x1, y1 = coords[0]
    x2, y2 = coords[1]
    
    # Calculating midpoint
    mid_x = (x1 + x2) / 2
    mid_y = (y1 + y2) / 2
    
    return (mid_x, mid_y)
# --------------

# ======================
def fn_get_index_of_bc_to_edit(str_boundary_to_edit, str_geom_path):
    # Edit the Internal boundary condition
    target_line = "BC Line Name=" + str_boundary_to_edit

    # Open the file and read its content
    with open(str_geom_path, 'r') as file:
        list_lines = file.readlines()

    # Find the index of the line starting with the target string
    index = None
    for i, line in enumerate(list_lines):
        if line.startswith(target_line):
            index = i
            break

    # Extract the following six lines into a list
    #list_boundary_lines = list_lines[index+2:index+7]
    return(index)
# ======================

# *****************************
def fn_build_internal_boundary_text(hdf_file_path, list_unique_indices_sorted, str_2d_area_name):
    
    gdf_firehose_cell = fn_get_gdf_of_cells_from_list(hdf_file_path, list_unique_indices_sorted, str_2d_area_name )

    # Create a shapley line inside the 'source' emmiter cell
    shp_line = fn_create_line_inside_polygon(gdf_firehose_cell.iloc[0]['geometry'])

    # Create a numpy array of the the shapel line coordinates
    coords = np.array(shp_line.coords)

    tup_mid_coords = fn_midpoint(coords)

    # Converting tuple to list
    list_mid_coords = list(tup_mid_coords)

    list_mid_coords_formatted = []
    for item in list_mid_coords:
        str_format = fn_number_round_digits(item, 16)
        list_mid_coords_formatted.append(str_format)

    # Format the boundarline coords as strings
    formatted_coords = fn_format_coords(coords, 16)

    first_pair = formatted_coords[0]
    last_pair = formatted_coords[-1]
    int_point_len = len(formatted_coords)

    # Join the formatted coordinates into a string
    str_start_point = 'BC Line Start Position= {} , {} \n'.format(*first_pair)
    str_mid_point = 'BC Line Middle Position= {} , {} \n'.format(*list_mid_coords_formatted)
    str_last_point = 'BC Line End Position= {} , {} \n'.format(*last_pair)
    str_line_arc = f'BC Line Arc= {int_point_len} \n'

    # Flatten the list of boundary condition points
    str_point_list = [item for sublist in formatted_coords for item in sublist]
    str_point_list += ' \n'

    # Join the elements into one continuous string
    str_line_points = ''.join(str_point_list)

    list_new_boundary_lines = [str_start_point,str_mid_point,str_last_point,str_line_arc,str_line_points]
    
    return(list_new_boundary_lines)
# *****************************

# -------------------
def fn_replace_boundary_lines_in_file(file_path, list_new_boundary_lines, index):
    # Read the contents of the file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Calculate the start and end indices for replacement
    start_index = index + 2
    end_index = index + 7

    # Perform the replacement
    lines[start_index:end_index] = list_new_boundary_lines

    # Write the modified contents back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)
# -------------------

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


# ---------------------------
def fn_create_hydrograph_lines(dict_flows):

    # get flow variables from the flow dictionary
    flt_lower_flow = dict_flows['lower_flow']
    flt_upper_flow = dict_flows['upper_flow']
    flt_delta_q = dict_flows['delta_q']
    int_hour_count = dict_flows['hour_count']
    
    list_flow_rows = []

    # get a list of flow steps
    list_flow_steps_int = fn_generate_step_integer(flt_lower_flow, flt_upper_flow, flt_delta_q)

    list_flows = []
    for item in list_flow_steps_int:
        for i in range(int_hour_count):
            list_flows.append(item)

    # Calculate the number of groups of 10
    int_groups_of_10 = len(list_flows) // 10

    # Calculate the remainder
    int_remainder = len(list_flows) % 10

    str_first_line = f"Flow Hydrograph= {len(list_flows)} \n"
    list_flow_rows.append(str_first_line)

    for i in range(int_groups_of_10):
        #print(i)
        list_current_row = list_flows[i*10:i*10+10]

        str_current_row = ""
        for item in list_current_row:
            # (8 characters with right justified with padding)
            str_flow = str(item).rjust(8)
            str_current_row += str_flow
        str_current_row += '\n'

        list_flow_rows.append(str_current_row)

    list_last_row = list_flows[int_groups_of_10*10:]
    str_current_row = ""
    for item in list_last_row:
        # (8 characters with right justified with padding)
        str_flow = str(item).rjust(8)
        str_current_row += str_flow
    str_current_row += '\n'
    list_flow_rows.append(str_current_row)

    return(list_flow_rows, len(list_flows))
# ---------------------------


# -------------------
def fn_copy_unsteady_flow_file(list_proj_title_files, int_u_to_copy, dict_flows, str_boundary_to_edit):
    # determine the unsteady flow files that are in the project

    list_unsteady_names = fn_list_filename_from_ras_prj(list_proj_title_files[0],"Unsteady File=")
    list_unsteady_int = fn_extract_numbers_from_strings(list_unsteady_names)

    # Splitting the path into folder and filename
    str_prj_folder, str_filename = os.path.split(list_proj_title_files[0])

    # Splitting the filename and its extension
    str_file_only, str_extension = os.path.splitext(str_filename)

    list_unsteady_fullpath = []

    for unsteady_item in list_unsteady_names:
        str_unsteady_file = str_file_only + "." + unsteady_item

        # Combine the folder path and filename
        str_unsteady_full_path = os.path.join(str_prj_folder, str_unsteady_file)

        list_unsteady_fullpath.append(str_unsteady_full_path)

    list_b_unsteady_exists = fn_list_of_file_exists(list_unsteady_fullpath)

    # Assuming all lists have the same length
    data = {
        'unsteady_int': list_unsteady_int,
        'unsteady_name': list_unsteady_names,
        'unsteady_path': list_unsteady_fullpath,
        'unsteady_exists': list_b_unsteady_exists,
    }

    df_unsteady = pd.DataFrame(data)

    # Filter the DataFrame
    df_filtered = df_unsteady[(df_unsteady['unsteady_int'] == int_u_to_copy) & (df_unsteady['unsteady_exists'])]

    # Check if any rows match the condition
    if not df_filtered.empty:
        # If there's at least one row matching the condition
        #print('Valid Unsteady Flow Match Found')
        print('---------------------------')

        # Determine the highest number in unsteady_int
        highest_u_int = df_unsteady['unsteady_int'].max()

        # Add one to the highest number
        next_u_int = highest_u_int + 1

        # Convert next_u_int to string with leading zero padding if necessary
        next_u_int_str = '{:02d}'.format(next_u_int)

        # copy geom file
        str_copy_from = df_filtered.iloc[0]['unsteady_path']
        str_copy_to_u = str_file_only + ".u" + next_u_int_str
        str_copy_to_u_full_path = os.path.join(str_prj_folder, str_copy_to_u)

        shutil.copy(str_copy_from, str_copy_to_u_full_path)
        print(f"Copied {str_copy_from} to {str_copy_to_u_full_path}")

        # ------------------
        # Add the new unsteady file to the HEC-RAS project

        # Read the contents of the file
        with open(list_proj_title_files[0], 'r') as file:
            lines = file.readlines()

        # Find the index of the last occurrence of a line starting with "Geom File="
        last_u_index = -1
        for i in range(len(lines)):
            if lines[i].strip().startswith("Unsteady File="):
                last_u_index = i

        # Insert a new line after the last occurrence of "Geom File="
        if last_u_index != -1:
            lines.insert(last_u_index + 1, "Unsteady File=u" + next_u_int_str + "\n")

        # Write the modified content back to the file
        with open(list_proj_title_files[0], 'w') as file:
            file.writelines(lines)
    
    # Determine where the boundary conditions begin
    if not os.path.exists(str_copy_to_u_full_path):
        print(f"File {str_copy_to_u_full_path} not found.")
    else:
        # Open the file and read its content
        with open(str_copy_to_u_full_path, "r") as file:
            # Read lines from the file
            lines = file.readlines()

            # Find indices and lines that start with "Boundary Location="
            matching_indices_and_lines = [(index, line.strip()) for index, line in enumerate(lines) if re.match(r'^Boundary Location=', line)]

            # Extract only the indices
            list_matching_indices = [index for index, _ in matching_indices_and_lines]

            # extract the lines text
            list_boundary_location_line = [value  for _, value in matching_indices_and_lines]

            # Extracting the seventh value 'Boundary Name' from each row
            list_boundary_names = [row.split(',')[7].strip() for row in list_boundary_location_line]

            # Create a df for the boundary condition rows
            data = {
                'boundary_int': list_matching_indices,
                'boundary_name': list_boundary_names
            }

            df_boundary = pd.DataFrame(data)
            
    with open(str_copy_to_u_full_path, "r") as file:
        # Read lines from the file
        lines = file.readlines()

        # Find indices and lines that start with "Boundary Location=" --- list of tuples
        list_tup_match = [(index, line.strip()) for index, line in enumerate(lines) if re.match(r'^Flow Hydrograph=', line)]
        
    # Create new columns
    df_boundary['is_hydro'] = False
    df_boundary['hydro_int_start'] = None
    df_boundary['hydro_count_int'] = None

    # Finding nearest boundary_int for each item in list_tup_match
    for item in list_tup_match:
        value_to_match = item[0]
        hydro_count_int = int(item[1].split('=')[-1])  # Parsing hydro_count_int from the second value in the tuple
        nearest_boundary_int = fn_find_nearest_boundary_int(value_to_match, df_boundary['boundary_int'])
        if nearest_boundary_int is not None:
            df_boundary.loc[df_boundary['boundary_int'] == nearest_boundary_int, 'is_hydro'] = True
            df_boundary.loc[df_boundary['boundary_int'] == nearest_boundary_int, 'hydro_int_start'] = value_to_match
            df_boundary.loc[df_boundary['boundary_int'] == nearest_boundary_int, 'hydro_count_int'] = hydro_count_int
            
    # determine the rows to delete from the original unsteady flow file
    df_inflow_rows = df_boundary[df_boundary['boundary_name'] == str_boundary_to_edit]

    # Check if any rows match the condition
    if not df_inflow_rows.empty:
        # If there's at least one row matching the condition
        print('---------------------------')
        #print('Valid Boundary Condition Found')

        # Grab the first row
        first_row = df_inflow_rows.iloc[0]

        # Access the "hydro_int_start" and "hydro_count_int" columns
        hydro_int_start = first_row['hydro_int_start']
        hydro_count_int = first_row['hydro_count_int']

        # Determine number of rows to delete

        # Calculate the number of groups of 10
        int_rows_to_delete = hydro_count_int // 10 + 2
    
    # get the unsteady flow rows and hours to simulate
    list_flow_rows, int_sim_time_hr = fn_create_hydrograph_lines(dict_flows)
    
    # Update the flow hydrograph for the selected boundary condition
    with open(str_copy_to_u_full_path, "r") as file:
        # Read lines from the file
        lines = file.readlines()

        int_end_index = hydro_int_start + int_rows_to_delete

        # Delete the rows from hydro_int_start to end_index
        del lines[hydro_int_start:int_end_index]

        # Insert lines from list_flow_rows at hydro_int_start index
        lines[hydro_int_start:hydro_int_start] = list_flow_rows

    # Write the modified lines back to the file
    with open(str_copy_to_u_full_path, "w") as file:
        file.writelines(lines)
        
    # Enforce a 1HOUR Interval ("Inverval=1HOUR") in the unsteady flow file

    # Read the file content
    with open(str_copy_to_u_full_path, "r") as file:
        lines = file.readlines()

    # Find indices and lines that start with "Boundary Location="
    list_tup_match = [(index, line.strip()) for index, line in enumerate(lines) if re.match(r'^Interval=', line)]

    # Replace lines in the list
    for index, _ in list_tup_match:
        lines[index] = "Interval=1HOUR\n"  # Replace the line with "Interval=1HOUR"

    # Write the modified lines back to the file
    with open(str_copy_to_u_full_path, "w") as file:
        file.writelines(lines)
        
    # --------------------- 
    # Change the name of the Unsteady Flow Title
    # Read the file content
    with open(str_copy_to_u_full_path, "r") as file:
        lines = file.readlines()

    str_run_name = dict_flows['run_name']
    
    # Find indices and lines that start with "Flow Title="
    list_tup_match = [(index, line.strip()) for index, line in enumerate(lines) if re.match(r'^Flow Title=', line)]

    # Replace lines in the list
    for index, _ in list_tup_match:
        lines[index] = "Flow Title=" + str_run_name + '_unsteady_flow' + "\n"

    # Write the modified lines back to the file
    with open(str_copy_to_u_full_path, "w") as file:
        file.writelines(lines)
        
    return(str_copy_to_u_full_path, int_sim_time_hr)
# -------------------


# ~~~~~~~~~~~~~~~~~~~~~~~~
def fn_copy_plan_file(list_proj_title_files, int_p_to_copy,
                      next_geom_int_str, next_u_int_str,
                      dict_flows):

    # determine the geometry files that are in the project
    list_plan_names = fn_list_filename_from_ras_prj(list_proj_title_files[0], "Plan File=")
    list_plan_int = fn_extract_numbers_from_strings(list_plan_names)

    # Splitting the path into folder and filename
    str_prj_folder, str_filename = os.path.split(list_proj_title_files[0])

    # Splitting the filename and its extension
    str_file_only, str_extension = os.path.splitext(str_filename)

    list_plan_fullpath = []

    for plan_item in list_plan_names:
        str_plan_file = str_file_only + "." + plan_item

        # Combine the folder path and filename
        str_plan_full_path = os.path.join(str_prj_folder, str_plan_file)

        list_plan_fullpath.append(str_plan_full_path)

    list_b_plan_exists = fn_list_of_file_exists(list_plan_fullpath)

    # Assuming all lists have the same length
    data = {
        'plan_int': list_plan_int,
        'plan_name': list_plan_names,
        'plan_path': list_plan_fullpath,
        'plan_exists': list_b_plan_exists,
    }

    df_plan = pd.DataFrame(data)

    # Filter the DataFrame
    df_filtered = df_plan[(df_plan['plan_int'] == int_p_to_copy) & (df_plan['plan_exists'])]

    # Check if any rows match the condition
    if not df_filtered.empty:
        # If there's at least one row matching the condition
        print('---------------------------')
        #print('Valid Plan Match Found')

        # Determine the highest number in plan_int
        highest_p_int = df_plan['plan_int'].max()

        # Add one to the highest number
        next_p_int = highest_p_int + 1

        # Convert next_p_int to string with leading zero padding if necessary
        next_p_int_str = '{:02d}'.format(next_p_int)

        # copy plan file
        str_copy_from = df_filtered.iloc[0]['plan_path']
        str_copy_to_p = str_file_only + ".p" + next_p_int_str
        str_copy_to_p_full_path = os.path.join(str_prj_folder, str_copy_to_p)

        shutil.copy(str_copy_from, str_copy_to_p_full_path)
        print(f"Copied {str_copy_from} to {str_copy_to_p_full_path}")

        # ------------------
        # Add the new plan file to the HEC-RAS project

        # Read the contents of the file
        with open(list_proj_title_files[0], 'r') as file:
            lines = file.readlines()

        # Find the index of the last occurrence of a line starting with "Plan File="
        last_p_index = -1
        for i in range(len(lines)):
            if lines[i].strip().startswith("Plan File="):
                last_p_index = i

        # Insert a new line after the last occurrence of "Plan File="
        if last_p_index != -1:
            lines.insert(last_p_index + 1, "Plan File=p" + next_p_int_str + "\n")

        # Write the modified content back to the file
        with open(list_proj_title_files[0], 'w') as file:
            file.writelines(lines)

    # ---------------------
    # Within the plan file, adjust the selected geometry file

    # Read the file content
    with open(str_copy_to_p_full_path, "r") as file:
        lines = file.readlines()

    # Find indices and lines that start with "Geom File="
    list_tup_match = [(index, line.strip()) for index, line in enumerate(lines) if re.match(r'^Geom File=', line)]

    # Replace lines in the list
    for index, _ in list_tup_match:
        lines[index] = "Geom File=g" + next_geom_int_str + "\n"

    # Write the modified lines back to the file
    with open(str_copy_to_p_full_path, "w") as file:
        file.writelines(lines)

    # ---------------------
    # Within the plan file, adjust the selected unsteady file

    # Read the file content
    with open(str_copy_to_p_full_path, "r") as file:
        lines = file.readlines()

    # Find indices and lines that start with "Flow File="
    list_tup_match = [(index, line.strip()) for index, line in enumerate(lines) if re.match(r'^Flow File=', line)]

    # Replace lines in the list
    for index, _ in list_tup_match:
        lines[index] = "Flow File=u" + next_u_int_str + "\n"

    # Write the modified lines back to the file
    with open(str_copy_to_p_full_path, "w") as file:
        file.writelines(lines)

    # ---------------------
    # Within the plan file, change the plan 

    # Read the file content
    with open(str_copy_to_p_full_path, "r") as file:
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
    with open(str_copy_to_p_full_path, "w") as file:
        file.writelines(lines)
        
    return(str_copy_to_p_full_path)
# ~~~~~~~~~~~~~~~~~~~~~~~~


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


# --------------
def fn_set_active_plan(list_proj_title_files, next_p_int_str):
    # Set the newly created plan to current in the project file.

    # Read the file content
    with open(list_proj_title_files[0], 'r') as file:
        lines = file.readlines()

    # Find indices and lines that start with "Geom File="
    list_tup_match = [(index, line.strip()) for index, line in enumerate(lines) if re.match(r'^Current Plan=', line)]

    # Replace lines in the list
    for index, _ in list_tup_match:
        lines[index] = "Current Plan=p" + next_p_int_str + "\n"  # Replace Simulation Date

    # Write the modified lines back to the file
    with open(list_proj_title_files[0], "w") as file:
        file.writelines(lines)
# --------------


# .........................................................
def fn_copy_hecras_2d_files(str_config_file_path,
                            str_model_hydrofabric_gpkg,
                            str_ras_path,
                            int_geom_to_copy, int_u_to_copy, int_p_to_copy,
                            flt_lower_flow, flt_upper_flow, flt_delta_q):
    
    # supress all warnings
    #warnings.filterwarnings("ignore", category=UserWarning )
    
    # ********************
    # ********************
    # TODO - 2024.05.11 - Mainstem Index needs to be provided
    # Temporary value for testing
    int_mainstem_idx = 41
    int_mainstem_idx = 26
    # ********************
    # ********************
    
    print(" ")
    print("+=================================================================+")
    print('|         COPY HEC-RAS 2D FILES TO CREATE "FIREHOSE" RUN          |')
    print("|                Created by Andy Carter, PE of                    |")
    print("|             Center for Water and the Environment                |")
    print("|                 University of Texas at Austin                   |")
    print("+-----------------------------------------------------------------+")

    
    print("  ---(c) INPUT GLOBAL CONFIGUTATION FILE: " + str_config_file_path)
    print("  ---(g) INPUT MODEL HYROFABRIC (GPKG): " + str_model_hydrofabric_gpkg)
    print("  ---(m) INPUT HEC-RAS FOLDER WITH FILES TO COPY: " + str_ras_path)
    print(" ")
    print("  ---(f) HEC-RAS FILES TO COPY: ")
    print("      -- GEOMETRY FILE:      " + str(int_geom_to_copy))
    print("      -- UNSTEADY FLOW FILE: " + str(int_u_to_copy))
    print("      -- PLAN FILE:          " + str(int_p_to_copy))
    print(" ")
    print("  ---(q) FLOW TO SIMULATE: ")
    print("      -- LOW FLOW:      " + str(flt_lower_flow))
    print("      -- HIGH FLOW:     " + str(flt_upper_flow))
    print("      -- FLOW STEP:     " + str(flt_delta_q))
    print("===================================================================")
    
    # --- Read variables from config.ini ---
    # Create a ConfigParser object
    config = configparser.ConfigParser()
    
    # Read the config.ini file
    config.read(str_config_file_path)
    
    # Get the variables in the [02_copy_geom] section
    if '02_copy_geom' in config:
        section = config['02_copy_geom']
        
        # Read variables in the section
        str_boundary_to_edit = section.get('str_boundary_to_edit', '')
        int_time_rolling_avg = section.getint('int_time_rolling_avg', 0)
        int_buffer_time = section.getint('int_buffer_time', 0)
        str_simulation_start_date = section.get('str_simulation_start_date', '')
    else:
        print("[02_copy_geom] section not found in the config file.")
    # --- ---
    
    # --- Read the geopackage ---
    gdf_area = gpd.read_file(str_model_hydrofabric_gpkg, layer='00_area_2d')
    gdf_mainstems = gpd.read_file(str_model_hydrofabric_gpkg, layer='03_flowpaths_stabilize')
    # --- ---
    
    # ---- example input: gdf_mainstems.iloc[17] ----
    dict_mainstem = gdf_mainstems.iloc[int_mainstem_idx].to_dict()
    int_idx_start_cell = int(dict_mainstem['idx_start_cell'])
    list_unique_indices_sorted = [int_idx_start_cell]
    
    dict_area_2d = gdf_area.iloc[0].to_dict()
    str_2d_area_name = dict_area_2d['area_2d_name']
    
    # Set the run time for stepped flows (this is the time for each constant flow)
    int_hour_count = int(dict_mainstem['travel_time_hr']) + 1 + int_time_rolling_avg + int_buffer_time
    
    # Add the constant flow to the dict_mainstem dictionary
    # this is so a run name can be created
    dict_mainstem['time_sim_hr'] = int_hour_count
    
    # create the run names (note: set flt_upper_flow to None for single flow simulation)
    str_run_name,str_short_name = fn_build_plan_names(flt_lower_flow, flt_upper_flow, flt_delta_q, dict_mainstem)
    
    # Convert the flow data to dictionary
    dict_flows = {
        "lower_flow": flt_lower_flow,
        "upper_flow": flt_upper_flow,
        "delta_q": flt_delta_q,
        "hour_count": int_hour_count,
        "run_name": str_run_name,
        "short_name": str_short_name
    }
    
    # Get a list of the HEC-RAS projects in a given directory
    list_proj_title_files = fn_list_of_ras_projects(str_ras_path)
    
    # Copy the geometry file
    str_geom_path, hdf_file_path = fn_copy_geom_file(list_proj_title_files, int_geom_to_copy, dict_flows)
    
    # TODO - 20240801 - need to copy gXX.hdf, cXX, xXX -- (XX is the source geometry model)
    
    # Get the index in the geometry file where the internal boundary condition needs an edit
    index = fn_get_index_of_bc_to_edit(str_boundary_to_edit, str_geom_path)
    
    # Create the text to insert as boundary condition 
    list_new_boundary_lines = fn_build_internal_boundary_text(hdf_file_path, list_unique_indices_sorted, str_2d_area_name)
    
    # replace boundary condition lines in file
    fn_replace_boundary_lines_in_file(str_geom_path, list_new_boundary_lines, index)
    
    # copy the unsteady flow file
    str_copy_to_u_full_path, int_sim_time_hr = fn_copy_unsteady_flow_file(list_proj_title_files, int_u_to_copy, dict_flows, str_boundary_to_edit)
    
    # --- inputs for plan file creation ---
    next_geom_int_str = os.path.splitext(str_geom_path)[1][2:]
    next_u_int_str = os.path.splitext(str_copy_to_u_full_path)[1][2:] 
    # --- ---
    
    # copy the plan flow file and modify
    str_copy_to_p_full_path = fn_copy_plan_file(list_proj_title_files, int_p_to_copy,
                                                next_geom_int_str, next_u_int_str,
                                                dict_flows)
    
    # Establish the flow simulation time range in the plan file
    fn_set_plan_simulation_date_range(str_copy_to_p_full_path,
                                      str_simulation_start_date,
                                      int_sim_time_hr)
    
    # --- input to set the active plan ---
    next_p_int_str = os.path.splitext(str_copy_to_p_full_path)[1][2:]
    # --- ---
    
    fn_set_active_plan(list_proj_title_files, next_p_int_str)
# .........................................................



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()
    
    parser = argparse.ArgumentParser(description='========= COPY HEC-RAS 2D FILES TO CREATE "FIREHOSE" RUN =========')
    
    parser.add_argument('-c',
                        dest = "str_config_file_path",
                        help=r'REQUIRED: Global configuration filepath Example: C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini',
                        required=False,
                        default=r'C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini',
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-g',
                        dest = "str_model_hydrofabric_gpkg",
                        help=r'REQUIRED: model specific hydrofabric (from step 01) Example: E:\sample_2d_output\BLE_LBSG_501_p02\model_hydrofabric.gpkg',
                        required=False,
                        default=r'E:\sample_2d_output\BLE_LBSG_501_p02\model_hydrofabric.gpkg',
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-m',
                        dest = "str_ras_path",
                        help=r'REQUIRED: folder containing HEC-RAS model to copy Example: E:\HECRAS_2D_12070205\base_model_20240414_copy',
                        required=False,
                        default=r'E:\HECRAS_2D_12070205\base_model_20240414_copy',
                        metavar='DIR',
                        type=str)
    
    parser.add_argument('-f',
                        dest = "tup_hecras_files_to_copy",
                        help=r'REQUIRED: three integers of files to copy (geom,unsteady,plan) Example: 1,1,1 :Three integers separated by commas',
                        required=False,
                        default=(1,1,1),
                        metavar='TUPLE',
                        type=tuple_of_ints)

    parser.add_argument('-q',
                        dest = "tup_flows",
                        help=r'REQUIRED: three flows (low, high, delta_q) Example: 500,2200,400 :Three integers separated by commas. High and delta_q can be None to simulate single flow',
                        required=False,
                        default=(500,None,None),
                        metavar='TUPLE',
                        type=tuple_of_ints)
    

    args = vars(parser.parse_args())
    
    str_config_file_path = args['str_config_file_path']
    str_model_hydrofabric_gpkg = args['str_model_hydrofabric_gpkg']
    str_ras_path = args['str_ras_path']
    
    int_geom_to_copy, int_u_to_copy, int_p_to_copy = args['tup_hecras_files_to_copy']
    flt_lower_flow, flt_upper_flow, flt_delta_q = args['tup_flows']
    

    fn_copy_hecras_2d_files(str_config_file_path,
                            str_model_hydrofabric_gpkg,
                            str_ras_path,
                            int_geom_to_copy, int_u_to_copy, int_p_to_copy,
                            flt_lower_flow, flt_upper_flow, flt_delta_q)
                            
    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = timedelta(seconds=flt_time_pass)
    
    print('Compute Time: ' + str(time_pass))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~