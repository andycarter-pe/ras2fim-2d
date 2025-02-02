# Script 04 - From the HEC-RAS output, computes time gradient and max WSEL
# for each step in the "firehose"
#
# This requires a HEC-RAS 2D model in an HDF plan containing at least one(1)
# 2D Flow Area.  This routine pulls the first found 2D Area.  It creates a 
# geopackage containing:
#  (00_area_2d) - polygon of the 2D area in NextGen projection
#  (01_stream_lines) - stream flow lines of the NextGen stream networks
#  (02_flow_points) - nexus points at the for each stream moved to the upstream
#   side of each stream flow line
#  (03_flowpaths_stabilize) - mainstem run on which a stabilizing low flow run will be
#   made to estimate stepped flow simulation times
#
# Created by: Andy Carter, PE
# Created - 2024.06.11

# ************************************************************
import h5py

from shapely.geometry import Point, LineString, MultiLineString, Polygon, MultiPolygon
from shapely.ops import linemerge
from shapely.ops import unary_union
from shapely.geometry import shape

import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.ops import unary_union
from shapely.geometry import shape
import re
import os
import sqlite3
import networkx as nx
import os
import configparser

import threading
import sys

import argparse
import time
import datetime
import warnings
# ************************************************************

# --- GLOBAL WARNING SUPRESSION ---
# Suppress PerformanceWarnings from GeoPandas
from pandas.errors import PerformanceWarning
warnings.filterwarnings("ignore", category=PerformanceWarning)

# Other warnings you may want to suppress
warnings.filterwarnings("ignore", category=FutureWarning, module="geopandas")
warnings.filterwarnings("ignore", category=RuntimeWarning)
# --- GLOBAL WARNING SUPRESSION ---

# ~~~~~~~~~~~~~~~~~~~~~~
# Define a function for the spinning cursor animation
def fn_spinning_cursor():
    while True:
        for cursor in '|/-\\':
            yield cursor
            
# Define a function to display the spinning cursor
def spinning_cursor(b_stop_spinner):
    spinner = fn_spinning_cursor()
    while not b_stop_spinner[0]:
        sys.stdout.write(next(spinner))
        sys.stdout.flush()
        time.sleep(0.1)
        sys.stdout.write('\b')
# ~~~~~~~~~~~~~~~~~~~~~~


# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist" % arg)
    else:
        # File exists so return the directory
        return arg
        return open(arg, 'r')  # return an open file handle
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


#--------------
def fn_get_plan_title_from_hdf(hdf_file_path):
    # Open the HDF5 file
    with h5py.File(hdf_file_path, 'r') as hdf_file:
        # Navigate to the 'Unsteady' group within the 'Results' group
        unsteady_group = hdf_file['/Results/Unsteady/']

        # Get the value of the 'PlanTitle' attribute
        plan_title = unsteady_group.attrs['Plan Title']

        str_plan_title = plan_title.decode('utf-8')[:-5]
        
    return(str_plan_title)
#--------------


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


# ~~~~~~~~~~~~~~~~~~~~~~~
def fn_parse_run_name(str_run_name):
    #str_run_name = '1884413_wb-2410249_wb-2410261_29-hr_14100-cfs_to_45000-cfs_step_2350-cfs'
    #str_run_name = '1884413_wb-2410249_wb-2410261_29-hr_14100-cfs'

    # Split the input string on '_' or '-'
    list_run_name = re.split(r'[_-]', str_run_name)

    # extract values from parsed string
    str_mainstem = list_run_name[0]
    str_start_stream = list_run_name[1] + "-" + list_run_name[2]
    str_end_stream = list_run_name[3] + "-" + list_run_name[4]
    int_hour_count = int(list_run_name[5])
    int_start_flow = int(list_run_name[7])

    if len(list_run_name) > 10:
        int_end_flow = int(list_run_name[10])
        int_delta_q = int(list_run_name[13])

        # determine the number of 'runs' that are int_hour_count long
        list_flow_steps_int = fn_generate_step_integer(int_start_flow, int_end_flow, int_delta_q)
        
        # dropping the first flow value from list of flows
        # the first profile has no precipitation - 2024.11.22
        list_flow_steps_int = list_flow_steps_int[1:]
    else:
        list_flow_steps_int = [int_start_flow]

    int_number_of_runs = len(list_flow_steps_int)

    # Create the dictionary
    dict_run_info = {
        'run_name': str_run_name,
        'mainstem': str_mainstem,
        'start_stream': str_start_stream,
        'end_stream': str_end_stream,
        'hour_count': int_hour_count,
        'flow_steps': list_flow_steps_int,
        'number_of_runs': int_number_of_runs
    }
    
    return dict_run_info
# ~~~~~~~~~~~~~~~~~~~~~~~


# ----------------
def fn_find_path_between_nodes(graph, start_node, end_node):
    # For Step #1
    visited = set()
    stack = [(start_node, [start_node])]

    while stack:
        current_node, path = stack.pop()
        if current_node == end_node:
            return path
        if current_node not in visited:
            visited.add(current_node)
            neighbors = list(graph.successors(current_node))
            for neighbor in neighbors:
                stack.append((neighbor, path + [neighbor]))
    return None
# ----------------


# ------------------
def fn_create_travel_path_shp(str_nextgen_gpkg, str_start_node, str_end_node, b_print_output):
    # For Step #1
    
    if b_print_output:
        print('Reading NextGEN hydrofabric...(~5 sec)...')    

    # read in the hydrofabric's nexus
    gdf_nexus = gpd.read_file(str_nextgen_gpkg, layer='nexus')

    # Read the nextgen hydrofabric
    gdf_flowpaths = gpd.read_file(str_nextgen_gpkg, layer='flowpaths')
    
    # Combine points and lines into a single GeoDataFrame - preperation for graph creation
    gdf_flow_network = gpd.GeoDataFrame(pd.concat([gdf_nexus,gdf_flowpaths], ignore_index=True),
                                        crs=gdf_flowpaths.crs)

    # Create a directed graph
    G = nx.DiGraph()

    # Add edges to the graph based on 'id' and 'toid' fields
    edges = gdf_flow_network[['id', 'toid']].values.tolist()
    G.add_edges_from(edges)

    try:
        # get a list of items between begining and end
        list_travel_path = fn_find_path_between_nodes(G, str_start_node, str_end_node)

        # get just the flowpaths from gdf_flowpath_with_attrib in list_travel_path
        gdf_streams_travel_path = gdf_flowpaths[gdf_flowpaths['id'].isin(list_travel_path)]

        # Assuming gdf_streams_travel_path is your GeoDataFrame
        gdf_streams_travel_path_sorted = gdf_streams_travel_path.sort_values(by='tot_drainage_areasqkm')

        # Assuming gdf_streams_travel_path_sorted is your GeoDataFrame
        geometries = gdf_streams_travel_path_sorted.geometry

        # Extract LineStrings from MultiLineStrings
        linestrings = []
        for geom in geometries:
            if isinstance(geom, LineString):
                linestrings.append(geom)
            elif isinstance(geom, MultiLineString):
                linestrings.extend(list(geom.geoms))

        # Merge the LineStrings into a single LineString
        shp_stream_travelpath = linemerge(linestrings)

        if b_print_output:
            print('Stream path found: ' + str(len(gdf_streams_travel_path)) + ' segments')
        
        return shp_stream_travelpath
    except:
        if b_print_output:
            print(f'Stream path between {str_start_node} and {str_end_node} not found')
        return(None)
# ------------------


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


# ------------
def fn_get_gradient_array_wet_cells(hdf_file_path,int_len_gradient):

    # Specify the HDF5 file path and group path
    group_path = '/Geometry/2D Flow Areas/'

    # Get names of HDF5 Group objects in the specified group
    list_group_names = fn_get_group_names(hdf_file_path, group_path)

    # determine all the cells that are 'wet'
    str_hdf_folder_results = '/Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/'
    str_hdf_folder_results = str_hdf_folder_results + list_group_names[0] + '/'

    str_cell_min_elev_path = group_path + list_group_names[0] + '/' + 'Cells Minimum Elevation'

    # Open the HDF file
    with h5py.File(hdf_file_path, 'r') as hdf_file:
        wsel_dataset_path = str_hdf_folder_results + 'Water Surface'
        wsel_data_per_cell = hdf_file[wsel_dataset_path][:]

        arr_min_elev_per_cell = hdf_file[str_cell_min_elev_path][:]

    # Transpose the array
    arr_wsel_data_per_cell_t = wsel_data_per_cell.T

    # Subtract arr_min_elev_per_cell from each row of arr_wsel_data_per_cell_t
    arr_depth_per_cell = arr_wsel_data_per_cell_t - arr_min_elev_per_cell[:, np.newaxis]

    # Identify nan values and replace them with zeros
    arr_depth_per_cell[np.isnan(arr_depth_per_cell)] = 0

    list_indices_per_column = []

    # Iterate over each column index
    for column_index in range(arr_depth_per_cell.shape[1]):
        # Get values in the current column
        column_values = arr_depth_per_cell[:, column_index]

        # Find indices where values are greater than 0
        indices = np.where(column_values > 0)[0]

        list_indices_per_column.append(indices)

    # Flatten the list of lists
    flattened_indices = [index for sublist in list_indices_per_column for index in sublist]

    # Get unique indices
    unique_indices = np.unique(flattened_indices)
    list_unique_indices_sorted = np.sort(unique_indices)

    # Filter the array to those cells that are wet
    arr_depth_wet_cells = arr_depth_per_cell[list_unique_indices_sorted]

    # Round all values to two decimal points
    arr_depth_wet_cells = np.round(arr_depth_wet_cells, 2)

    # Filter the rows
    arr_wsel_wet_cells = arr_wsel_data_per_cell_t[list_unique_indices_sorted]
    
    # Round all values to two decimal points
    arr_wsel_wet_cells = np.round(arr_wsel_wet_cells, 2)
    
    # Create a mask where arr_depth_wet_cells is zero
    mask_zero = (arr_depth_wet_cells == 0)

    # Apply this mask to arr_wsel_wet_cells to set corresponding elements to NaN
    arr_wsel_wet_cells[mask_zero] = np.nan
    # ------------
    
    return list_unique_indices_sorted,arr_depth_wet_cells,arr_wsel_wet_cells
# ------------


# *************************
def fn_create_run_stable_hour_and_max_wsel_v3(int_hour_count,
                                              flt_max_allowed_gradient,
                                              arr_wsel_wet_cells,
                                              int_len_gradient,
                                              list_unique_indices_sorted):
    
    
    # Suppress specific warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    warnings.filterwarnings("ignore", category=FutureWarning, module="geopandas")
    # Suppress NumPy warnings
    np.seterr(invalid='ignore')
    
    # For Step #3
    
    # Revised 2024.05.28 - Looking at rolling average of each window seperately
    
    # TODO - Logic Error - 2024.05.28 - For the first pass... the stable hour should
    # be the first hour that is less than the flt_max_allowed_gradient
    # --right of the max value is wrong as there is no surge
    
    # Determine the number of hours (columns) in this array
    num_columns = arr_wsel_wet_cells.shape[1]

    # Integer divide to get the number of runs
    int_number_of_runs = num_columns // int_hour_count

    if int_number_of_runs > 0:
        list_arr_hours_to_stability = []
        list_arr_max_wsel = []

        for i in range(int_number_of_runs):
            int_start = i * int_hour_count
            int_end = int_start + int_hour_count - 1

            #print(f'{int_start}:{int_end}')
            
            # determine if first run
            b_is_first_run = False
            
            if int_start == 0:
                b_is_first_run = True
            
            # ----- Max WSEL -----
            # Limit the range of wsel array for each constant firehose run
            arr_wsel_wet_cells_run = arr_wsel_wet_cells[:, int_start:int_end]
            
            # Check if arr_wsel_wet_cells_run is entirely NaN for all runs
            if np.isnan(arr_wsel_wet_cells_run).all():
                arr_max_wsel_run = np.full(arr_wsel_wet_cells_run.shape[0], np.nan)  # Set max values to NaN for all rows
            else:
                # Get the max WSEL for each run, ignoring NaNs
                arr_max_wsel_run = np.nanmax(arr_wsel_wet_cells_run, axis=1)
            
            # Ensure rows that were entirely NaN are handled
            all_nan_mask = np.isnan(arr_wsel_wet_cells_run).all(axis=1)
            arr_max_wsel_run[all_nan_mask] = np.nan
            
            # Append the result to the list
            list_arr_max_wsel.append(arr_max_wsel_run)
            # -----  -----
            
            '''
            # ----- Max WSEL -----
            # Limit the range of wsel array for each constant firehose run
            arr_wsel_wet_cells_run = arr_wsel_wet_cells[:, int_start:int_end]

            # To supress warnings of "RuntimeWarning: All-NaN slice encountered arr_max_wsel_run = np.nanmax(arr_wsel_wet_cells_run, axis=1)"
            np.seterr(invalid='ignore')
            
            # Get the max WSEL for each run, ignoring NaNs
            arr_max_wsel_run = np.nanmax(arr_wsel_wet_cells_run, axis=1)

            # Check for rows that were entirely NaN and set their max to NaN
            all_nan_mask = np.isnan(arr_wsel_wet_cells_run).all(axis=1)
            arr_max_wsel_run[all_nan_mask] = np.nan

            # Append the result to the list
            list_arr_max_wsel.append(arr_max_wsel_run)
            # -----  -----
            '''

            # Calculate the wsel gradient for each run
            arr_wsel_gradient_short = np.gradient(arr_wsel_wet_cells[:,int_start:int_end], axis=1)

            # gradients are set to absolute value
            arr_wsel_gradient_short_abs = np.absolute(arr_wsel_gradient_short)

            # Calculate the rolling average
            kernel = np.ones(int_len_gradient) / int_len_gradient
            arr_wsel_gradient_rolling_avg = np.apply_along_axis(lambda m: np.convolve(m, kernel, mode='valid'),
                                                                axis=1, arr=arr_wsel_gradient_short_abs)

            # padding the arr_wsel_gradient_rolling_avg back to original dimensions
            # this adds NaN back to each row infront of the rolling avg array
            arr_wsel_gradient_rolling_avg = np.pad(arr_wsel_gradient_rolling_avg, ((0, 0), (int_len_gradient - 1, 0)),
                                                     mode='constant', constant_values=np.nan)

            # Initialize arrays for results
            arr_max_positions = np.empty(arr_wsel_gradient_rolling_avg.shape[0], dtype=int)
            arr_max_values = np.empty(arr_wsel_gradient_rolling_avg.shape[0])
            revised_positions = np.empty(arr_wsel_gradient_rolling_avg.shape[0], dtype=int)
            
            # Iterate through each row
            for j in range(arr_wsel_gradient_rolling_avg.shape[0]):
                row = arr_wsel_gradient_rolling_avg[j]
                if np.isnan(row).all():
                    arr_max_positions[j] = -1  # Value indicating all NaNs
                    arr_max_values[j] = np.nan  # Value indicating all NaNs
                    revised_positions[j] = -1
                else:
                    if b_is_first_run:
                        # Find the position of the first value that is not nan
                        arr_max_positions[j] = np.argmax(~np.isnan(row))
                        
                        # set the value to that in position of arr_max_positions[j]
                        arr_max_values[j] = row[arr_max_positions[j]]
                    else:
                        # Find the position of the max gradient (pulse)
                        arr_max_positions[j] = np.nanargmax(row)
                        arr_max_values[j] = np.nanmax(row)

                    if arr_max_values[j] < flt_max_allowed_gradient:
                        revised_positions[j] = arr_max_positions[j]
                    else:
                        found = False
                        for k in range(arr_max_positions[j] + 1, len(row)):
                            if row[k] < flt_max_allowed_gradient:
                                revised_positions[j] = k
                                found = True
                                break
                        if not found:
                            revised_positions[j] = -1

            list_arr_hours_to_stability.append(revised_positions)

        # Combine lists into an array
        arr_hours_to_stability = np.array(list_arr_hours_to_stability)

        # Transpose array
        arr_hours_to_stability = arr_hours_to_stability.T

        # Subtract int_len_gradient from each value
        # TODO - Removed 2024.06.11 - MAC
        #arr_hours_to_stability -= int_len_gradient

        # Replace all values less than or equal to 0 with -1
        arr_hours_to_stability[arr_hours_to_stability < 0] = -1
        # +++++++++++++

        # Combine lists into an array
        arr_max_wsel = np.array(list_arr_max_wsel)

        # Transpose array
        arr_max_wsel = arr_max_wsel.T


    # Note: arr_positions_int_all_runs:  that a value of -1 means that during this run the cell was never 'stable'
    # Note: arr_max_wsel_all_runs:  a nan value means that during this 'run' the cell was never wet

    # ------ create the dataframe of the runs 'hours_to_stable' and 'max_wsel' for each run

    # Convert list_unique_indices_sorted to a DataFrame
    df = pd.DataFrame({'cell_idx': list_unique_indices_sorted})

    # Add columns for arr_hours_to_stability
    for j in range(arr_hours_to_stability.shape[1]):
        df[f'hours_to_stable_{j+1}'] = arr_hours_to_stability[:, j]

    # Add columns for arr_max_wsel
    for j in range(arr_max_wsel.shape[1]):
        df[f'wsel_max_{j+1}'] = arr_max_wsel[:, j]

    return df
# *************************


# ~~~~~~~~~~~~~~~~~~~
def fn_create_gdf_wet_cells(hdf_file_path, list_unique_indices_sorted):
    # For Step #4
    
    group_path = '/Geometry/2D Flow Areas/'
    
    # Get names of HDF5 Group objects in the specified group
    list_group_names = fn_get_group_names(hdf_file_path, group_path)
    
    str_hdf_folder_2darea = group_path + list_group_names[0] + '/'
    hdf5_file_path = hdf_file_path

    # Location of Face Point Coordinates in HDF5
    str_facepoint_coords = str_hdf_folder_2darea + 'FacePoints Coordinate'

    # Open the HDF5 file
    with h5py.File(hdf5_file_path, 'r') as hdf_file:
        # Extract X and Y coordinates
        x_coordinates = hdf_file[str_facepoint_coords][:, 0]
        y_coordinates = hdf_file[str_facepoint_coords][:, 1]

    # Create a pandas DataFrame
    df_facepoints = pd.DataFrame({'X': x_coordinates, 'Y': y_coordinates})

    # Location of Indices of face points making up the cells
    str_cells_facepoint_indexes = str_hdf_folder_2darea + 'Cells FacePoint Indexes'

    # Open the HDF5 file
    with h5py.File(hdf5_file_path, 'r') as hdf_file:
        # Extract face points coordinate data
        facepoints_data = hdf_file[str_cells_facepoint_indexes][:]

        # Extract the projection
        projection_wkt = hdf_file.attrs['Projection'].decode('utf-8')
        # TODO - 2024.11.22 **hardcoded**
        #projection_wkt = 'PROJCS["NAD_1983_StatePlane_Texas_Central_FIPS_4203_Feet",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["False_Easting",2296583.333],PARAMETER["False_Northing",9842500.0],PARAMETER["Central_Meridian",-100.333333333333],PARAMETER["Standard_Parallel_1",31.8833333333333],PARAMETER["Standard_Parallel_2",30.1166666666667],PARAMETER["Latitude_Of_Origin",29.6666666666667],UNIT["US survey foot",0.304800609601219]]'

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

    # create a new coloumn that contains the cell index
    gdf_cells['cell_idx'] = gdf_cells.index
    
    return gdf_cells
# ~~~~~~~~~~~~~~~~~~~


# ----------
def fn_find_nearest_flowpath_and_distance(centroid, gdf_mainstems):
    # For Step #5
    # Calculate the distances from the centroid to all mainstem lines
    distances = gdf_mainstems.distance(centroid)
    # Find the index of the minimum distance
    int_nearest_index = distances.idxmin()
    flt_nearest_distance = distances.min()
    
    str_nearest_id = gdf_mainstems.loc[int_nearest_index]['id']
    
    return str_nearest_id, flt_nearest_distance
# ----------

# ..........................
def fn_create_buffered_stream_polygons(str_model_hydrofabric_path,
                                       gdf_cells,
                                       shp_travel_path,
                                       flt_mainstem,
                                       gpkg_path,
                                       int_buffer_cells):
    # For Step #5
    # Note: gdf_stream_lines is assumed as CRS of EPSG:5070
    # shp_travel_path is also EPSG:5070
    
    str_stream_lines_layer = '01_stream_lines'
    
    # Load the specified layer
    # determined in step 01
    gdf_stream_lines = gpd.read_file(str_model_hydrofabric_path, layer=str_stream_lines_layer)
    
    # Filter the GeoDataFrame
    gdf_streams_filtered_to_mainstem = gdf_stream_lines.loc[gdf_stream_lines['mainstem'] == flt_mainstem]
    
    # Filter the GeoDataFrame
    gdf_streams_on_flowpath = gdf_streams_filtered_to_mainstem[gdf_streams_filtered_to_mainstem.geometry.intersects(shp_travel_path)]
    
    # Get a list of unique 'id' values along this mainstem's flowpath
    arr_unique_ids = gdf_streams_on_flowpath['id'].unique()
    
    # from NextGen hydrofabric - get all the streamms where the 'mainstem == flt_mainstem' Example: 1884894.0

    # Load the GeoPackage layers
    gdf_flowpaths = gpd.read_file(gpkg_path, layer='flowpaths')
    gdf_nexus = gpd.read_file(gpkg_path, layer='nexus')
    
    # find all the reaches that have a 'mainstem' == flt_mainstem
    gdf_selected_mainstem_streams = gdf_flowpaths[gdf_flowpaths['mainstem'] == flt_mainstem]

    # Set the crs of the mainstem to crs of cells
    gdf_selected_mainstem_streams = gdf_selected_mainstem_streams.to_crs(gdf_cells.crs)

    # find the most downstream mainstem
    gdf_largest_drainage_area = gdf_selected_mainstem_streams.loc[gdf_selected_mainstem_streams['tot_drainage_areasqkm'].idxmax()]
    str_toid_with_largest_drainage_area = gdf_largest_drainage_area['toid']

    
    gdf_selected_rows = gdf_nexus[gdf_nexus['id'] == str_toid_with_largest_drainage_area]

    if len(gdf_selected_rows) == 1:
        # only one downstream node found
        str_dwn_stream_id = gdf_selected_rows.iloc[0]['toid']

        # Finding the stream
        gdf_downstream_stream = gdf_flowpaths[gdf_flowpaths['id'] == str_dwn_stream_id]

        # Set the crs of the mainstem to crs of cells
        gdf_downstream_stream = gdf_downstream_stream.to_crs(gdf_cells.crs)

        # Concatenating DataFrames
        gdf_all_streams = pd.concat([gdf_downstream_stream, gdf_selected_mainstem_streams])
    else:
        gdf_all_streams = gdf_selected_mainstem_streams
        
    # Apply the function to each centroid in gdf_cells
    nearest_info = gdf_cells['geom_centroid'].apply(lambda x: fn_find_nearest_flowpath_and_distance(x, gdf_all_streams))

    gdf_cells['nearest_flowpath'] = nearest_info.apply(lambda x: x[0])
    gdf_cells['distance_to_nearest_flowpath'] = nearest_info.apply(lambda x: x[1])

    # List of unique streams in gdf_all_streams
    list_unique_stream_ids = gdf_all_streams['id'].unique().tolist()

    # remove the downstream stream
    # we don't want this polygons as it is on another river
    list_unique_stream_ids.remove(str_dwn_stream_id)

    # For each item in list_unique_stream_ids, create polygons of the disolved cells in gdf_cells.  
    # Disolve on 'nearest_flowpath'

    list_dissolved_polygons_gdf = []
    
    for stream_id in list_unique_stream_ids:
        # Filter cells associated with the current stream ID
        cells_for_stream = gdf_cells[gdf_cells['nearest_flowpath'] == stream_id]

        # Dissolve cells based on the 'nearest_flowpath' column
        dissolved = cells_for_stream.dissolve(by='nearest_flowpath')

        # Drop specified columns
        dissolved.drop(columns=['cell_idx', 'distance_to_nearest_flowpath'], inplace=True)

        # Append the dissolved polygon to the list
        list_dissolved_polygons_gdf.append(dissolved)
        
    # There can be streams along the mainsteam that have no geometry.  This is because the 
    # 'firehose' is downstream of these reach or the water never got to that reach.
    # Drop those reaches without geometry from list_dissolved_polygons_gdf

    # Remove reaches without geometry from list_dissolved_polygons_gdf
    list_dissolved_polygons_gdf_filtered = []

    for gdf_merged_poly in list_dissolved_polygons_gdf:
        if not gdf_merged_poly.empty:
            list_dissolved_polygons_gdf_filtered.append(gdf_merged_poly)

    list_dissolved_polygons_gdf = list_dissolved_polygons_gdf_filtered

    list_buffered_shp_polys = []
    list_mainstem = []

    for gdf_merged_poly in list_dissolved_polygons_gdf:

        str_mainstem = gdf_merged_poly.index[0]
        list_mainstem.append(str_mainstem)

        shp_polygon_to_check = gdf_merged_poly['geometry'].iloc[0]

        for i in range(int_buffer_cells):
            gdf_intersecting_cells = gdf_cells[gdf_cells.geometry.touches(shp_polygon_to_check)]
            merged_geometry = unary_union(gdf_intersecting_cells.geometry)
            merged_polygon = shp_polygon_to_check.union(merged_geometry)
            shp_polygon_to_check = merged_polygon

        list_buffered_shp_polys.append(shp_polygon_to_check)
        
    # The entire run's wet cells
    # Create a GeoDataFrame from the list of shapely objects and list of strings
    gdf_buffered_limits = gpd.GeoDataFrame({'flowpath': list_mainstem, 'geometry': list_buffered_shp_polys})

    # Set the crs of the newly created geodataframe
    gdf_buffered_limits.crs =  gdf_cells.crs
    
    # Convert arr_unique_ids to a set for faster lookup
    arr_unique_ids_set = set(arr_unique_ids)

    # Use isin to create a boolean mask
    mask = gdf_buffered_limits['flowpath'].isin(arr_unique_ids_set)

    # Apply the mask to select rows
    gdf_buffers_in_flowpath = gdf_buffered_limits[mask]

    # gdf_buffers_in_flowpath are all the polygons of the 'wet' cells that are
    # nearest to the streams within the simulated flowpath
    
    return gdf_buffers_in_flowpath, arr_unique_ids, gdf_all_streams
# ..........................


# ------------
def fn_merge_flowpath_buffers(gdf_buffers_in_flowpath):
    # For Step #7
    # merge gdf_buffers_in_flowpath into a single polygon

    # First, ensure that the 'geometry' column contains Polygon or MultiPolygon objects
    # This is usually done when reading the GeoDataFrame, but just in case
    gdf_buffers_in_flowpath['geometry'] = gdf_buffers_in_flowpath['geometry'].apply(lambda geom: geom if isinstance(geom, MultiPolygon) else MultiPolygon([geom]))

    # Then, merge all the polygons into a single MultiPolygon
    merged_geometry = gdf_buffers_in_flowpath['geometry'].unary_union

    # Create a new GeoDataFrame with the merged geometry
    gdf_merged_flooded_cells = gpd.GeoDataFrame(geometry=[merged_geometry])

    # Set the crs of the newly created geodataframe
    gdf_merged_flooded_cells.crs =  gdf_buffers_in_flowpath.crs
    
    return gdf_merged_flooded_cells
# ------------


# ~~~~~~~~~~~~~~~~~~
def fn_make_gdf_cells_with_runs(gdf_cells, 
                                gdf_merged_flooded_cells, 
                                df_stable_hr_max_wsel, 
                                int_len_gradient,
                                list_flow_steps,
                                str_run_name):

    # For Step #8
    # Clip gdf_cells using gdf_merged_flooded_cells geometry
    # This is to get only cells on the sstreams within the simulated flowpath
    gdf_cells_clipped = gpd.clip(gdf_cells, gdf_merged_flooded_cells.geometry)

    # Drop coloumns from gdf_cells_clipped
    gdf_cells_clipped_subset = gdf_cells_clipped.drop(columns=['geom_centroid', 'distance_to_nearest_flowpath'])

    # Perform left join
    gdf_cells_with_runs = pd.merge(gdf_cells_clipped_subset, df_stable_hr_max_wsel, on='cell_idx', how='left')

    # Correct the hours to stable... this is a moving average, so subtract this moving average to get to
    # correct travel time for the given flow.

    # for each coloumn in df_stable_hr_max_wsel that has a name that starts with
    # 'hours_to_stable', subract int_len_gradient.  If any of the revised values
    # in these coloumns are less than or equal to zero, set their value of -1

    # Loop through columns that start with 'hours_to_stable'
    for col in gdf_cells_with_runs.columns:
        if col.startswith('hours_to_stable'):
            # Subtract int_len_gradient from the column values
            #gdf_cells_with_runs[col] -= int_len_gradient

            # Set values less than or equal to zero to -1
            #gdf_cells_with_runs[col] = gdf_cells_with_runs[col].apply(lambda x: -1 if x <= 0 else x)
            
            # TODO - 2024.06.11 - Possible need to Turn-on subtraction - MAC
            pass

    # Create a coloumn for each flow step, set the entire contents of that cell to that flow
    # Iterate through the list and create columns
    for i, flow in enumerate(list_flow_steps, 1):
        col_name = f'flow_{i}'
        gdf_cells_with_runs[col_name] = flow

    # Add the run_name to every cell
    gdf_cells_with_runs['run_name'] = str_run_name

    # if gdf_cells_with_runs contains any geomety that is not Polygon, drop those rows
    gdf_cells_with_runs = gdf_cells_with_runs[gdf_cells_with_runs.geometry.geom_type == 'Polygon']
    
    return(gdf_cells_with_runs)
# ~~~~~~~~~~~~~~~~~~


# --------------
def fn_df_from_list_of_lists(list_of_lists_input, str_col_name):
    
    # For Step #9
    # Transpose the list of lists to make each sublist a column
    transposed_data = list(map(list, zip(*list_of_lists_input)))

    # Create column names
    num_columns = len(transposed_data[0])
    column_names = [f'{str_col_name}_{i+1}' for i in range(num_columns)]


    # Create the DataFrame
    df = pd.DataFrame(transposed_data, columns=column_names)
    
    return(df)
# --------------


# .......................
def fn_add_time_values_to_buffers_in_flowpath(gdf_cells_with_runs,
                                              gdf_buffers_in_flowpath,
                                              list_flow_steps,
                                              int_buffer_cells,
                                              str_run_name):

    # For Step #9
    
    # List comprehension to find columns that start with 'wsel_max_'
    wsel_columns = [col for col in gdf_cells_with_runs.columns if col.startswith('wsel_max_')]
    
    # For each gdf_buffers_in_flowpath

    # get the 'time to stable' for each buffered flowpath, the count of cells in the area
    # and the count of unstable cells in the area

    # Create empty lists to store results
    list_of_highest_hours = []
    list_of_lowest_hours = []
    list_of_cell_count = []
    list_of_negative_one_count = []

    # Iterate over each row in gdf_buffered_polygon

    for i in range(len(wsel_columns)):
        str_wsel_col_name = 'wsel_max_' + str(i+1)
        str_stable_hr_name = 'hours_to_stable_' + str(i+1)
        

        # reset each loop
        highest_hours = []
        lowest_hours = []
        cell_count = []
        negative_one_count = []

        for index, row in gdf_buffers_in_flowpath.iterrows():
            # Intersect the polygon geometry with gdf_cells_wsel
            intersected_cells = gdf_cells_with_runs[gdf_cells_with_runs.intersects(row['geometry'])]
            
            # Filter to all the rows in intersected_cells that have values in str_wsel_col_name that are not NaN
            # These are the 'wet' cells for the given run, for a given stream
            wet_cells_in_stream_poly = intersected_cells[intersected_cells[str_wsel_col_name].notna()]

            # If there are intersected 'wet' cells
            if not wet_cells_in_stream_poly.empty:
                # Get the 'time to stable' for each buffered flowpath
                highest_hours.append(wet_cells_in_stream_poly[str_stable_hr_name].max())
                lowest_hours.append(wet_cells_in_stream_poly[str_stable_hr_name].min())

                # Get the count of cells in the area
                cell_count.append(len(wet_cells_in_stream_poly))

                # Get the count of unstable cells in the area (where 'wsel_max' is -1)
                negative_one_count.append(len(wet_cells_in_stream_poly[wet_cells_in_stream_poly[str_stable_hr_name] == -1]))
            else:
                # If no intersection, set to None
                highest_hours.append(None)
                lowest_hours.append(None)
                cell_count.append(0)
                negative_one_count.append(0)
            
        # lists of lists (for each run)
        list_of_highest_hours.append(highest_hours)
        list_of_lowest_hours.append(lowest_hours)
        list_of_cell_count.append(cell_count)
        list_of_negative_one_count.append(negative_one_count)


    df_highest_hr = fn_df_from_list_of_lists(list_of_highest_hours, 'highest_hr_to_stable')
    df_lowest_hr = fn_df_from_list_of_lists(list_of_lowest_hours, 'lowest_hr_to_stable')
    df_cell_count = fn_df_from_list_of_lists(list_of_cell_count, 'cell_count')
    df_unstable_cell_count = fn_df_from_list_of_lists(list_of_negative_one_count, 'unstable_cell_count')

    # Concatenate the DataFrames
    df_combined = pd.concat([df_highest_hr, df_lowest_hr, df_cell_count, df_unstable_cell_count], axis=1)

    # add each run's flow
    for i, flow in enumerate(list_flow_steps, 1):
        col_name = f'flow_{i}'
        df_combined[col_name] = flow

    # Add values consistent to all runs
    df_combined['buffer_cell_count'] = int_buffer_cells
    df_combined['run_name'] = str_run_name

    list_flowpath = gdf_buffers_in_flowpath['flowpath'].tolist()
    df_combined['flowpath'] = list_flowpath

    # left join gdf_buffers_in_flowpath and df_combined on 'flowpath'
    gdf_buffers_in_flowpath_w_values = gdf_buffers_in_flowpath.merge(df_combined, on='flowpath', how='left')
    
    return(gdf_buffers_in_flowpath_w_values)
# .......................


# ~~~~~~~~~~~~~~~~~~
def fn_compute_stream_stats(gdf_all_streams,arr_unique_ids,
                            gdf_cells_with_runs,
                            list_flow_steps,
                            str_run_name):

    # For Step #10
    # Evaluating the stream statistics (% wet, % stable)

    # Select rows where 'id' column is in list_unique_stream_ids
    gdf_streams_in_simulation = gdf_all_streams[gdf_all_streams['id'].isin(arr_unique_ids)]

    gdf_streams_in_simulation_light = gdf_streams_in_simulation[['geometry','id', 'mainstem',
                                                                 'order','tot_drainage_areasqkm']].copy()

    # Calculte length of the entire stream reach
    gdf_streams_in_simulation_light['length'] = gdf_streams_in_simulation_light.geometry.length

    # List comprehension to find columns that start with 'wsel_max_'
    wsel_columns = [col for col in gdf_cells_with_runs.columns if col.startswith('wsel_max_')]

    list_df_runs = []

    # Iterate over each run
    for i in range(len(wsel_columns)):
        str_wsel_col_name = 'wsel_max_' + str(i+1)
        str_stable_hr_name = 'hours_to_stable_' + str(i+1)
        str_dist_wet = 'dist_wet_' + str(i+1)
        str_dist_stable = 'dist_stable_' + str(i+1)
        str_max_hr_stable = 'stream_cl_hr_to_stable_' + str(i+1)
        str_min_hr_stable = 'stream_cl_min_hr_to_stable_' + str(i+1)

        # Filter to all the rows in intersected_cells that have values in str_wsel_col_name that are not NaN
        # These are the 'wet' cells for the given run
        wet_cells_in_run = gdf_cells_with_runs[gdf_cells_with_runs[str_wsel_col_name].notna()]

        # Perform geometric overlay to compute intersection
        gdf_cell_intersection = gpd.overlay(gdf_streams_in_simulation_light, wet_cells_in_run, how='intersection')

        # Calculate length of each line
        gdf_cell_intersection['segment_length'] = gdf_cell_intersection.geometry.length

        # create pandas series that is max of stable hour of cells intersecting stream by unique 'id'
        ps_max_hr_to_stable_on_stream = gdf_cell_intersection.groupby('id')[str_stable_hr_name].max().rename(str_max_hr_stable)
        
        # ~~~~~~~~~
        # also getting the minimum hour
        # Filter out rows where str_stable_hr_name is -1
        gdf_filtered = gdf_cell_intersection[gdf_cell_intersection[str_stable_hr_name] != -1]
        
        # create pandas series that is min of stable hour of cells intersecting stream by unique 'id'
        ps_min_hr_to_stable_on_stream = gdf_filtered.groupby('id')[str_stable_hr_name].min().rename(str_min_hr_stable)
        # ~~~~~~~~~

        # create pandas series that is sum of segment_length by unique 'id'
        ps_sum_segment_length = gdf_cell_intersection.groupby('id')['segment_length'].sum().rename(str_dist_wet)

        # create pandas series that is sum of segment_length by unique 'id' where the stream is 'stable'
        filtered_df = gdf_cell_intersection[gdf_cell_intersection[str_stable_hr_name] > 0]
        ps_sum_segment_length_stable = filtered_df.groupby('id')['segment_length'].sum().rename(str_dist_stable)

        # Merge the Series into a DataFrame
        df_stats_this_run = pd.concat([ps_sum_segment_length, ps_sum_segment_length_stable,ps_max_hr_to_stable_on_stream, ps_min_hr_to_stable_on_stream], axis=1)

        list_df_runs.append(df_stats_this_run)

    # Combine the dataframes along the columns
    df_combined_stream_lengths = pd.concat(list_df_runs, axis=1)

    # Perform the left join
    gdf_stream_w_dist = gdf_streams_in_simulation_light.merge(df_combined_stream_lengths, on='id', how='left')

    # get the percentage of wet and stable for each run
    for i in range(len(wsel_columns)):
        str_dist_wet = 'dist_wet_' + str(i+1)
        str_dist_stable = 'dist_stable_' + str(i+1)
        str_perct_wet = 'perct_wet_' + str(i+1)
        str_perct_stable = 'perct_stable_' + str(i+1)

        # Compute percentage of wet length
        gdf_stream_w_dist[str_perct_wet] = round((gdf_stream_w_dist[str_dist_wet] / gdf_stream_w_dist['length']) * 100,1)

        # Compute percentage of stable length
        gdf_stream_w_dist[str_perct_stable] = round((gdf_stream_w_dist[str_dist_stable] / gdf_stream_w_dist['length']) * 100,1)

    # add each run's flow
    for i, flow in enumerate(list_flow_steps, 1):
        col_name = f'flow_{i}'
        gdf_stream_w_dist[col_name] = flow

    # Add the run_name to every stream
    gdf_stream_w_dist['run_name'] = str_run_name
    
    return(gdf_stream_w_dist)
# ~~~~~~~~~~~~~~~~~~


# ''''''''''''''''''''''''''''''
def fn_dict_of_ras_model(hdf_file_path):

    # Create a table of the hydraulic model attributes for hydraulic results geopackage

    # Open the HDF5 file
    with h5py.File(hdf_file_path, 'r') as hdf:
        # Navigate to the specified path
        plan_info = hdf['/Plan Data/Plan Information']

        # Get the attributes and convert values to strings
        dict_plan_attr = {attr: plan_info.attrs[attr].tobytes().decode('utf-8') for attr in plan_info.attrs}

        # Convert the dictionary to a DataFrame
        df1 = pd.DataFrame([dict_plan_attr])

    # Open the HDF file
    with h5py.File(hdf_file_path, 'r') as hdf:
        # Navigate to the specified path
        model_info = hdf['/']

        # Get the attributes and convert values to strings
        dict_model_attr = {attr: model_info.attrs[attr].tobytes().decode('utf-8') for attr in model_info.attrs}

        # Convert the dictionary to a DataFrame
        df2 = pd.DataFrame([dict_model_attr])

    # Combine the two dataframes
    df_model_info = pd.concat([df2, df1], axis=1)

    # Get the plan name -- remove '_plan' from the string
    df_model_info['run_name'] = df_model_info['Plan Name'].apply(lambda x: x[:-5])

    # Move 'run_name' to the first coloumn
    cols = df_model_info.columns.tolist()
    cols = ['run_name'] + [col for col in cols if col != 'run_name']
    df_model_info = df_model_info[cols]
    
    return (df_model_info)
# ''''''''''''''''''''''''''''''


#-------------
def fn_output_gpkg(str_output_path,
                   str_run_name,
                   gdf_stream_w_dist,
                   gdf_cells_with_runs,
                   gdf_buffers_in_flowpath_w_values,
                   df_model_info):
    
    if not os.path.exists(str_output_path):
        os.makedirs(str_output_path)

    str_gpkg_name = 'hydrualic_results_' + str_run_name + '.gpkg'
    str_output_gpkg_path = os.path.join(str_output_path, str_gpkg_name)

    # Write GeoDataFrames to GeoPackage as separate layers
    gdf_stream_w_dist.to_file(str_output_gpkg_path, layer='03_streams_ln', driver="GPKG")
    gdf_cells_with_runs.to_file(str_output_gpkg_path, layer='02_cells_wsel_ar', driver="GPKG")
    gdf_buffers_in_flowpath_w_values.to_file(str_output_gpkg_path, layer='01_flowpath_flooded_cells_ar', driver="GPKG")


    # ---- Add the HEC-RAS Table with sqlite ----
    # ---- No Spatial Data ----
    # Connect to the GeoPackage
    conn = sqlite3.connect(str_output_gpkg_path)

    # Add the DataFrame to the GeoPackage as a table
    df_model_info.to_sql(name='00_hec_ras_infor_tbl', con=conn, if_exists='replace', index=False)

    # Commit changes
    conn.commit()

    # Close connection
    conn.close()
#-------------


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


# ............................
def fn_compute_max_wsel_and_stable_hr(str_config_file_path,
                                      str_hdf_file_path,
                                      str_model_hydrofabric_path,
                                      str_output_dir,
                                      b_print_output):
    
    # supress all warnings
    warnings.filterwarnings("ignore", category=UserWarning )
    np.seterr(invalid='ignore')
    warnings.filterwarnings("ignore", category=FutureWarning, module="geopandas")
    
    if b_print_output:
        print(" ")
        print("+=================================================================+")
        print("|        DETERMINE TRAVEL TIME AND MAX WSEL FOR FLOW STEPS        |")
        print("|                Created by Andy Carter, PE of                    |")
        print("|             Center for Water and the Environment                |")
        print("|                 University of Texas at Austin                   |")
        print("+-----------------------------------------------------------------+")
    
        
        print("  ---(c) INPUT GLOBAL CONFIGUTATION FILE: " + str_config_file_path)
        print("  ---(i) INPUT HEC-RAS PLAN (HDF): " + str_hdf_file_path)
        print("  ---(g) INPUT HEC-RAS 2D AREA HYDROFABRIC (GPKG): " + str_model_hydrofabric_path)
        print("  ---(o) OUTPUT PATH: " + str_output_dir)
        print("  ---[r] PRINT OUTPUT: " + str(b_print_output))
        print("===================================================================")
    
    if not os.path.exists(str_output_dir):
        os.mkdir(str_output_dir)
        
    # --- Read variables from config.ini ---
    # Create a ConfigParser object
    config = configparser.ConfigParser()
    
    # Read the config.ini file
    config.read(str_config_file_path)
    
    # ----
    # Get the variables in the [global_input_files] section
    if 'global_input_files' in config:
        section = config['global_input_files']
        
        # Read variables in the section
        gpkg_path = section.get('str_nextgen_gpkg', '')
    else:
        print("[global_input_files] section not found in the config file.")
    # ----
    
    # ----
    # Get the variables in the [02_copy_geom] section
    if '02_copy_geom' in config:
        section = config['02_copy_geom']
        
        # Read variables in the section
        int_len_gradient = section.getint('int_time_rolling_avg', 0)
    else:
        print("[02_copy_geom] section not found in the config file.")
    # ----
    
    # ----
    # Get the variables in the [04_time_gradient] section
    if '04_time_gradient' in config:
        section = config['04_time_gradient']
        
        # Read variables in the section
        flt_max_allowed_gradient = section.getfloat('flt_max_allowed_gradient', 0)
        int_buffer_cells = section.getint('int_buffer_cells', 0)
    else:
        print("[04_time_gradient] section not found in the config file.")
    # ----
    
    # from output plan HDF, get the run name and parse values
    str_run_name = fn_get_plan_title_from_hdf(str_hdf_file_path)
    dict_run_info = fn_parse_run_name(str_run_name)
    
    # --- inputs from the run name (dict_run_info) ---
    int_hour_count = dict_run_info['hour_count']
    str_mainstem = str(dict_run_info['mainstem'])
    flt_mainstem = float(dict_run_info['mainstem'])
    list_flow_steps = dict_run_info['flow_steps']
    int_number_of_runs = dict_run_info['number_of_runs']
    
    str_start_node = dict_run_info['start_stream']
    str_end_node = dict_run_info['end_stream']
    # ---
    
    b_stop_spinner = [False]  # Use a list to allow modification from other threads
    spinner_thread = threading.Thread(target=spinning_cursor, args=(b_stop_spinner,))
    
    if b_print_output:
        spinner_thread.start()
    
    try:
        if b_print_output:
            print()
            #print('Step 1')
        shp_travel_path = fn_create_travel_path_shp(gpkg_path, str_start_node, str_end_node, b_print_output)
        
        if b_print_output:
            print()
            print('Determining wet cells...')
        list_unique_indices_sorted,arr_depth_wet_cells,arr_wsel_wet_cells = fn_get_gradient_array_wet_cells(str_hdf_file_path,
                                                                                                            int_len_gradient)
        #print()
        #print('Step 3')
        df_stable_hr_max_wsel = fn_create_run_stable_hour_and_max_wsel_v3(int_hour_count,
                                                                          flt_max_allowed_gradient,
                                                                          arr_wsel_wet_cells,
                                                                          int_len_gradient,
                                                                          list_unique_indices_sorted)
        
        if b_print_output:
            print()
            #print('Step 4')
            print('Extracting cell geometry...')
        gdf_cells = fn_create_gdf_wet_cells(str_hdf_file_path, list_unique_indices_sorted)
        

        #print('Step 5')
        # Compute the centroid for every cell's polygon
        gdf_cells['geom_centroid'] = gdf_cells.geometry.centroid
        
        
        if b_print_output:
            print()
            #print('Step 6')
            print('Grouping cells by NextGEN streams...')
        gdf_buffers_in_flowpath, arr_unique_ids, gdf_all_streams = fn_create_buffered_stream_polygons(str_model_hydrofabric_path,
                                                                                                      gdf_cells,
                                                                                                      shp_travel_path,
                                                                                                      flt_mainstem,
                                                                                                      gpkg_path,
                                                                                                      int_buffer_cells)
        
        
        if b_print_output:
            #print('Step 7')
            print()
            print('Making stream buffers (overlap)...')
        gdf_merged_flooded_cells = fn_merge_flowpath_buffers(gdf_buffers_in_flowpath)
        
        
        if b_print_output:
            #print()
            print('Step 8')
        gdf_cells_with_runs = fn_make_gdf_cells_with_runs(gdf_cells, 
                                                          gdf_merged_flooded_cells, 
                                                          df_stable_hr_max_wsel, 
                                                          int_len_gradient,
                                                          list_flow_steps,
                                                          str_run_name)
        
        if b_print_output:
            #print()
            print('Step 9')
        gdf_buffers_in_flowpath_w_values = fn_add_time_values_to_buffers_in_flowpath(gdf_cells_with_runs,
                                                                                     gdf_buffers_in_flowpath,
                                                                                     list_flow_steps,
                                                                                     int_buffer_cells,
                                                                                     str_run_name)
        
        if b_print_output:
            print()
            print('Evaluating stream statistics...')
        gdf_stream_w_dist = fn_compute_stream_stats(gdf_all_streams,
                                                    arr_unique_ids,
                                                    gdf_cells_with_runs,
                                                    list_flow_steps,
                                                    str_run_name)
        
        
        #print()
        #print('Step 11')
        df_model_info = fn_dict_of_ras_model(str_hdf_file_path)
        
        if b_print_output:
            print()
            print('Writing output geopackage...')
        fn_output_gpkg(str_output_dir,str_run_name,
                       gdf_stream_w_dist,gdf_cells_with_runs,gdf_buffers_in_flowpath_w_values,
                       df_model_info)
    finally:
        if b_print_output:
            b_stop_spinner[0] = True
            spinner_thread.join()
            
    if b_print_output:
        print()
# ............................


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()
    
    parser = argparse.ArgumentParser(description='======== DETERMINE TRAVEL TIME AND MAX WSEL FOR FLOW STEPS ========')
    
    parser.add_argument('-c',
                        dest = "str_config_file_path",
                        help=r'REQUIRED: Global configuration filepath Example: C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini',
                        required=True,
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-i',
                        dest = "str_hdf_file_path",
                        help=r'REQUIRED: path to plan HDF file (HDF) Example: E:\ras2fim_test_20250107\03_run_hdf\1884650_wb-2410416_wb-2410422_21-hr_7518-cfs_to_16893-cfs_step_500-cfs.p01.hdf',
                        required=True,
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-g',
                        dest = "str_model_hydrofabric_path",
                        help=r'REQUIRED: models 2D area hydrofabric (geopackage) Example: E:\ras2fim_test_20250107\01_model_hydrofabric\model_hydrofabric.gpkg',
                        required=True,
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-o',
                        dest = "str_output_dir",
                        help=r'REQUIRED: directory to write model hydrofabric geopackage Example: E:\ras2fim_test_20250107',
                        required=True,
                        metavar='DIR',
                        type=str)
    
    parser.add_argument('-r',
                        dest = "b_print_output",
                        help=r'OPTIONAL: Print output messages Default: True',
                        required=False,
                        default=True,
                        metavar='T/F',
                        type=fn_str_to_bool)
    

    args = vars(parser.parse_args())
    
    str_config_file_path = args['str_config_file_path']
    str_hdf_file_path = args['str_hdf_file_path']
    str_model_hydrofabric_path = args['str_model_hydrofabric_path']
    str_output_dir = args['str_output_dir']
    b_print_output = args['b_print_output']

    fn_compute_max_wsel_and_stable_hr(str_config_file_path,
                                      str_hdf_file_path,
                                      str_model_hydrofabric_path,
                                      str_output_dir,
                                      b_print_output)

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)
    
    print('Compute Time: ' + str(time_pass))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~