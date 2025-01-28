# Script 01 - Create a geopackage of the NextGen stream network from the HEC-RAS Model
#
# This requires as HEC-RAS 2D model in an HDF plan containing at least one(1)
# 2D Flow Area.  This routine pulls the first found 2D Area.  It creates a 
# geopackage containing:
#  (00_area_2d) - polygon of the 2D area in NextGen projection
#  (01_stream_lines) - stream flow lines of the NextGen stream networks
#  (02_flow_points) - nexus points at the for each stream moved to the upstream
#   side of each stream flow line
#  (03_flowpaths) - run on which a flow run will be simulated
#
# Created by: Andy Carter, PE
# Created - 2024.11.06
# Revised - 2025.01.04
# ************************************************************


# ************************************************************
import h5py

from shapely.geometry import Point, LineString, Polygon, MultiPolygon, box, mapping, shape
import geopandas as gpd
import pandas as pd
import networkx as nx
import numpy as np
import os
import configparser

from shapely.geometry import MultiLineString, LineString, Point
from shapely.ops import linemerge

import argparse
import time
import datetime
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
def fn_get_gdf_of_2d_area(hdf_file_path):

    # Specify the HDF5 file path and group path
    str_hdf_geom_path = '/Geometry/2D Flow Areas/'

    # Get names of HDF5 Group objects in the specified group
    list_group_names = fn_get_group_names(hdf_file_path, str_hdf_geom_path)

    b_has_2d_area = False

    if len(list_group_names) > 1:
        print('   -- Multiple 2D areas found -- Using first area:', list_group_names[0])
        b_has_2d_area = True
    elif len(list_group_names) == 0:
        print('Error: No 2D areas found')
    else:
        # Only one 2D area found
        b_has_2d_area = True

    if b_has_2d_area:
        str_perimeter_points = str_hdf_geom_path + list_group_names[0] + '/' + 'Perimeter'

        # Open the HDF file
        with h5py.File(hdf_file_path, 'r') as hdf_file:
            arr_perim_points = hdf_file[str_perimeter_points][:]

            # Extract the projection
            projection_wkt = hdf_file.attrs['Projection'].decode('utf-8')

            str_terrain_filename = hdf_file['/Geometry/'].attrs['Terrain Filename'].decode('utf-8')

        # Convert the array of perimeter points into a Polygon
        shp_polygon_geom = Polygon(arr_perim_points)

        # Create a GeoDataFrame
        gdf_2d_area_polygon = gpd.GeoDataFrame(index=[0], crs=projection_wkt, geometry=[shp_polygon_geom])
        
        gdf_2d_area_polygon['area_2d_name'] = list_group_names[0]
        gdf_2d_area_polygon['hdf_path'] = hdf_file_path
        gdf_2d_area_polygon['terrain_path'] = str_terrain_filename
        gdf_2d_area_polygon['prj_wkt_ras_model'] = projection_wkt
        
        return(gdf_2d_area_polygon)
    else:
        pass
        # return nothing as there is nothing to return
# ------------------------


# ------------------------
def fn_extract_flowpaths_gdf(gdf_2d_area_polygon, gpkg_nextgen_path):
    
    # this is throwing a ShapelyDeprecationWarning: Iteration over multi-part 
    # geometries is deprecated and will be removed in Shapely 2.0. Use the 
    # `geoms` property to access the constituent parts of a multi-part geometry.

    # From the 2d area, get the nextgen hydrofabric
    print('   -- Reading Nextgen Hydrofabric ~5 seconds...')

    # Read the nextgen hydrofabric
    gdf_flowpaths = gpd.read_file(gpkg_nextgen_path, layer='flowpaths')

    gdf_2d_area_polygon_nextgen = gdf_2d_area_polygon.to_crs(gdf_flowpaths.crs)

    # Extract the first polygon
    shp_first_poly = gdf_2d_area_polygon_nextgen.geometry.iloc[0]

    # Get the bounding box coordinates
    bbox = shp_first_poly.bounds

    # Create a bounding box geometry using the CRS of the GeoPackage
    bbox_geom = box(*bbox)
    
    # Filter lines within the bounding box
    gdf_ln_within_bbox = gdf_flowpaths[gdf_flowpaths.geometry.intersects(bbox_geom)]

    # Get just the stream lines that are within or intersect the 2d area
    gdf_ln_within_2d_area = gdf_ln_within_bbox[gdf_ln_within_bbox.geometry.within(shp_first_poly) | gdf_ln_within_bbox.geometry.intersects(shp_first_poly)]

    # Get a unique list of 'id' from gdf_ln_within_2d_area
    list_unique_ids = gdf_ln_within_2d_area['id'].unique().tolist()

    # Read the GeoPackage file to get the 'flowpath_attributes' table
    gdf_flowpath_attrib = gpd.read_file(gpkg_nextgen_path, layer='flowpath_attributes')

    # Select rows where 'id' is in list_unique_ids
    gdf_attrib_flowpath_2darea = gdf_flowpath_attrib[gdf_flowpath_attrib['id'].isin(list_unique_ids)]

    # Specify the columns to keep
    columns_to_keep = ['id', 'rl_gages', 'rl_NHDWaterbodyComID', 'So', 'ChSlp']

    # Drop all columns except the specified ones
    gdf_attrib_flowpath_2darea_subset = gdf_attrib_flowpath_2darea[columns_to_keep]

    # Perform left join on 'id' to add the selecte attributes
    gdf_flowpath_with_attrib = gdf_ln_within_2d_area.merge(gdf_attrib_flowpath_2darea_subset, on='id', how='left')

    # Check if each line in gdf_flowpath_with_attrib is entirely within the polygon
    gdf_flowpath_with_attrib['within_2darea'] = gdf_flowpath_with_attrib.geometry.within(gdf_2d_area_polygon_nextgen.iloc[0].geometry)
    
    return(gdf_flowpath_with_attrib)
# ------------------------


# ------------------
def fn_create_upstream_flowpath_points(gdf_flowpath_with_attrib, gdf_2d_area_polygon_nextgen):

    # get a point that is at the upstream end of every flowpath

    # Create an empty list to store the points
    list_upstream_points = []

    # Iterate over each geometry in the GeoDataFrame
    list_id = []

    for index,row in gdf_flowpath_with_attrib.iterrows():

        geom_items = row['geometry']
        list_id.append(row['id'])

        for geom in geom_items:
            # Check if the geometry is a LineString or MultiLineString
            if geom.geom_type == 'LineString':
                endpoint = Point(geom.coords[0])  # Get the beginning of the LineString
                list_upstream_points.append(endpoint)
            elif geom.geom_type == 'MultiLineString':
                for part in geom:
                    endpoint = Point(part.coords[0])  # Get the beginning of each part
                    list_upstream_points.append(endpoint)

    # Create a new GeoDataFrame from the points
    gdf_upstream_points = gpd.GeoDataFrame(geometry=list_upstream_points, crs=gdf_flowpath_with_attrib.crs)

    # Add the id list to the gdf
    gdf_upstream_points['id'] = list_id

    # Create a coloumn that states if a startpoint is within the 2d area
    gdf_upstream_points['within_2darea'] = gdf_upstream_points.geometry.within(gdf_2d_area_polygon_nextgen.iloc[0].geometry)
    
    # add attribute to points from the corresponding stream
    # Specify the columns to keep
    columns_to_keep = ['id', 'mainstem', 'tot_drainage_areasqkm', 'order']

    # Drop all columns except the specified ones
    gdf_flowpath_for_join = gdf_flowpath_with_attrib[columns_to_keep]

    # Perform left join on 'id' to add the selecte attributes
    gdf_upstream_points = gdf_upstream_points.merge(gdf_flowpath_for_join, on='id', how='left')
    
    return(gdf_upstream_points)
# ------------------


# +++++++++++++++++++++++++++++
def fn_create_cell_gdf(hdf5_file_path, str_hdf_folder_2darea):
    """
    Create a GeoDataFrame of cells from HDF5 data.

    Parameters:
    hdf5_file_path (str): The file path to the HDF5 file.
    str_hdf_folder_2darea (str): The folder containing 2D area data within the HDF5 file.

    Returns:
    GeoDataFrame: A GeoDataFrame containing polygons representing cells, 
                  constructed from the face point coordinates extracted from the HDF5 file.
    """
    # location of Face Point Coordiantes in HDF5
    str_facepoint_coords = str_hdf_folder_2darea + 'FacePoints Coordinate'

    # Open the HDF5 file
    with h5py.File(hdf5_file_path, 'r') as hdf_file:
        # Extract X and Y coordinates
        x_coordinates = hdf_file[str_facepoint_coords][:, 0]
        y_coordinates = hdf_file[str_facepoint_coords][:, 1]

    # Create a pandas DataFrame
    df_facepoints = pd.DataFrame({'X': x_coordinates, 'Y': y_coordinates})

    # location of Indecies of facepoints making up the cells
    str_cells_facepoint_indexes = str_hdf_folder_2darea + 'Cells FacePoint Indexes'

    # Open the HDF5 file
    with h5py.File(hdf5_file_path, 'r') as hdf_file:
        # Extract FacePoints Coordinate data
        facepoints_data = hdf_file[str_cells_facepoint_indexes][:]

        # Extract the projection
        projection_wkt = hdf_file.attrs['Projection'].decode('utf-8')

    # Create a pandas DataFrame from the array
    #df_cells_by_facepoints = pd.DataFrame(facepoints_data)

    # Create a GeoDataFrame to store the polygons
    geometry = []

    for row in facepoints_data:
        polygon_coords = []

        for idx in row:
            if idx != -1:
                x = df_facepoints.loc[idx, 'X']
                y = df_facepoints.loc[idx, 'Y']
                polygon_coords.append((x, y))
        # Connect to the first point to close the polygon
        polygon_coords.append(polygon_coords[0])
        geometry.append(Polygon(polygon_coords))

    # Create a GeoDataFrame
    gdf_cells = gpd.GeoDataFrame(geometry=geometry, columns=['geometry'], crs=projection_wkt)
    
    # Add a new column 'idx' to gdf_cells with index values
    gdf_cells['idx'] = gdf_cells.index

    return gdf_cells
# +++++++++++++++++++++++++++++


# ..........................
def fn_assign_mannings_per_stream(gdf_cells, gdf_flowpaths, hdf_file_path):
    
    # Specify the HDF5 file path and group path
    str_hdf_2darea_root_folder = '/Geometry/2D Flow Areas/'

    # Get names of HDF5 Group objects in the specified group
    list_group_names = fn_get_group_names(hdf_file_path, str_hdf_2darea_root_folder)

    str_hdf_folder_2darea = str_hdf_2darea_root_folder + list_group_names[0] + '/'
    
    # Get cell's Mannings at center point
    str_cell_center_mannings = str_hdf_folder_2darea + "Cells Center Manning's n"
    
    # Open the HDF5 file
    with h5py.File(hdf_file_path, 'r') as hdf_file:
        arr_cell_mannings = hdf_file[str_cell_center_mannings][:]

    # Reproject gdf_flowpaths to the CRS of gdf_cells
    gdf_flowpaths_local_crs = gdf_flowpaths.to_crs(gdf_cells.crs)

    # Keep only the 'geometry' and 'id' columns
    gdf_flowpaths_local_crs_lean = gdf_flowpaths_local_crs[['geometry', 'id']]
    
    # about 30 seconds on test grid

    print("   -- Determining stream thalweg Manning's 'n'...")

    # Intersect the cells and the stream centerlines
    gdf_stream_segments_per_cell = gpd.overlay(gdf_flowpaths_local_crs_lean, gdf_cells, how="intersection")

    # compute the length for each segment
    gdf_stream_segments_per_cell['length'] = gdf_stream_segments_per_cell['geometry'].length

    # Set Manning's coefficient for each segment based on 'idx' field
    gdf_stream_segments_per_cell['mannings'] = gdf_stream_segments_per_cell['idx'].map(lambda idx: arr_cell_mannings[int(idx)])

    # Calculate Manning's x length for each segment
    gdf_stream_segments_per_cell['mannings_x_len'] = gdf_stream_segments_per_cell['mannings'] * gdf_stream_segments_per_cell['length']

    # Summary values by stream segment

    # Create a pandas dataframe
    df = gdf_stream_segments_per_cell[['id', 'length', 'mannings_x_len']].copy()

    # Group by 'id' and sum the 'length' and 'mannings_x_len'
    df_stream_summary = df.groupby('id').agg({'length': 'sum', 'mannings_x_len': 'sum'}).reset_index()

    # Rename columns for clarity
    df_stream_summary.columns = ['id', 'total_length', 'total_mannings_x_len']

    # Compute reach averaged Mannings 'n' values (three decimal points)
    df_stream_summary['manning']= round(df_stream_summary['total_mannings_x_len'] / df_stream_summary['total_length'],3)

    # Filter down to just the values required
    df_stream_summary_lean = df_stream_summary[['id', 'manning']]
    
    # Perform the left join operation
    gdf_streams_w_mannings = pd.merge(gdf_flowpaths, df_stream_summary_lean, on='id', how='left')
    
    return(gdf_streams_w_mannings)
# ..........................


# ------------------
def fn_calculate_travel_time(row):
    
    # this is an estimate of "low flow travel time" - This assumes a hydraulic radius
    
    # *************
    #For a low flow travel time, this is the Rh that is being assumed (need a conservative estimate)
    # TODO - MAC - 2024.06.11 - Pass fn_calculate_travel_time into this function
    flt_assumed_hydraulic_radius = 1.0
    # *************
    
    flt_mannings = row['manning']
    flt_length_m = row['lengthkm'] * 1000
    flt_slope_m_per_m = row['So']
    flt_time_sec = (flt_length_m * flt_mannings) / ((flt_assumed_hydraulic_radius ** 0.667) * (flt_slope_m_per_m ** 0.5))
    flt_time_hr = round(flt_time_sec / 3600, 2)
    return flt_time_hr
# ------------------


# ,,,,,,,,,,,,,,,,,,,,,,,,
def fn_travel_time_and_peak_flow_est(gdf_streams,
                                     str_limiting_discharge_csv,
                                     flt_q_ratio,
                                     flt_assumed_hydraulic_radius):
    
 
    # Load the CSV file into a NumPy array
    data = np.genfromtxt(str_limiting_discharge_csv, delimiter=',')
    data = data[data[:,0].argsort()]

    gdf_streams['da_sq_mile'] = gdf_streams['tot_drainage_areasqkm'] * 0.386102

    # Perform linear interpolation for each row in df and round the result to the nearest integer
    gdf_streams['q_limiting'] = gdf_streams['da_sq_mile'].apply(lambda x: round(np.interp(x, data[:,0], data[:,1])))

    # Upper flow limit is ratio of limiting discharge == flt_q_ratio
    gdf_streams['q_upper_limit'] = round(gdf_streams['q_limiting'] * flt_q_ratio,0)

    gdf_streams['travel_time_hr'] = gdf_streams.apply(fn_calculate_travel_time, axis=1)
    
    return(gdf_streams)
# ,,,,,,,,,,,,,,,,,,,,,,,,


# ..........................
def fn_compute_lower_q_limit(gdf_points):
    # for each item in list_unique_mainstem_values, select the rows in gdf_points where gdf_points['mainstem'] == item.
    # sort assending by tot_drainage_areasqkm
    
    list_all_ids = []
    list_all_q_lower_limits = []

    # Get unique values in the 'mainstem' column
    list_unique_mainstem_values = gdf_points['mainstem'].unique().tolist()

    # Iterate over each unique 'mainstem' value
    for flt_mainstem_value in list_unique_mainstem_values:

        # Filter rows where 'mainstem' equals the current unique value
        gdf_filtered_points = gdf_points[gdf_points['mainstem'] == flt_mainstem_value]

        # Sort the filtered rows by 'tot_drainage_areasqkm' in ascending order
        gdf_filtered_points_sorted = gdf_filtered_points.sort_values(by='tot_drainage_areasqkm', ascending=True)

        # create a list of the 'id'
        list_ids = gdf_filtered_points_sorted['id'].tolist()

        # create a list of the 'q_upper_limit'
        list_q_upper_limits = gdf_filtered_points_sorted['q_upper_limit'].tolist()

        # add 0 as first value and remove the last (downstream) value
        list_q_lower_limits = [0] + list_q_upper_limits[:-1]

        # append the values of this list_ids to list_all_ids
        list_all_ids.extend(list_ids)

        # append the values of list_q_lower_limits to list_all_q_lower_limits
        list_all_q_lower_limits.extend(list_q_lower_limits)

    # Create a DataFrame from the lists
    df_lower_limit_flow = pd.DataFrame({
        'id': list_all_ids,
        'q_lower_limit': list_all_q_lower_limits
    })

    # Left join gdf_points with df_lower_limit_flow on 'id'
    gdf_points = gdf_points.merge(df_lower_limit_flow, on='id', how='left')
    
    return(gdf_points)
# ..........................


# ---------------
def fn_end_node_and_travel_time(gdf_points,gdf_streams):
    list_start_node = []
    list_end_node = []
    list_travel_time = []

    for index1,row1 in gdf_points.iterrows():

        str_start_node = row1['id']
        flt_mainstem = row1['mainstem']
        flt_drainage_area = row1['tot_drainage_areasqkm']

        # Filter streams where 'mainstem' equals the current unique value
        gdf_filtered_stream = gdf_streams[gdf_streams['mainstem'] == flt_mainstem]

        # Sort the filtered rows by 'tot_drainage_areasqkm' in ascending order
        gdf_filtered_stream_sorted = gdf_filtered_stream.sort_values(by='tot_drainage_areasqkm', ascending=True)

        # Drop rows where 'tot_drainage_areasqkm' < flt_drainage_area
        gdf_filtered_stream_sorted = gdf_filtered_stream_sorted[gdf_filtered_stream_sorted['tot_drainage_areasqkm'] >= flt_drainage_area]

        # Explore the sorted GeoDataFrame
        gdf_filtered_stream_sorted.explore()

        # Initialize the total travel time
        flt_total_time = 0

        for index, row in gdf_filtered_stream_sorted.iterrows():

            # Check if 'rl_NHDWaterbodyComID' is null
            if pd.isnull(row['rl_NHDWaterbodyComID']):
                flt_total_time += row['travel_time_hr']
                str_end_node = row['id']
            else:
                # Break the loop if 'rl_NHDWaterbodyComID' is not null
                break

        #if flt_total_time == 0: ... this is a waterbody
        list_start_node.append(str_start_node)
        list_end_node.append(str_end_node)
        list_travel_time.append(flt_total_time)

    # Round to two decimal places
    list_travel_time = [round(flt, 2) for flt in list_travel_time]

    # Create a DataFrame from the lists
    df_nodes_and_time = pd.DataFrame({
        'id': list_start_node,
        'id_end_node': list_end_node,
        'travel_time_hr': list_travel_time
    })

    # Left join gdf_points with df_lower_limit_flow on 'id'
    gdf_points = gdf_points.merge(df_nodes_and_time, on='id', how='left')
    
    return(gdf_points)
# ---------------


# ----------------------
def fn_remove_networks_without_outlets(gdf_streams,
                                       str_input_hydrofabric_gpkg, 
                                       str_hdf_file_path,
                                       gdf_2d_area_polygon):
    
    # Get the streams in the 2D area's network that flow though an outlet (Extrenal 2D boundary condition)
    print("   -- Evaluating stream network outlets...")

    # getting a list of nexus id's from the flowpaths 'toid'
    list_unique_nodes = gdf_streams['toid'].unique().tolist()

    # read in the hydrofabric's nexus
    gdf_nexus = gpd.read_file(str_input_hydrofabric_gpkg, layer='nexus')

    # filter to just nexus in list_unique_nodes
    gdf_nexus_on_streams = gdf_nexus[gdf_nexus['id'].isin(list_unique_nodes)]

    # Combine points and lines into a single GeoDataFrame - preperation for graph creation
    gdf_flow_network = gpd.GeoDataFrame(pd.concat([gdf_nexus_on_streams, gdf_streams], ignore_index=True),
                                        crs=gdf_streams.crs)

    # Create a directed graph
    G = nx.DiGraph()

    # Add edges to the graph based on 'id' and 'toid' fields
    edges = gdf_flow_network[['id', 'toid']].values.tolist()
    G.add_edges_from(edges)

    # Now, let's add attributes to the nodes using the additional columns
    for index, row in gdf_flow_network.iterrows():
        node_id = row['id']
        attributes = {'mainstem': row['mainstem'],
                    'tot_drainage_areasqkm': row['tot_drainage_areasqkm']}
        G.nodes[node_id].update(attributes)

    # Calculate weakly connected components
    weakly_connected_components = list(nx.weakly_connected_components(G))

    # Create a dictionary to map node ID to weakly connected component ID
    component_map = {}
    for idx, component in enumerate(weakly_connected_components):
        for node in component:
            component_map[node] = idx

    # Add a new column to gdf_streams to identify the weakly connected component
    gdf_streams['network_group'] = gdf_streams['toid'].map(component_map)

    # Calculate the number of weakly connected components
    num_components = nx.number_weakly_connected_components(G)
    #print("   Number of connected stream networks:", num_components)

    # list of connected network indecies...not yet checking outlets
    list_networks = gdf_streams['network_group'].unique().tolist()

    # ---- reading the boundary condition lines from HEC-RAS HDF ----
    # create a gdf of the boundary conditions
    str_hdf_boundary_conditions_path = '/Geometry/Boundary Condition Lines/'

    str_bc_attrib = str_hdf_boundary_conditions_path + 'Attributes'
    str_polyline_info = str_hdf_boundary_conditions_path + 'Polyline Info'
    str_polyline_parts = str_hdf_boundary_conditions_path + 'Polyline Parts'
    str_polyline_points = str_hdf_boundary_conditions_path + 'Polyline Points'

    # Open the HDF file
    with h5py.File(str_hdf_file_path, 'r') as hdf_file:
        arr_bc_attributes = hdf_file[str_bc_attrib][:]
        arr_polyline_info = hdf_file[str_polyline_info][:]
        arr_polyline_parts = hdf_file[str_polyline_parts][:]
        arr_polyline_points = hdf_file[str_polyline_points][:]


    if len(arr_bc_attributes) == len(arr_polyline_parts):
        #print('   Note: Only one part per boundary condition line')

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # create a shapely list of boundary condition polylines
        list_polyline_shp = [
            LineString(arr_polyline_points[row[0]:row[0] + row[1]])
            for row in arr_polyline_info]

        # Decode byte strings and convert to DataFrame from arr_bc_attributes
        df_bc_attrib = pd.DataFrame(arr_bc_attributes)

        # Decode byte strings in the DataFrame
        for col in df_bc_attrib.columns:
            if df_bc_attrib[col].dtype == object:
                df_bc_attrib[col] = df_bc_attrib[col].str.decode('utf-8')

        # Drop the "length" column from df_bc_attrib
        df_bc_attrib = df_bc_attrib.drop(columns=['Length'])   

        # Create a geodataframe of teh boundary condition lines
        gdf_bc_lines = gpd.GeoDataFrame(df_bc_attrib, geometry=list_polyline_shp, crs=gdf_2d_area_polygon.crs)

        # Reproject gdf_bc_lines to the CRS of gdf_streams
        gdf_bc_lines_nextgen = gdf_bc_lines.to_crs(gdf_streams.crs)

        str_2d_area_name = gdf_2d_area_polygon.iloc[0]['area_2d_name']

        # filter the gdf_bc_lines_nextgen to the SA-2D == str_2d_area_name and Type == "External"
        gdf_external_bc_on_area = gdf_bc_lines_nextgen[(gdf_bc_lines_nextgen['SA-2D'] == str_2d_area_name) & (gdf_bc_lines_nextgen['Type'] == 'External')]

        # determine what "weakly connected components" directed graphs flow through a external boundary condition

        # intersect gdf_streams with gdf_external_bc_on_area
        gdf_graphs_with_outflow = gpd.overlay(gdf_streams,
                                              gdf_external_bc_on_area,
                                              how='intersection',
                                              keep_geom_type=False)

        # get a unique list from network_group coloumn from gdf_graphs_with_outflow
        list_network_groups_w_outlets = gdf_graphs_with_outflow['network_group'].unique().tolist()

        # Determine list of networks without outlets
        # TODO - MAC - 20240911 - Error Here -- not determining netowrks correctly?
        list_networks_wo_outlets = [item for item in list_networks if item not in list_network_groups_w_outlets]

        #if len(list_networks_wo_outlets) > 0:
        #    print(f'Removing {len(list_networks_wo_outlets)} networks with no outlets')
        #    for item in list_networks_wo_outlets:
        #        # remove all streams where network has no outlet
        #        gdf_streams = gdf_streams[gdf_streams['network_group'] != item]
        #else:
        #    pass
        #    print('All streams networks have outlets')
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    else:
        print("Error: boundary condition lines contain multiple parts")

    return(gdf_streams)
# ----------------------


# ''''''''''''''''''''''''''''''
def fn_create_eval_travel_path(gdf_streams, gdf_points):
    
    # Filter columns for simplicity
    gdf_stream_simple = gdf_streams[['mainstem', 'geometry', 'rl_NHDWaterbodyComID', 'tot_drainage_areasqkm']].copy()

    # Convert MultiLineStrings to LineStrings
    gdf_stream_simple['geometry'] = gdf_stream_simple['geometry'].apply(
        lambda geom: geom if geom.geom_type == 'LineString' else geom.geoms[0] if geom.geom_type == 'MultiLineString' else None
    )

    # Remove rows where geometry could not be converted
    gdf_stream_simple = gdf_stream_simple[gdf_stream_simple['geometry'].notnull()]

    # Initialize the list to collect rows for the final GeoDataFrame
    eval_travel_path_data = []

    for index1, row1 in gdf_points.iterrows():
        str_start_node = row1['id']
        flt_mainstem = row1['mainstem']
        flt_drainage_area = row1['tot_drainage_areasqkm']

        # Filter streams by 'mainstem' and drainage area
        gdf_filtered_stream_sorted = (
            gdf_stream_simple[gdf_stream_simple['mainstem'] == flt_mainstem]
            .sort_values(by='tot_drainage_areasqkm')
            .query('tot_drainage_areasqkm >= @flt_drainage_area')
            .reset_index(drop=True)
        )

        b_found_waterbody = False

        # TODO -- 2024.11.05 -- need to remove waterbodies!!
        
        # Travel path stops at first waterbody found
        for index, row in gdf_filtered_stream_sorted.iterrows():
            if pd.notnull(row['rl_NHDWaterbodyComID']):
                b_found_waterbody = True
                break
                
        if b_found_waterbody:
            gdf_filtered_stream_sorted = gdf_filtered_stream_sorted.iloc[:index]

        combined_linestring = gdf_filtered_stream_sorted['geometry'].unary_union
        
        # on waterbodies... the combined_linestring might be empty (None)
        if combined_linestring and not combined_linestring.is_empty:
            if combined_linestring.geom_type == 'MultiLineString':
                combined_linestring = LineString([pt for line in combined_linestring.geoms for pt in line.coords])

            eval_travel_path_data.append({
                'id': str_start_node,
                'mainstem': flt_mainstem,
                'q_lower_limit': row1['q_lower_limit'],
                'q_upper_limit': row1['q_upper_limit'],
                'q_limiting': row1['q_limiting'],
                'travel_time_hr': row1['travel_time_hr'],
                'geometry': combined_linestring
            })

    # Create final GeoDataFrame
    gdf_eval_travel_path = gpd.GeoDataFrame(eval_travel_path_data, crs=gdf_stream_simple.crs)

    # Retain only rows in gdf_eval_travel_path where travel_time_hr not equal to 0
    gdf_eval_travel_path = gdf_eval_travel_path[gdf_eval_travel_path['travel_time_hr'] != 0]
    
    return(gdf_eval_travel_path)
# ''''''''''''''''''''''''''''''


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


# .........................................................
def fn_create_models_nextgen_gpkg(str_config_file_path,
                                  str_hdf_file_path,
                                  str_output_dir,
                                  b_print_output):
    
    
    # supress all warnings
    warnings.filterwarnings("ignore", category=UserWarning )
    
    print(" ")
    if b_print_output:
        print("+=================================================================+")
        print("|         EXTRACT NEXTGEN HYDROFABRIC FOR HEC-RAS 2D AREA         |")
        print("|                Created by Andy Carter, PE of                    |")
        print("|             Center for Water and the Environment                |")
        print("|                 University of Texas at Austin                   |")
        print("+-----------------------------------------------------------------+")
    
        
        print("  ---(c) INPUT GLOBAL CONFIGURATION FILE: " + str_config_file_path)
        print("  ---(i) INPUT HEC-RAS PLAN (HDF): " + str_hdf_file_path)
        print("  ---(o) OUTPUT PATH: " + str_output_dir)
        print("  ---[r] PRINT OUTPUT: " + str(b_print_output))
        print("===================================================================")
    else:
        print('Step 1: Extract NextGEN hydrofabric from HEC-RAS Model')
    
    if not os.path.exists(str_output_dir):
        os.mkdir(str_output_dir)
        
    # parse str_input_hydrofabric_gpkg and str_limiting_discharge_csv from config.ini
    
    # --- Read variables from config.ini ---
    # Create a ConfigParser object
    config = configparser.ConfigParser()
    
    # Read the config.ini file
    config.read(str_config_file_path)
    
    # Get the variables in the [02_copy_geom] section
    if 'global_input_files' in config:
        section = config['global_input_files']
        
        # Read variables in the section
        str_input_hydrofabric_gpkg = section.get('str_nextgen_gpkg', '')
        str_limiting_discharge_csv = section.get('str_limiting_discharge_csv', '')
    else:
        print("[global_input_files] section not found in the config file.")
        
        
    # Get the [01_nextgen_hydrofabric] variables... flt_q_ratio = 0.30 & flt_assumed_hydraulic_radius = 1.0
    # Check config file for values
    # ----
    # Get the variables in the [01_nextgen_hydrofabric] section
    if '01_nextgen_hydrofabric' in config:
        section = config['01_nextgen_hydrofabric']
        
        # Read variables in the section
        flt_q_ratio = section.getfloat('flt_q_ratio', 0)
        flt_assumed_hydraulic_radius = section.getfloat('flt_assumed_hydraulic_radius', 0)
    else:
        print("[01_nextgen_hydrofabric] section not found in the config file.")
    # ----
    
    gdf_2d_area_polygon = fn_get_gdf_of_2d_area(str_hdf_file_path)
    gdf_flowpath_with_attrib = fn_extract_flowpaths_gdf(gdf_2d_area_polygon, str_input_hydrofabric_gpkg)
    gdf_2d_area_polygon_nextgen = gdf_2d_area_polygon.to_crs(gdf_flowpath_with_attrib.crs)
    
    gdf_upstream_points = fn_create_upstream_flowpath_points(gdf_flowpath_with_attrib,
                                                             gdf_2d_area_polygon_nextgen)
    
    print('   -- Extracting HEC-RAS computational cell polygons...')

    # Specify the HDF5 file path and group path
    str_hdf_2darea_root_folder = '/Geometry/2D Flow Areas/'

    # Get names of HDF5 Group objects in the specified group
    list_group_names = fn_get_group_names(str_hdf_file_path, str_hdf_2darea_root_folder)

    str_hdf_folder_2darea = str_hdf_2darea_root_folder + list_group_names[0] + '/'
    gdf_cells = fn_create_cell_gdf(str_hdf_file_path, str_hdf_folder_2darea)
    
    # As the cells are already extracted, it is faster to do this now
    gdf_streams_w_mannings = fn_assign_mannings_per_stream(gdf_cells, gdf_flowpath_with_attrib, str_hdf_file_path)

    gdf_streams_w_traveltimes =fn_travel_time_and_peak_flow_est(gdf_streams_w_mannings,
                                                               str_limiting_discharge_csv,
                                                               flt_q_ratio,
                                                               flt_assumed_hydraulic_radius)
    
    # Perform a left join on 'id', adding 'q_limiting' and 'q_upper_limit' to gdf_points
    gdf_upstream_points = gdf_upstream_points.merge(
        gdf_streams_w_traveltimes[['id', 'q_limiting', 'q_upper_limit']],
        on='id',
        how='left')
    
    # get the lower limit flow for each point
    gdf_upstream_points = fn_compute_lower_q_limit(gdf_upstream_points)
    
    # get the downstream node of the run and the travel time to that node
    gdf_upstream_points = fn_end_node_and_travel_time(gdf_upstream_points, gdf_streams_w_traveltimes)
    
    # Filter out streams that don't have outlets (cross boundary conditions)
    gdf_streams_w_traveltimes = fn_remove_networks_without_outlets(gdf_streams_w_traveltimes,
                                                                   str_input_hydrofabric_gpkg, 
                                                                   str_hdf_file_path,
                                                                   gdf_2d_area_polygon)
    
    
    # Keep only rows in gdf_upstream_points where the 'id' exists in gdf_streams_w_traveltimes
    # this is to remove streams on networks without outlets
    gdf_upstream_points = gdf_upstream_points[gdf_upstream_points['id'].isin(gdf_streams_w_traveltimes['id'])]
    
    # Keep only rows where 'within_2darea' is True
    gdf_upstream_points = gdf_upstream_points[gdf_upstream_points['within_2darea'] == True]
    
    gdf_run_stream_paths = fn_create_eval_travel_path(gdf_streams_w_traveltimes, gdf_upstream_points)
    
    print('   -- Writing output geopackage...')
    
    str_folder_for_gpkg = os.path.join(str_output_dir, '01_model_hydrofabric')
    # If the folder does not exist, create it
    os.makedirs(str_folder_for_gpkg, exist_ok=True)
    
    # HARD CODED OPTPUT NAME OF GEOPACKAGE
    str_output_geopackage_name = 'model_hydrofabric.gpkg'
    str_gpkg_to_create = os.path.join(str_folder_for_gpkg, str_output_geopackage_name)
    
    gdf_2d_area_polygon_nextgen = gdf_2d_area_polygon.to_crs(gdf_flowpath_with_attrib.crs)

    # Write GeoDataFrames to GeoPackage as separate layers
    
    gdf_run_stream_paths.to_file(str_gpkg_to_create, layer='03_flowpaths', driver="GPKG")
    gdf_upstream_points.to_file(str_gpkg_to_create, layer='02_flow_points', driver="GPKG")
    gdf_streams_w_traveltimes.to_file(str_gpkg_to_create, layer='01_stream_lines', driver="GPKG")
    gdf_2d_area_polygon_nextgen.to_file(str_gpkg_to_create, layer='00_area_2d', driver="GPKG")
    
    #print('Complete')
    print("+-----------------------------------------------------------------+")
    
    return(gdf_2d_area_polygon,gdf_upstream_points, gdf_streams_w_traveltimes)
# .........................................................


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()
    
    parser = argparse.ArgumentParser(description='========= EXTRACT NEXTGEN HYDROFABRIC FOR HEC-RAS 2D AREA =========')
    
    parser.add_argument('-c',
                        dest = "str_config_file_path",
                        help=r'REQUIRED: Global configuration filepath Example: C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini',
                        required=True,
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-i',
                        dest = "str_hdf_file_path",
                        help=r'REQUIRED: path to plan HDF file (HDF) Example: E:\ras2fim_test_20250107\00_base_hecras\source_ras.g02.hdf',
                        required=True,
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-o',
                        dest = "str_output_dir",
                        help=r'REQUIRED: output directory of RAS2FIM-2D output Example: E:\ras2fim_test_20250107',
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
    str_output_dir = args['str_output_dir']
    b_print_output = args['b_print_output']

    fn_create_models_nextgen_gpkg(str_config_file_path,
                                  str_hdf_file_path,
                                  str_output_dir)

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)
    
    print('Compute Time: ' + str(time_pass))
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~