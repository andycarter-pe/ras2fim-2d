# Script 01 - Create a geopackage of the NextGen stream network from the HEC-RAS Model
#
# This requires as HEC-RAS 2D model in an HDF plan containing at least one(1)
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
# Created - 2024.05.08
# Last revised - 2024.05.11

# ************************************************************
import h5py

from shapely.geometry import Point, LineString, Polygon, MultiPolygon, box, mapping, shape
import geopandas as gpd
import pandas as pd
import networkx as nx
import numpy as np
import os
import configparser

from shapely.geometry import MultiLineString, LineString
from shapely.ops import linemerge

from shapely.ops import nearest_points
from shapely.geometry import Point

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
        print('Multiple 2D areas found -- Using first area:', list_group_names[0])
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
    print('Reading Nextgen Hydrofabric ~5 seconds...')

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


# ----------------------
def fn_remove_networks_without_outlets(gdf_flowpath_with_attrib,
                                        gpkg_nextgen_path, 
                                        hdf_file_path,
                                        gdf_2d_area_polygon):
    
    # Get the streams in the 2D area's network that flow though an outlet (Extrenal 2D boundary condition)
    print("   Evaluating stream network outlets...")

    # getting a list of nexus id's from the flowpaths 'toid'
    list_unique_nodes = gdf_flowpath_with_attrib['toid'].unique().tolist()

    # read in the hydrofabric's nexus
    gdf_nexus = gpd.read_file(gpkg_nextgen_path, layer='nexus')

    # filter to just nexus in list_unique_nodes
    gdf_nexus_on_streams = gdf_nexus[gdf_nexus['id'].isin(list_unique_nodes)]

    # Combine points and lines into a single GeoDataFrame - preperation for graph creation
    gdf_flow_network = gpd.GeoDataFrame(pd.concat([gdf_nexus_on_streams, gdf_flowpath_with_attrib], ignore_index=True),
                                        crs=gdf_flowpath_with_attrib.crs)

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

    # Add a new column to gdf_flowpath_with_attrib to identify the weakly connected component
    gdf_flowpath_with_attrib['network_group'] = gdf_flowpath_with_attrib['toid'].map(component_map)

    # Calculate the number of weakly connected components
    #num_components = nx.number_weakly_connected_components(G)
    #print("   Number of connected stream networks:", num_components)

    # list of connected network indecies...not yet checking outlets
    list_networks = gdf_flowpath_with_attrib['network_group'].unique().tolist()

    # ---- reading the boundary condition lines from HEC-RAS HDF ----
    # create a gdf of the boundary conditions
    str_hdf_boundary_conditions_path = '/Geometry/Boundary Condition Lines/'

    str_bc_attrib = str_hdf_boundary_conditions_path + 'Attributes'
    str_polyline_info = str_hdf_boundary_conditions_path + 'Polyline Info'
    str_polyline_parts = str_hdf_boundary_conditions_path + 'Polyline Parts'
    str_polyline_points = str_hdf_boundary_conditions_path + 'Polyline Points'

    # Open the HDF file
    with h5py.File(hdf_file_path, 'r') as hdf_file:
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

        # Reproject gdf_bc_lines to the CRS of gdf_flowpath_with_attrib
        gdf_bc_lines_nextgen = gdf_bc_lines.to_crs(gdf_flowpath_with_attrib.crs)

        str_2d_area_name = gdf_2d_area_polygon.iloc[0]['area_2d_name']

        # filter the gdf_bc_lines_nextgen to the SA-2D == str_2d_area_name and Type == "External"
        gdf_external_bc_on_area = gdf_bc_lines_nextgen[(gdf_bc_lines_nextgen['SA-2D'] == str_2d_area_name) & (gdf_bc_lines_nextgen['Type'] == 'External')]

        # determine what "weakly connected components" directed graphs flow through a external boundary condition

        # intersect gdf_flowpath_with_attrib with gdf_external_bc_on_area
        gdf_graphs_with_outflow = gpd.overlay(gdf_flowpath_with_attrib,
                                              gdf_external_bc_on_area,
                                              how='intersection',
                                              keep_geom_type=False)

        # get a unique list from network_group coloumn from gdf_graphs_with_outflow
        list_network_groups_w_outlets = gdf_graphs_with_outflow['network_group'].unique().tolist()

        # Determine list of networks without outlets
        list_networks_wo_outlets = [item for item in list_networks if item not in list_network_groups_w_outlets]

        if len(list_networks_wo_outlets) > 0:
            print(f'Removing {len(list_networks_wo_outlets)} networks with no outlets')
            for item in list_networks_wo_outlets:
                # remove all streams where network has no outlet
                gdf_flowpath_with_attrib = gdf_flowpath_with_attrib[gdf_flowpath_with_attrib['network_group'] != item]
        else:
            pass
            #print('All streams have outlets')
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    else:
        print("Error: boundary condition lines contain multiple parts")
        
    return(gdf_flowpath_with_attrib)
# ----------------------


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

    print("Determining stream thalweg Manning's 'n'...")

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


# ---------------
def fn_determine_starting_hecras_cells(hdf_file_path,
                                       gdf_upstream_points,
                                       gdf_flowpath_with_attrib):

    # get the HEC-RAS cell number for each upstream point

    print('Extracting HEC-RAS computational cell polygons...')

    # Specify the HDF5 file path and group path
    str_hdf_2darea_root_folder = '/Geometry/2D Flow Areas/'

    # Get names of HDF5 Group objects in the specified group
    list_group_names = fn_get_group_names(hdf_file_path, str_hdf_2darea_root_folder)

    str_hdf_folder_2darea = str_hdf_2darea_root_folder + list_group_names[0] + '/'
    gdf_cells = fn_create_cell_gdf(hdf_file_path, str_hdf_folder_2darea)

    print(f'  -- Number of cells in {list_group_names[0]}: {len(gdf_cells)}')

    # create points that are in the same projection as the HEC-RAS Model
    gdf_upstream_points_hecras_projection = gdf_upstream_points.to_crs(gdf_cells.crs)

    print(f'  -- Determining starting cell for {len(gdf_upstream_points_hecras_projection)} points...')

    # Create a spatial index for the polygon GeoDataFrame
    sindex = gdf_cells.sindex

    # Create an empty list to store the indices of cells that each point is inside
    cell_indices = []

    # Iterate over each point in the point GeoDataFrame
    for point in gdf_upstream_points_hecras_projection.geometry:
        # Find the index of the cell that contains the point
        possible_matches_index = list(sindex.intersection(point.bounds))
        possible_matches = gdf_cells.iloc[possible_matches_index]
        precise_matches = possible_matches[possible_matches.contains(point)]
        if len(precise_matches) > 0:
            # Append the index of the cell to the list
            cell_indices.append(precise_matches.index[0])
        else:
            # Point is not inside any cell
            cell_indices.append(None)

    # Add the cell indices as a new column to the point GeoDataFrame
    gdf_upstream_points['idx_start_cell'] = cell_indices

    # Filter the for only points within 2darea
    gdf_upstream_points_inside_area = gdf_upstream_points[(gdf_upstream_points['within_2darea'] == True)]
    
    # ----- Now, get the Mannings n for the thalweg for each stream -----
    # As the cells are already extracted, it is faster to do this now
    gdf_streams_w_mannings = fn_assign_mannings_per_stream(gdf_cells, gdf_flowpath_with_attrib, hdf_file_path)
    # ----------

    return(gdf_upstream_points_inside_area, gdf_streams_w_mannings)
# ---------------


# --------------------
def fn_find_most_downstream_node(graph, start_node):
    visited = set()
    stack = [start_node]

    while stack:
        current_node = stack.pop()
        if current_node not in visited:
            visited.add(current_node)
            neighbors = list(graph.successors(current_node))  # Get successors (outgoing edges)
            stack.extend(neighbors)

    return max(visited)  # Return the highest node ID found
# ----------------------


# ----------------
def fn_find_path_between_nodes(graph, start_node, end_node):
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


# ----------------
def fn_compute_travel_distance_per_point(gdf_flowpath_with_attrib, 
                                         gpkg_nextgen_path, 
                                         gdf_upstream_points):

    print('Compute travel distance per simulation...')
    
    # getting a list of nexus id's from the flowpaths 'toid'
    list_unique_nodes = gdf_flowpath_with_attrib['toid'].unique().tolist()

    # read in the hydrofabric's nexus
    gdf_nexus = gpd.read_file(gpkg_nextgen_path, layer='nexus')

    # filter to just nexus in list_unique_nodes
    gdf_nexus_on_streams = gdf_nexus[gdf_nexus['id'].isin(list_unique_nodes)]

    # Combine points and lines into a single GeoDataFrame - preperation for graph creation
    gdf_flow_network = gpd.GeoDataFrame(pd.concat([gdf_nexus_on_streams, gdf_flowpath_with_attrib], ignore_index=True),
                                        crs=gdf_flowpath_with_attrib.crs)

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

    # -- Some Graph Statistics ---
    print('Hydrofabric Statistics:')
    # Number of nodes
    num_nodes = G.number_of_nodes()
    print("   Number of nodes:", num_nodes)

    # Number of edges
    num_edges = G.number_of_edges()
    print("   Number of streams:", num_edges)

    # Average degree
    avg_degree = round(sum(dict(G.degree()).values()) / num_nodes, 2)
    print("   Average degree:", avg_degree)

    # Calculate the number of weakly connected components
    num_components = nx.number_weakly_connected_components(G)
    print("   Number of connected networks:", num_components)
    # -- End Some Graph Statistics ---
    
    # -- Determine travel distance for each node in gdf_upstream_points
    list_travel_dist_flt = []
    list_avg_slope = []

    for index, row in gdf_upstream_points.iterrows():

        str_start_node = row['id']
        str_downstream_node = fn_find_most_downstream_node(G, str_start_node)

        # Create a list of nodes (flowpath and nexus) down the travel path
        list_travel_path = fn_find_path_between_nodes(G, str_start_node, str_downstream_node)

        # get just the flowpaths from gdf_flowpath_with_attrib in list_travel_path
        gdf_streams_travel_path = gdf_flowpath_with_attrib[gdf_flowpath_with_attrib['id'].isin(list_travel_path)]

        # vertical fall per flowpath
        gdf_streams_travel_path = gdf_streams_travel_path.copy()
        gdf_streams_travel_path['delta_z'] = gdf_streams_travel_path['lengthkm'] * gdf_streams_travel_path['So']

        # Sum the values in the 'delta_z' column and rounded 
        flt_total_delta_z = round(gdf_streams_travel_path['delta_z'].sum(), 5)

        # Sum the values in the 'lengthkm' column and round to three decimal places
        flt_total_length_km = round(gdf_streams_travel_path['lengthkm'].sum(), 3)

        # Sum the values in the 'lengthkm' column and round to three decimal places
        flt_avg_slope = flt_total_delta_z / flt_total_length_km

        list_travel_dist_flt.append(flt_total_length_km)
        list_avg_slope.append(flt_avg_slope)

    gdf_upstream_points['flow_dist_km'] = list_travel_dist_flt
    gdf_upstream_points['slope_average'] = list_avg_slope
    
    return(gdf_upstream_points)
# ----------------


# ------------------
def fn_calculate_travel_time(row):
    
    # this is an estimate of "low flow travel time" - This assumes a hydraulic radius
    
    # *************
    #For a low flow travel time, this is the Rh that is being assumed (need a conservative estimate)
    flt_assumed_hydraulic_radius = 1.0
    # *************
    
    flt_mannings = row['manning']
    flt_length_m = row['lengthkm'] * 1000
    flt_slope_m_per_m = row['So']
    flt_time_sec = (flt_length_m * flt_mannings) / ((flt_assumed_hydraulic_radius ** 0.667) * (flt_slope_m_per_m ** 0.5))
    flt_time_hr = round(flt_time_sec / 3600, 2)
    return flt_time_hr
# ------------------

# --------
def fn_convert_to_list_of_linestrings(geometry):
    if isinstance(geometry, MultiLineString):
        result = []
        for part in geometry.geoms:
            result.append(LineString(part.coords))
        return result
    elif isinstance(geometry, LineString):
        return [geometry]
    else:
        return []
# --------


# ,,,,,,,,,,,,,,,,,,,,,,,,
def fn_travel_time_and_peak_flow_est(gdf_streams, gdf_flow_points, str_limiting_discharge_csv):
    
    # *************
    # Assumed ration of the limiting discharge for peak flow estimate
    flt_q_ratio = 0.30
    # *************
    
    # Load the CSV file into a NumPy array
    data = np.genfromtxt(str_limiting_discharge_csv, delimiter=',')
    data = data[data[:,0].argsort()]

    gdf_streams['da_sq_mile'] = gdf_streams['tot_drainage_areasqkm'] * 0.386102

    # Perform linear interpolation for each row in df and round the result to the nearest integer
    gdf_streams['q_limiting'] = gdf_streams['da_sq_mile'].apply(lambda x: round(np.interp(x, data[:,0], data[:,1])))

    # Upper flow limit is ratio of limiting discharge == flt_q_ratio
    gdf_streams['q_upper_limit'] = round(gdf_streams['q_limiting'] * flt_q_ratio,0)

    # Assuming gdf_streams is your GeoDataFrame
    gdf_streams['travel_time_hr'] = gdf_streams.apply(fn_calculate_travel_time, axis=1)
    
    # Dissolve the travel paths (flowpaths) by 'mainstem' attribute
    gdf_mainstems = gdf_streams.dissolve(by='mainstem')

    # Reset the index
    gdf_mainstems.reset_index(inplace=True)

    # Keep only the 'mainstem' and 'geometry' columns
    gdf_mainstems = gdf_mainstems[['mainstem', 'geometry']]

    # Assuming gdf_mainstems is a GeoDataFrame with geometry column containing LineString or MultiLineString
    for idx, row in gdf_mainstems.iterrows():
        geom = row['geometry']
        if isinstance(geom, MultiLineString):
            merged_line = linemerge(geom)
            if merged_line.is_empty:
                pass
            else:
                # Replace the MultiLineString with the merged LineString
                gdf_mainstems.at[idx, 'geometry'] = merged_line
                
    # Filter rows where rl_NHDWaterbodyComID is not null
    gdf_waterbody_flowpaths = gdf_streams[gdf_streams['rl_NHDWaterbodyComID'].notnull()]

    # Dissolve the travel paths by 'mainstem' attribute
    gdf_waterbody_flowpaths_disolve = gdf_waterbody_flowpaths.dissolve(by='mainstem')

    # Reset the index
    gdf_waterbody_flowpaths_disolve.reset_index(inplace=True)

    # Keep only the 'mainstem' and 'geometry' columns
    gdf_waterbody_flowpaths_disolve = gdf_waterbody_flowpaths_disolve[['mainstem', 'geometry']]

    # Assuming gdf_mainstems is a GeoDataFrame with geometry column containing LineString or MultiLineString
    for idx, row in gdf_waterbody_flowpaths_disolve.iterrows():
        geom = row['geometry']
        if isinstance(geom, MultiLineString):
            merged_line = linemerge(geom)
            if merged_line.is_empty:
                pass
            else:
                # Replace the MultiLineString with the merged LineString
                gdf_waterbody_flowpaths_disolve.at[idx, 'geometry'] = merged_line
                
    gdf_mainstems_revised = gdf_mainstems
    
    for index, row in gdf_waterbody_flowpaths_disolve.iterrows():
    
        # Extract the mainstem value
        mainstem_value = row['mainstem']

        # Extract the corresponding linestring from gdf_mainstems
        mainstem_linestring = gdf_mainstems[gdf_mainstems['mainstem'] == mainstem_value].geometry.iloc[0]

        # Clip the linestring from gdf_waterbody_flowpaths with the linestring from gdf_waterbody_flowpaths
        mainstem_differance = mainstem_linestring.difference(row.geometry)

        list_linestrings = fn_convert_to_list_of_linestrings(mainstem_differance)

        if len(list_linestrings) > 0:
            # Delete row in gdf_mainstems_revied where 'mainstem' = mainstem_value
            gdf_mainstems_revised = gdf_mainstems_revised[gdf_mainstems_revised['mainstem'] != mainstem_value]

            # Append items in list_linestrings to gdf_mainstems_revied as new rows
            for linestring in list_linestrings:
                gdf_mainstems_revised = pd.concat([gdf_mainstems_revised,
                                                   pd.DataFrame({'mainstem': [mainstem_value], 'geometry': [linestring]})],
                                                  ignore_index=True)

    gdf_mainstems_revised.crs = gdf_mainstems.crs
    
    # Determine the lines in gdf_streams that intersect the line in gdf_mainstems_revised
    # They must have the same 'mainstem' value
    # return a dataframe of the intersecting lines.

    # For each mainstem (with Waterbodies removed, determine the low flow travel time estiamte for that mainstem run)
    list_travel_time_hr = []

    for index, row in gdf_mainstems_revised.iterrows():

        # Get the value of 'mainstem' from the first row of gdf_mainstems_revised
        mainstem_value = row['mainstem']

        # Select all rows from gdf_streams where 'mainstem' matches mainstem_value
        selected_rows = gdf_streams[gdf_streams['mainstem'] == mainstem_value]

        # Get the geometry from gdf_mainstems_revised
        mainstems_geometry = row['geometry']

        # Check which rows in selected_rows have geometries that are completely covered by mainstems_geometry
        covered_rows = selected_rows[selected_rows.geometry.apply(lambda x: x.covered_by(mainstems_geometry))]

        # Sort covered_rows by 'da_sq_mile' in ascending order
        sorted_covered_rows = covered_rows.sort_values(by='da_sq_mile', ascending=True)

        # Sum the 'travel_time_hr' column in sorted_covered_rows
        total_travel_time_hr = round(sorted_covered_rows['travel_time_hr'].sum(),2)

        # create a list of the travel times
        list_travel_time_hr.append(total_travel_time_hr)
        
    gdf_mainstems_revised['travel_time_hr'] = list_travel_time_hr
    
    # ---- limiting dicharge per points ----
    gdf_flow_points['da_sq_mile'] = gdf_flow_points['tot_drainage_areasqkm'] * 0.386102

    # Perform linear interpolation for each row in df and round the result to the nearest integer
    gdf_flow_points['q_limiting'] = gdf_flow_points['da_sq_mile'].apply(lambda x: round(np.interp(x, data[:,0], data[:,1])))

    # Upper flow limit is ratio of limiting discharge == flt_q_ratio
    gdf_flow_points['q_upper_limit'] = round(gdf_flow_points['q_limiting'] * flt_q_ratio,0)
    # ---- ----

    # Add upper limit flow and starting cell index stabilizing run's flowpath
    list_start_cell_index = []
    list_q_upper_limit = []
    list_node_name = []

    # get the starting point's cell index for each stabilizing run
    for idx, row in gdf_mainstems_revised.iterrows():
        mainstem = row['mainstem']
        linestring = row['geometry']
        start_point = Point(linestring.coords[0])
        points_same_mainstem = gdf_flow_points[gdf_flow_points['mainstem'] == mainstem]['geometry']
        nearest_point_idx = None
        min_distance = float('inf')
        for point in points_same_mainstem:
            distance = start_point.distance(point)
            if distance < min_distance:
                min_distance = distance
                nearest_point = point
                nearest_point_idx = gdf_flow_points[gdf_flow_points['geometry'] == nearest_point].index[0]

        if nearest_point_idx != None:
            # append the starting cell index
            int_start_cell = int(gdf_flow_points.iloc[nearest_point_idx]['idx_start_cell'])
            list_start_cell_index.append(int_start_cell)

            # append the upper flow limit for that reach
            flt_q_upper_limit = gdf_flow_points.iloc[nearest_point_idx]['q_upper_limit']
            list_q_upper_limit.append(flt_q_upper_limit)
            
            # append the node name to list for that reach
            str_node_name = gdf_flow_points.iloc[nearest_point_idx]['id']
            list_node_name.append(str_node_name)
        else:
            list_start_cell_index.append(None)
            list_q_upper_limit.append(None)
            list_node_name.append(None)

    gdf_mainstems_revised['idx_start_cell'] = list_start_cell_index
    gdf_mainstems_revised['q_upper_limit'] = list_q_upper_limit
    gdf_mainstems_revised['id_start_node'] = list_node_name
    
    return(gdf_streams, gdf_flow_points, gdf_mainstems_revised)

# .........................................................
def fn_create_models_nextgen_gpkg(str_config_file_path,
                                  str_hdf_file_path,
                                  str_output_dir):
    
    
    # supress all warnings
    warnings.filterwarnings("ignore", category=UserWarning )
    
    print(" ")
    print("+=================================================================+")
    print("|         EXTRACT NEXTGEN HYDROFABRIC FOR HEC-RAS 2D AREA         |")
    print("|                Created by Andy Carter, PE of                    |")
    print("|             Center for Water and the Environment                |")
    print("|                 University of Texas at Austin                   |")
    print("+-----------------------------------------------------------------+")

    
    print("  ---(c) INPUT GLOBAL CONFIGUTATION FILE: " + str_config_file_path)
    print("  ---(i) INPUT HEC-RAS PLAN (HDF): " + str_hdf_file_path)
    print("  ---(o) OUTPUT PATH: " + str_output_dir)
    print("===================================================================")
    
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
        
        
    # TODO - 2024.05.13 - Get the [01_nextgen_hydrofabric] variables... flt_q_ratio = 0.30 & flt_assumed_hydraulic_radius = 1.0
    # instead of hard coding these inside of this file
    # --- ---
        
    gdf_2d_area_polygon = fn_get_gdf_of_2d_area(str_hdf_file_path)
    gdf_flowpath_with_attrib = fn_extract_flowpaths_gdf(gdf_2d_area_polygon, str_input_hydrofabric_gpkg)
    gdf_2d_area_polygon_nextgen = gdf_2d_area_polygon.to_crs(gdf_flowpath_with_attrib.crs)
    
    # append the gdf_flowpath_with_attrib to only networks with outlets
    gdf_flowpath_with_attrib = fn_remove_networks_without_outlets(gdf_flowpath_with_attrib,
                                                                  str_input_hydrofabric_gpkg,
                                                                  str_hdf_file_path,
                                                                  gdf_2d_area_polygon)
    
    gdf_upstream_points = fn_create_upstream_flowpath_points(gdf_flowpath_with_attrib,
                                                             gdf_2d_area_polygon_nextgen)
    
    gdf_upstream_points,gdf_flowpath_with_attrib = fn_determine_starting_hecras_cells(str_hdf_file_path,
                                                                                      gdf_upstream_points,
                                                                                      gdf_flowpath_with_attrib)
    
    gdf_upstream_points = fn_compute_travel_distance_per_point(gdf_flowpath_with_attrib,
                                                               str_input_hydrofabric_gpkg,
                                                               gdf_upstream_points)
    
    # Travel times, peak flow, remove water bodies, create mainstems... etc...
    gdf_streams,gdf_flow_points,gdf_mainstems_revised = fn_travel_time_and_peak_flow_est(gdf_flowpath_with_attrib,
                                                                                         gdf_upstream_points,
                                                                                         str_limiting_discharge_csv)
    
    
    print('Writing output geopackage...')
    
    str_output_geopackage_name = 'model_hydrofabric.gpkg'
    str_gpkg_to_create = os.path.join(str_output_dir, str_output_geopackage_name)
    
    gdf_2d_area_polygon_nextgen = gdf_2d_area_polygon.to_crs(gdf_flowpath_with_attrib.crs)

    # Write GeoDataFrames to GeoPackage as separate layers
    gdf_mainstems_revised.to_file(str_gpkg_to_create, layer='03_flowpaths_stabilize', driver="GPKG")
    gdf_flow_points.to_file(str_gpkg_to_create, layer='02_flow_points', driver="GPKG")
    gdf_streams.to_file(str_gpkg_to_create, layer='01_stream_lines', driver="GPKG")
    gdf_2d_area_polygon_nextgen.to_file(str_gpkg_to_create, layer='00_area_2d', driver="GPKG")
        
# .........................................................

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()
    
    parser = argparse.ArgumentParser(description='========= EXTRACT NEXTGEN HYDROFABRIC FOR 2D HEC-RAS AREA =========')
    
    parser.add_argument('-c',
                        dest = "str_config_file_path",
                        help=r'REQUIRED: Global configuration filepath Example: C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini',
                        required=False,
                        default=r'C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini',
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-i',
                        dest = "str_hdf_file_path",
                        help=r'REQUIRED: path to plan HDF file (polygons) Example: E:\HECRAS_2D_12070205\BLE_12070205_Engineering_Models\Engineering Models\Hydraulic Models\RAS_Submittal\LBSG_501\Input\BLE_LBSG_501.p02.hdf',
                        required=False,
                        default=r'E:\HECRAS_2D_12070205\BLE_12070205_Engineering_Models\Engineering Models\Hydraulic Models\RAS_Submittal\LBSG_501\Input\BLE_LBSG_501.p02.hdf',
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-o',
                        dest = "str_output_dir",
                        help=r'REQUIRED: directory to write model hydrofabric geopackage Example: E:\sample_2d_output\BLE_LBSG_501_p02',
                        required=False,
                        default=r'E:\sample_2d_output\BLE_LBSG_501_p02',
                        metavar='DIR',
                        type=str)
    

    args = vars(parser.parse_args())
    
    str_config_file_path = args['str_config_file_path']
    str_hdf_file_path = args['str_hdf_file_path']
    str_output_dir = args['str_output_dir']

    fn_create_models_nextgen_gpkg(str_config_file_path,
                                  str_hdf_file_path,
                                  str_output_dir)

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)
    
    print('Compute Time: ' + str(time_pass))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~