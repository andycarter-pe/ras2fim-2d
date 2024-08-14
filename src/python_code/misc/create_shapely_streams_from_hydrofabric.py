# Script XX - Using the provided NextGEN hydrofabric, create a shapely object
# for the stream path from the starting node to the ending node.  Note:  the
# returned shapely object is in EPSG:5070
#
# Created by: Andy Carter, PE
# Created - 2024.05.20

# ************************************************************

import argparse
import time
import datetime

import networkx as nx
import geopandas as gpd
import pandas as pd
import os

from shapely.geometry import MultiLineString, LineString
from shapely.ops import linemerge
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


# ------------------
def fn_create_travel_path_shp(str_nextgen_gpkg, str_start_node, str_end_node):

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

        print('Stream path found: ' + str(len(gdf_streams_travel_path)) + ' segments')
        
        return shp_stream_travelpath
    except:
        print(f'Stream path between {str_start_node} and {str_end_node} not found')
        return(None)
# ------------------


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()
    
    parser = argparse.ArgumentParser(description='========== CREATE SHAPELY BETWEEN STREAMS IN HYDROFABRIC =========')
    
    parser.add_argument('-i',
                        dest = "str_nextgen_gpkg_filepath",
                        help=r'REQUIRED: NextGEN hydrofabric geopackage Example: E:\ras2fim-2d\nextgen-test\nextgen_12.gpkg',
                        required=False,
                        default=r'E:\ras2fim-2d\nextgen-test\nextgen_12.gpkg',
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-s',
                        dest = "str_start_node",
                        help=r"REQUIRED: Starting stream 'id' in hydrofabric Example:wb-2410249",
                        required=False,
                        default="wb-2410249",
                        metavar='STRING',
                        type=str)
    
    parser.add_argument('-e',
                        dest = "str_end_node",
                        help=r"REQUIRED: Ending stream 'id' in hydrofabric Example:wb-2410262",
                        required=False,
                        default="wb-2410262",
                        metavar='STRING',
                        type=str)
    
    args = vars(parser.parse_args())
    
    str_nextgen_gpkg_filepath = args['str_nextgen_gpkg_filepath']
    str_start_node = args['str_start_node']
    str_end_node = args['str_end_node']
    
    shp_stream_travelpath = fn_create_travel_path_shp (str_nextgen_gpkg_filepath, str_start_node, str_end_node)

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)
    
    print('Compute Time: ' + str(time_pass))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~