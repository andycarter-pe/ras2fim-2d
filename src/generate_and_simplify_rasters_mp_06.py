# Script 06 - Step 05 created a flood rasters of depth and wsel for each 'run'
# over the entire 'mainstem' river.  Now we need to 'cut' these runs into rasters
# and polygons to each 'reach' in the mainstem.  The created rasters for WSEL 
# and depth are in compliance with InFRM reaster guidance.  The location of these
# rasters in addition to floodplain polygon limits are appended to the hydraulic
# model results geopackage.
#
# This requires as HEC-RAS 2D model from which the following inputs are provided:
#  (1) hydraulic model results geodatabase for a given mainstem
#  (2) detailed depth and wsel rasters for each flow along entire mainstem
#  (3) directory to write the output rasters
#
# Created by: Andy Carter, PE
# Created - 2024.07.10

# ************************************************************
import os
import geopandas as gpd
import pandas as pd
from shapely.geometry import box
import rasterio
import rioxarray as rxr
from rioxarray import merge

import numpy as np
import rasterio as rio

import tqdm
import multiprocessing as mp

import configparser
import argparse
import time
import datetime
import warnings
# ************************************************************


# ==============================
def fn_create_raster_products(int_flow_step,
                              str_path_large_dem_folder,
                              str_output_folder,
                              gdf_flowpath_limits,
                              str_output_crs,
                              flt_desired_res):
    
    # TODO - 2024.07.11 -- Need to supress the warnings
    # supress all warnings
    warnings.filterwarnings("ignore", category=UserWarning )

    # --- Create the floodplain limits polygons ---
    str_expected_floodplain_filename = '11_wsel_max_' + str(int_flow_step) + '_floodplain_ar.gpkg'
    str_expected_floodplain_filepath = os.path.join(str_path_large_dem_folder, str_expected_floodplain_filename)

    if os.path.exists(str_expected_floodplain_filepath):
        
        gdf_floodplain_ar = gpd.read_file(str_expected_floodplain_filepath)
        shp_floodplain_ar = gdf_floodplain_ar.iloc[0]['geometry']

        # Create an empty DataFrame
        data = []
        columns = ['run_name', 'flow_cfs', 'flowpath', 'clipped_geometry']
        gdf_results = pd.DataFrame(data, columns=columns)

        # Iterate through each row in gdf_flowpath_limits
        for index, row in gdf_flowpath_limits.iterrows():
            shp_flowpath_geom = row['geometry']

            # Clip shp_flowpath_geom to shp_floodplain_ar
            shp_clipped_geom = shp_flowpath_geom.intersection(shp_floodplain_ar)

            # Handle multi-part geometries
            if shp_clipped_geom.geom_type == 'MultiLineString':
                for geom in shp_clipped_geom.geoms:
                    gdf_results = pd.concat([gdf_results, pd.DataFrame([{
                        'run_name': row['run_name'],
                        'flow_cfs': row['flow_' + str(int_flow_step)],
                        'flowpath': row['flowpath'],
                        'clipped_geometry': geom
                    }])], ignore_index=True)
            else:
                gdf_results = pd.concat([gdf_results, pd.DataFrame([{
                    'run_name': row['run_name'],
                    'flow_cfs': row['flow_' + str(int_flow_step)],
                    'flowpath': row['flowpath'],
                    'clipped_geometry': shp_clipped_geom
                }])], ignore_index=True)

        gdf_results['geometry'] = gdf_results['clipped_geometry']

        # Drop the 'clipped_geometry' column
        gdf_results = gdf_results.drop(columns=['clipped_geometry'])

        # Convert DataFrame to GeoDataFrame
        gdf_results = gpd.GeoDataFrame(gdf_results, geometry='geometry')

        # Set the coordinate reference system of new geodataframe
        gdf_results.crs = gdf_floodplain_ar.crs
        
        # HARD CODED - For now, turning off depth compuations
        # ******
        b_compute_depth_raster = False
        # ******

        # ---- Now, create the depth raster for the entire run ----
        if b_compute_depth_raster:
            b_depth_raster_found = False
            
            str_expected_depth_filename = '09_wsel_max_' + str(int_flow_step) + '_depth_nudge.tif'
            str_expected_depth_filepath = os.path.join(str_path_large_dem_folder, str_expected_depth_filename)
            
            if os.path.exists(str_expected_depth_filepath):
            
                b_depth_raster_found = True
                
                #print(f'Converting Depth Raster: {str(int_flow_step)}')
                # convert the depth DEM to InFRM compliant rasters
                # Create a raster for the entire run's depth
                with rxr.open_rasterio(str_expected_depth_filepath) as xds_raw_dem:
                     # reproject the DEM to the requested CRS
                    xds_reproject = xds_raw_dem.rio.reproject(str_output_crs)
    
                    # change the depth values to integers representing 1/10th interval steps
                    # Example - a cell value of 25 = a depth of 2.5 units (feet or meters)
                    xds_reproject_scaled = ((xds_reproject * 10) + 0.5) // 1
    
                    # set the n/a cells to a value of 65535
                    xds_reproject_scaled = xds_reproject_scaled.fillna(65535)
    
                    # set the nodata value to 65535 - InFRM compliant
                    xds_reproject_scaled = xds_reproject_scaled.rio.set_nodata(65535)
    
                    # using the merge on a single raster to allow for user supplied
                    # raster resolution of the output
                    xds_depth_raster = rxr.merge.merge_arrays(xds_reproject_scaled,
                                                                     res=(flt_desired_res),
                                                                     nodata=(65535))
            else:
                print(f'Depth raster does not exist: {str_expected_depth_filepath}')
        # ---- ----

        # ---- Now, create the wsel raster for the entire run ----
        b_wsel_raster_found = False
        
        str_expected_wsel_filename = '10_wsel_max_' + str(int_flow_step) + '_wsel_nudge.tif'
        str_expected_wsel_filepath = os.path.join(str_path_large_dem_folder, str_expected_wsel_filename)
        
        if os.path.exists(str_expected_wsel_filepath):
        
            b_wsel_raster_found = True
            
            #print(f'Converting WSEL Raster: {str(int_flow_step)}')
            # convert the depth DEM to InFRM compliant rasters
            # Create a raster for the entire run's depth
            
            with rxr.open_rasterio(str_expected_wsel_filepath) as xds_raw_dem:
                 # reproject the DEM to the requested CRS
                xds_reproject = xds_raw_dem.rio.reproject(str_output_crs)

                # change the depth values to integers representing 1/10th interval steps
                xds_reproject_scaled = xds_reproject

                # set the n/a cells to a value of 65535
                xds_reproject_scaled = xds_reproject_scaled.fillna(-3.4028235e+38)

                # set the nodata value to 65535 - InFRM compliant
                xds_reproject_scaled = xds_reproject_scaled.rio.set_nodata(-3.4028235e+38)

                # using the merge on a single raster to allow for user supplied
                # raster resolution of the output
                xds_wsel_raster = rxr.merge.merge_arrays(xds_reproject_scaled,
                                                                 res=(flt_desired_res),
                                                                 nodata=(-3.4028235e+38))
        else:
            print(f'WSEL raster does not exist: {str_expected_depth_filepath}')
        # ---- ----
        
        # === Clip the depth raster to the polygons in gdf_results ===
        if b_compute_depth_raster:
            if b_depth_raster_found:
                
                
                #print('Clipping Depth Rasters')
                str_depth_folder = os.path.join(str_output_folder, '00_depth_rasters')
                
                if not os.path.exists(str_depth_folder):
                    # Create the directory
                    os.makedirs(str_depth_folder)
                
                # Reproject gdf_results to str_output_crs
                gdf_flowpath_limits_3857 = gdf_results.to_crs(str_output_crs)
    
                list_depth_rasters = []
                
                for index, row in gdf_flowpath_limits_3857.iterrows():
                    geom = row['geometry']
                    
                    # ----- create clipped depth rasters-----
                    xds_clipped_depth = xds_depth_raster.rio.clip([geom])
                    str_depth_name = 'depth_' + row['flowpath'] + '_' + str(row['flow_cfs']) + '.tif'
                    str_depth_filepath = os.path.join(str_depth_folder, str_depth_name)
                
                    # compress and write out depth raster as unsigned 16 bit integer
                    xds_clipped_depth.rio.to_raster(str_depth_filepath,compress='lzw',dtype="uint16")
                    
                    list_depth_rasters.append(str_depth_filepath)
                    
                gdf_results['depth_raster_path'] = list_depth_rasters
        # ====
        
        # ++++ Clip the wsel raster to the polygons in gdf_results ++++
        if b_wsel_raster_found:
            
            #print('Clipping WSEL Rasters')
            
            str_wsel_folder = os.path.join(str_output_folder, '01_wsel_rasters')
            
            if not os.path.exists(str_wsel_folder):
                # Create the directory
                os.makedirs(str_wsel_folder)
            
            # Reproject gdf_results to str_output_crs
            gdf_flowpath_limits_3857 = gdf_results.to_crs(str_output_crs)

            list_wsel_rasters = []
            
            for index, row in gdf_flowpath_limits_3857.iterrows():
                geom = row['geometry']
                
                if geom is None or geom.is_empty:
                    pass
                    # this segment has no geometry
                else:
                    # ----- create clpped wsel rasters -----
                    xds_clipped_wsel = xds_wsel_raster.rio.clip([geom])
                    str_wsel_name = 'wsel_' + row['flowpath'] + '_' + str(row['flow_cfs']) + '.tif'
                    str_wsel_filepath = os.path.join(str_wsel_folder, str_wsel_name)
    
                    # compress and write out wsel raster as signed 32 bit integer
                    xds_clipped_wsel.rio.to_raster(str_wsel_filepath,compress='lzw',dtype="float32")
    
                    list_wsel_rasters.append(str_wsel_filepath)
                
            gdf_results['wsel_raster_path'] = list_wsel_rasters
        # ++++
        
        # Append the gdf_results for this given flow to an overall gdf
        return(gdf_results)

    else:
        print(f'Floodplain geopackage does not exist: {str_expected_floodplain_filepath}')
# ==============================


# -------------
def fn_process_flow_step(dict_of_parameters):
    
    # function used for the multiprocessing -- returns a gdf
    
    return fn_create_raster_products(dict_of_parameters['int_flow_step'],
                                     dict_of_parameters['str_path_large_dem_folder'],
                                     dict_of_parameters['str_output_folder'],
                                     dict_of_parameters['gdf_flowpath_limits'],
                                     dict_of_parameters['str_output_crs'],
                                     dict_of_parameters['flt_desired_res'])
# -------------


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


# ++++++++++++++++++++++++++++
def fn_extract_and_simplify_rasters(str_config_file_path,
                                    str_flood_limits_gpkg,
                                    str_path_large_dem_folder,
                                    str_output_dir,
                                    b_print_output):
    
    # supress all warnings
    warnings.filterwarnings("ignore", category=UserWarning )
    
    if b_print_output:
        print(" ")
        print("+=================================================================+")
        print("|        EXTRACT AND SIMPLIFY RASTERS FROM MAINSTEM RUNS          |")
        print("|                Created by Andy Carter, PE of                    |")
        print("|             Center for Water and the Environment                |")
        print("|                 University of Texas at Austin                   |")
        print("+-----------------------------------------------------------------+")
    
        
        print("  ---(c) INPUT GLOBAL CONFIGUTATION FILE: " + str_config_file_path)
        print("  ---(i) INPUT HYDRAULIC RESULTS (GPKG): " + str_flood_limits_gpkg)
        print("  ---(d) INPUT FOLDER OF MAINSTEM RASTERS: " + str_path_large_dem_folder)
        print("  ---(o) OUTPUT FOLDER: " + str_output_dir)
        print("  ---[r] PRINT OUTPUT: " + str(b_print_output))
        print("===================================================================")

    
    # --- Read variables from config.ini ---
    # Create a ConfigParser object
    config = configparser.ConfigParser()
    
    # Read the config.ini file
    config.read(str_config_file_path)
    
    # Get the variables in the [02_copy_geom] section
    if '06_output_products' in config:
        section = config['06_output_products']
        
        # Read variables in the section
        str_output_crs = section.get('str_output_crs', '')
        flt_desired_res = section.getfloat('flt_desired_res', 0)
    else:
        print("[06_output_products] section not found in the config file.")
        
    
    # -------
    
    
    # Load the geopackage to get flooded limits of flowpath
    gdf_flowpath_limits = gpd.read_file(str_flood_limits_gpkg, layer='01_flowpath_flooded_cells_ar')
    
    # reproject gdf_flowpath_limits to str_output_crs == 'EPSG:3857'
    gdf_flowpath_limits_reproj = gdf_flowpath_limits.to_crs(str_output_crs)
    
    # List of columns that start with 'flow_'
    list_flow_columns = [col for col in gdf_flowpath_limits_reproj.columns if col.startswith('flow_')]
    
    # Prepare list of dictionaries for parameters
    list_of_dict = [{'int_flow_step': i,
                     'str_path_large_dem_folder': str_path_large_dem_folder,
                     'str_output_folder': str_output_dir,
                     'gdf_flowpath_limits': gdf_flowpath_limits_reproj,
                     'str_output_crs': str_output_crs,
                     'flt_desired_res': flt_desired_res}
                    for i in range(1, len(list_flow_columns) + 1)]
    
    # **********
    # HARED CODED MAX NUMBER OF PROCESSORS
    int_max_processors = 16
    # **********
    
    num_processors = (mp.cpu_count() - 1)
    
    if num_processors > int_max_processors:
        num_processors = int_max_processors
    
    # --- multiprocessing funtion ---
    l = len(list_flow_columns)
    p = mp.Pool(processes = num_processors)
    
    if b_print_output:
        # run tqdm on the multiprocessing step
        list_returned_gdf = list(tqdm.tqdm(p.imap(fn_process_flow_step, list_of_dict),
                                            total = l,
                                            desc='Create rasters',
                                            bar_format = "{desc}:({n_fmt}/{total_fmt})|{bar}| {percentage:.1f}%",
                                            ncols=65))
    else:
        # run the multiprocessing step without the tqdm status bar
        list_returned_gdf = list(p.imap(fn_process_flow_step, list_of_dict))
            
    p.close()
    p.join()
    
    # TODO - merge the geodataframes and append / create table in str_flood_limits_gpkg
# ++++++++++++++++++++++++++++

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()
    
    parser = argparse.ArgumentParser(description='======== EXTRACT AND SIMPLIFY RASTERS FROM MAINSTEM RUNS ========')
    
    parser.add_argument('-c',
                        dest = "str_config_file_path",
                        help=r'REQUIRED: Global configuration filepath Example: C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini',
                        required=True,
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-i',
                        dest = "str_flood_limits_gpkg",
                        help=r'REQUIRED: hydraulic results geopackage (GPKG) Example: E:\sample_2d_output\BLE_LBSG_501_p02\hydrualic_results_1884413_wb-2410249_wb-2410261_30-hr_500-cfs_to_3000-cfs_step_600-cfs.gpkg',
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-d',
                        dest = "str_path_large_dem_folder",
                        help=r'REQUIRED: path to folder containing mainstem WSEL and Depth grids (TIF) Example: E:\ras2fim_test_20250107\processing\05_large_terrain',
                        required=True,
                        metavar='DIR',
                        type=str)
    
    parser.add_argument('-o',
                        dest = "str_output_dir",
                        help=r'REQUIRED: directory to write simplified rasters Example: E:\ras2fim_test_20250105\processing\06_simple_rasters',
                        required=False,
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
    str_flood_limits_gpkg = args['str_flood_limits_gpkg']
    str_path_large_dem_folder = args['str_path_large_dem_folder']
    str_output_dir = args['str_output_dir']
    b_print_output = args['b_print_output']

    fn_extract_and_simplify_rasters(str_config_file_path,
                                    str_flood_limits_gpkg,
                                    str_path_large_dem_folder,
                                    str_output_dir,
                                    b_print_output)

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)
    
    print('Compute Time: ' + str(time_pass))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~