# Script 07 - From the geoTIFFs created in step 6, convert these into a raster
# stack that is a netCDF file.  This netCDF will have the terrain and all the 
# water surface elevations pixel aligned.
#
# NOTE: this uses the 'nccopy' installed on the machine to compress the netcdf
#
# Created by: Andy Carter, PE
# Created - 2025.01.08

# ************************************************************
import os
import geopandas as gpd
import pandas as pd

import multiprocessing as mp
import tqdm
from functools import partial

import re
import fnmatch 
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import calculate_default_transform, reproject, Resampling
import numpy as np
import rioxarray
import xarray as xr
from datetime import datetime
from shapely.geometry import box
import subprocess

import configparser
import argparse
import time
import warnings
# ************************************************************


#----------------
def fn_fetch_q_upper_limit(gdf_model_stream_lines, str_input_id):
    # Assuming gdf_model_stream_lines is your GeoDataFrame and str_input_id is the input id to search for
    filtered_row = gdf_model_stream_lines[gdf_model_stream_lines['id'] == str_input_id]

    if not filtered_row.empty:
        # Get 'q_upper_limit' as integer
        int_q_upper_limit = int(filtered_row['q_upper_limit'].iloc[0])
    else:
        int_q_upper_limit = 1000000000
        
    return(int_q_upper_limit)
#----------------


# ...................................
def fn_generate_terrain_raster(str_ref_filepath, str_terrain_dem_path, str_type, str_id):
    """
    Generates a terrain raster clipped and aligned to a reference binary raster.

    Parameters:
    -----------
    str_ref_filepath : str
        Path to the binary raster that defines the bounding box for clipping the terrain raster.

    str_terrain_dem_path : str
        Path to the terrain Digital Elevation Model (DEM) used in the HEC-RAS model.

    str_type : str
        The type of data being processed, typically representing a specific attribute (e.g., water surface elevation).
        This string is used in the output filename.
        Example: 'wsel'

    str_id : str
        Unique identifier for the stream reach. This identifier will be used in the output filename.
        Example: 'wb-123456'

    Output File:
    ------------
    The output filename format is `{str_type}_{str_id}_0_align.tif`
    Example: 'wsel_wb-123456_0_align.tif'
    """
    
    # Load the terrain DEM and binary raster using context management
    with rioxarray.open_rasterio(str_terrain_dem_path) as xr_terrain_raster, \
         rioxarray.open_rasterio(str_ref_filepath) as xr_binary_raster:

        # Get bounds of the binary raster
        tup_bounds = xr_binary_raster.rio.bounds()

        # Create a geometry box for clipping
        geometry = box(tup_bounds[0], tup_bounds[1], tup_bounds[2], tup_bounds[3])

        # Clip the terrain to the limits of the binary raster
        xr_clipped_raster = xr_terrain_raster.rio.clip([geometry], crs=xr_binary_raster.rio.crs)

        # Identify the no_data value from the terrain raster
        terrain_nodata = xr_terrain_raster.rio.nodata

        # Replace the no_data values with np.nan and ensure the data type is float32
        xr_clipped_raster = xr_clipped_raster.where(xr_clipped_raster != terrain_nodata, np.nan).astype('float32')

        # Get the transformation and shape of the binary raster
        binary_transform = xr_binary_raster.rio.transform()
        tup_binary_shape = (xr_binary_raster.shape[1], xr_binary_raster.shape[2])  # Height, Width

        # Prepare an empty array for the aligned raster
        arr_aligned = np.empty((xr_clipped_raster.shape[0], *tup_binary_shape), dtype='float32')

        # Reproject the clipped raster to align with the binary raster
        for i in range(xr_clipped_raster.shape[0]):  # Iterate over each band
            reproject(
                source=xr_clipped_raster.data[i],
                destination=arr_aligned[i],
                src_transform=xr_clipped_raster.rio.transform(),
                src_crs=xr_clipped_raster.rio.crs,
                dst_transform=binary_transform,
                dst_crs=xr_binary_raster.rio.crs,
                resampling=Resampling.nearest
            )

    # Convert aligned array to xarray DataArray
    xr_aligned_raster = xr.DataArray(
        arr_aligned, 
        coords=xr_binary_raster.coords, 
        dims=xr_binary_raster.dims  # This retains the original dimensions
    )

    # Define the parent directory and build the output filename
    str_parent_directory = os.path.dirname(str_ref_filepath)
    str_filename = f"{str_type}_{str_id}_0_align.tif"
    str_filename = f"terrain_{str_id}_0_align.tif"
    str_outpath = os.path.join(str_parent_directory, str_filename)

    # Save the aligned raster to a new file with LZW compression
    xr_aligned_raster.rio.to_raster(str_outpath, driver="GTiff", compress='LZW', nodata=np.nan)

    #print(f"Clipped terrain raster saved to: {str_outpath}")
    
    return str_outpath
# ...................................


# --------------
def fn_replace_nodata_with_negative_infinity(str_raster_path):
    """
    Replace NaN values in the raster with -3.4028234663852886e+38 and save over the original raster.

    Parameters:
        str_raster_path (str): The path to the raster file to modify.
    """
    # NoData value to use
    nodata_value = -3.4028234663852886e+38

    # Open the original raster file
    with rasterio.open(str_raster_path) as src:
        # Read the data
        data = src.read(1)  # Read the first band

        # Replace NaN values with the specified NoData value
        data[np.isnan(data)] = nodata_value  # Replace NaNs with NoData value

        # Update metadata
        meta = src.meta.copy()
        meta.update({
            'dtype': 'float32',  # Ensure the data type is float32
            'nodata': nodata_value,  # Set the new NoData value
        })

    # Save the modified data back to the original raster file
    with rasterio.open(str_raster_path, 'w', **meta) as dst:
        dst.write(data, 1)  # Write the modified data to the first band
# --------------


# ---------------------
def fn_align_rasters_to_ref(list_grouped_filepaths, str_ref_filepath, no_data_value):
    list_align_rasters = []
    
    for filepath in list_grouped_filepaths:
        dir_path, filename = os.path.split(filepath)
        dir_path += "_align"
        os.makedirs(dir_path, exist_ok=True)  # Create the output directory if it doesn't exist
        str_output_aligned = os.path.join(dir_path, filename[:-4] + "_align.tif")
        
        # Open the reference raster
        with rasterio.open(str_ref_filepath) as ref:
            ref_meta = ref.meta.copy()
        
        # Open the input raster
        with rasterio.open(filepath) as src:
            # Update the metadata to match the reference
            ref_meta.update({
                'height': ref.height,
                'width': ref.width,
                'crs': ref.crs,
                'transform': ref.transform,
                'nodata': no_data_value,
                'compress': 'lzw'  # Change this to your desired compression method
            })
            
            # Perform the reproject
            with rasterio.open(str_output_aligned, 'w', **ref_meta) as dst:
                for i in range(1, src.count + 1):
                    reproject(
                        source=rasterio.band(src, i),
                        destination=rasterio.band(dst, i),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=ref_meta['transform'],
                        dst_crs=ref_meta['crs'],
                        resampling=Resampling.bilinear,
                        nodata=no_data_value
                    )
        list_align_rasters.append(str_output_aligned)
    return(list_align_rasters)
# ---------------------


# -----------------
def fn_compress_netcdf_nccopy(str_netcdf_filepath):
    
    str_output_netcdf_filepath = str_netcdf_filepath[:-7] + ".nc"
    
    # Construct the nccopy command
    # Note:  uses the nccopy command which is installed in this environment
    command = ['nccopy', '-d', '5', str_netcdf_filepath, str_output_netcdf_filepath]
    
    # Call the nccopy command using subprocess
    try:
        subprocess.run(command, check=True)
        #print(f"Compressed file created: {str_output_netcdf_filepath}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
        
    # Delete the uncompressed wsel
    if os.path.exists(str_netcdf_filepath):
        os.remove(str_netcdf_filepath)
# -----------------


# ----------------
# Function to extract the integer from the filename
def fn_extract_integer(file_path):
    # Use a regular expression to find the integer in the filename
    match = re.search(r'_(\d+)_align\.tif$', file_path)
    if match:
        return int(match.group(1))  # Return the extracted integer
    return None  # Return None if no integer is found
# ----------------


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def fn_create_netcdf(list_align_rasters, str_type, str_id, str_input_dir, str_terrain_filepath):
    
    # *******
    str_flow_units = 'cfs' # hard coded
    str_raster_vertical_units = 'feet' # hard coded
    # *******
    
    # Sort raster files based on the extracted integers
    sorted_tif_files = sorted(list_align_rasters, key=fn_extract_integer)
    
    # Parse out the flow values
    list_file_names = [os.path.basename(file_path) for file_path in sorted_tif_files]
    
    # Regular expression to extract the flow values
    list_flows_str = [re.search(r'_(\d+)_align', filename).group(1) for filename in list_file_names]
    
    # Convert the extracted flow values to integers
    list_flows_int = list(map(int, list_flows_str))
    
    # Initialize a list to hold the datasets
    da_list = []
    
    # Loop through each GeoTIFF and load it
    for tiff_file in sorted_tif_files:
        # Load the GeoTIFF as an xarray DataArray
        da = rioxarray.open_rasterio(tiff_file)

        # Assign a name to the DataArray (e.g., 'depth')
        da.name = str_type  # or use any other name that reflects the data

        # Use squeeze to remove unnecessary dimensions (like the 'band' dimension if it's 1)
        da = da.squeeze()  # Ensures the DataArray stays 3D (y, x) instead of 4D
        
        # Expand dims for 'flow' and prepare for concatenation
        da_list.append(da.expand_dims(dim='flow'))
        
    # Load the terrain raster as a DataArray
    terrain_da = rioxarray.open_rasterio(str_terrain_filepath)
    terrain_da = terrain_da.squeeze()
    terrain_da = terrain_da.rename("terrain")
        
    # Concatenate all datasets along the 'flow' dimension and assign flow values as coordinates
    da_xr = xr.concat(da_list, dim='flow')
    
    # Add the terrain data to the DataArray
    da_xr['terrain'] = terrain_da
    
    da_xr = da_xr.assign_coords(flow=("flow", list_flows_int))
    
    da_xr.attrs['units'] = str_raster_vertical_units 
    da_xr.coords['flow'].attrs['units'] = str_flow_units
    da_xr.coords['terrain'].attrs['units'] = 'feet'
    
    # output folder is to be one folder up from list_align_rasters[0]
    str_parent_dir = os.path.dirname(list_align_rasters[0])
    str_parent_dir = os.path.dirname(str_parent_dir)
    str_netcdf_folder = os.path.join(str_parent_dir, '02_' + str_type + "_nc")
    
    # Create the netcdf folder
    os.makedirs(str_netcdf_folder, exist_ok=True)
    str_netcdf_filename = f"{str_type}_{str_id}_big.nc"
    
    # location to store the netcdf file
    str_netcdf_filepath = os.path.join(str_netcdf_folder, str_netcdf_filename)
    
    # Save the combined dataset as a NetCDF file
    da_xr.to_netcdf(str_netcdf_filepath)
    
    # ----------
    # Add some global attributes to the file
    # Load the existing NetCDF file
    ds = xr.load_dataset(str_netcdf_filepath)

    # Add a global attribute
    ds.attrs['00_stream_id'] = str_id
    ds.attrs['01_type'] = str_type

    ds.attrs['02_description'] = 'RAS2FIM-2D Raster Stack'
    ds.attrs['03_author'] = 'Created by Andy Carter, PE'
    ds.attrs['04_organization'] = 'Center for Water and the Environment'
    ds.attrs['05_institution'] = 'University of Texas at Austin'
    ds.attrs['06_history'] = f"Created on {datetime.now().strftime('%Y-%m-%d %H:%M')}"

    ds.attrs['07_raster_folder'] = str_input_dir

    # Save the updated dataset back to the NetCDF file
    ds.to_netcdf(str_netcdf_filepath)
    # ----------

    #Compress the original netCDF file -- "has a _big.nc" suffix
    fn_compress_netcdf_nccopy(str_netcdf_filepath)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ''''''''''''''''''''''''''''''''''''''
def fn_align_raster_groups(list_grouped_filepaths, str_type, str_id, str_input_dir, str_terrain_dem_path):
    #print(f"type: {str_type}, ID: {str_id} , Len: {len(list_grouped_filepaths)}")
    
    # Initialize variables to find the largest raster dimensions
    max_left = float('inf')
    max_right = float('-inf')
    max_top = float('-inf')
    max_bottom = float('inf')

    rasters_info = []

    # Initialize variables to find the largest raster dimensions
    max_width = 0
    max_height = 0
    raster_data = []

    # List to store raster resolutions
    list_raster_resolutions = []

    # Loop through the raster files to determine the max dimensions and load the data
    for filepath in list_grouped_filepaths:
        try:
            with rasterio.open(filepath) as src:
                transform = src.transform
                width = src.width
                height = src.height
                max_width = max(max_width, width)
                max_height = max(max_height, height)
                bounds = src.bounds
                
                crs = src.crs # this will grab the crs 
                no_data_value = src.nodata

                # Read the raster data and append it to the raster_data list
                data = src.read(1)  # Read the first band
                raster_data.append(data)

                # Update the overall bounds
                max_left = min(max_left, bounds.left)
                max_right = max(max_right, bounds.right)
                max_top = max(max_top, bounds.top)
                max_bottom = min(max_bottom, bounds.bottom)

                # Store raster info for later processing
                rasters_info.append((filepath, width, height, transform))

                # Get the resolution and add to the list (pixel size in x and y directions)
                x_res = transform[0]  # Pixel width
                y_res = -transform[4]  # Pixel height (negative because y-coordinates decrease upwards)
                list_raster_resolutions.append((x_res, y_res))
        except Exception as e:
            print(f"Error reading {filepath}: {e}")

    num_rasters = len(raster_data)

    #print(f"{num_rasters} x {max_width} x {max_height}")
    #print(f"Max Left: {max_left}, Max Right: {max_right}, Max Top: {max_top}, Max Bottom: {max_bottom}")

    list_grouped_filepaths[0]

    # Open the raster file and get the data type
    with rasterio.open(list_grouped_filepaths[0]) as src:
        raster_dtype = src.dtypes[0]  # Get the dtype of the first band

    #print(f"data type: {raster_dtype}")
    
    # get the parent folder of the rasters
    # Get the parent directory of str_input_dir
    str_parent_dir = os.path.dirname(list_grouped_filepaths[0])
    str_parent_dir += '_align' # hard coding the folder suffix
    
    os.makedirs(str_parent_dir, exist_ok=True)
    
    str_ref_raster_filename = f"{str_type}_{str_id}_bbox_ref.tif"
    str_ref_filepath = os.path.join(str_parent_dir, str_ref_raster_filename)
    
    # Create an example array (random values or specific data)
    data = np.random.randint(0, 1, size=(height, width), dtype=np.uint8)
    
    # Define the GeoTransform parameters using from_origin
    transform = from_origin(
        max_left,  # top-left x
        max_top,    # top-left y
        list_raster_resolutions[0][0],          # pixel size in the x-direction
        list_raster_resolutions[0][0]           # pixel size in the y-direction (absolute value)
    )

    # Define the metadata for the GeoTIFF
    metadata = {
        'driver': 'GTiff',
        'count': 1,
        'dtype': raster_dtype,
        'width': max_width,
        'height': max_height,
        'crs': crs,
        'transform': transform,
        'compress': 'lzw',
    }
    
    # Create the GeoTIFF
    with rasterio.open(str_ref_filepath , 'w', **metadata) as dst:
        dst.write(data, 1)

    #print(f"Limits reference GeoTIFF created at: {str_ref_filepath}")
    #print(crs)
    
    # From str_ref_filepath, get a terrain raster and add to wsel_raster_align
    # as a wsel with a flow of zero
    str_terrain_filepath = fn_generate_terrain_raster(str_ref_filepath, str_terrain_dem_path, str_type, str_id)
    
    # change the terrain raster NaN to negative infinity
    fn_replace_nodata_with_negative_infinity(str_terrain_filepath)
    
    # Transform (padding) all the rasters in list_grouped_filepaths to limits of str_ref_filepath
    list_align_rasters = fn_align_rasters_to_ref(list_grouped_filepaths, str_ref_filepath, no_data_value)
    
    # Delete the reference raster (str_ref_filepath)
    #if os.path.exists(str_ref_filepath):
    #    os.remove(str_ref_filepath)
    
    # from list_align_rasters create a netCDF of the aligned raster data
    fn_create_netcdf(list_align_rasters,str_type, str_id, str_input_dir, str_terrain_filepath)
    
    #print('----------------------')
# ''''''''''''''''''''''''''''''''''''''


# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist" % arg)
    else:
        # File exists so return the directory
        return arg
        return open(arg, 'r')  # return an open file handle
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


# -------------------------
# Function to process a single ID, with extra arguments
def fn_process_single_id(str_id, df_valid_rows, str_type, str_simple_raster_path, str_terrain_dem_path):
    list_full_path = []
    
    # Get list of filenames where 'id' == str_id from df_valid_rows
    list_filenames = df_valid_rows[df_valid_rows['id'] == str_id]['filename'].tolist()
    
    # Generate the full paths
    for filename in list_filenames:
        str_full_path = os.path.join(str_simple_raster_path, filename)
        list_full_path.append(str_full_path)
    
    # Call the processing function
    fn_align_raster_groups(list_full_path, str_type, str_id, str_simple_raster_path, str_terrain_dem_path)
# -------------------------


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
def fn_prepare_netcdf_files(str_config_file_path,
                            str_terrain_filename,
                            str_output_dir,
                            b_print_output):
    
    # supress all warnings
    warnings.filterwarnings("ignore", category=UserWarning )
    
    print(" ")
    if b_print_output:
        print("+=================================================================+")
        print("|               PREPARE NETCDF FROM WSEL RASTERS                  |")
        print("|                Created by Andy Carter, PE of                    |")
        print("|             Center for Water and the Environment                |")
        print("|                 University of Texas at Austin                   |")
        print("+-----------------------------------------------------------------+")
    
        
        print("  ---(c) INPUT GLOBAL CONFIGUTATION FILE: " + str_config_file_path)
        print("  ---(t) TERRAIN FILENAME (TIF): " + str_terrain_filename)
        print("  ---(o) OUTPUT FOLDER: " + str_output_dir)
        print("  ---[r] PRINT OUTPUT: " + str(b_print_output))
        print("===================================================================")
    else:
        print('Step 7: Prepare NetCDF form WSEL rasters')
    
    
    str_hydraulic_results_gpkg_path = os.path.join(str_output_dir, '04_wsel_and_hr')
    
    # To determine the peak flow for a given stream reach
    str_model_hydrofabric = os.path.join(str_output_dir, '01_model_hydrofabric', 'model_hydrofabric.gpkg')
    gdf_model_stream_lines = gpd.read_file(str_model_hydrofabric, layer='01_stream_lines')
    
    str_simple_raster_path = os.path.join(str_output_dir, '06_simple_rasters', '01_wsel_rasters')
    
    str_terrain_dem_path = os.path.join(str_output_dir, '02_model_copies', 'source_terrain', str_terrain_filename)
    
    str_type = 'wsel'
    
    # --- Read variables from config.ini ---
    # Create a ConfigParser object
    config = configparser.ConfigParser()
    
    # Read the config.ini file
    config.read(str_config_file_path)
    
    # Get the variables in the [07_prep_netcdf] section
    if '07_prep_netcdf' in config:
        section = config['07_prep_netcdf']
    
        # Read variables in the section
        flt_allowed_varianvce = section.getfloat('flt_allowed_varianvce', 0)
    else:
        print("[07_prep_netcdf] section not found in the config file.")
    # ---
    
    int_grids = 0
    
    # Initialize an empty pandas DataFrame
    df_valid_rows = pd.DataFrame(columns=['mainstem', 'id', 'run_name', 'flow', 'perct_wet', 'perct_stable', 'filename'])
    
    # Loop over each file in the directory
    for filename in os.listdir(str_hydraulic_results_gpkg_path):
        if filename.endswith('.gpkg'):
            # Construct the full file path
            filepath = os.path.join(str_hydraulic_results_gpkg_path, filename)
            
            # Load the '03_streams_ln' layer
            gdf = gpd.read_file(filepath, layer='03_streams_ln')
    
            # Count fields starting with "flow"
            flow_fields = [col for col in gdf.columns if col.startswith('flow')]
            int_num_flow_fields = len(flow_fields)
            
            for index, row in gdf.iterrows():
                #print('--------')
                #print(f"{row['mainstem']} : {row['id']} : {row['run_name']}")
                
                # determine the peak flow of row['id']
                int_q_upper_limit = fn_fetch_q_upper_limit(gdf_model_stream_lines, row['id'])
                #print(f'Stream upper flow: {int_q_upper_limit}')
    
                for i in range(1,int_num_flow_fields + 1):
                    str_flow_field = 'flow_' + str(i)
                    str_perct_wet_field = 'perct_wet_' + str(i)
                    str_perct_stable_field = 'perct_stable_' + str(i)
                    
                    # Determine if this result is a backwater from a downstream reach
                    if row[str_flow_field] <=  int_q_upper_limit:
                        #print(f"    {row[str_flow_field]} : {row[str_perct_wet_field]} : {row[str_perct_stable_field]}")
                        
                        try:
                            flt_perct_wet = float(row[str_perct_wet_field])
                            flt_perct_stable = float(row[str_perct_stable_field])
                            flt_delta = flt_perct_wet - flt_perct_stable
                            if flt_delta <= flt_allowed_varianvce:
                                int_grids += 1
                                
                                str_expected_filename = 'wsel_' + row['id'] + '_' + str(row[str_flow_field]) + '.tif'
                                
                                # add row to new/existing dataframe
                                df_valid_rows = pd.concat([
                                    df_valid_rows,
                                    pd.DataFrame([{
                                        'mainstem': row['mainstem'],
                                        'id': row['id'],
                                        'run_name': row['run_name'],
                                        'flow': row[str_flow_field],
                                        'perct_wet': flt_perct_wet,
                                        'perct_stable': flt_perct_stable,
                                        'filename': str_expected_filename
                                    }])
                                ], ignore_index=True)
                            
                        except:
                            # Could not compute... 'None value likely'
                            pass
                    else:
                        #print(f"    **BACKWATER: {row[str_flow_field]} : {row[str_perct_wet_field]} : {row[str_perct_stable_field]}")
                        pass
                    
    print(f'   -- Valid Rasters Expected: {int_grids}')
    
    # Check if each filename exists in the directory
    df_valid_rows['file_exists'] = df_valid_rows['filename'].apply(
        lambda x: os.path.isfile(os.path.join(str_simple_raster_path, x))
    )
    
    # Count rows where file_exists is False
    int_count_found = df_valid_rows['file_exists'].value_counts().get(True, 0)
    
    print(f"   -- Matching Rasters Found: {int_count_found} of {int_grids}")
    
    # filter to only rows where 'file_exists' == True
    df_valid_rows = df_valid_rows[df_valid_rows['file_exists'] == True]
    
    # Get a list of unique 'id' values
    list_unique_ids = df_valid_rows['id'].unique().tolist()
    
    # **********
    # HARD-CODED MAX NUMBER OF PROCESSORS
    int_max_processors = 16
    # **********

    # Determine the number of processors to use
    num_processors = mp.cpu_count() - 1
    if num_processors > int_max_processors:
        num_processors = int_max_processors

    # Create a partial function to include additional arguments
    partial_process_single_id = partial(
        fn_process_single_id,
        df_valid_rows=df_valid_rows,
        str_type=str_type,
        str_simple_raster_path=str_simple_raster_path,
        str_terrain_dem_path=str_terrain_dem_path
    )

    # Create a multiprocessing pool
    with mp.Pool(processes=num_processors) as pool:
        # Use tqdm to show progress
        list(tqdm.tqdm(
            pool.imap(partial_process_single_id, list_unique_ids),
            total=len(list_unique_ids),
            desc='   -- Create netCDFs',
            bar_format="{desc}:({n_fmt}/{total_fmt})|{bar}| {percentage:.1f}%",
            ncols=75
        ))
        
    print("+-----------------------------------------------------------------+")
# ++++++++++++++++++++++++++++


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()
    
    parser = argparse.ArgumentParser(description='======== PREPARE NETCDF FROM WSEL RASTERS ========')
    
    parser.add_argument('-c',
                        dest = "str_config_file_path",
                        help=r'REQUIRED: Global configuration filepath Example: C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini',
                        required=True,
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-t',
                        dest = "str_terrain_filename",
                        help=r'REQUIRED: terrain raster filename in 02_model_copies\source_terrain Example: 120702050405.tif',
                        metavar='FILE',
                        type=str)
    
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
    str_terrain_filename = args['str_terrain_filename']
    str_output_dir = args['str_output_dir']
    b_print_output = args['b_print_output']

    fn_prepare_netcdf_files(str_config_file_path,
                            str_terrain_filename,
                            str_output_dir,
                            b_print_output)

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)
    
    print('Compute Time: ' + str(time_pass))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
