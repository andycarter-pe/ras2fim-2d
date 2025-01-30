# This is the main orchestration script for "ras2fim-2d".  It attempts to convert
# 2D HEC-RAS models into a set of flood inundation mapping (fim)
# library of rasters (netCDF) for a corresponding National Water Model (NWM) 
# NextGEN stream segments within the HEC-RAS model's 2D area.
#
# This is built to run on a Windows machine and requires that HEC-RAS v6.1
# be installed prior to execution.  Additionally, 'nccopy' needs to be
# installed on the host machine
#
# Created by: Andy Carter, PE
# 2025.01.08

# ************************************************************
import argparse
import time
import datetime
import warnings
import os
import shutil
import rioxarray as rxr
import multiprocessing as mp
from functools import partial
import tqdm

# Import modules
from copy_source_folder_00 import fn_copy_folder
from create_model_hydrofabric_gpkg_01 import fn_create_models_nextgen_gpkg
from spawn_hecras_copies_02 import fn_spawn_hecras_copies
##from run_hecras_windows_03 import fn_run_hec_ras_models
from compute_max_wsel_and_stable_hr_04 import fn_compute_max_wsel_and_stable_hr
from generate_depth_wsel_rasters_mp_05 import fn_generate_depth_wsel_rasters
from generate_and_simplify_rasters_mp_06 import fn_extract_and_simplify_rasters
from prep_for_netcdf_mp_07 import fn_prepare_netcdf_files
from filter_netcdf_wsel_mp_08 import fn_filter_netcdf
# ************************************************************
# Additional comments if needed

# ************************************************************

# -------------
def fn_parse_tuple(value):
    # Remove any spaces and parentheses, then split by commas
    return tuple(map(int, value.strip('()').split(',')))
# -------------


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


# ..................
def fn_copy_ras_hdf_plans(str_model_path_02, str_hdf_path_03):
    # Walk through the source directory
    for root, dirs, files in os.walk(str_model_path_02):
        for file in files:
            if file.endswith(".p01.hdf"):
                # Full path of the source file
                source_file_path = os.path.join(root, file)
                
                # Copy the file to the destination directory
                shutil.copy(source_file_path, str_hdf_path_03)
# ..................        


# -------------------
def fn_determine_tif_in_dir_with_largest_bbox(str_terrain_folder):
    
    tif_files = []
    
    # Walk through the directory to find all .tif files
    for root, dirs, files in os.walk(str_terrain_folder):
        for file in files:
            if file.endswith('.tif'):
                tif_files.append(os.path.join(root, file))

    largest_area = 0
    largest_tif = None

    # Process each GeoTIFF
    for tif_file in tif_files:
        # Open the GeoTIFF
        data = rxr.open_rasterio(tif_file, masked=True)

        # Get the bounding box (minx, miny, maxx, maxy)
        bounds = data.rio.bounds()
        minx, miny, maxx, maxy = bounds

        # Calculate the area of the bounding box
        area = (maxx - minx) * (maxy - miny)

        # Determine the largest bounding box
        if area > largest_area:
            largest_area = area
            largest_tif = tif_file

    # Output the result
    return (largest_tif)
# -------------------

# ------------
def fn_process_hdf(args):
    str_config_file_path, hdf_filepath, hydrofabric_path, output_path, b_print_output = args
    fn_compute_max_wsel_and_stable_hr(str_config_file_path, hdf_filepath, hydrofabric_path, output_path, b_print_output)
# ------------

# +++++++++++++++++++++++++++++
def fn_ras2fim_2d(str_source_folder, 
                  tup_hecras_files_to_copy, 
                  str_config_file_path, 
                  str_output_dir,
                  tup_step_range):
    
    
    # To supress the printing of Step output text
    b_print_output = False
    
    # supress all warnings
    warnings.filterwarnings("ignore", category=UserWarning )
    
    print(" ")
    print("+=================================================================+")
    print("|                       RUN RAS2FIM-2D                            |")
    print("|                Created by Andy Carter, PE of                    |")
    print("|             Center for Water and the Environment                |")
    print("|                 University of Texas at Austin                   |")
    print("+-----------------------------------------------------------------+")

    print("  ---(i) SOURCE HEC-RAS FOLDER: " + str_source_folder)
    print("  ---(o) OUTPUT FOLDER: " + str_output_dir)
    print("  ---(c) INPUT GLOBAL CONFIGUTATION FILE: " + str_config_file_path)
    print("  ---(f) FILES TO COPY (geom, unsteady flow, plan): " + ', '.join(map(str, tup_hecras_files_to_copy)))
    print("  ---[s] Optional: Steps to run (start, stop): "+ ', '.join(map(str, tup_step_range)))
    print("===================================================================")
    
    int_start_step = tup_step_range[0]
    int_end_step = tup_step_range[1]
    
    
    # **********
    # HARDCODED MAX NUMBER OF PROCESSORS
    int_max_processors = 16
    # **********

    # Determine number of processors to use
    num_processors = min(mp.cpu_count() - 1, int_max_processors)
    

    # === Step 00
    if int_start_step == 0:
        fn_copy_folder(str_source_folder, str_output_dir, b_print_output)
        
    
    str_copied_model = os.path.join(str_output_dir, '00_base_hecras')
    str_hdf_file_path = fn_determine_ras_hdf(str_copied_model, tup_hecras_files_to_copy[2], 'p')
    
    # === Step 01
    if int_start_step <= 1 and int_end_step >= 1:
        fn_create_models_nextgen_gpkg(str_config_file_path,
                                      str_hdf_file_path,
                                      str_output_dir,
                                      b_print_output)
    
    str_model_hydrofabric_path_01 = os.path.join(str_output_dir,'01_model_hydrofabric', 'model_hydrofabric.gpkg')
    
    # === Step 02
    if int_start_step <= 2 and int_end_step >= 2:
        fn_spawn_hecras_copies(str_config_file_path,
                               str_output_dir,
                               tup_hecras_files_to_copy,
                               b_print_output)
        
    
    # === Step 03
    if int_start_step <= 3 and int_end_step >= 3:
        #fn_run_hec_ras_models(str_config_file_path, str_output_dir)
        print('  ** This step is removed on the Linux container **')
    
    
    
    str_model_path_02 = os.path.join(str_output_dir,'02_model_copies')
    str_hdf_path_03 = os.path.join(str_output_dir,'03_run_hdf')
    
    # Make the subfolders for processing
    os.makedirs(str_hdf_path_03, exist_ok=True)
    
    # Copy all the results "p01.hdf" from HEC-RAS to a new folder
    fn_copy_ras_hdf_plans(str_model_path_02, str_hdf_path_03)
    
    # ---- Preperation for steps 04, 05 and 06
    
    str_output_step_04 = os.path.join(str_output_dir,'04_wsel_and_hr')
    str_output_step_05 = os.path.join(str_output_dir,'05_large_terrain')
    str_output_step_06 = os.path.join(str_output_dir,'06_simple_rasters')
    
    # Make the subfolders for processing
    os.makedirs(str_output_step_04, exist_ok=True)
    
    # Directory containing source_terrain
    str_terrain_folder = os.path.join(str_output_dir,'02_model_copies', 'source_terrain')
    str_terrain_path = fn_determine_tif_in_dir_with_largest_bbox(str_terrain_folder)
    
    # List all .hdf files in the specified directory
    list_hdf_files = [f for f in os.listdir(str_hdf_path_03) if f.endswith('.hdf')]
    list_hdf_filepath = [os.path.join(str_hdf_path_03,f) for f in list_hdf_files]
    
    # remove '.p01.hdf' -- seven chars
    list_gpkg_files = ['hydrualic_results_' + f[:-8] + '.gpkg' for f in list_hdf_files]
    list_gpkg_filepath = [os.path.join(str_output_step_04,f) for f in list_gpkg_files]
    
    list_output_05_foldername = [f[:-8] for f in list_hdf_files]
    list_output_05_fullpath = [os.path.join(str_output_step_05,f) for f in list_output_05_foldername]
    
    # Create each step 05 folder if it does not already exist
    for folder in list_output_05_fullpath:
        os.makedirs(folder, exist_ok=True)
    
    # Prepare arguments for multiprocessing
    args = [(str_config_file_path, hdf_filepath, str_model_hydrofabric_path_01, str_output_step_04, b_print_output) for hdf_filepath in list_hdf_filepath]
    
    # === Step 04
    if int_start_step <= 4 and int_end_step >= 4:
        # Use multiprocessing with tqdm
        with mp.Pool(processes=num_processors) as pool:
            list(
                tqdm.tqdm(
                    pool.imap(fn_process_hdf, args),
                    total=len(args),
                    desc="Step 4: Create Model Geopackage",
                    bar_format="{desc}:({n_fmt}/{total_fmt})|{bar}| {percentage:.1f}%",
                    ncols=75
                )
            )
    
    # === Step 05
    if int_start_step <= 5 and int_end_step >= 5:
        print(" ")
        print('Step 5: Generate depth and WSEL rasters')
        
        # Note:  Step 05 script is multithreaded on a 'per-ru'n basis
        # Therefore, this step is running in serial with each 'run' running in parallel
        
        # Initialize tqdm progress bar
        with tqdm.tqdm(total=len(list_hdf_filepath),
                       desc='   -- Generating Rasters',
                       bar_format="{desc}:({n_fmt}/{total_fmt})|{bar}| {percentage:.1f}%",
                       ncols=75) as pbar:
            for int_count, item in enumerate(list_hdf_filepath):
                fn_generate_depth_wsel_rasters(
                    str_config_file_path,
                    list_hdf_filepath[int_count],
                    str_terrain_path,
                    list_gpkg_filepath[int_count],
                    list_output_05_fullpath[int_count],
                    b_print_output
                )
                pbar.update(1)  # Update progress bar after each iteration
        print("+-----------------------------------------------------------------+")
    
    
    # === Step 06
    int_count = 0
    if int_start_step <= 6 and int_end_step >= 6:
        print(" ")
        print('Step 6: Simplify rasters')
        
        # Note:  Step 05 script is multithreaded on a 'per-ru'n basis
        # Therefore, this step is running in serial with each 'run' running in parallel
        
        # Initialize tqdm progress bar
        with tqdm.tqdm(total=len(list_hdf_filepath),
                       desc='   -- Simplify Rasters',
                       bar_format="{desc}:({n_fmt}/{total_fmt})|{bar}| {percentage:.1f}%",
                       ncols=75) as pbar:
            for int_count, item in enumerate(list_hdf_filepath):
                fn_extract_and_simplify_rasters(
                    str_config_file_path,
                    list_gpkg_filepath[int_count],
                    list_output_05_fullpath[int_count],
                    str_output_step_06,
                    b_print_output
                )
                pbar.update(1)  # Update progress bar after each iteration
        print("+-----------------------------------------------------------------+")
    
    # === Step 07
    str_terrain_filename = os.path.basename(str_terrain_path)
    
    if int_start_step <= 7 and int_end_step >= 7:
        fn_prepare_netcdf_files(str_config_file_path,
                                str_terrain_filename,
                                str_output_dir,
                                b_print_output)
    
    # === Step 08
    if int_start_step <= 8 and int_end_step >= 8:
        fn_filter_netcdf(str_config_file_path,
                         str_output_dir,
                         b_print_output)
# +++++++++++++++++++++++++++++


# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist" % arg)
    else:
        # File exists so return the directory
        return arg
        return open(arg, 'r')  # return an open file handle
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# example string to run model:  C:\Users\civil\ras2fim-2d\src\20250107_python>
  #python ras2fim-2d.py 
  # -i E:\taylor_source_hecras 
  # -o E:\ras2fim2d_test_20250108 
  # -c C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini 
  # -f "(1,1,1)" 
  # -s "(6,6)"
  
  # python ras2fim-2d.py -i E:\taylor_source_hecras -o E:\ras2fim2d_test_20250108 -c C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini -f "(1,1,1)" -s "(6,6)"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()
    
    parser = argparse.ArgumentParser(description='========= RUN RAS2FIM-2D =========')
    
    parser.add_argument('-i',
                        dest = "str_source_folder",
                        help=r'REQUIRED: Source folder to copy',
                        required=True,
                        metavar='DIR',
                        type=str)
    
    parser.add_argument('-o',
                        dest = "str_output_dir",
                        help=r'REQUIRED: Destination folder',
                        required=True,
                        metavar='DIR',
                        type=str)
    
    parser.add_argument('-c',
                        dest = "str_config_file_path",
                        help=r'REQUIRED: Global configuration filepath Example: C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini',
                        required=True,
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-f',
                        dest = "tup_hecras_files_to_copy",
                        help=r'REQUIRED: Files to copy from source HEC-RAS (geom, unsteady flow, plan) Example: "(1,1,1)"',
                        required=True,
                        metavar='TUPLE',
                        type=fn_parse_tuple)
    
    parser.add_argument('-s',
                        dest = "tup_step_range",
                        help=r'OPTIONAL: Steps to run (start, end) Default: "(0,8)"',
                        required=False,
                        metavar='TUPLE',
                        default='(0,8)',
                        type=fn_parse_tuple)
    

    args = vars(parser.parse_args())
    
    str_source_folder = args['str_source_folder']
    str_output_dir = args['str_output_dir']
    str_config_file_path = args['str_config_file_path']
    tup_hecras_files_to_copy = args['tup_hecras_files_to_copy']
    tup_step_range = args['tup_step_range']
    
    fn_ras2fim_2d(str_source_folder, 
                  tup_hecras_files_to_copy, 
                  str_config_file_path, 
                  str_output_dir,
                  tup_step_range)

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)
    
    print('Compute Time: ' + str(time_pass))
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~