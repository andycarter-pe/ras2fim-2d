# Script 04 - Run HEC-RAS 2D models
# on a windows machine
#
# Created by: Andy Carter, PE
# Revised - 2025.01.07
# ************************************************************


# ************************************************************
import win32com.client
import os
#import sys
#import shutil

import argparse
import time
import datetime
import warnings
import configparser
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


# ..........................................
def fn_run_hec_ras_models(str_config_file_path, str_output_folder):
    
    # supress all warnings
    warnings.filterwarnings("ignore", category=UserWarning )
    
    print(" ")
    print("+=================================================================+")
    print("|             RUN HEC-RAS MODELS (WINDOWS VERSION)                |")
    print("|                Created by Andy Carter, PE of                    |")
    print("|             Center for Water and the Environment                |")
    print("|                 University of Texas at Austin                   |")
    print("+-----------------------------------------------------------------+")

    
    print("  ---(c) INPUT GLOBAL CONFIGURATION FILE: " + str_config_file_path)
    print("  ---(o) OUTPUT PATH: " + str_output_folder)
    
    print("===================================================================")
    
    str_copied_model = os.path.join(str_output_folder, '02_model_copies')
    
    # --- Read variables from config.ini ---
    # Create a ConfigParser object
    config = configparser.ConfigParser()
    
    # Read the config.ini file
    config.read(str_config_file_path)
    
    # Get the variables in the [03_run_hec_ras] section
    if '03_run_hec_ras' in config:
        section = config['03_run_hec_ras']
    
        # Read variables in the section
        str_hecras_controller = section.get('str_hecras_controller', '')
    else:
        print("[03_run_hec_ras] section not found in the config file.")
    # --- ---
    
    # ----- Get a list of the HEC-RAS projects
    # Initialize a list to store the full paths of all .prj files
    list_prj_files = []
    
    # Walk through the directory and its subdirectories
    for root, dirs, files in os.walk(str_copied_model):
        for filename in files:  # Renamed to `filename` to avoid conflict
            if filename.endswith('.prj'):
                str_filepath = os.path.join(root, filename)
                
                # Attempt to open the file and read the first 11 characters
                try:
                    with open(str_filepath, 'r') as file_obj:
                        str_first_chars = file_obj.read(11)
                        
                        # Check for 'Proj Title=' at the beginning of the file
                        if str_first_chars == 'Proj Title=':
                            list_prj_files.append(str_filepath)
                except UnicodeDecodeError:
                    # Skip binary files or files that can't be read as text
                    print(f"Skipped binary or unreadable file: {str_filepath}")
    # -----
    
    
    # --- Run all the models
    int_count = 0
    
    for str_prj_path in list_prj_files:
        
        int_count += 1
        print(f'Model {int_count} of {len(list_prj_files)}')
        
        # Record start time
        start_time = time.time()
        
        hec = win32com.client.Dispatch(str_hecras_controller)
        hec.ShowRas()
        hec.Project_Open(str_prj_path)
        
        # Blocking mode set to True, so we can not check file creation while running
        is_blocking_mode = True
        
        # to be populated: number and list of messages, blocking mode
        NMsg, TabMsg, block = None, None, is_blocking_mode
        
        # computations of the current plan
        v1, NMsg, TabMsg, v2 = hec.Compute_CurrentPlan(NMsg, TabMsg, is_blocking_mode)
        
        hec.QuitRas()  # Close HEC-RAS
        
        # Calculate and print runtime for this iteration
        end_time = time.time()
        run_time = end_time - start_time
        print(f'Runtime for model {int_count}: {run_time:.2f} seconds\n')
        print('------------------')
# --- 
# ..........................................


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()
    
    parser = argparse.ArgumentParser(description='========= RUN HEC-RAS MODELS (WINDOWS MACHINE) =========')
    
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
    

    args = vars(parser.parse_args())
    
    str_config_file_path = args['str_config_file_path']
    str_output_dir = args['str_output_dir']
    
    fn_run_hec_ras_models(str_config_file_path, str_output_dir)

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)
    
    print('Compute Time: ' + str(time_pass))
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~