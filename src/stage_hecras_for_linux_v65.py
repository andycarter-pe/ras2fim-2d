# Purpose: Script automates the preparation of HEC-RAS project files for 
# unsteady flow simulations on Linux-based platforms, such as Docker or HPC 
# environments. It identifies valid HEC-RAS .prj files, generates necessary 
# output files, and organizes them for efficient transfer and execution. 
# Designed for use on a Windows machine with HEC-RAS installed, the script 
# streamlines workflows for large-scale simulations and parallel processing.
#
#
# Created by: Andy Carter, PE
# Revised - 2025.01.16
# ************************************************************

# ************************************************************
import win32com.client
import os

from tqdm import tqdm
import multiprocessing as mp

import argparse
import time
import datetime
import warnings
import configparser
import shutil
import re
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


# --------------
def fn_is_binary_string(file_path):
    textchars = bytearray({7, 8, 9, 10, 12, 13, 27} | set(range(0x20, 0x100)) - {0x7F})
    with open(file_path, "rb") as file:
        return bool(file.read(1024).translate(None, textchars))
# --------------


# -----------------------
def fn_find_prj_files(directory):
    # os.swalk a directory and get a list of the files
    # with a .prj extension
    prj_files = []
    
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".prj"):  # Check if the file has a .prj extension
                prj_files.append(os.path.join(root, file))
    
    return prj_files
# -----------------------

# --------------
def fn_filter_hecras_prj_files(str_filepath_for_prj):
    # for a list of .prj files, determine which are
    # HEC-RAS .prj files
    
    str_check = "Current Plan"  # Validation string for HEC-RAS .prj files
    list_valid_prj_files = []
    
    for file_path in str_filepath_for_prj:
        if fn_is_binary_string(file_path):
            # Don't want the binary files
            # Likely they are projection files
            pass
        else:
            # Open file as text and check for HEC-RAS identifier
            with open(file_path, "r") as file:
                if any(str_check in line for line in file):
                    list_valid_prj_files.append(file_path)
    
    return list_valid_prj_files
# --------------


# ~~~~~~~~~~~~~~~~~~~~~~~~~~

def fn_copy_files(str_prj_filepath, str_output_directory, list_needed_files):
    # Parse filename without extension
    str_filename_without_ext = os.path.splitext(os.path.basename(str_prj_filepath))[0]

    # Create the output folder
    str_target_folder = os.path.join(str_output_directory, str_filename_without_ext)
    os.makedirs(str_target_folder, exist_ok=True)
    #print(f"Created directory: {str_target_folder}")

    # Copy needed files to the created folder
    for file in list_needed_files:
        if os.path.exists(file):  # Ensure the file exists before copying
            shutil.copy(file, str_target_folder)
            #print(f"Copied: {file} -> {target_folder}")
        else:
            print(f"File not found, skipping: {file}")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~        



# ------------------
def fn_run_hecras_v66(str_hecras_controller, str_prj_filepath, list_needed_files, str_output_directory):
    
    # This is a fuction that is used to create the needed files from
    # a HEC-RAS so that an unsteady run can be performed on a Linux
    # Container (TACC, Docker, AWS, etc...)
    
    try:
        # Launch HEC-RAS and open project
        hec = win32com.client.Dispatch(str_hecras_controller)
        hec.Project_Open(str_prj_filepath)

        # Start computation in blocking mode
        # to allow monitoring of output
        is_blocking_mode = False
        NMsg, TabMsg = None, None

        #print("Starting HEC-RAS computation...")
        hec.Compute_CurrentPlan(NMsg, TabMsg, is_blocking_mode)
        
        while True:
            # Monitor files
            #print('---------')
            all_files_exist = all(os.path.exists(file) for file in list_needed_files)
            
            if all_files_exist:
                # if there is a line in list_needed_files[3] that looks like 
                #Message  | Computation engine: {randon text ... anything} not found.
                b_match_fouund = fn_find_engine_message(list_needed_files[3])
                
                if b_match_fouund:
                    #print("Files are ready, stopping HEC-RAS...")
                    hec.QuitRas()  # stop computation
                    
                    #print(str_output_directory)
                    fn_copy_files(str_prj_filepath, str_output_directory, list_needed_files)
                    
                    break
                
        # Wait before checking again
        time.sleep(1.0)
        
    except Exception as e:
        print(f"Error: {e}")
    finally:
        # Clean up COM properly
        try:
            if hec:
                hec.QuitRas()
        except Exception as e:
            print(f"Error quitting HEC-RAS: {e}")
# ------------------


# ....................
def fn_find_engine_message(file_path):
    """
    Reads a text file and determines if there's a line that matches the pattern:
    "Message  | Computation engine:" + random text + "not found."
    
    Args:
        file_path (str): Path to the text file.
    
    Returns:
        bool: True if a matching line is found, False otherwise.
    """
    pattern = r"Message\s+\|\s+Computation engine:.*?not found\."

    try:
        with open(file_path, 'r') as file:
            for line in file:
                if re.search(pattern, line):
                    return True
        return False
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return False
    except Exception as e:
        print(f"An error occurred: {e}")
        return False
# ....................


# ---------------
def fn_toggle_ras_unsteady_exe(str_original_path, str_renamed_path):

    try:
        if os.path.exists(str_original_path):
            # Rename
            os.rename(str_original_path, str_renamed_path)
            print(f"Renamed to: {str_renamed_path}")
        elif os.path.exists(str_renamed_path):
            # Rename back to original
            os.rename(str_renamed_path, str_original_path)
            print(f"Renamed to: {str_original_path}")
        else:
            print("Neither version of the file was found.")
    except PermissionError:
        print("Permission denied: Run the script with administrator privileges.")
    except Exception as e:
        print(f"An error occurred: {e}")
# ---------------


# ------------------
def fn_process_flow_step(args):
    """Process a single HEC-RAS project path with its controller."""
    prj_path, hecras_controller, str_output_directory = args
    list_needed_files = [
        prj_path[:-4] + '.b01',
        prj_path[:-4] + '.p01.tmp.hdf',
        prj_path[:-4] + '.x01',
        prj_path[:-4] + '.p01.computeMsgs.txt'
    ]
    fn_run_hecras_v66(hecras_controller, prj_path, list_needed_files, str_output_directory)
# ------------------


# ++++++++++++++++++++++++++++
def fn_prepare_hecras_for_linux(str_config_file_path,
                                str_serach_directory,
                                str_output_directory,
                                b_print_output):
    

    print(" ")
    if b_print_output:
        print("+=================================================================+")
        print("|           PREPARE HEC-RAS RUNS FOR LINUX SIMULATION             |")
        print("+-----------------------------------------------------------------+")
        print("|          THIS SCRIPT MUST BE RUN ON A WINDOWS MACHINE           |")
        print("|                   WITH HEC-RAS INSTALLED                        |")
        print("+-----------------------------------------------------------------+")
        print("|                Created by Andy Carter, PE of                    |")
        print("|             Center for Water and the Environment                |")
        print("|                 University of Texas at Austin                   |")
        print("+-----------------------------------------------------------------+")
    
        print("  ---(c) INPUT GLOBAL CONFIGURATION FILE: " + str_config_file_path)
        print("  ---(i) INPUT DIRECTORY OF HEC-RAS RUNS: " + str_serach_directory)
        print("  ---(o) OUTPUT DIRECTORY OF LINUX STAGED FILES: " + str_output_directory)
        print("  ---[r] PRINT OUTPUT: " + str(b_print_output))
        print("===================================================================")
    else:
        print('Step X: Prepare HEC-RAS for Linux')


    os.makedirs(str_output_directory, exist_ok=True)
    
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
    
    # ************* HARD CODED OVERRIDE
    str_hecras_controller = 'RAS66.HECRASController'
    # *************
    
    original_exe_path = r'C:\Program Files (x86)\HEC\HEC-RAS\6.6\x64\RasUnsteady.exe'
    renamed_exe_path = r'C:\Program Files (x86)\HEC\HEC-RAS\6.6\x64\x_RasUnsteady.exe'
    
    # Change the RasUnsteady.exe to x_RasUnsteady.exe
    fn_toggle_ras_unsteady_exe(original_exe_path, renamed_exe_path)
    

    # **********
    # HARED CODED MAX NUMBER OF PROCESSORS
    int_max_processors = 16
    # **********
    
    num_processors = (mp.cpu_count() - 1)
    
    if num_processors > int_max_processors:
        num_processors = int_max_processors
    
    # --- multiprocessing funtion ---
    
    list_prj_files = fn_find_prj_files(str_serach_directory)
    list_prj_path = fn_filter_hecras_prj_files(list_prj_files)
    
    l = len(list_prj_path)
    
    if l > 0:
        # Prepare arguments for multiprocessing
        list_args = [(prj_path, str_hecras_controller, str_output_directory) for prj_path in list_prj_path]
    
        with mp.Pool(processes=num_processors) as pool:
            # Wrap the multiprocessing call with tqdm for progress tracking
            for _ in tqdm(pool.imap(fn_process_flow_step, list_args),
                           total=l,
                           desc='Processing HEC-RAS Projects',
                           bar_format="{desc}:({n_fmt}/{total_fmt})|{bar}| {percentage:.1f}%",
                           ncols=75):
                pass
            
    # Change the x_RasUnsteady.exe to RasUnsteady.exe
    fn_toggle_ras_unsteady_exe(original_exe_path, renamed_exe_path)
# ++++++++++++++++++++++++++++


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()
    
    parser = argparse.ArgumentParser(description='======== PREPARE HEC-RAS RUNS FOR LINUX SIMULATION ========')
    
    parser.add_argument('-c',
                        dest = "str_config_file_path",
                        help=r'REQUIRED: Global configuration filepath Example: C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini',
                        required=False,
                        default=r'C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini',
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-i',
                        dest = "str_input_folder",
                        help=r'REQUIRED: folder containing HEC-RAS models Example: E:\ras2fim2d_test_20250108\02_model_copies',
                        required=False,
                        default=r'E:\ras2fim2d_test_20250116\02_model_copies',
                        metavar='DIR',
                        type=str)
    
    parser.add_argument('-o',
                        dest = "str_output_dir",
                        help=r'REQUIRED: directory to write Linux staged models Example: E:\ras2fim2d_test_20250116\xx-stage',
                        required=False,
                        default=r'E:\ras2fim2d_test_20250116\xx-stage',
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
    str_input_folder = args['str_input_folder']
    str_output_dir = args['str_output_dir']
    b_print_output = args['b_print_output']
    
    fn_prepare_hecras_for_linux(str_config_file_path,
                                str_input_folder,
                                str_output_dir,
                                b_print_output)

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)
    
    print('Compute Time: ' + str(time_pass))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~