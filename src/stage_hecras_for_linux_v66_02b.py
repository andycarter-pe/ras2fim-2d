# Script 2b - Stage HEC-RAS for Linux Runs
# Purpose: Script automates the preparation of HEC-RAS project files for
# unsteady flow simulations on Linux-based platforms, such as Docker or HPC
# environments. It identifies valid HEC-RAS .prj files, generates necessary
# output files, and organizes them for transfer and execution.
# Designed for use on a Windows machine with HEC-RAS installed, the script
# streamlines workflows for large-scale simulations and parallel processing.
#
# Created by: Andy Carter, PE
# Revised - 2026.09.04
# ************************************************************

# ************************************************************
import win32com.client
import pythoncom
import os

from tqdm import tqdm
import multiprocessing as mp

import argparse
import time
import datetime
import configparser
import shutil
import re
# ************************************************************


# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist" % arg)
    else:
        # File exists so return the path
        return arg
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
    # Walk a directory and get a list of the files with a .prj extension
    prj_files = []

    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".prj"):  # Check if the file has a .prj extension
                prj_files.append(os.path.join(root, file))

    return prj_files
# -----------------------


# --------------
def fn_filter_hecras_prj_files(list_filepath_for_prj):
    # For a list of .prj files, determine which are HEC-RAS .prj files
    str_check = "Current Plan"  # Validation string for HEC-RAS .prj files
    list_valid_prj_files = []

    for file_path in list_filepath_for_prj:
        if fn_is_binary_string(file_path):
            # Skip binary files (likely projection files, etc.)
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

    # Copy needed files to the created folder
    for file in list_needed_files:
        if os.path.exists(file):  # Ensure the file exists before copying
            shutil.copy(file, str_target_folder)
        else:
            print(f"File not found, skipping: {file}")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~


# ....................
def fn_find_engine_message(file_path):
    """
    Reads a text file and determines if there's a line that matches the pattern:
    "Message  | Computation engine:" + random text + "not found."

    This is the signal that HEC-RAS has finished writing the preprocessed files
    and has hit the (intentionally) renamed RasUnsteady.exe.

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
        return False
    except Exception as e:
        print(f"An error occurred reading {file_path}: {e}")
        return False
# ....................


# ------------------
def fn_run_hecras_v66(str_hecras_controller,
                      str_prj_filepath,
                      list_needed_files,
                      str_output_directory,
                      flt_timeout_sec):
    """
    Create the needed preprocessed files from a HEC-RAS project so that an
    unsteady run can be performed on a Linux container (TACC, Docker, AWS...).

    Returns a (bool_success, str_message) tuple so the parent can report status.
    """
    hec = None
    try:
        # Launch HEC-RAS and open project
        hec = win32com.client.Dispatch(str_hecras_controller)
        hec.Project_Open(str_prj_filepath)

        # Start computation in non-blocking mode so we can monitor output files
        is_blocking_mode = False
        NMsg, TabMsg = None, None

        hec.Compute_CurrentPlan(NMsg, TabMsg, is_blocking_mode)

        flt_start = time.time()

        while True:
            all_files_exist = all(os.path.exists(f) for f in list_needed_files)

            if all_files_exist:
                # Wait for the "engine not found" message in the computeMsgs file,
                # which confirms preprocessing wrote everything we need.
                if fn_find_engine_message(list_needed_files[3]):
                    hec.QuitRas()
                    fn_copy_files(str_prj_filepath, str_output_directory, list_needed_files)
                    return (True, f"OK: {os.path.basename(str_prj_filepath)}")

            # Timeout guard so a stuck run cannot block forever
            if time.time() - flt_start > flt_timeout_sec:
                try:
                    hec.QuitRas()
                except Exception:
                    pass
                return (False, f"TIMEOUT after {flt_timeout_sec:.0f}s: "
                               f"{os.path.basename(str_prj_filepath)}")

            # Actually wait between checks (this was the main bug: the sleep was
            # outside the loop, so this spun at 100% CPU and starved HEC-RAS).
            time.sleep(1.0)

    except Exception as e:
        return (False, f"ERROR on {os.path.basename(str_prj_filepath)}: {e}")
    finally:
        # Clean up COM properly
        try:
            if hec is not None:
                hec.QuitRas()
        except Exception:
            pass
# ------------------


# ---------------
def fn_toggle_ras_unsteady_exe(str_original_path, str_renamed_path):
    try:
        if os.path.exists(str_original_path):
            os.rename(str_original_path, str_renamed_path)
            print(f"Renamed to: {str_renamed_path}")
        elif os.path.exists(str_renamed_path):
            os.rename(str_renamed_path, str_original_path)
            print(f"Renamed to: {str_original_path}")
        else:
            print("Neither version of the RasUnsteady exe was found.")
    except PermissionError:
        print("Permission denied: Run the script with administrator privileges.")
    except Exception as e:
        print(f"An error occurred renaming the exe: {e}")
# ---------------


# ------------------
def fn_process_flow_step(args):
    """Process a single HEC-RAS project path with its controller (worker)."""
    prj_path, hecras_controller, str_output_directory, flt_timeout_sec = args

    list_needed_files = [
        prj_path[:-4] + '.b01',
        prj_path[:-4] + '.p01.tmp.hdf',
        prj_path[:-4] + '.x01',
        prj_path[:-4] + '.p01.computeMsgs.txt'
    ]

    # COM must be initialized in each worker process
    pythoncom.CoInitialize()
    try:
        return fn_run_hecras_v66(hecras_controller,
                                 prj_path,
                                 list_needed_files,
                                 str_output_directory,
                                 flt_timeout_sec)
    finally:
        pythoncom.CoUninitialize()
# ------------------


# ++++++++++++++++++++++++++++
def fn_prepare_hecras_for_linux(str_config_file_path,
                                str_search_directory,
                                str_output_directory,
                                int_processes,
                                flt_timeout_sec,
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
        print("  ---(i) INPUT DIRECTORY OF HEC-RAS RUNS: " + str_search_directory)
        print("  ---(o) OUTPUT DIRECTORY OF LINUX STAGED FILES: " + str_output_directory)
        print("  ---(p) WORKER PROCESSES: " + str(int_processes))
        print("  ---(t) PER-PROJECT TIMEOUT (s): " + str(flt_timeout_sec))
        print("  ---[r] PRINT OUTPUT: " + str(b_print_output))
        print("===================================================================")
    else:
        print('Step 2b: Prepare HEC-RAS for Linux')

    os.makedirs(str_output_directory, exist_ok=True)

    # --- Read variables from config.ini ---
    config = configparser.ConfigParser()
    config.read(str_config_file_path)

    str_hecras_controller = ''
    if '03_run_hec_ras' in config:
        section = config['03_run_hec_ras']
        str_hecras_controller = section.get('str_hecras_controller', '')
    else:
        print("[03_run_hec_ras] section not found in the config file.")
    # --- ---

    # ************* HARD CODED OVERRIDE
    str_hecras_controller = 'RAS66.HECRASController'
    # *************

    original_exe_path = r'C:\Program Files (x86)\HEC\HEC-RAS\6.6\x64\RasUnsteady.exe'
    renamed_exe_path = r'C:\Program Files (x86)\HEC\HEC-RAS\6.6\x64\x_RasUnsteady.exe'

    # Find and validate the HEC-RAS projects up front
    list_prj_files = fn_find_prj_files(str_search_directory)
    list_prj_path = fn_filter_hecras_prj_files(list_prj_files)
    l = len(list_prj_path)

    if l == 0:
        print("No valid HEC-RAS .prj files found. Nothing to do.")
        return

    # Clamp the worker count to something sane
    num_processors = max(1, min(int_processes, mp.cpu_count(), l))

    list_results = []

    # Rename RasUnsteady.exe -> x_RasUnsteady.exe. Wrap the whole compute in
    # try/finally so we ALWAYS rename it back, even on Ctrl-C or error.
    fn_toggle_ras_unsteady_exe(original_exe_path, renamed_exe_path)
    try:
        if num_processors == 1:
            # Drive HEC-RAS COM from the MAIN process. The HECRASController is an
            # STA COM server and is unreliable when launched inside a spawned
            # multiprocessing worker (even a single one) -- it typically hangs at
            # Dispatch / Compute_CurrentPlan with no error until the timeout trips.
            # The reliable pattern is one project at a time, in-process.
            pythoncom.CoInitialize()
            try:
                for prj_path in tqdm(list_prj_path,
                                     total=l,
                                     desc='Processing HEC-RAS Projects',
                                     bar_format="{desc}:({n_fmt}/{total_fmt})|{bar}| {percentage:.1f}%",
                                     ncols=75):
                    list_needed_files = [
                        prj_path[:-4] + '.b01',
                        prj_path[:-4] + '.p01.tmp.hdf',
                        prj_path[:-4] + '.x01',
                        prj_path[:-4] + '.p01.computeMsgs.txt'
                    ]
                    list_results.append(
                        fn_run_hecras_v66(str_hecras_controller,
                                          prj_path,
                                          list_needed_files,
                                          str_output_directory,
                                          flt_timeout_sec))
            finally:
                pythoncom.CoUninitialize()
        else:
            list_args = [
                (prj_path, str_hecras_controller, str_output_directory, flt_timeout_sec)
                for prj_path in list_prj_path
            ]

            with mp.Pool(processes=num_processors) as pool:
                for result in tqdm(pool.imap(fn_process_flow_step, list_args),
                                   total=l,
                                   desc='Processing HEC-RAS Projects',
                                   bar_format="{desc}:({n_fmt}/{total_fmt})|{bar}| {percentage:.1f}%",
                                   ncols=75):
                    list_results.append(result)
    finally:
        # Always restore the original exe name
        fn_toggle_ras_unsteady_exe(original_exe_path, renamed_exe_path)

    # --- Report ---
    list_failures = [msg for ok, msg in list_results if not ok]
    print("\n-----------------------------------------------------------------")
    print(f"Completed {l - len(list_failures)} / {l} projects successfully.")
    if list_failures:
        print("The following projects did not complete:")
        for msg in list_failures:
            print("  - " + msg)
    print("-----------------------------------------------------------------")
# ++++++++++++++++++++++++++++


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()

    parser = argparse.ArgumentParser(
        description='======== PREPARE HEC-RAS RUNS FOR LINUX SIMULATION ========')

    parser.add_argument('-c',
                        dest="str_config_file_path",
                        help=r'REQUIRED: Global configuration filepath '
                             r'Example: C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini',
                        required=False,
                        default=r'C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini',
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument('-i',
                        dest="str_input_folder",
                        help=r'REQUIRED: folder containing HEC-RAS models '
                             r'Example: E:\ras2fim2d_test_20250108\02_model_copies',
                        required=False,
                        default=r'E:\ras2fim2d_test_20250116\02_model_copies',
                        metavar='DIR',
                        type=str)

    parser.add_argument('-o',
                        dest="str_output_dir",
                        help=r'REQUIRED: directory to write Linux staged models '
                             r'Example: E:\ras2fim2d_test_20250116\xx-stage',
                        required=False,
                        default=r'E:\ras2fim2d_test_20250116\xx-stage',
                        metavar='DIR',
                        type=str)

    parser.add_argument('-p',
                        dest="int_processes",
                        help=r'OPTIONAL: number of parallel HEC-RAS workers. '
                             r'HEC-RAS COM does not parallelize well; start low. Default: 4',
                        required=False,
                        default=4,
                        metavar='INT',
                        type=int)

    parser.add_argument('-t',
                        dest="flt_timeout_sec",
                        help=r'OPTIONAL: per-project timeout in seconds. Default: 900',
                        required=False,
                        default=900.0,
                        metavar='SEC',
                        type=float)

    parser.add_argument('-r',
                        dest="b_print_output",
                        help=r'OPTIONAL: Print output messages Default: True',
                        required=False,
                        default=True,
                        metavar='T/F',
                        type=fn_str_to_bool)

    args = vars(parser.parse_args())

    fn_prepare_hecras_for_linux(args['str_config_file_path'],
                                args['str_input_folder'],
                                args['str_output_dir'],
                                args['int_processes'],
                                args['flt_timeout_sec'],
                                args['b_print_output'])

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)

    print('Compute Time: ' + str(time_pass))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~