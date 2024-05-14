# Script 03 - Using the HEC-RAS controller, run the HEC-RAS 'current plan' in
# the provided project.
#
# Created by: Andy Carter, PE
# Created - 2024.05.14

# ************************************************************

import win32com.client
import os
import sys

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


# ~~~~~~~~~~~~~~~~~~~~~~
# Define a function for the spinning cursor animation
def fn_spinning_cursor():
    while True:
        for cursor in '|/-\\':
            yield cursor
# ~~~~~~~~~~~~~~~~~~~~~~

# -------
def fn_filepath_expected(str_input_filepath, 
                         str_input_extension, 
                         str_output_extension):

    # From a given input path of a HEC-RAS file, determine if the run has created
    # the required Linux run inputs
    
    # Split the file path into directory, base filename, and extension
    directory, filename = os.path.split(str_input_filepath)
    base_filename, extension = os.path.splitext(filename)

    # Example: Replace 'gXX' with 'bXX' in the extension
    mod_extension = extension.replace(str_input_extension, str_output_extension)

    # Construct the file path
    str_ammended_filepath = os.path.join(directory, base_filename + mod_extension)
    
    return(str_ammended_filepath)
# -------


# .........................................................
def fn_run_hecras_model(str_config_file_path,
                        str_ras_projectpath):
    
    # supress all warnings
    #warnings.filterwarnings("ignore", category=UserWarning )
  
    print(" ")
    print("+=================================================================+")
    print('|             RUN HEC-RAS 2D MODEL (WINDOWS MACHINE)              |')
    print("|                Created by Andy Carter, PE of                    |")
    print("|             Center for Water and the Environment                |")
    print("|                 University of Texas at Austin                   |")
    print("+-----------------------------------------------------------------+")

    
    print("  ---(c) INPUT GLOBAL CONFIGUTATION FILE: " + str_config_file_path)
    print("  ---(p) INPUT MODEL HYROFABRIC (GPKG): " + str_ras_projectpath)
    print("===================================================================")
    
    # --- Read variables from config.ini ---
    # Create a ConfigParser object
    config = configparser.ConfigParser()
    
    # Read the config.ini file
    config.read(str_config_file_path)
    
    # Get the variables in the [02_copy_geom] section
    if '03_run_hec_ras' in config:
        section = config['03_run_hec_ras']
        
        # Read variables in the section
        str_hecras_controller = section.get('str_hecras_controller', '')
    else:
        print("[03_run_hec_ras] section not found in the config file.")
    # --- ---
    
    # TODO 2024.05.14 - the config file is not getting the controller string correctly 
    # Currently hard coding... need to fix
    
    hec = win32com.client.Dispatch('RAS610.HECRASController')
    
    hec.ShowRas()
    #hec.Compute_HideComputationWindow()
    
    print('Loading HEC-RAS project...')
    
    # Run the active HEC-RAS project 
    hec.Project_Open(str_ras_projectpath)   # opening HEC-RAS
    
    # ---------
    # TODO 2024.05.14 - This dictionary definition does not work with 6.5
    # Currently using 6.1.0 to match Linux version.
    
    # dictionary of this current HEC-RAS run
    dict_hec_info = {
        "HEC-RAS_Version": hec.HECRASVersion(),
        "Project_Path": hec.CurrentProjectFile(),
        "Project_Title": hec.CurrentProjectTitle(),
        "Plan_Path": hec.CurrentPlanFile(),
        "Geometry_Path": hec.CurrentGeomFile(),
        "Geometry_HDF_Path": hec.CurrentGeomHDFFile(),
        "Unsteady_File_Path": hec.CurrentUnsteadyFile()
    }
    
    # name of 2d areas in active geometry
    tup_2dareas = hec.Geometry_Get2DFlowAreas()
    
    # Extracting the string part, stripping whitespace, and converting it to a list of strings
    list_of_strings_2darea = [item.strip() for sublist in tup_2dareas if isinstance(sublist, tuple) for item in sublist]
    
    # Add the list_of_strings_2darea to the dict_hec_info under a suitable key
    dict_hec_info["2D_Flow_Area_Names"] = list_of_strings_2darea
    # ---------
    
    
    # Determine the list of files that are to be created from this run
    # This could be in preperation for Linux simulation.
    
    # HEC-RAS should create a matching cXX for the geometry file gXX
    str_expected_c_file = fn_filepath_expected(hec.CurrentGeomFile(), 'g', 'c')
    
    # HEC-RAS should create a matching bXX for the geometry file pXX
    str_expected_b_file = fn_filepath_expected(hec.CurrentPlanFile(), 'p', 'b')
    
    # HDF Version of the plan
    str_plan_hdf_filepath = hec.CurrentPlanFile() + ".tmp.hdf"
    
    # create a list of the expected files
    list_expected_files = [str_expected_c_file, str_expected_b_file, str_plan_hdf_filepath]
    # ---------
    
    # -- Running the HEC-RAS Model --
    
    # Blocking mode set to false, so we can check file creation while running
    is_blocking_mode = False
    
    # option to kill HEC-RAS
    b_kill_hecras = False
    
    # to be populated: number and list of messages, blocking mode
    NMsg, TabMsg, block = None, None, is_blocking_mode
    
    print('Loading current plan ...')
    
    # computations of the current plan
    v1, NMsg, TabMsg, v2 = hec.Compute_CurrentPlan(NMsg, TabMsg, is_blocking_mode)
    
    start_time = time.time()
    
    print('Starting Simulation (please wait)...')
    
    spinner = fn_spinning_cursor()
    
    while not all(os.path.exists(file_path) for file_path in list_expected_files) or b_kill_hecras:
        
        # Print the spinning cursor
        sys.stdout.write(next(spinner))
        sys.stdout.flush()
        time.sleep(0.1)
        sys.stdout.write('\b')
        
        list_files_exists = [os.path.exists(file_path) for file_path in list_expected_files]
    
        if all(list_files_exists):
            
            # All the files exist
            b_kill_hecras = True
            
            # kill the run
            hec.Compute_Cancel()
            break  # Exit the loop once all files exist
    
    print('COMPLETE')
    
    hec.QuitRas()  # close HEC-RAS
        
# .........................................................


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()
    
    parser = argparse.ArgumentParser(description='============ RUN A HEC-RAS 2D MODEL (WINDOWS MACHINE) ============')
    
    parser.add_argument('-c',
                        dest = "str_config_file_path",
                        help=r'REQUIRED: Global configuration filepath Example: C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini',
                        required=False,
                        default=r'C:\Users\civil\ras2fim-2d\src\python_code\config_global.ini',
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    parser.add_argument('-p',
                        dest = "str_ras_projectpath",
                        help=r'REQUIRED: Path to HEC-RAS project (.prj) Example: E:\HECRAS_2D_12070205\base_model_20240414_copy\BLE_LBSG_501.prj',
                        required=False,
                        default=r'E:\HECRAS_2D_12070205\base_model_20240414_copy\BLE_LBSG_501.prj',
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))
    
    args = vars(parser.parse_args())
    
    str_config_file_path = args['str_config_file_path']
    str_ras_projectpath = args['str_ras_projectpath']
    
    fn_run_hecras_model(str_config_file_path,
                        str_ras_projectpath)

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)
    
    print('Compute Time: ' + str(time_pass))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~