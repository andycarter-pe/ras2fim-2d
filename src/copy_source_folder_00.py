# Script 00 - Copy folder from source to destination '00_base_hecras' folder
#
# Created by: Andy Carter, PE
# Created - 2025.01.07
# ************************************************************

# ************************************************************
import shutil
import os
import argparse
import time
import datetime
import warnings
# ************************************************************

# +++++++++++++++++++++++++++++++
def fn_copy_all_files(str_input_folder, str_output_folder):
    # Check if the input folder exists
    if not os.path.exists(str_input_folder):
        raise ValueError(f"Input folder {str_input_folder} does not exist.")
    
    # Create the base output folder if it does not exist
    if not os.path.exists(str_output_folder):
        os.makedirs(str_output_folder)
    
    # Define the subdirectory for '00_base_hecras'
    str_output_folder = os.path.join(str_output_folder, '00_base_hecras')

    # Check if the folder exists; if yes, delete its contents
    if os.path.exists(str_output_folder):
        # Remove all contents of the folder
        for item in os.listdir(str_output_folder):
            item_path = os.path.join(str_output_folder, item)
            if os.path.isdir(item_path):
                shutil.rmtree(item_path)  # Remove directories
            else:
                os.remove(item_path)  # Remove files
    else:
        # Create the folder if it doesn't exist
        os.makedirs(str_output_folder)

    # Walk through the input folder
    for item in os.listdir(str_input_folder):
        src = os.path.join(str_input_folder, item)
        dest = os.path.join(str_output_folder, item)

        if os.path.isdir(src):
            # Recursively copy subdirectories
            shutil.copytree(src, dest)
        else:
            # Copy files
            shutil.copy2(src, dest)

    print(f"   -- All contents from {str_input_folder} have been copied to {str_output_folder}")
# +++++++++++++++++++++++++++++++


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


# ---------------------------------------------
def fn_copy_folder(str_source_folder, str_output_dir, b_print_output):
    
    # supress all warnings
    warnings.filterwarnings("ignore", category=UserWarning )
    
    print(" ")
    if b_print_output:
        print("+=================================================================+")
        print("|              COPY SOURCE HEC-RAS TO DESTINATION                 |")
        print("|                Created by Andy Carter, PE of                    |")
        print("|             Center for Water and the Environment                |")
        print("|                 University of Texas at Austin                   |")
        print("+-----------------------------------------------------------------+")
    
        
        print("  ---(i) INPUT SOURCE FOLDER TO COPY: " + str_source_folder)
        print("  ---(o) OUTPUT DESTINATION FOLDER: " + str_output_dir)
        print("  ---[r] PRINT OUTPUT: " + str(b_print_output))
        print("===================================================================")
    else:
        print('Step 0: Copy Source HEC-RAS to destination folder')
    
    
    fn_copy_all_files(str_source_folder, str_output_dir)
    print("+-----------------------------------------------------------------+")
# --------------------------------------------- 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()
    
    parser = argparse.ArgumentParser(description='========= COPY SOURCE HEC-RAS TO DESTINATION =========')
    
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
    
    parser.add_argument('-r',
                        dest = "b_print_output",
                        help=r'OPTIONAL: Print output messages Default: True',
                        required=False,
                        default=True,
                        metavar='T/F',
                        type=fn_str_to_bool)

    args = vars(parser.parse_args())
    
    str_source_folder = args['str_source_folder']
    str_output_dir = args['str_output_dir']
    b_print_output = args['b_print_output']

    fn_copy_all_files(str_source_folder, str_output_dir, b_print_output)

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)
    
    print('Compute Time: ' + str(time_pass))
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~