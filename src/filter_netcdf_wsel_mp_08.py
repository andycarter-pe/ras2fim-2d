# Script 08 - Filter the created netCDF.
#
# Created by: Andy Carter, PE
# Created - 2025.01.08

# ************************************************************
import xarray as xr
import numpy as np
import os
from datetime import datetime
import matplotlib.pyplot as plt
import configparser
import multiprocessing as mp
import tqdm
from functools import partial

import argparse
import time
import warnings
# ************************************************************


# -----------------
def fn_retain_indices_within_difference(arr, threshold):
    # Initialize the list to store indices
    indices_to_retain = [0]  # Always retain the first element
    
    # Loop through the array starting from the second element
    for i in range(1, len(arr)):
        # Check if the difference with the previous value is within the threshold
        if abs(arr[i] - arr[indices_to_retain[-1]]) <= threshold:
            # Only retain the index if difference is within threshold and it's the start of a new group
            continue
        else:
            i = i - 1 # to get values that are below the threshold
            indices_to_retain.append(i)
    
    return indices_to_retain
# -----------------


# ----
def fn_interpolate_x(x1, y1, x2, y2, y3):
    if x1 == x2:
        if y1 == y2 and y1 == y3:
            return x1 # or x2, they are the same
        return None  # Vertical line, cannot interpolate x from y unless it's the same point

    if (y3 < min(y1, y2) or y3 > max(y1, y2)) and y1 != y2:
        raise ValueError("y3 is outside the range of y1 and y2.")

    x3 = x1 + (y3 - y1) * (x2 - x1) / (y2 - y1)
    return x3
# ----

# ----------------
def fn_filter_netcdf_layers(str_netcdf_path, str_out_dir, flt_threshold):
    
    str_filename = os.path.split(str_netcdf_path)[1]

    # Between each layer in ds['wsel'] determine the average difference 
    # between the values, ignoring NAN cells in either of the
    # two compared layers -- this is to remove the "unstable" layers

    # This uses only the cells that are 'wet' for all the flows in the data cube
    # For example: only cells that are wet for the lowest flow are evaluated
    # across the entire cube

    # Load the NetCDF file
    ds = xr.open_dataset(str_netcdf_path)

    # Access the 'wsel' DataArray directly
    wsel = ds['wsel'].values  # Convert to NumPy array for faster computation

    # Create a mask of finite values
    valid_mask = np.isfinite(wsel)

    # Calculate differences between consecutive layers
    differences = np.diff(wsel, axis=0)  # Shape: (flow-1, y, x)

    # Mask invalid cells in consecutive layers
    valid_pairs = valid_mask[:-1, :, :] & valid_mask[1:, :, :]
    masked_differences = np.where(valid_pairs, differences, np.nan)

    # Compute the average difference for each pair of consecutive layers
    arr_average_differences = np.nanmean(masked_differences, axis=(1, 2))

    # Get a list of indices where the average difference is negative
    # This is likely a balooning due to model instability
    negative_diff_indices = np.where(arr_average_differences < 0)[0]

    # Create a boolean mask for the layers to keep
    num_flows = ds['wsel'].shape[0]  # Total number of flow layers
    keep_mask = np.ones(num_flows, dtype=bool)  # Start with all True
    keep_mask[negative_diff_indices] = False  # Set indices to False for removal

    # Remove the unwanted layers
    ds_filtered = ds.isel(flow=keep_mask)

    #print(f"Layers dropped: {len(negative_diff_indices)}")

    # from valid cells, subtract the ds_filtered['terrain'] value to create a depth value
    # then get the average depth value per layer

    # Access the 'wsel' and 'terrain' DataArrays
    wsel_filtered = ds_filtered['wsel'].values  # Convert to NumPy array
    terrain = ds_filtered['terrain'].values  # Convert to NumPy array

    # Create a mask where all layers have non-NaN values at each cell
    valid_mask_all_layers = np.all(np.isfinite(wsel_filtered), axis=0)  # Shape: (y, x)

    # Mask the wsel array to keep only valid cells
    valid_wsel = np.where(valid_mask_all_layers, wsel_filtered, np.nan)

    # Subtract terrain values to calculate depth, ensuring terrain is applied only to valid cells
    depth = np.where(valid_mask_all_layers, valid_wsel - terrain, np.nan)

    # Compute the average depth per layer
    average_depth_per_layer = np.nanmean(depth, axis=(1, 2))  # Shape: (flow,)

    # Extract the 'flow' values from the dataset
    arr_flow_values = ds_filtered['flow'].values  # Convert to a NumPy array

    # Insert 0 at the beginning of both arrays using np.insert()
    arr_flow_values = np.insert(arr_flow_values, 0, 0)
    average_depth_per_layer = np.insert(average_depth_per_layer, 0, 0)

    # -----
    # Filter the layers to only those layers where the 
    # average change (depth) between layers is just less than flt_threshold
    list_indices = fn_retain_indices_within_difference(average_depth_per_layer, flt_threshold)

    # add the last value
    list_indices.append(len(average_depth_per_layer)-1)

    # Remove duplicates in list_indices
    list_unique_indices = list(dict.fromkeys(list_indices))
    
    #print(list_unique_indices)
    
    # ~~~~~~ Create revised netCDFs ~~~~~~~
    list_layers_to_keep = list_unique_indices[1:] # drop the zero layer
    # create a new netCDF file that retains only those wsel layers with these indecies
    
    # Subtract 1 from each item in list_layers_to_keep
    list_layers_to_keep = [index - 1 for index in list_layers_to_keep]
    
    # Select the layers to keep using `isel`
    ds_final_filtered = ds_filtered.isel(flow=list_layers_to_keep)
    
    # Specify output file path
    # Create the directory if it doesn't exist
    if not os.path.exists(str_out_dir):
        os.makedirs(str_out_dir)
        
    output_filename = f"filter_{os.path.basename(str_netcdf_path)}"
    output_path = os.path.join(str_out_dir, output_filename)
    
    # ~~~~~~ End Create revised netCDFs ~~~~~~~
    
    # Filter the average_depth_per_layer array using the indices
    arr_filtered_average_depth = average_depth_per_layer[list_unique_indices]

    # Filter the average_depth_per_layer array using the indices
    arr_filtered_flows = arr_flow_values[list_unique_indices]
    
    # Determine the depths and flows needed that are missing for
    # Desired vertical gradation

    # Calculate the difference between consecutive values
    arr_differences = np.diff(arr_filtered_average_depth)

    list_tup_needed = []

    list_tup_start_index_val = []
    for index, value in enumerate(arr_differences):
        if value > flt_threshold:
            list_tup_start_index_val.append((index, value))

    for tup_item in list_tup_start_index_val:
        flt_flow_start = arr_filtered_flows[tup_item[0]]

        # Ensure you don't go beyond the bounds of the array
        if tup_item[0] + 1 < len(arr_filtered_flows):
            flt_flow_end = arr_filtered_flows[tup_item[0] + 1]
            flt_start_depth = arr_filtered_average_depth[tup_item[0]]
            flt_end_depth = arr_filtered_average_depth[tup_item[0] + 1]

            # number of values to interpolate
            int_needed_steps = int(tup_item[1] // flt_threshold)

            for i in range(int_needed_steps):
                flt_y3 = flt_start_depth + (i+1) * flt_threshold

                # interpolate values
                flt_flow_needed = fn_interpolate_x(flt_flow_start, flt_start_depth,
                                                   flt_flow_end, flt_end_depth,
                                                   flt_y3)

                list_tup_needed.append((round(flt_flow_needed,0), round(flt_y3,3)))
        else:
            print(f"Index {tup_item[0] + 1} is out of range. Cannot calculate flt_flow_end.")
    
    # Extract x and y values from list_tup_needed
    if list_tup_needed:  # Check if the list is not empty
        needed_flows, needed_depths = zip(*list_tup_needed)
    else:
        needed_flows = []
        needed_depths = []
        
    '''
    # ---- Print a graphic ---
    # ------------------------
    
    # TODO - 2025.01.08 - Save the rating curve graphic
    
    # Create figure and axis objects
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot the original data with smaller markers
    ax.plot(arr_flow_values, average_depth_per_layer, marker='+', markersize=6, linestyle='-', color='b', label='RAS2FIM-2D')

    # Plot the original data with smaller markers
    ax.plot(arr_filtered_flows, arr_filtered_average_depth, marker='o', markersize=8, linestyle=' ', color='r', label='Retained Layers')

    # Plot the interpolated data
    ax.plot(needed_flows, needed_depths, marker='o', markersize=6, linestyle=' ', color='g', label='Needed Values')

    # Add labels and title
    ax.set_xlabel('Flow (cfs)')
    ax.set_ylabel('Depth (ft)')
    ax.set_title(str_filename)
    ax.legend()

    # Add a grid for better readability
    ax.grid(True)

    # Set background color of the entire figure
    fig.patch.set_facecolor('#d3d3d3')  # Set gray background color for the entire figure

    # Show the plot
    plt.show()
    # ---- End Print a graphic ---
    # ----------------------------
    '''
    
    # Append the additional filter data to the netCDF
    # as global attributes
    ds_final_filtered.attrs['08_filter_history'] = f"Filtered on {datetime.now().strftime('%Y-%m-%d %H:%M')}"
    ds_final_filtered.attrs['09_vertical_filter'] = flt_threshold
    ds_final_filtered.attrs['10_still_needed_values'] = list(needed_flows)  # Converts it to a list
    
    # Save the new dataset to a NetCDF file
    ds_final_filtered.to_netcdf(output_path)
    
    # Close the datasets
    ds.close()
    ds_filtered.close()
    ds_final_filtered.close()
    
    '''
    print(f"Vertical gradation: {flt_threshold}")
    print(f"  Points available: {len(wsel)}")
    print(f"  Points dropped for stability issues: {len(negative_diff_indices)}")
    print(f"  Points retained: {len(arr_filtered_flows) -1 }")
    print(f"  Additional points needed: {len(needed_flows)}")
    '''
# ----------------


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
def fn_filter_netcdf(str_config_file_path, str_output_dir,b_print_output):
    
    # supress all warnings
    warnings.filterwarnings("ignore", category=UserWarning )
    
    print(" ")
    if b_print_output:
        print("+=================================================================+")
        print("|                   FILTER NETCDF WSEL FILES                      |")
        print("|                Created by Andy Carter, PE of                    |")
        print("|             Center for Water and the Environment                |")
        print("|                 University of Texas at Austin                   |")
        print("+-----------------------------------------------------------------+")
    
        
        print("  ---(c) INPUT GLOBAL CONFIGUTATION FILE: " + str_config_file_path)
        print("  ---(o) OUTPUT FOLDER: " + str_output_dir)
        print("  ---[r] PRINT OUTPUT: " + str(b_print_output))
        print("===================================================================")
    else:
        print('Step 8: Filter NetCDF WSEL files')
    
    str_directory = os.path.join(str_output_dir, '06_simple_rasters', '02_wsel_nc')
    str_out_dir = os.path.join(str_output_dir, '06_simple_rasters', '03_wsel_nc_filtered')
    
    # --- Read variables from config.ini ---
    # Create a ConfigParser object
    config = configparser.ConfigParser()
    
    # Read the config.ini file
    config.read(str_config_file_path)
    
    # Get the variables in the [08_filter_netcdf] section
    if '08_filter_netcdf' in config:
        section = config['08_filter_netcdf']
    
        # Read variables in the section
        flt_threshold = section.getfloat('flt_threshold', 0)
    else:
        print("[08_filter_netcdf] section not found in the config file.")
    # ---
    
    list_nc_files = [os.path.join(str_directory, file) for file in os.listdir(str_directory) if file.endswith(".nc")]
    
    # Prepare arguments for multiprocessing
    args = [(str_netcdf_path, str_out_dir, flt_threshold) for str_netcdf_path in list_nc_files]
    
    # **********
    # HARED CODED MAX NUMBER OF PROCESSORS
    int_max_processors = 16
    # **********
    
    num_processors = (mp.cpu_count() - 1)
    
    if num_processors > int_max_processors:
        num_processors = int_max_processors
    
    # Use multiprocessing with tqdm
    with mp.Pool(processes=num_processors) as pool:
        # Use tqdm to show progress
        list(
            tqdm.tqdm(
                pool.imap(
                    partial(fn_filter_netcdf_layers, str_out_dir=str_out_dir, flt_threshold=flt_threshold),
                    [str_netcdf_path for str_netcdf_path, _, _ in args]
                ),
                total=len(list_nc_files),
                desc='   -- Filter NetCDFs',
                bar_format="{desc}:({n_fmt}/{total_fmt})|{bar}| {percentage:.1f}%",
                ncols=75
            )
        )
    print("+-----------------------------------------------------------------+")
# ++++++++++++++++++++++++++++

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()
    
    parser = argparse.ArgumentParser(description='======== FILTER NETCDF WSEL FILES ========')
    
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
    
    parser.add_argument('-r',
                        dest = "b_print_output",
                        help=r'OPTIONAL: Print output messages Default: True',
                        required=False,
                        default=True,
                        metavar='T/F',
                        type=fn_str_to_bool)
    
    args = vars(parser.parse_args())
    
    str_config_file_path = args['str_config_file_path']
    str_output_dir = args['str_output_dir']
    b_print_output = args['b_print_output']

    fn_filter_netcdf(str_config_file_path,
                     str_output_dir,
                     b_print_output)

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)
    
    print('Compute Time: ' + str(time_pass))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~