<!-- Contributed by CLB Engineering for Andy Carter, PE, 2026.10.27. Distilled from
copy_source_folder_00.py, determine_fim_per_reach_04b.py and filter_netcdf_wsel_mp_08.py. -->

# ras2fim-2d step template

Fill the `<...>` placeholders. Keep the spacing, rule lines and blank lines exactly as shown;
they copy Andy's files. Delete blocks marked `[CONFIG]` or `[MP]` when the step doesn't need them.

| Placeholder | Example |
|---|---|
| `<NN>` | `09`, `05b` |
| `<Purpose>` | `Check the filtered netCDF global attributes` |
| `<TITLE>` | `CHECK NETCDF ATTRIBUTES`, centered in the 65-character banner field |
| `<Title>` | `Check NetCDF attributes` (the one-line `Step 9:` message; no leading zero) |
| `<entry>` | `fn_check_netcdf_attributes` (no step number, like `fn_filter_netcdf`) |
| `<NN_section>` | `09_check_netcdf` |
| `<NN_folder>` | `09_netcdf_check` |
| `<YYYY.MM.DD>` | today, e.g. `2026.10.27` |

Banner line: `"|" + "<TITLE>".center(65) + "|"` gives the exact width of the `+====+` rule.
Argument lines: `(x)` for required inputs, `[x]` for optional ones; the letter matches the flag
(`(c)` config, `(i)` input, `(o)` output, `[r]` print output).

```python
# Script <NN> - <Purpose>
#
# <What it reads (which NN_ folder or file) and what it writes (which NN_ folder).>
#
# Created by: Andy Carter, PE
# Created - <YYYY.MM.DD>
# ************************************************************

# ************************************************************
import os
import configparser                 # [CONFIG]
import multiprocessing as mp        # [MP]
import tqdm                         # [MP]
from functools import partial       # [MP]

import argparse
import time
import datetime
import warnings
# ************************************************************


# ----------------
def fn_<helper>(str_input_path, str_out_dir, flt_<setting>):
    # one job: plain procedural code, prefixed names (str_, int_, flt_, b_, list_, dict_,
    # tup_, gdf_, df_, ds_, arr_)
    pass
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
def <entry>(str_config_file_path,
            str_output_dir,
            b_print_output):
    
    # supress all warnings
    warnings.filterwarnings("ignore", category=UserWarning )
    
    print(" ")
    if b_print_output:
        print("+=================================================================+")
        print("|                       <TITLE>                                   |")
        print("|                Created by Andy Carter, PE of                    |")
        print("|             Center for Water and the Environment                |")
        print("|                 University of Texas at Austin                   |")
        print("+-----------------------------------------------------------------+")
    
        
        print("  ---(c) INPUT GLOBAL CONFIGURATION FILE: " + str_config_file_path)
        print("  ---(o) OUTPUT FOLDER: " + str_output_dir)
        print("  ---[r] PRINT OUTPUT: " + str(b_print_output))
        print("===================================================================")
    else:
        print('Step <N>: <Title>')
    
    # [CONFIG] --- Read variables from config.ini ---
    # Create a ConfigParser object
    config = configparser.ConfigParser()
    
    # Read the config.ini file
    config.read(str_config_file_path)
    
    # Get the variables in the [<NN_section>] section
    if '<NN_section>' in config:
        section = config['<NN_section>']
    
        # Read variables in the section
        flt_<setting> = section.getfloat('flt_<setting>', 0)
    else:
        print("[<NN_section>] section not found in the config file.")
    # ---
    
    str_input_dir = os.path.join(str_output_dir, '<previous NN_ folder>')
    str_out_dir = os.path.join(str_output_dir, '<NN_folder>')
    os.makedirs(str_out_dir, exist_ok=True)
    
    list_input_files = [os.path.join(str_input_dir, f) for f in os.listdir(str_input_dir) if f.endswith(".<ext>")]
    
    # [MP] **********
    # HARD CODED MAX NUMBER OF PROCESSORS
    int_max_processors = 16
    # **********
    
    num_processors = (mp.cpu_count() - 1)
    
    if num_processors > int_max_processors:
        num_processors = int_max_processors
    
    # Use multiprocessing with tqdm
    with mp.Pool(processes=num_processors) as pool:
        list(
            tqdm.tqdm(
                pool.imap(
                    partial(fn_<helper>, str_out_dir=str_out_dir, flt_<setting>=flt_<setting>),
                    list_input_files
                ),
                total=len(list_input_files),
                desc='   -- <Short Desc>',
                bar_format="{desc}:({n_fmt}/{total_fmt})|{bar}| {percentage:.1f}%",
                ncols=75
            )
        )
    # [MP] end -- a serial step calls fn_<helper> in a plain for loop instead
    
    print("+-----------------------------------------------------------------+")
# ++++++++++++++++++++++++++++


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()
    
    parser = argparse.ArgumentParser(description='========= <TITLE> =========')
    
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

    <entry>(str_config_file_path,
            str_output_dir,
            b_print_output)

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)
    
    print('Compute Time: ' + str(time_pass))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
```

Notes taken from his files:
- Comment rules vary by function (`# ----`, `# ++++`, `# ....`, `# ^^^^`); the same rule opens
  and closes each function, two blank lines between functions. `# ++++` usually frames the entry
  function, `# ^^^^` the copied `is_valid_file`.
- The orchestrator calls every step with `b_print_output = False`, so the `else:` one-liner is
  what a full run prints. Standalone runs default `-r` to True and print the banner.
- `is_valid_file` and `fn_str_to_bool` are copied into each script on purpose. Don't import them
  from another module.
- Each `-x` flag: `dest` with a type prefix, `help=r'REQUIRED: ...'` or
  `r'OPTIONAL: ... Default: <value>'`, a `metavar` (`FILE`, `DIR`, `T/F`, `TUPLE`).
