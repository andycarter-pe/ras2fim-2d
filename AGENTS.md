<!-- Contributed by CLB Engineering for Andy Carter, PE, 2026.10.27. Andy owns these rules; edit freely. -->
# Agent rules for ras2fim-2d
Author and maintainer: Andy Carter, PE (Center for Water and the Environment, University of Texas at Austin).
These conventions are deliberate. Match the surrounding code; do not "modernize" it.
If you cannot tell whether something is a convention or a mistake, ask.

## Pipeline
- `src/ras2fim-2d.py` runs steps 0-8. `-s "(start,end)"` picks the range (default `"(0,8)"`); each step is
  called inside `if int_start_step <= N and int_end_step >= N:`.
- Step scripts are `src/<verb>_<noun>_NN.py`. Sub-steps add a letter (`02b`, `02c`, `04b`); `_mp_` marks
  multiprocessing steps (`generate_depth_wsel_rasters_mp_05.py`). The orchestrator imports one entry
  function per step (`from filter_netcdf_wsel_mp_08 import fn_filter_netcdf`).
- Steps 02b/02c need Windows + HEC-RAS 6.6 (imported under `if platform.system() == "Windows":`). Step 03 runs
  HEC-RAS 6.5 Docker workers from `02b_prep_for_ras/run_docker_windows_parallel.bat`, outside `ras2fim-2d.py`.
- Steps hand off through numbered folders under the output dir (`-o`): `00_base_hecras`, `01_model_hydrofabric`,
  `02_model_copies`, `02b_prep_for_ras`, `03_run_hdf`, `04_wsel_and_hr`, `05_large_terrain`, `06_simple_rasters/`
  (`01_wsel_rasters`, `01_wsel_rasters_align`, `02_wsel_nc`, `03_wsel_nc_filtered`). Keep names and layout.
- Every step also runs on its own (entry function with `b_print_output`, plus a `__main__` argparse block).
- A new version of a step is a new file beside the old one (`generate_depth_wsel_rasters_mp_05_v2.py`), with the
  same entry function name and arguments, switched on by changing one import in `ras2fim-2d.py`.
  Keep the old file. Do not edit, archive or delete it unless asked.

## Code style
- File header, then the imports between two `# ****` lines (60 `*`), as in `filter_netcdf_wsel_mp_08.py`:
  ```
  # Script 08 - Filter the created netCDF.
  #
  # Created by: Andy Carter, PE
  # Created - 2025.01.08
  # Revised - 2026.09.07
  ```
  For a meaningful change (not whitespace) append `# Revised - YYYY.MM.DD -- <what changed>`.
  Never remove or edit existing Created/Revised lines or the `Created by` credit.
- Open and close every function with the same comment rule: `# ----------------` above `def` and again after
  the last line. Copy a rule style already in the file (`# ----`, `# ++++`, `# ....`); `is_valid_file` uses
  `# ^^^^`, and the `__main__` block is wrapped in `#~~~~` lines.
- Prefix new names by type: `str_`, `int_`, `flt_`, `b_`, `list_`, `dict_`, `tup_`, `gdf_`, `df_`, `ds_`, `arr_`.
  Functions are `fn_snake_case`. Do not rename existing unprefixed names.
- Plain procedural functions with `#` comments. No classes, type hints, dataclasses, `logging` or shared
  utility modules, and do not convert existing string concatenation to f-strings. The helpers
  `is_valid_file(parser, arg)` and `fn_str_to_bool(value)` are copied into each script; keep copying them.
  Leave Andy's commented-out code in place.
- Output is `print`. With `b_print_output` True, print the banner exactly in this form, then one line per input:
  ```
  print("+=================================================================+")
  print("|                   FILTER NETCDF WSEL FILES                      |")
  print("|                Created by Andy Carter, PE of                    |")
  print("|             Center for Water and the Environment                |")
  print("|                 University of Texas at Austin                   |")
  print("+-----------------------------------------------------------------+")
  print("  ---(o) OUTPUT FOLDER: " + str_output_dir)     # (x) required flag, [x] optional
  print("  ---[r] PRINT OUTPUT: " + str(b_print_output))
  print("===================================================================")
  ```
  Otherwise print one line, `print('Step 8: Filter NetCDF WSEL files')`. End a step with
  `print("+-----------------------------------------------------------------+")`.
- Progress bars: `desc='   -- Filter NetCDFs'`, `bar_format="{desc}:({n_fmt}/{total_fmt})|{bar}| {percentage:.1f}%"`,
  `ncols=75`. Pools use `int_max_processors = 16` and `mp.cpu_count() - 1`. `__main__` starts with
  `flt_start_run = time.time()` and ends with `print('Compute Time: ' + str(time_pass))`.
- argparse: short flags (`-c` config, `-i` input, `-o` output dir, `-r` print output), one argument per call:
  `dest = "str_config_file_path"`, `help=r'REQUIRED: Global configuration filepath Example: C:\...\config_global.ini'`,
  `required=True`, `metavar='FILE'` / `'DIR'` / `'TUPLE'` / `'T/F'`, `type=lambda x: is_valid_file(parser, x)`.
  Optional arguments: `help=r'OPTIONAL: Print output messages Default: True'`, `required=False`, `default=...`.

## Configs
- Two variants, both LF: `src/config_global.ini` (Docker paths, `/global_input/...`) and
  `src/config_global_windows.ini` (`E:\ras2fim-2d-windows-resources\...`). A key change goes into both.
- Sections: `[global_input_files]`, then one per step: `[01_nextgen_hydrofabric]`, `[02_copy_geom]`,
  `[03_run_hec_ras]`, `[04_time_gradient]`, `[05_raster_engine]`, `[06_output_products]`, `[07_prep_netcdf]`,
  `[08_filter_netcdf]`. Each is framed by its own rule (`# ......`, `# >>>>>`, `# ,,,,,`); keep it.
- Keys use the variable prefixes (`flt_q_ratio`, `str_hecras_controller`) with a `#` comment above saying what
  and units. To change a value, comment out the old line and add the new one below it; never delete old values.
- Read keys the way the step already does: `if '08_filter_netcdf' in config:` / `section = config['08_filter_netcdf']`
  / `flt_threshold = section.getfloat('flt_threshold', 0)`. A new optional key gets a fallback that keeps today's
  behaviour, so older config files still run. Preserve comments, blank lines and line endings.

## Do not
- Rename keys, functions, files, folders, gpkg layers or NetCDF attributes, even misspelled ones. They are
  interfaces. Known examples: `hydrualic_results_<run>.gpkg`, `flt_allowed_varianvce`, the banner text
  `INPUT GLOBAL CONFIGUTATION FILE`, and `spawn_hecras_copies_emit_02.py` whose header says "Script 03".
  New NetCDF global attributes take the next number after `10_still_needed_values`.
- Fix a suspected bug you were not asked about (e.g. step 05 checks `'04_time_gradient' in config` before reading
  `[06_output_products]`). Report it instead.
- Reformat, reorder imports, or tidy code you were not asked to change.
- Add dependencies, or change `environment_ras2fim2d.yml`, the Dockerfile or the Python version (3.8.12), without
  asking. New code must run on Python 3.8.
- Commit or push. Stage the change and propose a short lowercase message ("revise step 04"); Andy commits.

## Checking work
- There is no test suite. `python -m py_compile` each edited file, then run the steps on sample data
  (`s3://rasfim-2d-sample`) and look at the numbered outputs and NetCDF attributes.
- Write test runs outside the repo (e.g. `E:\ras2fim2d_test_YYYYMMDD`). `sample_output/` is committed
  reference output; do not overwrite it unless asked.
- Steps 02b, 02c and 03 need Windows and HEC-RAS. Report exactly what you ran, what you saw, and what you could not run.

## Tools in this repo
- Project skills, if any, live in `.claude/skills/<name>/SKILL.md` (type `/` to list them). They follow this file.
