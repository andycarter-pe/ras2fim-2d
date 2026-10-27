---
name: new-step
description: Create a new numbered ras2fim-2d pipeline step (src/<verb>_<noun>_NN.py, sub-steps like 04b, _mp_ for multiprocessing) or a new _vN version of an existing step, in Andy Carter's house style, then (after asking) wire it into src/ras2fim-2d.py and the -s "(start,end)" step range and add its [NN_section] to both config_global*.ini files. Use when asked to add, create, copy forward, or version a step script.
argument-hint: "<verb_noun> <NN> [purpose]  |  v<N> <existing_step.py> [what changes]"
metadata:
  contributed-by: "CLB Engineering for Andy Carter, PE, 2026.10.27"
---

# /new-step for ras2fim-2d

Follow AGENTS.md if present. These rules apply even without it:
- Match the existing step files; don't "modernize". No classes, type hints, logging module,
  dataclasses or new dependencies unless Andy asks.
- Never renumber, rename or edit an existing step, key or output folder (even misspelled ones like
  `hydrualic_results_` or `CONFIGUTATION`); they are interfaces.
- Leave every line you didn't change byte-for-byte as it was, including trailing spaces on blank
  lines (many of his blank lines carry indentation). In new code, indent blank lines inside a
  block the way the surrounding file does.
- **Ask before editing `src/ras2fim-2d.py` or any config file.** Show the exact lines first.
- Don't commit. Show `git status` and propose a short lowercase commit message; Andy commits.

Two modes. Arguments: `$ARGUMENTS`. If they don't make the mode, name, number and purpose clear,
ask once, in one message, before writing anything.

## The pipeline today (read the real files, this is a map)

| NN | File | Entry function | Writes to `<output>/` |
|---|---|---|---|
| 00 | `copy_source_folder_00.py` | `fn_copy_folder` | `00_base_hecras` |
| 01 | `create_model_hydrofabric_gpkg_01.py` | `fn_create_models_nextgen_gpkg` | `01_model_hydrofabric` |
| 02 | `spawn_hecras_copies_emit_02.py` | `fn_spawn_hecras_copies` | `02_model_copies` |
| 02b | `stage_hecras_for_linux_v66_02b.py` (Windows only) | `fn_prepare_hecras_for_linux` | `02b_prep_for_ras` |
| 02c | `create_hec_ras_scripts_for_linux_02c.py` (Windows only) | `fn_hec_ras_scripts_for_linux` | `02b_prep_for_ras` |
| 03 | (removed; HEC-RAS runs in the Linux container) | | `03_run_hdf` |
| 04 | `compute_max_wsel_and_stable_hr_04.py` | `fn_compute_max_wsel_and_stable_hr` | `04_wsel_and_hr` |
| 04b | `determine_fim_per_reach_04b.py` | `fn_determine_fims_per_reach` | `04_wsel_and_hr/raster_write_plan.csv` |
| 05 | `generate_depth_wsel_rasters_mp_05.py` | `fn_generate_depth_wsel_rasters` | `05_large_terrain` |
| 06 | `generate_and_simplify_rasters_mp_06.py` | `fn_extract_and_simplify_rasters` | `06_simple_rasters` |
| 07 | `prep_for_netcdf_mp_07.py` | `fn_prepare_netcdf_files` | `06_simple_rasters/02_wsel_nc` |
| 08 | `filter_netcdf_wsel_mp_08.py` | `fn_filter_netcdf` | `06_simple_rasters/03_wsel_nc_filtered` |

Entry functions in this repo are `fn_<descriptive_name>` **without** the step number. Their
arguments are prefixed paths, usually `str_config_file_path` first, then inputs, then
`str_output_dir`, and always `b_print_output` last.

## Mode A: a new step

1. **Name it.** `src/<verb>_<noun>_NN.py`; add `_mp_` before the number if it uses
   multiprocessing (`filter_netcdf_wsel_mp_08.py`). After step 08 use `09`. Between two steps use
   a letter sub-step of the earlier one (`05b`), never a renumber. Confirm the name, entry function,
   inputs, output folder (`NN_<name>` under the output directory) and config keys with Andy.
2. **Write it from [step-template.md](step-template.md),** reading the closest existing step as
   well (`filter_netcdf_wsel_mp_08.py` for an `_mp_` step, `copy_source_folder_00.py` or
   `determine_fim_per_reach_04b.py` for a serial one). It must have:
   the header block (`# Script NN - ...`, `Created by: Andy Carter, PE`, `Created - YYYY.MM.DD`);
   the `# ****` import block; a comment rule above and below every function; the banner and
   `---(x) LABEL:` lines inside `if b_print_output:` with `print('Step 9: ...')` in the `else:`
   (no leading zero, as in `'Step 8:'`, `'Step 2b:'`);
   the tqdm `bar_format="{desc}:({n_fmt}/{total_fmt})|{bar}| {percentage:.1f}%"`, `ncols=75` for
   any progress bar; copied `is_valid_file` / `fn_str_to_bool`; the `__main__` argparse block with
   short flags, prefixed `dest`, `REQUIRED:` / `OPTIONAL: ... Default:` help; and the closing
   `Compute Time:` print.
3. Use LF line endings and no tabs, like the other files in `src/`.
4. `python -m py_compile src/<new_file>.py`. Run `python src/<new_file>.py -h` too if the
   `ras2fim2d` env is active.
5. Wire it in (below), after asking.

## Mode B: a new version (`_vN`) of an existing step

1. Pick the name: the old stem plus `_v2` (or the next free N): `generate_depth_wsel_rasters_mp_05.py`
   becomes `generate_depth_wsel_rasters_mp_05_v2.py`. The `_v66_` in
   `stage_hecras_for_linux_v66_02b.py` is the HEC-RAS version, not a revision.
2. Copy forward with `cp` (keeps bytes and line endings). **Never edit, move or delete the old
   file.**
3. In the new file's header, keep `Created by` and the created date, and add a line under the
   last `Created` / `Revised` line: `# Revised - YYYY.MM.DD -- <what changed>`.
4. Make the requested change in the new file only. Keep the entry function name so the switch is
   one import line.
5. `python -m py_compile` the new file; `git status` must show the old file unmodified.
6. Offer the switch, and edit only if Andy says yes:
   ```
   -from generate_depth_wsel_rasters_mp_05 import fn_generate_depth_wsel_rasters
   +from generate_depth_wsel_rasters_mp_05_v2 import fn_generate_depth_wsel_rasters
   ```

## Wiring into src/ras2fim-2d.py (ask first, show the diff)

- **Import:** add `from <module> import <entry>` under `# Import modules`, in step order. A step
  that needs `win32com` goes inside the `if platform.system() == "Windows":` block like 02b/02c.
- **Step range:** `-s "(start,end)"` is parsed by `fn_parse_tuple` into integers, and each step
  runs inside `if int_start_step <= NN and int_end_step >= NN:` under a `# === Step NN` comment
  (step 00 alone uses `== 0`). Add the new block in order. A letter sub-step has no range number
  of its own: call it inside its parent's block, as 02b/02c run inside Step 02 and 04b inside
  Step 05.
- If the entry function takes a path the orchestrator builds, define it before the block as the
  file does for steps 04-06 (`str_output_step_04 = os.path.join(str_output_dir,'04_wsel_and_hr')`).
- Pass `b_print_output` (the orchestrator sets it to `False`).
- A step numbered above 08 is outside the default range: the `-s` argument has
  `default='(0,8)'` and help `Default: "(0,8)"`, and the README says steps `0`-`8`. Ask whether to
  raise the default; don't edit the README unless asked.
- Add `# Revised - YYYY.MM.DD -- add step NN <name>` under the header's `# Updated -` line.
- `python -m py_compile src/ras2fim-2d.py`.

## Config section (ask first)

If the step reads settings, add `[NN_<name>]` to **both** `src/config_global.ini` and
`src/config_global_windows.ini`, in step order. If the `/config-sync` skill is present, use its
`-t section ... -b -` action. Otherwise use exact-text edits and match the neighbours: its own rule line
above and below (`# .................`), the header, a blank line, a `#` comment above each
prefixed key, and two blank lines between sections. Both files are LF and end without a final
newline; keep it that way. Ask for per-file values when a value is a path. Check after:
```
git diff -- src/*.ini
python -c "import configparser,sys;[print(f, dict(p['<NN_section>'])) for f in sys.argv[1:] for p in [configparser.ConfigParser()] if p.read(f)]" src/config_global.ini src/config_global_windows.ini
```

## Finish

Report: files created or changed, the `py_compile` results, the diff of the orchestrator and
configs, what you did not run (most steps need the sample HEC-RAS data or Windows HEC-RAS), and a
proposed commit message such as `add step 09 check netcdf attributes`.
