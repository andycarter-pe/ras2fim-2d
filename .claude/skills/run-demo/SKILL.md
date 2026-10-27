---
name: run-demo
description: Run the ras2fim-2d sample HEC-RAS 2D model through a chosen step range, in the published Docker image or the local conda env, then check every numbered output folder, read the NetCDF attributes and return a short plain summary (steps run, folders, attributes, compute time, anything missing). Use when Andy asks to run the demo or sample data, or to check that a change still runs. Never edits code.
argument-hint: "docker|docker-src|local \"(start,end)\" [previous-run-folder]"
context: fork
agent: general-purpose
disallowed-tools: Edit Write NotebookEdit
metadata:
  contributed-by: "CLB Engineering for Andy Carter, PE, 2026.10.27"
---

# /run-demo for ras2fim-2d

Run Andy's sample model the way README.md says, check the outputs and report back. You run in
your own context. Keep the long output here and return only the summary in section 7.

Follow AGENTS.md if present. These rules apply even without it:
- Never edit, create or delete anything in the repository: no code, config, README or Dockerfile
  changes. Write only inside the run folder you create.
- Never write into the repo's `sample_output/`. It is the reference copy.
- Do not commit, push, or remove Docker images. Remove any container you started.
- Report exactly what ran. If a step could not run, name it and say why. Never call a step
  passed unless you ran it and checked its folder.

## 1. Read the arguments
Arguments: `$ARGUMENTS`
- Word 1 is the run mode. Andy chooses it every time; there is no default.
  - `docker`: the image named in README.md (today `civileng127/ras2fim2d:v01`), using the code
    baked into the image.
  - `docker-src`: the same image with this working copy's `src/` mounted over
    `/ras2fim-2d/src`, so it runs the code on disk now.
  - `local`: the `ras2fim2d` conda env from `environment_ras2fim2d.yml` on this machine.
- Word 2 is the step range in the same form as `-s`, e.g. `"(0,2)"` or `"(4,8)"`.
- Word 3 (optional) is a previous run folder to compare against.

If the mode or the step range is missing or not one of the above, run nothing. Return only
this question, with the table from section 2:
"Run in Docker (`docker`, or `docker-src` to test your working copy) or in your local conda env
(`local`)? Which steps, as \"(start,end)\"? On Linux, use (0,2) or (4,8). Example:
`/run-demo docker \"(4,8)\"`."

## 2. What each step needs
| Step | Script | Output in the run folder | Runs where |
|---|---|---|---|
| 0 | `copy_source_folder_00.py` | `00_base_hecras/` | Docker or Windows |
| 1 | `create_model_hydrofabric_gpkg_01.py` | `01_model_hydrofabric/model_hydrofabric.gpkg` | Docker or Windows |
| 2 | `spawn_hecras_copies_emit_02.py` | `02_model_copies/` (one folder per run, plus `source_terrain/`, `projection/`) | Docker or Windows |
| 2b, 2c | `stage_hecras_for_linux_v66_02b.py`, `create_hec_ras_scripts_for_linux_02c.py` | `02b_prep_for_ras/` and `run_docker_windows_parallel.bat` | **Windows only**: needs HEC-RAS 6.6 (`RAS66.HECRASController`) and pywin32 |
| 3 | none | prints `** This step is removed on the Linux container **` | n/a |
| HEC-RAS | `02b_prep_for_ras\run_docker_windows_parallel.bat` (4 HEC-RAS 6.5 workers; 02c writes `civileng127/ras_v65:v0`) | `*.p01.tmp.hdf` in `02b_prep_for_ras/` | Windows batch. This skill never runs it |
| 4 | `compute_max_wsel_and_stable_hr_04.py` | copies `*.p01.tmp.hdf` to `03_run_hdf/*.p01.hdf` (overwrites), writes `04_wsel_and_hr/hydrualic_results_*.gpkg` | Docker or Windows |
| 5 | `determine_fim_per_reach_04b.py`, then `generate_depth_wsel_rasters_mp_05.py` | `04_wsel_and_hr/raster_write_plan.csv`, `05_large_terrain/<run>/` | Docker or Windows |
| 6 | `generate_and_simplify_rasters_mp_06.py` | `06_simple_rasters/01_wsel_rasters/`, `01_wsel_rasters_align/` | Docker or Windows |
| 7 | `prep_for_netcdf_mp_07.py` (uses `nccopy`) | `06_simple_rasters/02_wsel_nc/wsel_<wb>.nc` | Docker or Windows (`nccopy` on PATH) |
| 8 | `filter_netcdf_wsel_mp_08.py` | `06_simple_rasters/03_wsel_nc_filtered/filter_wsel_<wb>.nc` | Docker or Windows |

Facts to keep straight (from `src/ras2fim-2d.py`):
- Steps 2b and 2c are imported and run only when `platform.system() == "Windows"`. In Docker,
  step 2 makes `02_model_copies/` and nothing else, so no `.bat` file appears. Say so; it is
  expected, not a failure.
- Steps 4 to 8 need the HEC-RAS results, `01_model_hydrofabric/model_hydrofabric.gpkg` and
  `02_model_copies/source_terrain/`. `ras2fim-2d.py` lists `03_run_hdf/*.hdf` before step 4
  copies `*.p01.tmp.hdf` in from `02b_prep_for_ras/`, so `03_run_hdf/` must already hold the
  `.p01.hdf` for a (4,8) run.
- Steps 2 and 4 each build the flow list from `[02_copy_geom]` in the config you pass. Docker uses
  `config_global.ini` (`flt_min_q = 100`); the sample results were made with
  `config_global_windows.ini` (`flt_min_q = 50`). Flows in a Docker run of (4,8) on the sample
  therefore differ from `sample_output/`. That is the config, not a code fault; say so.

## 3. Get the sample data and make the run folder
- Input model: `s3://rasfim-2d-sample/sample_input` (public; README "Sample Data"). Full
  reference output, 00 through 06, including the HEC-RAS results:
  `s3://rasfim-2d-sample/sample_output`. Download with
  `aws s3 sync s3://rasfim-2d-sample/<prefix> <dest> --no-sign-request`, or with the Python
  fallback in [reference.md](reference.md). Reuse an earlier download if Andy has one.
- The repo's `sample_output/` is a partial copy of that run: `01_model_hydrofabric`,
  `03_run_hdf` (the HEC-RAS result), `04_wsel_and_hr`, `06_simple_rasters` and
  `filter_wsel_wb-2427466.nc`. It has no `00`, `02`, `02b` or `05`, so it cannot seed steps
  4 to 8 by itself. Use its NetCDF files as the reference for the final products.
- Make a new run folder outside the repo: `ras2fim2d_demo_YYYYMMDD_HHMM` next to the repo, unless
  Andy named one. For a range starting at 4, copy `00_base_hecras`, `01_model_hydrofabric`,
  `02_model_copies`, `02b_prep_for_ras` and `03_run_hdf` from the S3 `sample_output` into it,
  and leave out `04` to `06` so the run rebuilds them.
- Global inputs: Docker uses `/global_input` baked into the image (`config_global.ini`). Local
  uses `[global_input_files]` in `src/config_global_windows.ini`
  (`E:\ras2fim-2d-windows-resources\...`). If those files are missing, stop and say so.

## 4. Run
Record the code that runs: `git rev-parse --short HEAD` and `git status --short` in the repo
and, for Docker, `docker run --rm <image> git -C /ras2fim-2d log -1 --format="%h %ad %s" --date=short`.
If the image is not present, `docker pull` it first.

Docker (README Phase 1 and Phase 4, without `-it` because there is no terminal here):
```
docker run --rm --name run-demo-ras2fim2d \
  -v <model_input>:/model_input -v <run_folder>:/model_output \
  civileng127/ras2fim2d:v01 \
  python ras2fim-2d.py -i /model_input -o /model_output -c config_global.ini -f "(1,1,1)" -s "(start,end)"
```
- `docker-src`: add `-e PYTHONDONTWRITEBYTECODE=1 -v <repo>/src:/ras2fim-2d/src:ro` before the
  image name. Read-only, because the repo tracks `src/__pycache__/*.pyc` and Python would
  rewrite them.
- In Git Bash on Windows, put `MSYS_NO_PATHCONV=1` in front so `/model_input` is not rewritten.

Local (README Phase 2 form, on Windows), with `PYTHONDONTWRITEBYTECODE=1` set for the same reason:
```
conda run -n ras2fim2d python <repo>\src\ras2fim-2d.py -i <model_input> -o <run_folder> -c <repo>\src\config_global_windows.ini -f "(1,1,1)" -s "(start,end)"
```

Start the command in the background with stdout and stderr going to
`<run_folder>/run_demo_log.txt`, and time it. Wait for it to finish. Do not read the whole log
into context: tqdm redraws lines with `\r`. Use `tr '\r' '\n' < log | grep -E ...` for
`Step `, `Compute Time:`, `Traceback`, `Error`, `removed on the Linux container`. If the log has
not grown for 20 minutes, stop the container and report it as hung.

After the run, `git status --short` in the repo must match what you recorded before. If it does
not, list the changed files in the summary. Do not restore them yourself.

## 5. Check the outputs
- For every step in the range, confirm its folder from the table exists and is not empty. Give
  the file count and total size. List anything missing.
- For each `.nc` in `02_wsel_nc/` and `03_wsel_nc_filtered/`, read the global attributes
  (`00_stream_id` to `07_raster_folder`, and from step 8 `08_filter_history`,
  `09_vertical_filter`, `10_still_needed_values`), the dimensions and the list of flows. Use the
  snippet in [reference.md](reference.md), run in the same image or env.
- If the repo's `sample_output/` has the same file, compare dimensions, flows and the
  non-date attributes, and say whether they match. Expect different flows from a Docker run
  (see the config note in section 2). For a like-for-like check, compare against an earlier
  run made with the same mode and config.

## 6. Compare and clean up
- If Andy gave a previous run folder and the andy-tools plugin is installed, hand both folders to
  its `results-compare` agent. If it is not installed, say the comparison was skipped and that
  `andy-tools:results-compare` can do it.
- `docker ps -a --filter name=run-demo-` must come back empty. Remove any container you
  started. Keep the run folder; it is the result.

## 7. Return this summary (plain text, short)
```
/run-demo ras2fim-2d -- <mode>, steps (<start>,<end>)
Code:     repo <hash> (<clean|uncommitted changes>); image <name:tag> code <hash date>
Input:    <model_input>     Output: <run_folder>
Ran:      steps <list>      Not run: <steps and why>
Folders:  <folder> <n files, size> OK|MISSING   (one line each)
NetCDF:   <file>: <stream_id>, <n> flows (<min>-<max> cfs), 09_vertical_filter=<v>; matches sample_output: yes|no|n/a
Time:     Compute Time <from log>; wall <measured>
Problems: <tracebacks, empty folders, warnings worth a look, or "none">
Compare:  <results-compare verdict, or "skipped">
Log:      <run_folder>/run_demo_log.txt
```
