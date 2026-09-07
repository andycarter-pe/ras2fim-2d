# Script 02c - Create scripts to run models assuming that we will be using
# a linux based docker container named "civileng127/ras_v65:v0" to run 
# these models... for TACC ... assumes apptainer SIF at /home1/08140/acarter/ras_v65_v0.sif
#
# Created by: Andy Carter, PE
# Created - 2026.09.07

# ************************************************************
import os
import argparse
import time
import datetime
import configparser
import shutil
import re
from pathlib import Path
import posixpath
import stat
# ************************************************************


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


# ----------------
def fn_build_commands(str_input_folder,
                      str_output_file_win_serial,
                      str_image,
                      str_mount_target):
    
    base = Path(str_input_folder)
    lines = []
    for hdf in sorted(base.rglob("*.p01.tmp.hdf")):
        host_path = hdf.parent          # folder to mount
        filename = hdf.name             # <name>.p01.tmp.hdf
        cmd = (
            f'docker run --rm -i -v "{host_path}:{str_mount_target}" {str_image} '
            f'/bin/bash -c "cd {str_mount_target} && RasUnsteady {filename} x01"'
        )
        lines.append(cmd)
    return lines
# ----------------


# ----------------
def fn_build_windows_serial_bat(str_input_folder,
                                str_output_file_win_serial,
                                str_image,
                                str_mount_target):
    
    lines = fn_build_commands(str_input_folder,str_output_file_win_serial,str_image,str_mount_target)
    header = ["@echo off", "setlocal"]
    body = []
    for i, cmd in enumerate(lines, 1):
        body.append(f'echo [%date% %time%] Running {i} of {len(lines)}')
        body.append(cmd)
        body.append("if errorlevel 1 echo   ^>^> FAILED, continuing...")
        body.append("")
    footer = ["echo Done.", "pause"]

    # .bat files want CRLF line endings
    with open(str_output_file_win_serial, "w", encoding="utf-8", newline="\r\n") as f:
        f.write("\n".join(header + [""] + body + footer) + "\n")
    print(f"Wrote {len(lines)} command(s) to {str_output_file_win_serial}")
# ----------------


# ----------------
def fn_find_parallel_jobs(str_input_folder):
    # walk the tree and return (host_path, filename) for each model
    base = Path(str_input_folder)
    return [(str(hdf.parent), hdf.name)
            for hdf in sorted(base.rglob("*.p01.tmp.hdf"))]
# ----------------


# ----------------
def fn_docker_command_parallel(str_host_path,
                               str_filename,
                               str_image,
                               str_mount_target,
                               str_docker_limits):
    
    str_limits = f"{str_docker_limits} " if str_docker_limits else ""
    # dropped -i (interactive): pointless for unattended parallel runs and it can
    # wedge a backgrounded container waiting on a stdin that never comes.
    return (
        f'docker run --rm {str_limits}-v "{str_host_path}:{str_mount_target}" {str_image} '
        f'/bin/bash -c "cd {str_mount_target} && RasUnsteady {str_filename} x01"'
    )
# ----------------


# ----------------
def fn_build_windows_parallel_bat(str_input_folder,
                                  str_output_file_win_serial,
                                  str_image,
                                  str_mount_target,
                                  int_max_parallel,
                                  int_poll_seconds,
                                  str_docker_limits):
    
    # parallel .bat lives alongside the serial one, in the same folder
    str_output_file_win_parallel = os.path.join(str_input_folder,
                                                'run_docker_windows_parallel.bat')

    list_jobs = fn_find_parallel_jobs(str_input_folder)
    int_total = len(list_jobs)

    list_lines = []
    a = list_lines.append

    a("@echo off")
    a("setlocal EnableDelayedExpansion")
    a("")
    a("REM ===== config (evaluated in BOTH the launcher and each worker) =====")
    a(f'set "MAX_PARALLEL={int_max_parallel}"')
    a(f'set "POLL={int_poll_seconds}"')
    a('set "WORKDIR=%~dp0.par"')
    a('set "LOCKDIR=%WORKDIR%\\locks"')
    a('set "LOGDIR=%WORKDIR%\\logs"')
    a('set "RESULTS=%WORKDIR%\\results.txt"')
    a("")
    a("REM Re-entry: workers relaunch this same .bat with a :runN argument.")
    a('if not "%~1"=="" goto %~1')
    a("")
    a("REM ===== launcher =====")
    a('if not exist "%LOCKDIR%" mkdir "%LOCKDIR%"')
    a('if not exist "%LOGDIR%" mkdir "%LOGDIR%"')
    a('del /q "%LOCKDIR%\\*.lock" >nul 2>&1')
    a(f'set "TOTAL={int_total}"')
    a('>"%RESULTS%" echo Run started %date% %time%  (max %MAX_PARALLEL% at a time, %TOTAL% jobs)')
    a("")
    for i in range(1, int_total + 1):
        a(f"call :launch {i}")
    a("")
    a("echo.")
    a("echo All %TOTAL% jobs launched. Waiting for running containers to finish...")
    a(":drain")
    a("call :count running")
    a("if !running! gtr 0 (")
    a("    timeout /t %POLL% /nobreak >nul")
    a("    goto drain")
    a(")")
    a("echo.")
    a("echo Done. Summary:")
    a('type "%RESULTS%"')
    a("pause")
    a("exit /b")
    a("")
    a("REM ===== helpers =====")
    a(":count")
    a('set "%~1=0"')
    a('for /f %%C in (\'dir /b /a-d "%LOCKDIR%\\*.lock" 2^>nul ^| find /c /v ""\') do set "%~1=%%C"')
    a("exit /b")
    a("")
    a(":launch")
    a('set "i=%~1"')
    a(":throttle")
    a("call :count running")
    a("if !running! geq %MAX_PARALLEL% (")
    a("    timeout /t %POLL% /nobreak >nul")
    a("    goto throttle")
    a(")")
    a('copy /y nul "%LOCKDIR%\\!i!.lock" >nul')
    a("echo [%date% %time%] launching job !i! of %TOTAL%")
    a('start "" /b cmd /c ""%~f0" :run!i!"')
    a("exit /b")
    a("")
    a("REM ===== workers (one block per model) =====")
    for i, (str_host_path, str_filename) in enumerate(list_jobs, 1):
        str_cmd = fn_docker_command_parallel(str_host_path, str_filename,
                                             str_image, str_mount_target, str_docker_limits)
        a(f":run{i}")
        a(f'{str_cmd} > "%LOGDIR%\\run{i}.log" 2>&1')
        a("if errorlevel 1 (")
        a(f'    >>"%RESULTS%" echo [%date% %time%] job {i} FAILED  {str_filename}  ^(see logs\\run{i}.log^)')
        a(") else (")
        a(f'    >>"%RESULTS%" echo [%date% %time%] job {i} OK      {str_filename}')
        a(")")
        a(f'del /q "%LOCKDIR%\\{i}.lock" >nul 2>&1')
        a("exit /b")
        a("")

    # .bat files want CRLF line endings
    with open(str_output_file_win_parallel, "w", encoding="utf-8", newline="\r\n") as f:
        f.write("\n".join(list_lines) + "\n")
    print(f"Wrote {int_total} job(s), max {int_max_parallel} in parallel, "
          f"to {str_output_file_win_parallel}")
# ----------------


# ----------------
def fn_find_tacc_jobs(str_local_base, str_tacc_base):
    # Walk the LOCAL tree; return (tacc_host_dir, filename) for each model.
    base = Path(str_local_base)
    if not base.is_dir():
        raise SystemExit(f"Local base dir does not exist on this machine: {str_local_base}")
    str_tacc_base = str_tacc_base.rstrip("/")
    list_jobs = []
    for hdf in sorted(base.rglob("*.p01.tmp.hdf")):
        str_rel = hdf.relative_to(base).as_posix()      # forward slashes, TACC-side
        str_tacc_full = f"{str_tacc_base}/{str_rel}"
        list_jobs.append((posixpath.dirname(str_tacc_full),
                          posixpath.basename(str_tacc_full)))
    return list_jobs
# ----------------


# ----------------
def fn_bash_single_quote(str_s):
    # Safely single-quote a string for bash.
    return "'" + str_s.replace("'", "'\\''") + "'"
# ----------------


# ----------------
def fn_build_tacc_parallel_sh(str_input_folder,
                              str_output_file_tacc_parallel,
                              str_tacc_base_dir,
                              str_sif,
                              str_mount_target,
                              str_slurm_account,
                              str_slurm_partition,
                              str_slurm_time,
                              str_slurm_jobname,
                              int_default_jobs,
                              int_threads_per_run):

    list_jobs = fn_find_tacc_jobs(str_input_folder, str_tacc_base_dir)
    int_total = len(list_jobs)
    str_hosts = " ".join(fn_bash_single_quote(h) for h, _ in list_jobs)
    str_files = " ".join(fn_bash_single_quote(f) for _, f in list_jobs)
    str_jobs_default = "" if int_default_jobs is None else str(int_default_jobs)

    str_script = f"""#!/bin/bash
#SBATCH -J {str_slurm_jobname}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A {str_slurm_account}
#SBATCH -p {str_slurm_partition}
#SBATCH -t {str_slurm_time}
#SBATCH -o {str_slurm_jobname}.%j.out
#SBATCH -e {str_slurm_jobname}.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL # Email notifications
#SBATCH --mail-user=andy.carter@austin.utexas.edu
#
# Run RasUnsteady models concurrently on ONE node via Apptainer.
#   sbatch run_ras_tacc_parallel.sh     # batch submit
#   idev -N 1 -n 1 -t 2:00:00           # ...then: ./run_ras_tacc_parallel.sh
# Override concurrency:  JOBS=8 ./run_ras_tacc_parallel.sh
#
# DO NOT run this directly on a login node -- it will hammer the shared host.

set -u -o pipefail

# ---- config ---------------------------------------------------------------
SIF={fn_bash_single_quote(str_sif)}
MOUNT_TARGET={fn_bash_single_quote(str_mount_target)}
THREADS_PER_RUN={int_threads_per_run}

CORES="${{SLURM_CPUS_ON_NODE:-$(nproc)}}"
JOBS="${{JOBS:-{str_jobs_default or '$(( CORES / THREADS_PER_RUN ))'}}}"
(( JOBS < 1 )) && JOBS=1

LOGROOT="${{LOGROOT:-$PWD/ras_run_logs}}"
LOGDIR="$LOGROOT/logs"
STATUSDIR="$LOGROOT/status"
RESULTS="$LOGROOT/results.txt"
mkdir -p "$LOGDIR" "$STATUSDIR"

# ---- preflight ------------------------------------------------------------
if ! command -v apptainer >/dev/null 2>&1; then
    module load tacc-apptainer 2>/dev/null || module load apptainer 2>/dev/null || true
fi
command -v apptainer >/dev/null 2>&1 || {{ echo "ERROR: apptainer not found (module load tacc-apptainer)"; exit 1; }}
[[ -f "$SIF" ]] || {{ echo "ERROR: SIF not found: $SIF"; exit 1; }}

TOTAL={int_total}
echo "Node cores: $CORES   Concurrency (JOBS): $JOBS   Models: $TOTAL"
echo "Logs: $LOGROOT"

# ---- model list (index-aligned) -------------------------------------------
HOSTS=( {str_hosts} )
FILES=( {str_files} )

# ---- one model ------------------------------------------------------------
run_one() {{
    local idx="$1" host="$2" file="$3"
    local log="$LOGDIR/run_${{idx}}.log"
    local status="$STATUSDIR/run_${{idx}}.status"
    local ts; ts=$(date '+%F %T')
    if apptainer exec -B "$host:$MOUNT_TARGET" "$SIF" \\
         /bin/bash -c 'cd "$1" && export OMP_NUM_THREADS="$2" && RasUnsteady "$3" x01' \\
         _ "$MOUNT_TARGET" "$THREADS_PER_RUN" "$file" > "$log" 2>&1
    then
        printf '[%s] job %-4s OK             %s\\n' "$ts" "$idx" "$file" > "$status"
    else
        local rc=$?
        printf '[%s] job %-4s FAILED(rc=%s)  %s  (see %s)\\n' "$ts" "$idx" "$rc" "$file" "$log" > "$status"
    fi
}}

# ---- bounded job pool -----------------------------------------------------
for i in "${{!HOSTS[@]}}"; do
    idx=$(( i + 1 ))
    while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n; done
    echo "[$(date '+%F %T')] launching job $idx of $TOTAL"
    run_one "$idx" "${{HOSTS[$i]}}" "${{FILES[$i]}}" &
done
wait

# ---- summary --------------------------------------------------------------
{{
    echo "Run finished $(date '+%F %T')"
    for k in $(seq 1 "$TOTAL"); do cat "$STATUSDIR/run_${{k}}.status" 2>/dev/null; done
}} | tee "$RESULTS"

ok=$(grep -c ' OK ' "$RESULTS" || true)
fail=$(grep -c 'FAILED' "$RESULTS" || true)
echo "Summary: $ok OK, $fail FAILED  (details in $RESULTS)"
"""

    # LF endings are essential -- CRLF makes bash choke with '\r' errors on TACC.
    with open(str_output_file_tacc_parallel, "w", encoding="utf-8", newline="\n") as f:
        f.write(str_script)
    try:  # harmless/no-op on Windows; real effect only matters after copy to TACC
        os.chmod(str_output_file_tacc_parallel, os.stat(str_output_file_tacc_parallel).st_mode | stat.S_IEXEC)
    except OSError:
        pass
    print(f"Wrote {int_total} model(s) to {str_output_file_tacc_parallel}")
    print("Copy it to TACC, then: chmod +x run_ras_tacc_parallel.sh && sbatch run_ras_tacc_parallel.sh")
# ----------------


# ++++++++++++++++++++++++++++
def fn_hec_ras_scripts_for_linux(str_input_folder,
                                 b_print_output):

    print(" ")
    if b_print_output:
        print("+=================================================================+")
        print("|  CREATE SCRIPTS TO RUN HEC_RAS MODELS LINUX DOCKER CONTAINER    |")
        print("+-----------------------------------------------------------------+")
        print("|                Created by Andy Carter, PE of                    |")
        print("|             Center for Water and the Environment                |")
        print("|                 University of Texas at Austin                   |")
        print("+-----------------------------------------------------------------+")
        print("  ---(i) INPUT DIRECTORY OF HEC-RAS RUNS: " + str_input_folder)
        print("  ---[r] PRINT OUTPUT: " + str(b_print_output))
        print("===================================================================")
    else:
        print('Step 2c: Scripts to run HEC-RAS Linux Docker Container')
        

    str_output_file_win_serial = os.path.join(str_input_folder, 'run_docker_windows_serial.bat')
    str_output_file_win_parallel = os.path.join(str_input_folder, 'run_docker_windows_parallel.bat')
    str_output_file_tacc_parallel = os.path.join(str_input_folder, 'run_ras_tacc_parallel.sh')
    
    str_image = "civileng127/ras_v65:v0"
    str_mount_target = "/ras/Linux_RAS_v65/mac-test"
    
    # ************************************************************
    # for the windows parallel batch script
    # Parallel-run tunables (Windows Docker Desktop)
    #   MAX_PARALLEL: containers running at once. Good start:
    #     (CPUs given to Docker Desktop) / (threads each RasUnsteady uses).
    #   POLL_SECONDS: seconds between checks while throttling / draining.
    #   DOCKER_LIMITS: optional per-container caps, e.g. "--cpus=2 --memory=8g"
    #                  (leave "" to disable).
    int_max_parallel = 4
    int_poll_seconds = 5
    str_docker_limits = ""
    # ************************************************************


    # ************************************************************
    # TACC parallel shell script
    str_tacc_base_dir =  "/scratch/08140/acarter/ras2fim_2d_output_20260906/02b_prep_for_ras"
    str_sif = "/home1/08140/acarter/ras_v65_v0.sif"
    
    str_slurm_account = "BCS25094"
    str_slurm_partition = "icx"
    str_slurm_time = "48:00:00"
    str_slurm_jobname = "ras_par"
    
    # Concurrency. Leave DEFAULT_JOBS = None to auto-pick (cores / THREADS_PER_RUN).
    # Override at runtime:  JOBS=8 ./run_ras_parallel.sh
    int_default_jobs = None
    int_threads_per_run = 1
    # ************************************************************

    fn_build_windows_serial_bat(str_input_folder, str_output_file_win_serial,
                                str_image, str_mount_target)
    
    fn_build_windows_parallel_bat(str_input_folder, str_output_file_win_parallel,
                                  str_image, str_mount_target,
                                  int_max_parallel, int_poll_seconds, str_docker_limits)
    
    fn_build_tacc_parallel_sh(str_input_folder, str_output_file_tacc_parallel,
                              str_tacc_base_dir, str_sif, str_mount_target,
                              str_slurm_account, str_slurm_partition,
                              str_slurm_time, str_slurm_jobname,
                              int_default_jobs, int_threads_per_run)
# ++++++++++++++++++++++++++++


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()

    parser = argparse.ArgumentParser(
        description='===== CREATE SCRIPTS TO RUN HEC-RAS LINUX CONTAINERS =====')

    parser.add_argument('-i',
                        dest="str_input_folder",
                        help=r'REQUIRED: folder containing HEC-RAS models '
                             r'Example:E:\ras2fim_2d_output_20260906\02b_prep_for_ras',
                        required=False,
                        default=r'E:\ras2fim_2d_output_20260906\02b_prep_for_ras',
                        metavar='DIR',
                        type=str)

    parser.add_argument('-r',
                        dest="b_print_output",
                        help=r'OPTIONAL: Print output messages Default: True',
                        required=False,
                        default=True,
                        metavar='T/F',
                        type=fn_str_to_bool)

    args = vars(parser.parse_args())

    fn_hec_ras_scripts_for_linux(args['str_input_folder'],
                                 args['b_print_output'])

    flt_end_run = time.time()
    flt_time_pass = (flt_end_run - flt_start_run) // 1
    time_pass = datetime.timedelta(seconds=flt_time_pass)

    print('Compute Time: ' + str(time_pass))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~