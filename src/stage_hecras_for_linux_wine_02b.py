# Script 02b - Stage HEC-RAS for Linux Runs -- preprocess on a Linux Docker host (Wine)
# Purpose: Optional alternative to stage_hecras_for_linux_v66_02b.py. Instead of
# running HEC-RAS 6.6 through the HECRASController on the Windows machine, each
# generated HEC-RAS project in 02_model_copies is preprocessed (geometry and
# 2D property tables, then the .p01.tmp.hdf) inside a Linux/Wine container on
# a Docker host reached over SSH. Windows remains the controller and both hosts
# see the same shared work folder. Nothing is staged to 02b_prep_for_ras until
# every project has finished and its output files have been checked.
#
# Used only when [03_run_hec_ras] b_prepare_on_linux_wine = True in the config.
# Container: containers/hecras-prepare/ (see its README.md)
#
# Part of ras2fim-2d -- Created by: Andy Carter, PE
# This file contributed by: CLB Engineering Corporation
# Created - 2026.09.09
# Revised - 2026.10.27 -- ras2fim-2d naming (fn_, prefixes), banner, standalone
#                          argparse block; behavior unchanged
# ************************************************************

# ************************************************************
import configparser
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shlex
import shutil
import subprocess
import uuid

from tqdm import tqdm

import argparse
import time
import datetime
# ************************************************************

# The container image must be pinned by its immutable SHA-256 digest
str_image_digest_pattern = r"(?:sha256:|[^\s]+@sha256:)[0-9a-f]{64}"

# Files handed to the existing compute container (same as the v66 step)
tup_required_suffixes = (".p01.tmp.hdf", ".b01", ".x01")


# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
def is_valid_file(parser, arg):
    if not Path(arg).exists():
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


# ----------------
def fn_load_wine_settings(str_config_file_path):
    # Read and check the [03_run_hec_ras] keys for the Linux/Wine preprocessing
    config = configparser.ConfigParser(interpolation=None)
    if not config.read(str_config_file_path):
        raise ValueError("Cannot read configuration: " + str(str_config_file_path))
    if "03_run_hec_ras" not in config:
        raise ValueError("Missing [03_run_hec_ras] configuration section")

    section = config["03_run_hec_ras"]
    dict_settings = {}
    for str_name in ("linux_host",
                     "windows_share",
                     "linux_share",
                     "prepare_image",
                     "prepare_version"):
        str_value = section.get("str_" + str_name, "").strip()
        if not str_value:
            raise ValueError("Missing [03_run_hec_ras] str_" + str_name)
        dict_settings[str_name] = str_value

    if not re.fullmatch(str_image_digest_pattern, dict_settings["prepare_image"]):
        raise ValueError("str_prepare_image must be pinned by SHA-256")
    if not re.fullmatch(r"[0-9]+\.[0-9]+(?:\.[0-9]+)?", dict_settings["prepare_version"]):
        raise ValueError("str_prepare_version must be a dotted HEC-RAS version")

    # Blank = use the HEC-RAS runtime bundled in the image
    dict_settings["wine_profile"] = section.get("str_wine_profile", "").strip()
    for str_name in ("linux_share", "wine_profile"):
        if str_name == "wine_profile" and not dict_settings[str_name]:
            continue
        if not PurePosixPath(dict_settings[str_name]).is_absolute():
            raise ValueError("str_" + str_name + " must be an absolute Linux path")

    dict_settings["windows_share"] = Path(dict_settings["windows_share"]).resolve(strict=True)
    dict_settings["cores_per_job"] = section.getint("int_prepare_cores_per_job", fallback=2)
    dict_settings["memory_gb"] = section.getint("int_prepare_memory_gb", fallback=6)
    dict_settings["use_ntsync"] = section.getboolean("b_use_ntsync", fallback=True)
    if not 1 <= dict_settings["cores_per_job"] <= 8:
        raise ValueError("int_prepare_cores_per_job must be an integer from 1 to 8")
    if dict_settings["memory_gb"] < 1:
        raise ValueError("Preparation memory limit must be positive")
    return dict_settings
# ----------------


# ----------------
def fn_windows_to_linux_path(dict_settings, str_path):
    # Convert a path under str_windows_share to the same path under str_linux_share
    path_local = Path(str_path).resolve(strict=True)
    try:
        path_relative = path_local.relative_to(dict_settings["windows_share"])
    except ValueError as exc:
        raise ValueError(str(path_local) + " is outside str_windows_share") from exc
    str_remote = str(PurePosixPath(dict_settings["linux_share"]) / path_relative.as_posix())
    if "," in str_remote or "\n" in str_remote:
        raise ValueError("Docker bind-mount paths cannot contain commas or newlines")
    return str_remote
# ----------------


# ----------------
def fn_run_docker_over_ssh(dict_settings, list_arguments, path_log, int_timeout_sec):
    # Run one 'docker ...' command on the Linux host; write its output to a log
    list_command = [
        "ssh",
        "-o", "BatchMode=yes",
        "-o", "ConnectTimeout=15",
        "-o", "ServerAliveInterval=30",
        "-o", "ServerAliveCountMax=3",
        dict_settings["linux_host"],
        shlex.join(["docker"] + list(list_arguments)),
    ]
    path_log.parent.mkdir(parents=True, exist_ok=True)
    with path_log.open("w", encoding="utf-8") as log:
        result = subprocess.run(list_command,
                                stdout=log,
                                stderr=subprocess.STDOUT,
                                timeout=int_timeout_sec)
    if result.returncode:
        raise RuntimeError("Container/SSH exit " + str(result.returncode)
                           + "; inspect " + str(path_log))
# ----------------


# ----------------
def fn_find_generated_projects(str_search_directory):
    # Find the step 02 HEC-RAS projects; each must use plan p01 -> g01 / u01
    path_search = Path(str_search_directory).resolve(strict=True)
    list_projects = []
    set_names = set()
    for path_prj in sorted(path_search.rglob("*.prj")):
        str_contents = path_prj.read_text(encoding="utf-8", errors="replace")
        if "Proj Title=" not in str_contents:
            # not a HEC-RAS project file (e.g. a projection .prj)
            continue
        if path_prj.stem in set_names:
            raise ValueError("Duplicate generated project name: " + path_prj.stem)
        if not re.search(r"(?m)^Current Plan=p01\s*$", str_contents):
            raise ValueError("Expected generated Current Plan=p01: " + str(path_prj))
        if re.findall(r"(?m)^Plan File=(\S+)\s*$", str_contents) != ["p01"]:
            raise ValueError("Expected one generated p01 plan: " + str(path_prj))

        path_plan = path_prj.with_suffix(".p01")
        str_plan = path_plan.read_text(encoding="utf-8", errors="replace")
        for str_reference in ("Geom File=g01", "Flow File=u01"):
            if not re.search(r"(?m)^" + re.escape(str_reference) + r"\s*$", str_plan):
                raise ValueError("Generated p01 lacks " + str_reference + ": " + str(path_plan))
        set_names.add(path_prj.stem)
        list_projects.append(path_prj)
    if not list_projects:
        raise ValueError("No generated HEC-RAS projects found in " + str(path_search))
    return list_projects
# ----------------


# ----------------
def fn_build_container_command(dict_settings, path_prj, flt_timeout_sec, str_run_id):
    # Arguments for 'docker run' of one preprocessing job (no network, read-only image)
    path_model_root = path_prj.parent.parent
    list_arguments = [
        "run", "--rm",
        "--name", str_run_id,
        "--network", "none",
        "--read-only",
        "--user", "1000:1000",
        "--cap-drop", "ALL",
        "--security-opt", "no-new-privileges",
        "--cpus", str(dict_settings["cores_per_job"]),
        "--memory", str(dict_settings["memory_gb"]) + "g",
        "--memory-swap", str(dict_settings["memory_gb"]) + "g",
        "--tmpfs", "/tmp:rw,nosuid,nodev,size=512m,mode=1777",
        "--env", "RAS2FIM_JOB_ROOT=/job/project",
        "--mount", "type=bind,src=" + fn_windows_to_linux_path(dict_settings, path_prj.parent) + ",dst=/job/project",
        "--mount", "type=bind,src=" + fn_windows_to_linux_path(dict_settings, path_model_root / "source_terrain") + ",dst=/job/source_terrain,readonly",
        "--mount", "type=bind,src=" + fn_windows_to_linux_path(dict_settings, path_model_root / "projection") + ",dst=/job/projection,readonly",
        "--mount", "type=volume,dst=/run/ras-job",
    ]
    if dict_settings.get("wine_profile"):
        list_arguments.extend([
            "--mount", "type=bind,src=" + dict_settings["wine_profile"] + ",dst=/runtime/wine-seed,readonly",
        ])
    if dict_settings["use_ntsync"]:
        list_arguments.extend(["--device", "/dev/ntsync"])
    list_arguments.extend([
        dict_settings["prepare_image"],
        "prepare",
        "--project", "/job/project/" + path_prj.name,
        "--plan", "01",
        "--num-cores", str(dict_settings["cores_per_job"]),
        "--timeout", str(int(flt_timeout_sec)),
        "--run-id", str_run_id,
        "--replace-generated",
    ])
    return list_arguments
# ----------------


# ----------------
def fn_file_sha256(path_file):
    digest = hashlib.sha256()
    with path_file.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
# ----------------


# ----------------
def fn_validate_receipt(dict_settings, path_prj, str_run_id):
    # Check the container's receipt (prepare.json) and re-hash the three output files
    path_receipt = path_prj.parent / ".ras-commander" / "runs" / str_run_id / "prepare.json"
    dict_receipt = json.loads(path_receipt.read_text(encoding="utf-8"))
    dict_expected = {
        "schema": "ras-commander-job/v1",
        "run_id": str_run_id,
        "command": "prepare",
        "status": "succeeded",
        "project": path_prj.name,
        "plan": "01",
        "geometry": "01",
    }
    if any(dict_receipt.get(key) != value for key, value in dict_expected.items()):
        raise ValueError("Receipt does not identify this p01/g01 preparation: " + str(path_receipt))

    dict_runtime = dict_receipt.get("runtime", {})
    if dict_runtime.get("kind") != "wine" or dict_runtime.get("hec_ras_version") != dict_settings["prepare_version"]:
        raise ValueError("Preparation receipt has the wrong HEC-RAS runtime")
    int_num_cores = dict_receipt.get("arguments", {}).get("num_cores")
    if isinstance(int_num_cores, bool) or not isinstance(int_num_cores, int) or int_num_cores != dict_settings["cores_per_job"]:
        raise ValueError("Preparation receipt has the wrong HEC-RAS core count")
    dict_result = dict_receipt.get("result", {})
    if (dict_result.get("timed_out") is not False
            or dict_result.get("full_result_copied") is not False
            or dict_result.get("signal_source") not in {"bco", "owned_process_artifacts"}):
        raise ValueError("Preparation receipt does not contain qualified result evidence")

    set_expected_paths = {path_prj.stem + str_suffix for str_suffix in tup_required_suffixes}
    list_artifacts = dict_receipt.get("artifacts", [])
    if len(list_artifacts) != 3 or {item.get("path") for item in list_artifacts} != set_expected_paths:
        raise ValueError("Preparation receipt does not name the required artifacts")
    for dict_artifact in list_artifacts:
        path_artifact = path_prj.parent / dict_artifact["path"]
        if (not path_artifact.is_file()
                or path_artifact.stat().st_size <= 0
                or path_artifact.stat().st_size != dict_artifact.get("size_bytes")
                or fn_file_sha256(path_artifact) != dict_artifact.get("sha256")):
            raise ValueError("Missing or changed preparation artifact: " + str(path_artifact))
# ----------------


# ----------------
def fn_run_one_preparation(dict_settings, path_prj, flt_timeout_sec, str_run_id):
    # Preprocess one project in a container, then check its receipt
    path_log = path_prj.parent / ".ras-commander" / "launcher" / (str_run_id + ".log")
    fn_run_docker_over_ssh(dict_settings,
                           fn_build_container_command(dict_settings, path_prj, flt_timeout_sec, str_run_id),
                           path_log,
                           int(flt_timeout_sec) + 300)
    fn_validate_receipt(dict_settings, path_prj, str_run_id)
# ----------------


# ----------------
def fn_stage_projects(list_projects, str_output_directory):
    # Copy .p01.tmp.hdf, .b01 and .x01 for each project to 02b_prep_for_ras/<project>
    path_output = Path(str_output_directory)
    for path_prj in list_projects:
        path_target = path_output / path_prj.stem
        path_target.mkdir(parents=True, exist_ok=True)
        for str_suffix in tup_required_suffixes:
            path_source = path_prj.parent / (path_prj.stem + str_suffix)
            shutil.copy2(path_source, path_target / path_source.name)
# ----------------


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
        print("|        -- PREPROCESS ON A LINUX DOCKER HOST (WINE) --           |")
        print("+-----------------------------------------------------------------+")
        print("|                Created by Andy Carter, PE of                    |")
        print("|             Center for Water and the Environment                |")
        print("|                 University of Texas at Austin                   |")
        print("|       Linux/Wine step contributed by CLB Engineering            |")
        print("+-----------------------------------------------------------------+")
        print("  ---(c) INPUT GLOBAL CONFIGURATION FILE: " + str(str_config_file_path))
        print("  ---(i) INPUT DIRECTORY OF HEC-RAS RUNS: " + str(str_search_directory))
        print("  ---(o) OUTPUT DIRECTORY OF LINUX STAGED FILES: " + str(str_output_directory))
        print("  ---(p) PARALLEL CONTAINERS: " + str(int_processes))
        print("  ---(t) PER-PROJECT TIMEOUT (s): " + str(flt_timeout_sec))
        print("  ---[r] PRINT OUTPUT: " + str(b_print_output))
        print("===================================================================")
    else:
        print('Step 2b: Prepare HEC-RAS for Linux (Linux/Wine host)')

    dict_settings = fn_load_wine_settings(str_config_file_path)
    list_projects = fn_find_generated_projects(str_search_directory)
    l = len(list_projects)
    int_workers = max(1, min(int(int_processes), l))
    path_output = Path(str_output_directory)
    path_output.mkdir(parents=True, exist_ok=True)

    str_attempt = "ras2fim-prep-" + uuid.uuid4().hex[:12]
    list_failures = []
    with ThreadPoolExecutor(max_workers=int_workers) as pool:
        dict_futures = {
            pool.submit(fn_run_one_preparation,
                        dict_settings,
                        path_prj,
                        flt_timeout_sec,
                        str_attempt + "-" + str(int_index).zfill(4)): path_prj
            for int_index, path_prj in enumerate(list_projects, start=1)
        }
        for future in tqdm(as_completed(dict_futures),
                           total=l,
                           desc='Processing HEC-RAS Projects',
                           bar_format="{desc}:({n_fmt}/{total_fmt})|{bar}| {percentage:.1f}%",
                           ncols=75):
            path_prj = dict_futures[future]
            try:
                future.result()
            except Exception as exc:
                list_failures.append(path_prj.stem + ": " + str(exc))

    # --- Report ---
    print("\n-----------------------------------------------------------------")
    print(f"Completed {l - len(list_failures)} / {l} projects successfully.")
    if list_failures:
        print("The following projects did not complete:")
        for str_msg in sorted(list_failures):
            print("  - " + str_msg)
    print("-----------------------------------------------------------------")

    if list_failures:
        raise RuntimeError("Step 2b preparation batch failed:\n" + "\n".join(sorted(list_failures)))

    # Nothing is staged until every receipt and hash is valid
    fn_stage_projects(list_projects, path_output)
    print("Step 2b: all " + str(l) + " projects prepared and staged for the existing compute container")
# ++++++++++++++++++++++++++++


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    flt_start_run = time.time()

    parser = argparse.ArgumentParser(
        description='======== PREPARE HEC-RAS RUNS FOR LINUX SIMULATION (LINUX/WINE HOST) ========')

    parser.add_argument('-c',
                        dest="str_config_file_path",
                        help=r'REQUIRED: Global configuration filepath with the [03_run_hec_ras] Linux/Wine keys '
                             r'Example: C:\Users\civil\ras2fim-2d\src\config_global_windows.ini',
                        required=True,
                        metavar='FILE',
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument('-i',
                        dest="str_input_folder",
                        help=r'REQUIRED: folder containing HEC-RAS models (on the shared work folder) '
                             r'Example: E:\ras2fim2d_test_20250108\02_model_copies',
                        required=True,
                        metavar='DIR',
                        type=str)

    parser.add_argument('-o',
                        dest="str_output_dir",
                        help=r'REQUIRED: directory to write Linux staged models '
                             r'Example: E:\ras2fim2d_test_20250108\02b_prep_for_ras',
                        required=True,
                        metavar='DIR',
                        type=str)

    parser.add_argument('-p',
                        dest="int_processes",
                        help=r'OPTIONAL: number of preprocessing containers to run at once. Default: 4',
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
