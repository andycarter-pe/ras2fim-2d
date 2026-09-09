# Script 02b - Preprocess HEC-RAS projects in parallel on a Linux Docker host.
# Windows remains the controller and both hosts use the same shared work folder.

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


_IMAGE_DIGEST = re.compile(r"(?:sha256:|[^\s]+@sha256:)[0-9a-f]{64}")
_REQUIRED_SUFFIXES = (".p01.tmp.hdf", ".b01", ".x01")


def _load_settings(str_config_file_path):
    config = configparser.ConfigParser(interpolation=None)
    if not config.read(str_config_file_path):
        raise ValueError("Cannot read configuration: " + str(str_config_file_path))
    if "03_run_hec_ras" not in config:
        raise ValueError("Missing [03_run_hec_ras] configuration section")

    section = config["03_run_hec_ras"]
    settings = {}
    for name in (
        "linux_host",
        "windows_share",
        "linux_share",
        "prepare_image",
        "prepare_version",
        "wine_profile",
    ):
        value = section.get("str_" + name, "").strip()
        if not value:
            raise ValueError("Missing [03_run_hec_ras] str_" + name)
        settings[name] = value

    if not _IMAGE_DIGEST.fullmatch(settings["prepare_image"]):
        raise ValueError("str_prepare_image must be pinned by SHA-256")
    if not re.fullmatch(r"[0-9]+\.[0-9]+(?:\.[0-9]+)?", settings["prepare_version"]):
        raise ValueError("str_prepare_version must be a dotted HEC-RAS version")
    for name in ("linux_share", "wine_profile"):
        if not PurePosixPath(settings[name]).is_absolute():
            raise ValueError("str_" + name + " must be an absolute Linux path")

    settings["windows_share"] = Path(settings["windows_share"]).resolve(strict=True)
    settings["cores_per_job"] = section.getint("int_prepare_cores_per_job", fallback=2)
    settings["memory_gb"] = section.getint("int_prepare_memory_gb", fallback=6)
    settings["use_ntsync"] = section.getboolean("b_use_ntsync", fallback=True)
    if settings["cores_per_job"] < 1 or settings["memory_gb"] < 1:
        raise ValueError("Preparation CPU and memory limits must be positive")
    return settings


def _remote_path(settings, path):
    local_path = Path(path).resolve(strict=True)
    try:
        relative = local_path.relative_to(settings["windows_share"])
    except ValueError as exc:
        raise ValueError(str(local_path) + " is outside str_windows_share") from exc
    remote = str(PurePosixPath(settings["linux_share"]) / relative.as_posix())
    if "," in remote or "\n" in remote:
        raise ValueError("Docker bind-mount paths cannot contain commas or newlines")
    return remote


def _ssh(settings, arguments, log_path, timeout):
    command = [
        "ssh",
        "-o", "BatchMode=yes",
        "-o", "ConnectTimeout=15",
        "-o", "ServerAliveInterval=30",
        "-o", "ServerAliveCountMax=3",
        settings["linux_host"],
        shlex.join(["docker"] + list(arguments)),
    ]
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        result = subprocess.run(
            command, stdout=log, stderr=subprocess.STDOUT, timeout=timeout
        )
    if result.returncode:
        raise RuntimeError(
            "Container/SSH exit " + str(result.returncode) + "; inspect " + str(log_path)
        )


def _find_projects(str_search_directory):
    search = Path(str_search_directory).resolve(strict=True)
    projects = []
    names = set()
    for project in sorted(search.rglob("*.prj")):
        contents = project.read_text(encoding="utf-8", errors="replace")
        if "Proj Title=" not in contents:
            continue
        if project.stem in names:
            raise ValueError("Duplicate generated project name: " + project.stem)
        if not re.search(r"(?m)^Current Plan=p01\s*$", contents):
            raise ValueError("Expected generated Current Plan=p01: " + str(project))
        if re.findall(r"(?m)^Plan File=(\S+)\s*$", contents) != ["p01"]:
            raise ValueError("Expected one generated p01 plan: " + str(project))

        plan_path = project.with_suffix(".p01")
        plan = plan_path.read_text(encoding="utf-8", errors="replace")
        for reference in ("Geom File=g01", "Flow File=u01"):
            if not re.search(r"(?m)^" + re.escape(reference) + r"\s*$", plan):
                raise ValueError("Generated p01 lacks " + reference + ": " + str(plan_path))
        names.add(project.stem)
        projects.append(project)
    if not projects:
        raise ValueError("No generated HEC-RAS projects found in " + str(search))
    return projects


def _container_command(settings, project, timeout_seconds, run_id):
    model_root = project.parent.parent
    arguments = [
        "run", "--rm",
        "--name", run_id,
        "--network", "none",
        "--read-only",
        "--user", "1000:1000",
        "--cap-drop", "ALL",
        "--security-opt", "no-new-privileges",
        "--cpus", str(settings["cores_per_job"]),
        "--memory", str(settings["memory_gb"]) + "g",
        "--memory-swap", str(settings["memory_gb"]) + "g",
        "--tmpfs", "/tmp:rw,nosuid,nodev,size=512m,mode=1777",
        "--env", "RAS2FIM_JOB_ROOT=/job/project",
        "--mount", "type=bind,src=" + _remote_path(settings, project.parent) + ",dst=/job/project",
        "--mount", "type=bind,src=" + _remote_path(settings, model_root / "source_terrain") + ",dst=/job/source_terrain,readonly",
        "--mount", "type=bind,src=" + _remote_path(settings, model_root / "projection") + ",dst=/job/projection,readonly",
        "--mount", "type=bind,src=" + settings["wine_profile"] + ",dst=/runtime/wine-seed,readonly",
        "--mount", "type=volume,dst=/run/ras-job",
    ]
    if settings["use_ntsync"]:
        arguments.extend(["--device", "/dev/ntsync"])
    arguments.extend([
        settings["prepare_image"],
        "prepare",
        "--project", "/job/project/" + project.name,
        "--plan", "01",
        "--timeout", str(int(timeout_seconds)),
        "--run-id", run_id,
        "--replace-generated",
    ])
    return arguments


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _validate_receipt(settings, project, run_id):
    receipt_path = project.parent / ".ras-commander" / "runs" / run_id / "prepare.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    expected = {
        "schema": "ras-commander-job/v1",
        "run_id": run_id,
        "command": "prepare",
        "status": "succeeded",
        "project": project.name,
        "plan": "01",
        "geometry": "01",
    }
    if any(receipt.get(key) != value for key, value in expected.items()):
        raise ValueError("Receipt does not identify this p01/g01 preparation: " + str(receipt_path))

    runtime = receipt.get("runtime", {})
    if runtime.get("kind") != "wine" or runtime.get("hec_ras_version") != settings["prepare_version"]:
        raise ValueError("Preparation receipt has the wrong HEC-RAS runtime")
    result = receipt.get("result", {})
    if (
        result.get("timed_out") is not False
        or result.get("full_result_copied") is not False
        or result.get("signal_source") not in {"bco", "owned_process_artifacts"}
    ):
        raise ValueError("Preparation receipt does not contain qualified result evidence")

    expected_paths = {project.stem + suffix for suffix in _REQUIRED_SUFFIXES}
    artifacts = receipt.get("artifacts", [])
    if len(artifacts) != 3 or {item.get("path") for item in artifacts} != expected_paths:
        raise ValueError("Preparation receipt does not name the required artifacts")
    for artifact in artifacts:
        path = project.parent / artifact["path"]
        if (
            not path.is_file()
            or path.stat().st_size <= 0
            or path.stat().st_size != artifact.get("size_bytes")
            or _sha256(path) != artifact.get("sha256")
        ):
            raise ValueError("Missing or changed preparation artifact: " + str(path))


def _run_preparation(settings, project, timeout_seconds, run_id):
    log_path = project.parent / ".ras-commander" / "launcher" / (run_id + ".log")
    print("Step 2b: prepare " + project.stem, flush=True)
    _ssh(
        settings,
        _container_command(settings, project, timeout_seconds, run_id),
        log_path,
        int(timeout_seconds) + 300,
    )
    _validate_receipt(settings, project, run_id)


def _stage_projects(projects, str_output_directory):
    output = Path(str_output_directory)
    for project in projects:
        target = output / project.stem
        target.mkdir(parents=True, exist_ok=True)
        for suffix in _REQUIRED_SUFFIXES:
            source = project.parent / (project.stem + suffix)
            shutil.copy2(source, target / source.name)


def fn_prepare_hecras_for_linux(
    str_config_file_path,
    str_search_directory,
    str_output_directory,
    int_processes,
    flt_timeout_sec,
    b_print_output,
):
    """Preprocess the complete batch, then stage it for the compute workers."""
    settings = _load_settings(str_config_file_path)
    projects = _find_projects(str_search_directory)
    workers = max(1, min(int(int_processes), len(projects)))
    output = Path(str_output_directory)
    output.mkdir(parents=True, exist_ok=True)

    if b_print_output:
        print(" ")
        print("+=================================================================+")
        print("|      PREPARE HEC-RAS RUNS ON THE LINUX DOCKER HOST             |")
        print("+-----------------------------------------------------------------+")
        print("  ---(i) INPUT DIRECTORY OF HEC-RAS RUNS: " + str(str_search_directory))
        print("  ---(o) OUTPUT DIRECTORY OF LINUX STAGED FILES: " + str(output))
        print("  ---(p) PARALLEL CONTAINERS: " + str(workers))
        print("  ---(t) PER-PROJECT TIMEOUT (s): " + str(flt_timeout_sec))
        print("===================================================================")
    else:
        print("Step 2b: Prepare HEC-RAS for Linux")

    attempt = "ras2fim-prep-" + uuid.uuid4().hex[:12]
    failures = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(
                _run_preparation,
                settings,
                project,
                flt_timeout_sec,
                attempt + "-" + str(index).zfill(4),
            ): project
            for index, project in enumerate(projects, start=1)
        }
        for future in as_completed(futures):
            project = futures[future]
            try:
                future.result()
            except Exception as exc:
                failures.append(project.stem + ": " + str(exc))

    if failures:
        raise RuntimeError("Step 2b preparation batch failed:\n" + "\n".join(sorted(failures)))

    # Phase barrier: nothing is staged until every receipt and hash is valid.
    _stage_projects(projects, output)
    print(
        "Step 2b: all " + str(len(projects))
        + " projects prepared and staged for the existing compute container",
        flush=True,
    )
