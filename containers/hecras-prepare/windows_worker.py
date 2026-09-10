"""Windows Python worker launched by the Linux controller through Wine."""

import argparse
import hashlib
import importlib.metadata
import json
import os
import tempfile
from pathlib import Path


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_ras_commander(wheel_path, expected_sha256):
    """Verify the retained wheel and the package imported by Windows Python."""
    actual_sha256 = sha256_file(wheel_path)
    if actual_sha256 != expected_sha256:
        raise RuntimeError("retained ras-commander wheel SHA-256 does not match")
    distribution = importlib.metadata.distribution("ras-commander")
    direct_url_text = distribution.read_text("direct_url.json")
    if not direct_url_text:
        raise RuntimeError("ras-commander installation has no wheel provenance")
    archive = json.loads(direct_url_text).get("archive_info")
    if not isinstance(archive, dict):
        raise RuntimeError("ras-commander was not installed from a wheel archive")
    hashes = archive.get("hashes")
    installed_sha256 = hashes.get("sha256") if isinstance(hashes, dict) else None
    if installed_sha256 is None:
        legacy_hash = archive.get("hash")
        if isinstance(legacy_hash, str) and legacy_hash.startswith("sha256="):
            installed_sha256 = legacy_hash.split("=", 1)[1]
    if installed_sha256 != expected_sha256:
        raise RuntimeError("installed ras-commander does not match the retained wheel")
    return {
        "ras_commander_wheel_sha256": actual_sha256,
        "ras_commander_distribution_version": distribution.version,
    }


def write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix="." + path.name + ".", dir=path.parent)
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def run_worker(project, plan, ras_executable, ras_commander_wheel,
               expected_ras_commander_wheel_sha256, timeout, replace_generated):
    """Prepare one plan with the public ras-commander preprocessing API."""
    from model_checks import preflight, normalize_inputs, validate_outputs

    provenance = verify_ras_commander(
        Path(ras_commander_wheel), expected_ras_commander_wheel_sha256
    )
    inputs, geometry_hdf, baseline, dependencies = preflight(project, plan)
    normalized = normalize_inputs(inputs)
    from ras_commander import (
        GeomPreprocessor,
        RasPlan,
        RasPreprocess,
        RasPrj,
        init_ras_project,
    )

    ras_object = RasPrj()
    init_ras_project(
        project,
        ras_version=ras_executable,
        ras_object=ras_object,
        load_results_summary=False,
        load_hdf_metadata=False,
        hide_intro=True,
        accept_tcu=True,
    )
    plan_path = RasPlan.get_plan_path(plan, ras_object=ras_object)
    if plan_path is None:
        raise RuntimeError("Selected plan " + plan + " could not be resolved")
    GeomPreprocessor.clear_geompre_files(plan_path, ras_object=ras_object)
    RasPlan.update_run_flags(
        plan_path, geometry_preprocessor=True, ras_object=ras_object
    )
    result = RasPreprocess.preprocess_plan(
        plan,
        ras_object=ras_object,
        max_wait=timeout,
        clear_existing=replace_generated,
        fix_line_endings=True,
    )
    if not result:
        raise RuntimeError(result.error or "HEC-RAS preprocessing failed")
    validation = validate_outputs(geometry_hdf, result.tmp_hdf_path, baseline)
    return {
        "success": bool(result),
        "plan": result.plan_number,
        "geometry": result.geometry_number,
        "tmp_hdf_path": str(result.tmp_hdf_path) if result.tmp_hdf_path else None,
        "b_file_path": str(result.b_file_path) if result.b_file_path else None,
        "x_file_path": str(result.x_file_path) if result.x_file_path else None,
        "elapsed_seconds": result.elapsed_seconds,
        "signal_source": result.signal_source,
        "full_result_copied": result.full_result_copied,
        "timed_out": result.timed_out,
        "error": result.error,
        "input_preparation": {
            "line_endings": "CRLF",
            "normalized_files": normalized,
            "dependencies": dependencies,
        },
        "hdf_validation": validation,
        **provenance,
    }


def build_parser():
    parser = argparse.ArgumentParser(description="ras-commander Wine worker")
    parser.add_argument("--project", required=True)
    parser.add_argument("--plan", required=True)
    parser.add_argument("--ras-executable", required=True)
    parser.add_argument("--ras-commander-wheel", required=True)
    parser.add_argument("--expected-ras-commander-wheel-sha256", required=True)
    parser.add_argument("--timeout", required=True, type=int)
    parser.add_argument("--result", required=True, type=Path)
    parser.add_argument("--replace-generated", action="store_true")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        payload = run_worker(
            args.project,
            args.plan,
            args.ras_executable,
            args.ras_commander_wheel,
            args.expected_ras_commander_wheel_sha256,
            args.timeout,
            args.replace_generated,
        )
    except Exception as exc:
        payload = {
            "success": False,
            "plan": args.plan,
            "geometry": None,
            "tmp_hdf_path": None,
            "b_file_path": None,
            "x_file_path": None,
            "elapsed_seconds": 0.0,
            "signal_source": None,
            "full_result_copied": False,
            "timed_out": False,
            "error": type(exc).__name__ + ": " + str(exc),
            "ras_commander_wheel_sha256": None,
            "ras_commander_distribution_version": None,
        }
    write_json(args.result, payload)
    return 0 if payload["success"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
