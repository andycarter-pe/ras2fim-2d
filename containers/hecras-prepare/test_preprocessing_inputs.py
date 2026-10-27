#!/usr/bin/env python3
"""Run disposable LF/CRLF/mixed and missing-mount regression cases with Docker.

Run this same script from Windows or Linux. Supply a spawned sample project
whose parent folder has sibling source_terrain and projection directories.
All model copies and evidence must be kept outside Git.
"""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time


def snapshot(folder):
    return {str(p.relative_to(folder)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in folder.rglob("*") if p.is_file() and ".ras-commander" not in p.parts}


def run_cases(image, source, work, container_user=None):
    source = Path(source).resolve(strict=True)
    work = Path(work).resolve()
    repo = Path(__file__).resolve().parents[2]
    if work == repo or repo in work.parents or work == source.parent or source.parent in work.parents:
        raise ValueError("Use a new external work directory, separate from the source model")
    work.mkdir(parents=True, exist_ok=False)
    inspected = subprocess.run(["docker", "image", "inspect", image], check=True,
                               capture_output=True, text=True)
    image_id = json.loads(inspected.stdout)[0]["Id"]
    (work / "image.json").write_text(inspected.stdout, encoding="utf-8")
    records = []
    for case in ("lf", "crlf", "mixed", "missing-dependencies"):
        folder = work / case / "project"
        shutil.copytree(source.parent, folder, copy_function=shutil.copyfile,
                        ignore=shutil.ignore_patterns(".ras-commander"))
        # Disposable copies must be writable by the image's default UID 1000.
        folder.chmod(0o777)
        for path in folder.rglob("*"):
            path.chmod(0o777 if path.is_dir() else 0o666)
            if path.suffix.lower() in (".prj", ".p01", ".g01", ".u01"):
                lines = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n").split(b"\n")
                data = b"\n".join(lines)
                if case == "crlf":
                    data = b"\r\n".join(lines)
                elif case == "mixed":
                    data = b"".join(line + (b"\r\n" if n % 2 else b"\n")
                                    for n, line in enumerate(lines[:-1])) + lines[-1]
                path.write_bytes(data)
        before = snapshot(folder)
        cmd = ["docker", "run", "--rm"]
        if container_user is not None:
            cmd += ["--user", container_user]
        elif case != "crlf":
            cmd += ["--user", "root"]
        cmd += ["--mount", "type=bind,src=" + str(folder) + ",dst=/job"]
        if case != "missing-dependencies":
            for dependency in ("projection", "source_terrain"):
                cmd += ["--mount", "type=bind,src=" + str(source.parent.parent / dependency)
                        + ",dst=/" + dependency + ",readonly"]
        cmd += [image_id, "prepare", "--project", "/job/" + source.name, "--plan", "01",
                "--replace-generated", "--timeout", "180", "--run-id", case]
        (folder.parent / "command.json").write_text(json.dumps(cmd, indent=2), encoding="utf-8")
        started = time.monotonic()
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=330)
        receipt = json.loads((folder / ".ras-commander/runs" / case / "prepare.json").read_text())
        record = {"case": case, "exit_code": result.returncode,
                  "seconds": round(time.monotonic() - started, 2), "receipt": receipt,
                  "stdout": result.stdout, "stderr": result.stderr}
        if case == "missing-dependencies":
            record["passed"] = (result.returncode != 0 and receipt["status"] == "failed"
                                and "Missing or empty model input" in receipt["error"]["message"]
                                and before == snapshot(folder))
        else:
            validation = receipt.get("result", {}).get("hdf_validation", {})
            record["passed"] = (result.returncode == 0 and receipt["status"] == "succeeded"
                                and bool(validation.get("geometry"))
                                and bool(validation.get("temporary_plan")))
        records.append(record)
        (work / "results.json").write_text(json.dumps(records, indent=2), encoding="utf-8")
        print(case, "PASS" if record["passed"] else "FAIL", flush=True)
    return all(record["passed"] for record in records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", required=True)
    parser.add_argument("--source-project", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--container-user", help="Docker user for every case; use root for Windows drive mounts")
    args = parser.parse_args()
    raise SystemExit(0 if run_cases(args.image, args.source_project, args.work_dir,
                                  args.container_user) else 1)
