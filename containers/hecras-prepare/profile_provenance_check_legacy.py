"""Metadata/import check only. Never launch HEC-RAS or open a model."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import sys
from pathlib import Path

EXPECTED_WHEEL = "dc0c2f9baa9db66afee01e34eb00d7a04f53e7b82c3340c786b2e9aef1795932"
EXPECTED_RAS = "b23bd359f47e2a869a5b931a98b461c43978d86c6250894e88a09b21f8aae99d"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


wheel = Path(r"C:\ras2fim-runtime\ras_commander-0.99.2-py3-none-any.whl")
ras_exe = Path(r"C:\Program Files (x86)\HEC\HEC-RAS\6.5\Ras.exe")
assert sha256(wheel) == EXPECTED_WHEEL, "Retained wheel mismatch"
assert sha256(ras_exe) == EXPECTED_RAS, "Ras.exe mismatch"
distribution = importlib.metadata.distribution("ras-commander")
assert distribution.version == "0.99.2", "Distribution version mismatch"
direct_url = json.loads(distribution.read_text("direct_url.json") or "{}")
assert not direct_url.get("dir_info", {}).get("editable", False), "Editable installation remains"
archive_info = direct_url.get("archive_info", {})
installed_sha = archive_info.get("hashes", {}).get("sha256")
if installed_sha is None:
    installed_sha = archive_info.get("hash", "").removeprefix("sha256=")
assert installed_sha == EXPECTED_WHEEL, "Installed wheel provenance mismatch"
import ras_commander
from ras_commander import RasCmdr, RasPlan, init_ras_project

module_path = Path(ras_commander.__file__).resolve()
module_path.relative_to(Path(r"C:\Python311\Lib\site-packages").resolve())
assert ras_commander.__version__ == "0.99.2", "Imported version mismatch"
print(json.dumps({
    "check": "windows-runtime-provenance-and-import",
    "success": True,
    "python_version": sys.version,
    "python_executable": sys.executable,
    "ras_commander_version": distribution.version,
    "ras_commander_module": str(module_path),
    "direct_url": direct_url.get("url"),
    "installed_archive_sha256": installed_sha,
    "retained_wheel_sha256": EXPECTED_WHEEL,
    "ras_exe_sha256": EXPECTED_RAS,
    "api_imports": ["RasCmdr", "RasPlan", "init_ras_project"],
    "hec_ras_executed": False,
}, indent=2))
