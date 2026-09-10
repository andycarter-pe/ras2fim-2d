import ast
import importlib.util
from pathlib import Path

import h5py
import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "model_checks", ROOT / "containers/hecras-prepare/model_checks.py")
checks = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checks)


def make_hdf(path, cells=2, hydraulic=True):
    with h5py.File(path, "w") as hdf:
        area = hdf.create_group("Geometry/2D Flow Areas/Area")
        area.create_dataset("Cells Center Coordinate", data=np.zeros((cells, 2)))
        if hydraulic:
            for key, width in [("Cells Volume Elevation", 2), ("Faces Area Elevation", 4)]:
                area.create_dataset(key + " Info", data=np.array([[0, 2]] * cells))
                area.create_dataset(key + " Values", data=np.ones((2, width)))
    return path


@pytest.mark.parametrize("raw", [b"a\nb\n", b"a\r\nb\r\n", b"a\rb\r", b"a\r\nb\n"])
def test_input_normalization_is_byte_preserving_and_idempotent(tmp_path, raw):
    path = tmp_path / "Model.g01"
    path.write_bytes(b"\xe9" + raw)
    checks.normalize_inputs([path])
    assert path.read_bytes() == b"\xe9a\r\nb\r\n"
    assert checks.normalize_inputs([path]) == []


def test_missing_terrain_raster_fails_before_text_is_modified(tmp_path):
    project = tmp_path / "Model.prj"
    project.write_bytes(b"Proj Title=Test\n")
    project.with_suffix(".p01").write_text("Geom File=g01\nFlow File=u01\n")
    project.with_suffix(".g01").write_bytes(b"Geom Title=Test\n")
    project.with_suffix(".u01").write_text("Flow Title=Test\n")
    project.with_suffix(".rasmap").write_text(
        '<RASMapper><RASProjectionFilename Filename="projection.prj"/>'
        '<Terrains><Layer Type="TerrainLayer" Filename="Terrain.hdf"/></Terrains></RASMapper>')
    (tmp_path / "projection.prj").write_text("projection")
    make_hdf(project.with_suffix(".g01.hdf"))
    with h5py.File(tmp_path / "Terrain.hdf", "w") as hdf:
        hdf.create_group("Terrain/raster").attrs["File"] = np.bytes_("missing.tif")
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    with pytest.raises(ValueError, match="missing.tif"):
        checks.preflight(project, "01")
    assert before == {p.name: p.read_bytes() for p in tmp_path.iterdir()}


def test_successful_output_requires_2d_mesh_in_both_files(tmp_path):
    geometry = make_hdf(tmp_path / "Model.g01.hdf")
    tmp = make_hdf(tmp_path / "Model.p01.tmp.hdf")
    baseline = checks.mesh_summary(geometry)
    checks.validate_outputs(geometry, tmp, baseline)
    with h5py.File(tmp, "a") as hdf:
        del hdf["Geometry/2D Flow Areas"]
        hdf.create_group("Results")  # Results presence does not prove preprocessing.
    with pytest.raises(ValueError, match="missing Geometry/2D Flow Areas"):
        checks.validate_outputs(geometry, tmp, baseline)


def test_empty_table_rows_are_valid_but_invalid_offsets_are_rejected(tmp_path):
    path = make_hdf(tmp_path / "geometry.hdf")
    with h5py.File(path, "a") as hdf:
        hdf["Geometry/2D Flow Areas/Area/Cells Volume Elevation Info"][1] = [2, 0]
    checks.mesh_summary(path, hydraulic=True)
    with h5py.File(path, "a") as hdf:
        hdf["Geometry/2D Flow Areas/Area/Cells Volume Elevation Info"][1] = [2, 1]
    with pytest.raises(ValueError, match="table values"):
        checks.mesh_summary(path, hydraulic=True)


def test_lost_cells_and_missing_hydraulic_tables_are_rejected(tmp_path):
    path = make_hdf(tmp_path / "geometry.hdf")
    with pytest.raises(ValueError, match="cell counts changed"):
        checks.mesh_summary(path, {"Area": {"cells": 3}}, hydraulic=True)
    with h5py.File(path, "a") as hdf:
        del hdf["Geometry/2D Flow Areas/Area/Faces Area Elevation Values"]
    with pytest.raises(ValueError, match="Faces Area Elevation"):
        checks.mesh_summary(path, hydraulic=True)


def test_final_spawn_geometry_write_retains_crlf(tmp_path):
    # Exercise the actual writer without importing unrelated GIS dependencies.
    source = (ROOT / "src/spawn_hecras_copies_emit_02.py").read_text(encoding="utf-8")
    tree = ast.parse(source)
    function = next(node for node in tree.body if isinstance(node, ast.FunctionDef)
                    and node.name == "fn_replace_boundary_lines_in_file")
    scope = {}
    exec(compile(ast.Module(body=[function], type_ignores=[]), "spawn_writer", "exec"), scope)
    path = tmp_path / "Model.g01"
    path.write_bytes(b"a\r\nb\r\nc\r\nd\r\ne\r\nf\r\ng\r\n")
    scope[function.name](path, ["replacement\n"], 0)
    assert path.read_bytes() == b"a\r\nb\r\nreplacement\r\n"
