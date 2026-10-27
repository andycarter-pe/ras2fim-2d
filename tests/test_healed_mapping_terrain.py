import ast
import json
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import rasterio
from rasterio.transform import Affine
from rasterio.warp import transform_bounds

import healed_mapping_terrain as terrain

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def inputs(tmp_path):
    base = tmp_path / "model.tif"
    data = np.full((64, 64), 1000, "float32")
    data[:, 28:32] = 990
    data[2, 2] = -9999
    transform = Affine(3.2808333333333337, 0, 1815452, 0, -3.2808333333333337, 13934762)
    with rasterio.open(
        base,
        "w",
        driver="GTiff",
        width=64,
        height=64,
        count=1,
        dtype="float32",
        nodata=-9999,
        crs="EPSG:2278",
        transform=transform,
    ) as out:
        out.write(data, 1)
    with rasterio.open(base) as src:
        west, south, east, north = transform_bounds(src.crs, "EPSG:3081", *src.bounds)
    library = tmp_path / "rem"
    tile = library / "15_SAT/dsm/Section-1.tif"
    tile.parent.mkdir(parents=True)
    rem = np.full((40, 40), -9999, "float32")
    rem[10:18, :] = 1000.1 / 3.2808333333333337
    with rasterio.open(
        tile,
        "w",
        driver="GTiff",
        width=40,
        height=40,
        count=1,
        dtype="float32",
        nodata=-9999,
        crs="EPSG:3081",
        transform=Affine(3, 0, np.floor(west) - 12, 0, -3, np.ceil(north) + 12),
    ) as out:
        out.write(rem, 1)
    config = tmp_path / "config.ini"
    config.write_text(
        "[05_raster_engine]\nstr_mapping_terrain=healed\nstr_rem_tiles=rem\n",
        encoding="utf8",
    )
    return config, base, tile, tmp_path / "cache"


def real_builder():
    return terrain._load_builder()


def test_disabled_does_not_load_optional_backend(tmp_path, monkeypatch):
    config = tmp_path / "config.ini"
    config.write_text("[05_raster_engine]\n", encoding="utf8")
    monkeypatch.setattr(
        terrain, "_load_builder", lambda: pytest.fail("unexpected backend import")
    )
    assert terrain.mapping_terrain(config, "model.tif", tmp_path / "cache") == (
        "model.tif",
        None,
    )
    assert not (tmp_path / "cache").exists()


def test_direct_builder_raise_only_cache_and_consistent_mapping(inputs, monkeypatch):
    real_builder()
    config, base, tile, cache = inputs
    before = base.read_bytes()
    path, record = terrain.mapping_terrain(config, base, cache)
    assert record["mapping_only"] and record["approval_state"] == "not_approved"
    assert base.read_bytes() == before
    with rasterio.open(path) as healed, rasterio.open(base) as model:
        a, b = healed.read(1), model.read(1)
        valid = b != model.nodata
        assert np.all(a[valid] >= b[valid]) and np.any(a[valid] > b[valid])
        assert a[2, 2] == -9999
        assert (healed.crs, healed.transform, healed.shape) == (
            model.crs,
            model.transform,
            model.shape,
        )
    qa = json.loads((Path(path).parent / "qa.json").read_text())
    assert qa["datum_check"]["status"] == "pass"
    assert qa["sampling"]["backend"] == "terrio"
    build, select = terrain._load_builder()
    monkeypatch.setattr(
        terrain,
        "_load_builder",
        lambda: (lambda *a, **k: pytest.fail("cache rebuilt"), select),
    )
    assert terrain.mapping_terrain(config, base, cache, allow_build=False) == (
        path,
        record,
    )
    step5, step6 = cache.parent / "run", cache.parent / "simple"
    terrain.bind_mapping_terrain(step5, record)
    terrain.bind_mapping_terrain(step6, terrain.step5_mapping_record(step5))
    terrain.check_netcdf_mapping_inputs(step6, record)
    with pytest.raises(terrain.MappingTerrainError, match="different|unrecorded"):
        terrain.check_netcdf_mapping_inputs(step6, None)
    changed = dict(record, input_sha256="other")
    with pytest.raises(terrain.MappingTerrainError, match="separate output"):
        terrain.bind_mapping_terrain(step5, changed)
    with pytest.raises(terrain.MappingTerrainError, match="separate output"):
        terrain.bind_mapping_terrain(step5, None)


def test_no_road_core_api_produces_unchanged_terrain(inputs):
    real_builder()
    config, base, tile, cache = inputs
    with rasterio.open(tile, "r+") as dst:
        dst.write(np.full(dst.shape, -9999, "float32"), 1)
    path, record = terrain.mapping_terrain(config, base, cache)
    with rasterio.open(path) as healed, rasterio.open(base) as model:
        np.testing.assert_array_equal(healed.read(1), model.read(1))
    assert (
        json.loads((Path(path).parent / "qa.json").read_text())["datum_check"]["status"]
        == "not_evaluable"
    )


def test_failed_datum_is_retained_but_never_reused(inputs):
    real_builder()
    config, base, tile, cache = inputs
    with rasterio.open(tile, "r+") as dst:
        values = dst.read(1)
        values[values != -9999] += 3
        dst.write(values, 1)
    for _ in range(2):
        with pytest.raises(terrain.MappingTerrainError, match="QA failed"):
            terrain.mapping_terrain(config, base, cache)
    failures = list(cache.glob("failed-*"))
    assert len(failures) == 2
    assert all(not (folder / "mapping_terrain.json").exists() for folder in failures)
    assert all(
        json.loads((folder / "qa.json").read_text())["datum_check"]["status"] == "fail"
        for folder in failures
    )


def test_cached_artifact_changes_fail_closed(inputs):
    real_builder()
    config, base, tile, cache = inputs
    path, record = terrain.mapping_terrain(config, base, cache)
    with rasterio.open(path, "r+") as dst:
        values = dst.read(1)
        values[3, 3] += 1
        dst.write(values, 1)
    with pytest.raises(terrain.MappingTerrainError, match="changed"):
        terrain.mapping_terrain(config, base, cache)


def test_step7_does_not_build_a_missing_variant(inputs):
    real_builder()
    config, base, tile, cache = inputs
    with pytest.raises(terrain.MappingTerrainError, match="Run Step 5"):
        terrain.mapping_terrain(config, base, cache, allow_build=False)
    assert not cache.exists()


@pytest.mark.parametrize(
    "setting", ["str_terrain_vertical_unit=m", "flt_datum_tolerance_ft=nan"]
)
def test_invalid_units_or_tolerance_rejected(inputs, setting):
    real_builder()
    config, base, tile, cache = inputs
    with config.open("a") as out:
        out.write(setting + "\n")
    with pytest.raises(terrain.MappingTerrainError):
        terrain.mapping_terrain(config, base, cache)


def test_existing_untracked_products_and_mismatched_step6_are_rejected(tmp_path):
    output = tmp_path / "products"
    output.mkdir()
    (output / "old.tif").write_bytes(b"legacy")
    record = {"input_sha256": "selected"}
    with pytest.raises(
        terrain.MappingTerrainError, match="no healed-terrain provenance"
    ):
        terrain.bind_mapping_terrain(output, record)
    with pytest.raises(terrain.MappingTerrainError, match="unrecorded"):
        terrain.check_netcdf_mapping_inputs(output, record)
    terrain.bind_mapping_terrain(output, None)
    terrain.check_netcdf_mapping_inputs(output, None)


def _hook(filename, function, stop_at):
    tree = ast.parse((ROOT / "src" / filename).read_text(encoding="utf8"))
    node = next(
        n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == function
    )
    # Execute the real entrypoint's prefix through its first heavy processing operation.
    prefix = []
    for statement in node.body:
        if stop_at in ast.unparse(statement):
            break
        prefix.append(statement)
    node.body = prefix
    # The v2 hooks are imported in the file's import block; rerun those imports
    # here so that monkeypatched helpers are the ones bound by the entrypoint.
    hooks = [
        n
        for n in tree.body
        if isinstance(n, ast.ImportFrom)
        and n.module in ("mapping_engine", "healed_mapping_terrain")
    ]
    module = ast.fix_missing_locations(ast.Module(body=hooks + [node], type_ignores=[]))
    namespace = {
        "os": os,
        "warnings": SimpleNamespace(filterwarnings=lambda *a, **k: None),
    }
    exec(compile(module, filename, "exec"), namespace)
    return namespace[function], namespace


def test_entrypoint_hooks_wire_one_surface_through_steps_5_6_7(tmp_path, monkeypatch):
    import mapping_engine

    monkeypatch.setattr(mapping_engine, "mapping_engine", lambda path: "legacy")
    monkeypatch.setattr(mapping_engine, "bind_mapping_engine", lambda *args: None)
    calls = []
    record = {"input_sha256": "selected"}

    def select(config, base, cache, **kwargs):
        calls.append(("select", str(base), str(cache), kwargs))
        return "healed.tif", record

    monkeypatch.setattr(terrain, "mapping_terrain", select)
    monkeypatch.setattr(
        terrain,
        "bind_mapping_terrain",
        lambda out, rec: calls.append(("bind", out, rec)),
    )
    fn, ns = _hook(
        "generate_depth_wsel_rasters_mp_05_v2.py",
        "fn_generate_depth_wsel_rasters",
        "gpd.read_file",
    )
    fn(
        "config",
        "plan",
        "model.tif",
        "hydraulics",
        str(tmp_path / "05_large_terrain/run"),
        False,
    )
    assert calls[0][1] == "model.tif" and calls[1][2] == record
    monkeypatch.setattr(terrain, "step5_mapping_record", lambda path: record)
    fn, ns = _hook(
        "generate_and_simplify_rasters_mp_06_v2.py",
        "fn_extract_and_simplify_rasters",
        "configparser.ConfigParser",
    )
    fn("config", "limits", "step5", "step6", False)
    assert calls[-1] == ("bind", "step6", record)
    monkeypatch.setattr(
        terrain,
        "check_netcdf_mapping_inputs",
        lambda out, rec: calls.append(("check", out, rec)),
    )
    fn, ns = _hook(
        "prep_for_netcdf_mp_07_v2.py", "fn_prepare_netcdf_files", "gpd.read_file"
    )
    fn("config", "model.tif", str(tmp_path), False)
    assert calls[-2][3] == {"allow_build": False}
    assert calls[-1][2] == record


def test_active_runtime_files_select_python_314():
    assert "python=3.14" in (ROOT / "Dockerfile").read_text()
    assert "python=3.14" in (ROOT / "environment_ras2fim2d.yml").read_text()
    for path in (ROOT / "src").glob("*.py"):
        ast.parse(
            path.read_text(encoding="utf8"), filename=str(path), feature_version=(3, 14)
        )


def test_netcdf_metadata_distinguishes_healed_from_model(tmp_path):
    simple = tmp_path / "simple"
    simple.mkdir()
    # Default legacy / model adds no attributes, so v2 NetCDF files match v1.
    assert terrain.netcdf_mapping_metadata(simple / "01_wsel_rasters") == {}
    terrain.bind_mapping_terrain(simple, {"input_sha256": "selected"})
    metadata = terrain.netcdf_mapping_metadata(simple / "01_wsel_rasters")
    assert metadata["11_terrain_kind"] == "healed_mapping_only"
    assert metadata["12_mapping_engine"] == "legacy"
    assert metadata["13_mapping_terrain_sha256"] == "selected"
    assert metadata["14_mapping_terrain_approval"] == "not_approved"
