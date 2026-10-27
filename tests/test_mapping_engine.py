"""Adapter boundaries plus a real Rust mapping smoke test on synthetic cells."""

import json

import geopandas as gpd
import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box

import mapping_engine as engine


def test_selection_defaults_and_native_mode(tmp_path):
    config = tmp_path / "config.ini"
    config.write_text("[05_raster_engine]\n")
    assert engine.mapping_engine(config) == "legacy"
    config.write_text("[05_raster_engine]\nstr_mapping_engine=rasmapper_sloping\n")
    assert engine.mapping_engine(config) == "rasmapper_sloping"
    config.write_text("[05_raster_engine]\nstr_mapping_engine=unknown\n")
    with pytest.raises(ValueError, match="str_mapping_engine"):
        engine.mapping_engine(config)
    # Legacy writes no marker; existing unmarked rasters still block a method change.
    engine.bind_mapping_engine(tmp_path / "output", "legacy")
    assert not (tmp_path / "output" / "mapping_engine.json").exists()
    (tmp_path / "output").mkdir()
    (tmp_path / "output" / "09_wsel.tif").write_bytes(b"")
    with pytest.raises(ValueError, match="no engine provenance"):
        engine.bind_mapping_engine(tmp_path / "output", "fast_step05")
    # A marked non-legacy folder rejects legacy.
    engine.bind_mapping_engine(tmp_path / "marked", "fast_step05")
    with pytest.raises(ValueError, match="engine changed"):
        engine.bind_mapping_engine(tmp_path / "marked", "legacy")


def test_real_rust_mapping_preserves_names_and_uses_selected_ground(
    tmp_path, monkeypatch
):
    from ras_commander import RasFim

    config = tmp_path / "config.ini"
    config.write_text(
        "[05_raster_engine]\nint_mapping_threads=1\n"
        "[06_output_products]\nstr_output_crs=EPSG:3081\nflt_desired_res=1\n"
    )
    cells = gpd.GeoDataFrame(
        {"cell_idx": [0, 1, 2, 3], "wsel_max_1": [301.0] * 4},
        geometry=[
            box(x, y, x + 10, y + 10)
            for x, y in [
                (1000000, 1000000),
                (1000010, 1000000),
                (1000000, 1000010),
                (1000010, 1000010),
            ]
        ],
        crs="EPSG:3081",
    )
    centers = cells[["cell_idx", "geometry"]].copy()
    centers.geometry = centers.centroid

    # Exercise the actual adapter and Rust arithmetic; replace only HDF/GPKG ingestion.
    def from_inputs(gpkg, hdf, terrain, output, **kwargs):
        assert gpkg == "results.gpkg" and hdf == "plan.hdf"
        return RasFim.map_profiles(cells, terrain, output, centers=centers, **kwargs)

    monkeypatch.setattr(RasFim, "map_ras2fim_results", from_inputs)
    outputs = []
    for height, label in [(300.0, "model"), (300.5, "healed")]:
        terrain = tmp_path / (label + ".tif")
        with rasterio.open(
            terrain,
            "w",
            driver="GTiff",
            width=40,
            height=40,
            count=1,
            dtype="float32",
            crs="EPSG:3081",
            nodata=-9999,
            transform=from_origin(999990, 1000030, 1, 1),
        ) as dst:
            dst.write(np.full((40, 40), height, dtype="float32"), 1)
        out = tmp_path / label
        info = engine.map_fast_step05(config, "plan.hdf", terrain, "results.gpkg", out)
        assert info["profiles"] == 1
        assert (out / "11_wsel_max_1_floodplain_ar.gpkg").is_file()
        manifest = json.loads((out / "fim_manifest.json").read_text())
        assert manifest["profiles"][0]["status"] == "complete"
        assert ".fast-step05-" not in (out / "fim_manifest.json").read_text()
        with rasterio.open(out / "09_wsel_max_1_depth_nudge.tif") as depth:
            values = depth.read(1)
            assert np.isfinite(values).any()
            np.testing.assert_allclose(values[np.isfinite(values)], 301 - height)
            outputs.append(np.isfinite(values))
    np.testing.assert_array_equal(*outputs)
    with pytest.raises(ValueError, match="already exist"):
        engine.map_fast_step05(config, "plan.hdf", terrain, "results.gpkg", out)


def test_catalog_packaging_rejects_a_different_engine(tmp_path):
    config = tmp_path / "config.ini"
    config.write_text("[05_raster_engine]\nstr_mapping_engine=rasmapper_sloping\n")
    folder = tmp_path / "step6"
    with pytest.raises(ValueError, match="Step6 mapping engine differs"):
        engine.check_catalog_engine(config, folder)
    engine.bind_mapping_engine(folder, "rasmapper_sloping")
    engine.check_catalog_engine(config, folder)
    assert engine.recorded_mapping_engine(folder) == "rasmapper_sloping"
