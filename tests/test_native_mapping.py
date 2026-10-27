"""Preserve Step4's window semantics, complete identity and catalog delivery CRS."""

import json
from pathlib import Path

import geopandas as gpd
import h5py
import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box

import native_mapping as native


@pytest.fixture
def catalog(tmp_path, monkeypatch):
    from ras_commander import HdfMesh
    from compute_max_wsel_and_stable_hr_04 import (
        fn_get_gradient_array_wet_cells,
        fn_create_run_stable_hour_and_max_wsel_v3,
    )

    hdf = tmp_path / "plan.hdf"
    values = np.full((8, 5), 100.0, dtype="float32")
    values[:, :3] = [[101.234, 102.0, 103.0]]
    values[2, 0] = 101.456
    values[3, 0] = 500.0  # Step4's existing exclusive window stop omits this sample.
    values[4:7, 0] = 104.321
    values[:, 4] = 110.123  # Ghost cell is supplied from this same window, not Max.
    with h5py.File(hdf, "w") as file:
        file.create_dataset(
            "Geometry/2D Flow Areas/M/Cells Minimum Elevation", data=np.full(5, 100.0)
        )
        file.create_dataset(
            "Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/M/Water Surface",
            data=values,
        )
        file["Results/Unsteady"].attrs["Plan Title"] = np.bytes_(
            "123_wb-1_wb-2_4-hr_100-cfs_to_200-cfs_plan"
        )
    keys = pd.DataFrame({"mesh_name": ["M"] * 5, "cell_id": [4, 1, 3, 0, 2]})
    monkeypatch.setattr(HdfMesh, "get_mesh_area_names", lambda path: ["M"])
    monkeypatch.setattr(
        HdfMesh, "get_mesh_cell_points", lambda *args, **kwargs: keys.copy()
    )
    ids, _, rounded = fn_get_gradient_array_wet_cells(hdf, 1)
    reference = fn_create_run_stable_hour_and_max_wsel_v3(4, 1.0, rounded, 1, ids)
    cells = gpd.GeoDataFrame(
        reference.set_index("cell_idx").loc[[2, 0]].reset_index(),
        geometry=[
            box(1000000, 1000000, 1000010, 1000010),
            box(1000010, 1000000, 1000020, 1000010),
        ],
        crs="EPSG:3081",
    )
    gpkg = tmp_path / "results.gpkg"
    cells.to_file(gpkg, layer="02_cells_wsel_ar")
    config = tmp_path / "config.ini"
    config.write_text(
        "[02_copy_geom]\nint_time_rolling_avg=1\n"
        "[05_raster_engine]\nstr_mapping_engine=rasmapper_sloping\nint_mapping_threads=1\n"
        "[06_output_products]\nstr_output_crs=EPSG:3857\nflt_desired_res=3\n"
    )
    return config, hdf, gpkg, reference, cells


def test_complete_profiles_match_existing_step4_including_rounding_ghosts_and_dry(
    catalog,
):
    config, hdf, gpkg, reference, _ = catalog
    table, _, record = native.prepare_catalog_profiles(config, hdf, gpkg)
    assert table.cell_id.tolist() == list(range(5))
    assert table.loc[3, ["wsel_max_1", "wsel_max_2"]].isna().all()
    assert table.loc[4, "wsel_max_1"] == pytest.approx(110.12)
    assert table.loc[0, "wsel_max_1"] == pytest.approx(101.46)
    assert table.loc[0, "wsel_max_2"] == pytest.approx(104.32)
    for profile in ("wsel_max_1", "wsel_max_2"):
        np.testing.assert_array_equal(
            table.set_index("cell_id").loc[reference.cell_idx, profile],
            reference[profile],
        )
    assert record["windows"][0]["sample_stop_exclusive"] == 3
    assert record["never_wet_cells_explicitly_dry"] == 1
    assert table.attrs["fim_source_hdf"] == str(hdf.resolve())


def test_stale_catalog_and_multiple_meshes_rejected(catalog, monkeypatch):
    from ras_commander import HdfMesh

    config, hdf, gpkg, _, cells = catalog
    cells.loc[0, "wsel_max_1"] += 1
    cells.to_file(gpkg, layer="02_cells_wsel_ar")
    with pytest.raises(ValueError, match="disagrees"):
        native.prepare_catalog_profiles(config, hdf, gpkg)
    monkeypatch.setattr(HdfMesh, "get_mesh_area_names", lambda path: ["M", "Other"])
    with pytest.raises(ValueError, match="one mesh"):
        native.prepare_catalog_profiles(config, hdf, gpkg)


def test_native_dispatch_passes_complete_profiles_and_delivers_correct_crs(
    catalog, monkeypatch, tmp_path
):
    from ras_commander import RasFim
    from mapping_engine import map_catalog

    config, hdf, gpkg, _, _ = catalog
    terrain = tmp_path / "terrain.tif"
    with rasterio.open(
        terrain,
        "w",
        driver="GTiff",
        count=1,
        dtype="float32",
        width=20,
        height=10,
        crs="EPSG:3081",
        nodata=-9999,
        transform=from_origin(1000000, 1000010, 1, 1),
    ) as dst:
        dst.write(np.full((10, 20), 100.0, "float32"), 1)

    def render(source, ground, folder, **kwargs):
        assert source == hdf and ground == terrain
        assert "profiles" not in kwargs
        supplied = kwargs["wsel_profiles"]
        assert len(supplied) == 5 and np.isfinite(supplied.loc[4, "wsel_max_1"])
        assert len(kwargs["polygon"]) == 2
        folder.mkdir()
        rows = []
        for profile in ("wsel_max_1", "wsel_max_2"):
            paths = {}
            for prefix, suffix, elevation in (
                ("09", "depth", 1.0),
                ("10", "wsel", 101.0),
            ):
                path = folder / f"{prefix}_{profile}_{suffix}_nudge.tif"
                with rasterio.open(ground) as src:
                    metadata = src.profile
                with rasterio.open(path, "w", **metadata) as dst:
                    dst.write(np.full((10, 20), elevation, "float32"), 1)
                paths["depth_tif" if prefix == "09" else "wse_tif"] = str(path)
            rows.append({"profile": profile, "status": "complete", **paths})
        result = pd.DataFrame(rows)
        result.attrs["profile_source"] = {"kind": "caller-cell-wsel/1"}
        (folder / "fim_manifest.json").write_text(json.dumps({"profiles": rows}))
        return result

    monkeypatch.setattr(RasFim, "map_sloping_profiles", render)
    before = hdf.read_bytes(), terrain.read_bytes()
    output = tmp_path / "mapped"
    result = map_catalog(config, hdf, terrain, gpkg, output)
    assert result == {"profiles": 2, "engine": "rasmapper_sloping"}
    assert (hdf.read_bytes(), terrain.read_bytes()) == before
    record = json.loads((output / "mapping_catalog.json").read_text())
    assert record["native_profile_source"]["kind"] == "caller-cell-wsel/1"
    assert ".native-catalog-" not in (output / "mapping_catalog.json").read_text()
    for profile in ("wsel_max_1", "wsel_max_2"):
        with rasterio.open(output / f"09_{profile}_depth_nudge.tif") as depth:
            assert depth.crs.to_epsg() == 3857 and depth.res == (3.0, 3.0)
            assert np.isfinite(depth.read(1)).any()
        with rasterio.open(
            output / "native" / f"09_{profile}_depth_nudge.tif"
        ) as depth:
            assert depth.crs.to_epsg() == 3081 and depth.res == (1.0, 1.0)
        limits = gpd.read_file(output / f"11_{profile}_floodplain_ar.gpkg")
        assert limits.crs.to_epsg() == 3857
        assert Path(record["products"][0]["depth_tif"]).is_file()
