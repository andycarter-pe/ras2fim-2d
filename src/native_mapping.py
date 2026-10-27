"""Map existing discharge windows through RasCommander's native sloping engine.

No hydraulic inputs are modified. Elevations retain the source model units/datum.
Native rasters remain on the selected terrain grid; delivery copies use the catalog
CRS/resolution. The consumer's existing first-mesh analysis limits this adapter to
one mesh. Complete topology/WSEL values remain available outside the output mask.
"""

import shutil
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from healed_mapping_terrain import MappingTerrainError, _file, _settings, _write

COMMANDER_COMMIT = "777c7746ceb35fc3b1c29d3718935c9013ae009c"


def prepare_catalog_profiles(config_path, hdf_path, results_gpkg):
    """Return complete named cell WSEL profiles, output mask, and window provenance.

    Reuse Step4's rounded WSEL and rounded-depth dry mask before spatial clipping.
    Cells proved never wet by that calculation are explicitly dry, including ghosts.
    Window endpoints intentionally match Step4's existing exclusive slice end.
    Reject mismatched/stale retained-cell results rather than mixing model runs.
    """
    import geopandas as gpd
    from ras_commander import HdfMesh
    from compute_max_wsel_and_stable_hr_04 import (
        fn_get_gradient_array_wet_cells,
        fn_get_plan_title_from_hdf,
        fn_parse_run_name,
    )

    names = HdfMesh.get_mesh_area_names(hdf_path)
    if len(names) != 1:
        raise MappingTerrainError(
            "Native catalog mapping requires one mesh: existing Step4 analyzes only the first mesh"
        )
    inputs = {"hydraulic_hdf": _file(hdf_path), "catalog_results": _file(results_gpkg)}
    cells = gpd.read_file(results_gpkg, layer="02_cells_wsel_ar")
    if cells.empty or cells.crs is None or cells.cell_idx.duplicated().any():
        raise MappingTerrainError(
            "Catalog cells must be nonempty, unique and have a CRS"
        )
    if not pd.api.types.is_integer_dtype(cells.cell_idx):
        raise MappingTerrainError("Catalog cell_idx must contain integer HDF cell IDs")
    config = _settings(config_path)
    gradient = config.getint("02_copy_geom", "int_time_rolling_avg")
    run_name = fn_get_plan_title_from_hdf(hdf_path)
    count = fn_parse_run_name(run_name)["hour_count"]
    if count < 2:
        raise MappingTerrainError("Catalog windows require at least two stored samples")
    wet_ids, _, rounded_wsel = fn_get_gradient_array_wet_cells(hdf_path, gradient)
    keys = HdfMesh.get_mesh_cell_points(
        hdf_path, include_ghost_cells=True, include_geometry=False
    )[["mesh_name", "cell_id"]].copy()
    if (
        keys.empty
        or keys.mesh_name.nunique() != 1
        or keys.mesh_name.iloc[0] != names[0]
    ):
        raise MappingTerrainError(
            "Complete cell identity does not match the analyzed mesh"
        )
    if keys.duplicated().any() or not set(wet_ids).issubset(set(keys.cell_id)):
        raise MappingTerrainError(
            "Wet-cell identities do not match the complete HDF mesh"
        )
    keys = keys.set_index("cell_id", drop=False).sort_index()
    windows = []
    profiles = [f"wsel_max_{i + 1}" for i in range(rounded_wsel.shape[1] // count)]
    stored = {c for c in cells.columns if c.startswith("wsel")}
    if not profiles or stored != set(profiles):
        raise MappingTerrainError(
            "Catalog profiles do not match the HDF discharge windows; rerun Step4"
        )
    if not set(cells.cell_idx).issubset(set(keys.cell_id)):
        raise MappingTerrainError("Catalog contains unknown HDF cell IDs")
    for i, profile in enumerate(profiles):
        start, stop = i * count, (i + 1) * count - 1
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            maxima = np.nanmax(rounded_wsel[:, start:stop], axis=1)
        keys[profile] = np.full(len(keys), np.nan, dtype="float32")
        keys.loc[wet_ids, profile] = maxima
        expected = keys.loc[cells.cell_idx, profile].to_numpy(dtype="float32")
        supplied = cells[profile].to_numpy(dtype="float32")
        if not np.array_equal(expected, supplied, equal_nan=True):
            raise MappingTerrainError(
                f"{profile} disagrees with Step4 results; rerun Step4 for this HDF"
            )
        windows.append(
            {"profile": profile, "sample_start": start, "sample_stop_exclusive": stop}
        )
    keys = keys.reset_index(drop=True)
    keys.attrs.update(fim_source_hdf=str(Path(hdf_path).resolve()), crs=cells.crs)
    provenance = {
        "schema": "ras2fim-2d-catalog-window-profiles/v1",
        "run_name": run_name,
        "mesh": names[0],
        "complete_cells": len(keys),
        "never_wet_cells_explicitly_dry": len(keys) - len(wet_ids),
        "policy": "Step4 two-decimal WSEL; rounded-depth zero is dry; existing exclusive window end",
        "windows": windows,
        **inputs,
        "mapping_terrain_usage": "independent TIFF; source hydraulic geometry/properties retained",
    }
    if inputs != {
        "hydraulic_hdf": _file(hdf_path),
        "catalog_results": _file(results_gpkg),
    }:
        raise MappingTerrainError(
            "Catalog source inputs changed during profile preparation"
        )
    return keys, cells, provenance


def _deliver_raster(source, target, crs, resolution):
    """Nearest delivery warp after native reconstruction; cover the aligned output grid."""
    import rasterio
    from rasterio.warp import (
        Resampling,
        aligned_target,
        calculate_default_transform,
        reproject,
    )

    with rasterio.open(source) as src:
        transform, width, height = calculate_default_transform(
            src.crs, crs, src.width, src.height, *src.bounds, resolution=resolution
        )
        transform, width, height = aligned_target(transform, width, height, resolution)
        metadata = src.profile.copy()
        if not metadata.get("tiled"):
            metadata.pop("blockxsize", None)
            metadata.pop("blockysize", None)
        metadata.update(
            crs=crs,
            transform=transform,
            width=width,
            height=height,
            dtype="float32",
            nodata=np.nan,
            compress="lzw",
        )
        with rasterio.open(target, "w", **metadata) as dst:
            reproject(
                rasterio.band(src, 1),
                rasterio.band(dst, 1),
                src_nodata=src.nodata,
                dst_nodata=np.nan,
                resampling=Resampling.nearest,
                num_threads=1,
            )


def _deliver_profile(row, folder, crs, resolution):
    """Keep existing 09/10/11 names and floodplain depth limits for Steps6-7."""
    import geopandas as gpd
    import rasterio
    from rasterio.features import shapes
    from shapely import union_all
    from shapely.geometry import shape

    paths = {}
    for column in ("depth_tif", "wse_tif"):
        target = folder / Path(row[column]).name
        _deliver_raster(row[column], target, crs, resolution)
        paths[column] = str(target)
    with rasterio.open(paths["depth_tif"]) as src:
        depth = src.read(1)
        wet = (depth > 0) & (depth <= 10000)
        geometries = [
            shape(g)
            for g, _ in shapes(wet.astype("uint8"), mask=wet, transform=src.transform)
        ]
        if geometries:
            target = folder / f"11_{row['profile']}_floodplain_ar.gpkg"
            gpd.GeoDataFrame(geometry=[union_all(geometries)], crs=src.crs).to_file(
                target
            )
            paths["floodplain_gpkg"] = str(target)
    return {"profile": row["profile"], **paths}


def map_native_catalog(config_path, hdf_path, terrain_path, results_gpkg, output_dir):
    """Map complete catalog WSEL profiles, retaining native and delivery rasters/QA.

    ``terrain_path`` is the independently selected original or healed TIFF. Writes
    only beneath output_dir. Native interpolation is plain sloping, with no extra
    shallow/face weighting. Delivery resampling is recorded separately from native
    reconstruction. Failure propagates; no selector or legacy fallback is allowed.
    """
    from ras_commander import RasFim

    config = _settings(config_path)
    threads = config.getint("05_raster_engine", "int_mapping_threads", fallback=8)
    working = config.getint(
        "05_raster_engine", "int_mapping_max_working_mib", fallback=1024
    )
    inflight = config.getint(
        "05_raster_engine", "int_mapping_inflight_profiles", fallback=1
    )
    if not 1 <= threads <= 64 or working <= 0 or not 1 <= inflight <= 64:
        raise MappingTerrainError(
            "Invalid native thread, memory or in-flight profile budget"
        )
    crs = config.get("06_output_products", "str_output_crs")
    resolution = config.getfloat("06_output_products", "flt_desired_res")
    if not np.isfinite(resolution) or resolution <= 0:
        raise MappingTerrainError(
            "Catalog output resolution must be positive and finite"
        )
    folder = Path(output_dir).resolve()
    folder.mkdir(parents=True, exist_ok=True)
    if (
        (folder / "native").exists()
        or any(folder.glob("09_*.tif"))
        or any(folder.glob("10_*.tif"))
        or any(folder.glob("11_*.gpkg"))
    ):
        raise MappingTerrainError(
            "Mapped products already exist; use a new output folder"
        )
    profiles, cells, provenance = prepare_catalog_profiles(
        config_path, hdf_path, results_gpkg
    )
    provenance.update(
        mapping_terrain=_file(terrain_path),
        commander_commit=COMMANDER_COMMIT,
        engine="rasmapper_sloping",
        render_mode="sloping",
        delivery={
            "crs": crs,
            "resolution": resolution,
            "resampling": "nearest",
            "grid_alignment": "aligned_target",
        },
    )
    with tempfile.TemporaryDirectory(
        prefix=".native-catalog-", dir=folder.parent
    ) as scratch:
        staged = Path(scratch)
        native = staged / "native"
        manifest = RasFim.map_sloping_profiles(
            hdf_path,
            terrain_path,
            native,
            wsel_profiles=profiles,
            polygon=cells,
            render_mode="sloping",
            threads=threads,
            max_working_mib=working,
            max_inflight_profiles=inflight,
        )
        if (
            manifest.empty
            or not manifest.status.isin(("complete", "no_flooded_pixels")).all()
        ):
            raise MappingTerrainError(
                "Native mapping did not produce every catalog profile"
            )
        products = [
            _deliver_profile(row, staged, crs, resolution)
            for _, row in manifest.iterrows()
        ]
        provenance.update(
            native_profile_source=manifest.attrs.get("profile_source"),
            products=products,
        )
        if (
            provenance["hydraulic_hdf"] != _file(hdf_path)
            or provenance["catalog_results"] != _file(results_gpkg)
            or provenance["mapping_terrain"] != _file(terrain_path)
        ):
            raise MappingTerrainError(
                "Mapping inputs changed during native catalog rendering"
            )
        _write(staged / "mapping_catalog.json", provenance)
        for marker in (native / "fim_manifest.json", staged / "mapping_catalog.json"):
            text = marker.read_text(encoding="utf8")
            marker.write_text(
                text.replace(
                    str(staged).replace("\\", "\\\\"), str(folder).replace("\\", "\\\\")
                ),
                encoding="utf8",
            )
        for path in staged.iterdir():
            shutil.move(str(path), folder / path.name)
    return {"profiles": len(manifest), "engine": "rasmapper_sloping"}
