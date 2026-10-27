"""Optional public RasCommander adapter for the existing Step05 mapping recipe.

All heights retain the input WSEL/terrain units. CRS and delivery resolution come
from [06_output_products]. Source HDF, GeoPackage and terrain are read-only.
"""

import json
import shutil
import tempfile
from importlib.metadata import version
from pathlib import Path

from healed_mapping_terrain import MappingTerrainError, _settings


def mapping_engine(config_path):
    """Read the independently configured catalog mapping method."""
    engine = (
        _settings(config_path)
        .get("05_raster_engine", "str_mapping_engine", fallback="legacy")
        .strip()
    )
    if engine not in ("legacy", "fast_step05", "rasmapper_sloping"):
        raise MappingTerrainError(
            "str_mapping_engine must be legacy, fast_step05 or rasmapper_sloping"
        )
    return engine


def bind_mapping_engine(output_dir, engine):
    """Prevent changing mapping methods in an existing Step05 product folder."""
    folder = Path(output_dir).resolve()
    marker = folder / "mapping_engine.json"
    if marker.exists():
        if json.loads(marker.read_text())["engine"] != engine:
            raise MappingTerrainError("Mapping engine changed; use a new output folder")
        return
    # Legacy writes no marker, so default v2 output folders are identical to v1.
    # An unmarked folder already reads as legacy (recorded_mapping_engine).
    if engine == "legacy":
        return
    if engine != "legacy" and any(folder.rglob("*.tif")):
        raise MappingTerrainError(
            "Existing products have no engine provenance; use a new output folder"
        )
    folder.mkdir(parents=True, exist_ok=True)
    marker.write_text(
        json.dumps({"schema": "ras2fim-2d-mapping-engine/v1", "engine": engine}) + "\n"
    )


def recorded_mapping_engine(output_dir):
    """Read the method bound to a raster folder; older unmarked products are legacy."""
    marker = Path(output_dir) / "mapping_engine.json"
    return json.loads(marker.read_text())["engine"] if marker.exists() else "legacy"


def check_catalog_engine(config_path, step6_dir):
    """Reject packaging rasters from another configured mapping method."""
    if mapping_engine(config_path) != recorded_mapping_engine(step6_dir):
        raise MappingTerrainError(
            "Step6 mapping engine differs; rerun Steps5-7 in a new output folder"
        )


def map_fast_step05(config_path, hdf_path, terrain_path, results_gpkg, output_dir):
    """Map existing catalog WSEL columns with the pinned public Rust Step05 engine.

    Writes the existing 09/10 TIFFs, 11 floodplain GeoPackages and engine manifest.
    Rejects existing mapped products. Failure propagates; no legacy fallback.
    """
    try:
        from ras_commander import RasFim
    except ImportError as exc:
        raise MappingTerrainError(
            "Install requirements-mapping-engine.txt for fast_step05"
        ) from exc
    config = _settings(config_path)
    threads = config.getint("05_raster_engine", "int_mapping_threads", fallback=8)
    if not 1 <= threads <= 64:
        raise MappingTerrainError("int_mapping_threads must be between 1 and 64")
    folder = Path(output_dir).resolve()
    folder.mkdir(parents=True, exist_ok=True)
    if (
        any(folder.glob("09_*.tif"))
        or any(folder.glob("10_*.tif"))
        or any(folder.glob("11_*.gpkg"))
    ):
        raise MappingTerrainError(
            "Mapped products already exist; use a new output folder"
        )
    with tempfile.TemporaryDirectory(
        prefix=".fast-step05-", dir=folder.parent
    ) as scratch:
        products = Path(scratch) / "products"
        manifest = RasFim.map_ras2fim_results(
            results_gpkg,
            hdf_path,
            terrain_path,
            products,
            output_crs=config.get("06_output_products", "str_output_crs"),
            output_resolution=config.getfloat("06_output_products", "flt_desired_res"),
            threads=threads,
            grid_anchor="ras2fim",
            coplanar="ignore",
        )
        # Product paths are rewritten after staging so the public manifest remains usable.
        metadata_path = products / "fim_manifest.json"
        metadata = metadata_path.read_text(encoding="utf8")
        metadata_path.write_text(
            metadata.replace(
                str(products).replace("\\", "\\\\"),
                str(folder.resolve()).replace("\\", "\\\\"),
            ),
            encoding="utf8",
        )
        for path in products.iterdir():
            shutil.move(str(path), folder / path.name)
    return {
        "profiles": len(manifest),
        "ras_commander": version("ras-commander"),
        "rust_engine": version("ras-commander-fast-map"),
    }


def map_catalog(config_path, hdf_path, terrain_path, results_gpkg, output_dir):
    """Dispatch an explicitly selected accelerated method; never change methods on failure."""
    selected = mapping_engine(config_path)
    if selected == "rasmapper_sloping":
        from native_mapping import map_native_catalog

        return map_native_catalog(
            config_path, hdf_path, terrain_path, results_gpkg, output_dir
        )
    if selected == "fast_step05":
        return map_fast_step05(
            config_path, hdf_path, terrain_path, results_gpkg, output_dir
        )
    raise MappingTerrainError(
        "Legacy mapping is handled by the existing Step05 entrypoint"
    )
