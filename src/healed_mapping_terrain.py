"""Mapping-only integration with the bundled BSD-licensed Terrio terrain healer."""

import configparser
import hashlib
import json
import shutil
import tempfile
from functools import lru_cache
from pathlib import Path


class MappingTerrainError(ValueError):
    """Requested healing failed qualification; never silently fall back."""


def _json(path):
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise MappingTerrainError(f"Expected a JSON object: {path}")
    return value


def _write(path, value):
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _stat(path):
    if not path.is_file() or path.is_symlink():
        raise MappingTerrainError(f"Missing or linked terrain input: {path}")
    s = path.stat()
    return (s.st_size, s.st_mtime_ns, s.st_ctime_ns, s.st_ino)


@lru_cache(maxsize=1024)
def _hashed(path_name, fingerprint):
    digest = hashlib.sha256()
    path = Path(path_name)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    if _stat(path) != fingerprint:
        raise MappingTerrainError(f"Input changed while hashing: {path}")
    return digest.hexdigest()


def _file(path):
    path = Path(path).resolve()
    fingerprint = _stat(path)
    return {"sha256": _hashed(str(path), fingerprint), "bytes": fingerprint[0]}


def _settings(config_path):
    config = configparser.ConfigParser(interpolation=None)
    if not config.read(str(config_path)):
        raise MappingTerrainError(f"Cannot read configuration: {config_path}")
    return config


def _load_builder():
    try:
        from rem_tiles import select_tiles
        from terrain_healer import SCHEMA, build_healed, configure_terrio

        configure_terrio()
        from terrio import __version__
    except ImportError as exc:
        raise MappingTerrainError(
            "Healing requires Python 3.14 and requirements-healed-terrain.txt"
        ) from exc
    if SCHEMA != "ras2fim-2d-healed-terrain-qa/v1":
        raise MappingTerrainError(
            "Require the bundled raise-only Terrio implementation"
        )
    if __version__ != "0.5.0":
        raise MappingTerrainError("Require the qualified Terrio 0.5.0 release")
    return build_healed, select_tiles


def mapping_terrain(config_path, model_terrain, cache_dir, *, allow_build=True):
    """Build/reuse a mapping-only healed TIFF, returning its path and provenance.

    [05_raster_engine] str_mapping_terrain defaults to model. When enabled, str_rem_tiles
    names a local district/dsm section library, resolved relative to the config file.
    The input is the native model DEM with terrain modifications already rasterized.
    Heights must be ft or ftUS to match this workflow's NetCDF feet contract. No input
    or HEC-RAS project is modified. All new rasters/QA live in cache_dir/<input hash>.
    Step 7 uses allow_build=False and requires the Step 5 cached variant to exist.
    """
    config = _settings(config_path)
    mode = config.get(
        "05_raster_engine", "str_mapping_terrain", fallback="model"
    ).strip()
    if mode not in ("model", "healed"):
        raise MappingTerrainError("str_mapping_terrain must be model or healed")
    if mode == "model":
        return str(model_terrain), None
    try:
        from terrio import rio

        build_healed, select_tiles = _load_builder()
        source_name = config.get(
            "05_raster_engine", "str_rem_tiles", fallback=""
        ).strip()
        if not source_name:
            raise MappingTerrainError(
                "Set str_rem_tiles to the local maintenance-section DSM library"
            )
        source = Path(source_name).expanduser()
        if not source.is_absolute():
            source = Path(config_path).resolve().parent / source
        if not source.is_dir() or not any(source.glob("*/dsm/*.tif")):
            raise MappingTerrainError(
                "str_rem_tiles must contain <district>/dsm/*.tif section DSMs"
            )
        unit = config.get(
            "05_raster_engine", "str_terrain_vertical_unit", fallback="ftUS"
        )
        if unit not in ("ft", "ftUS"):
            raise MappingTerrainError(
                "The existing mapping/NetCDF workflow requires terrain in feet"
            )
        tolerance = config.getfloat(
            "05_raster_engine", "flt_datum_tolerance_ft", fallback=0.5
        )
        if not 0 <= tolerance < float("inf"):
            raise MappingTerrainError(
                "flt_datum_tolerance_ft must be finite and nonnegative"
            )
        with rio.open(model_terrain) as base:
            paths, considered = select_tiles(source, tuple(base.bounds), base.crs)
        inputs = {
            "method": "ras2fim-2d-healed-terrain/v1",
            "model_terrain": _file(model_terrain),
            "tiles": [
                {"section": p.relative_to(source).as_posix(), **_file(p)} for p in paths
            ],
            "vertical_unit": unit,
            "datum_tolerance_ft": tolerance,
            "terrio": "0.5.0",
        }
        digest = hashlib.sha256(
            json.dumps(
                inputs, sort_keys=True, separators=(",", ":"), allow_nan=False
            ).encode()
        ).hexdigest()
        cache = Path(cache_dir)
        folder = cache / digest
        record_path = folder / "mapping_terrain.json"
        if folder.exists():
            record = _json(record_path)
            if record.get("inputs") != inputs or record.get("input_sha256") != digest:
                raise MappingTerrainError(
                    "Cached mapping terrain has inconsistent provenance"
                )
            for name in (
                "healed.tif",
                "difference.tif",
                "validity.tif",
                "qa.json",
                "review.html",
            ):
                if record.get("outputs", {}).get(name) != _file(folder / name):
                    raise MappingTerrainError(f"Cached mapping terrain changed: {name}")
            _check_qa(_json(folder / "qa.json"))
            return str(folder / "healed.tif"), record
        if not allow_build:
            raise MappingTerrainError(
                "Run Step 5 with this healing configuration before Step 7"
            )
        cache.mkdir(parents=True, exist_ok=True)
        scratch = Path(tempfile.mkdtemp(prefix=".healing-", dir=cache))
        try:
            qa = build_healed(
                Path(model_terrain),
                paths,
                scratch,
                vertical_unit=unit,
                datum_tolerance_ft=tolerance,
            )
            qa["sections_considered"] = considered
            _write(scratch / "qa.json", qa)
            from terrain_review import write_review

            write_review(model_terrain, scratch, qa)
            try:
                _check_qa(qa)
            except MappingTerrainError:
                # Retain signed datum diagnostics, but no reusable success marker.
                scratch.rename(
                    cache
                    / (
                        "failed-"
                        + digest
                        + "-"
                        + scratch.name.removeprefix(".healing-")
                    )
                )
                raise
            if inputs["model_terrain"] != _file(model_terrain) or any(
                row != {"section": p.relative_to(source).as_posix(), **_file(p)}
                for row, p in zip(inputs["tiles"], paths)
            ):
                raise MappingTerrainError("Terrain inputs changed during healing")
            record = {
                "schema": "ras2fim-2d-mapping-terrain/v1",
                "kind": "healed",
                "input_sha256": digest,
                "inputs": inputs,
                "outputs": {
                    name: _file(scratch / name)
                    for name in (
                        "healed.tif",
                        "difference.tif",
                        "validity.tif",
                        "qa.json",
                        "review.html",
                    )
                },
                "mapping_only": True,
                "approval_state": "not_approved",
            }
            _write(scratch / "mapping_terrain.json", record)
            scratch.rename(folder)
        finally:
            if scratch.exists():
                if scratch.resolve().parent != cache.resolve():
                    raise MappingTerrainError(
                        "Refusing cleanup outside the terrain cache"
                    )
                shutil.rmtree(scratch)
        return str(folder / "healed.tif"), record
    except MappingTerrainError:
        raise
    except (OSError, ValueError, TypeError, KeyError, ImportError) as exc:
        raise MappingTerrainError(
            f"Cannot prepare healed mapping terrain: {exc}"
        ) from exc


def _check_qa(qa):
    if (
        qa.get("schema") != "ras2fim-2d-healed-terrain-qa/v1"
        or qa.get("mapping_only") is not True
        or qa.get("sampling", {}).get("backend") != "terrio"
        or qa.get("datum_check", {}).get("status") not in ("pass", "not_evaluable")
        or qa.get("counts", {}).get("lowered_cells_over_0p1_ft") != 0
    ):
        raise MappingTerrainError(
            "Healed-terrain QA failed; inspect its retained signed datum diagnostics"
        )


def _same_selection(record, expected):
    if record is None or expected is None:
        return record is expected
    return record.get("input_sha256") == expected.get("input_sha256")


def bind_mapping_terrain(output_dir, record):
    """Pin a mapping step's inputs; reject mixing terrain variants in a raster folder."""
    folder = Path(output_dir)
    marker = folder / "mapping_terrain.json"
    if marker.exists():
        if not _same_selection(_json(marker), record):
            raise MappingTerrainError(
                "Mapping terrain changed; use a separate output folder"
            )
    elif record is not None:
        if any(folder.rglob("*.tif")) or any(folder.rglob("*.nc")):
            raise MappingTerrainError(
                "Existing rasters have no healed-terrain provenance; use a new output folder"
            )
        folder.mkdir(parents=True, exist_ok=True)
        _write(marker, record)


def step5_mapping_record(step5_dir):
    """Read one run's mapping provenance; None denotes the default model terrain."""
    marker = Path(step5_dir) / "mapping_terrain.json"
    return _json(marker) if marker.exists() else None


def check_netcdf_mapping_inputs(step6_dir, record):
    """Require Step 7 terrain to agree with the Step 6 raster provenance."""
    marker = Path(step6_dir) / "mapping_terrain.json"
    previous = _json(marker) if marker.exists() else None
    if not _same_selection(previous, record):
        raise MappingTerrainError(
            "Step 6 rasters use different or unrecorded mapping terrain; rerun Steps 5-7 in a new output folder"
        )


def netcdf_mapping_metadata(simple_wsel_dir):
    """Global NetCDF attributes describing the terrain embedded by Step 7.

    None are added for the default legacy engine with model terrain, so default
    NetCDF files match v1. Otherwise the keys continue the existing numbering
    after Step 8's 10_still_needed_values.
    """
    from mapping_engine import recorded_mapping_engine

    engine = recorded_mapping_engine(Path(simple_wsel_dir).parent)
    record = step5_mapping_record(Path(simple_wsel_dir).parent)
    if record is None:
        if engine == "legacy":
            return {}
        return {"11_terrain_kind": "model", "12_mapping_engine": engine}
    return {
        "11_terrain_kind": "healed_mapping_only",
        "12_mapping_engine": engine,
        "13_mapping_terrain_sha256": record["input_sha256"],
        "14_mapping_terrain_approval": "not_approved",
    }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Prepare mapping terrain without running hydraulics"
    )
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--terrain", required=True, type=Path)
    parser.add_argument("--cache", required=True, type=Path)
    args = parser.parse_args()
    path, record = mapping_terrain(args.config, args.terrain, args.cache)
    print(json.dumps({"terrain": path, "provenance": record}, indent=2))
