"""Raise-only road restoration for mapping terrain; does not modify hydraulic inputs.

Copyright (c) 2026, CLB Engineering Corporation.
Distributed under the BSD-3-Clause license; see NOTICE.
Adapted with permission from FIM Commander revision c0aefd016a66de22e19258cf654bfd65e7d0890d.
REM elevations are metres; vertical datum is assumed NAVD88 and checked diagnostically.
"""

from __future__ import annotations

import os
import time
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np

NODATA = -9999.0
#: US survey feet per metre (1 m = 3937/1200 ftUS).
FT_US_PER_M = 1.0 / 0.3048006096012192
_FT_PER_M = 1.0 / 0.3048
_SCALE = {"ftUS": FT_US_PER_M, "ft": _FT_PER_M, "m": 1.0}
SCHEMA = "ras2fim-2d-healed-terrain-qa/v1"
#: One lattice tolerance for every grid comparison of the healed terrain, in cells: RAS Mapper
#: writes a cell size that differs from the crop by about 1e-14 ft, which is not a different grid.
LATTICE_TOLERANCE = 1e-6
#: Block edge in cells for the sampling loop (a multiple of the 256-cell TIFF tile).
BLOCK = 1024
#: The corridor core starts this many REM cells (3 m) inside the footprint edge.
CORE_REM_CELLS = 2
REM_CELL_M = 3.0

#: Mask values written to ``validity.tif``.
MASK_BASE = 0
MASK_REM = 1
MASK_REM_FILLS_NODATA = 2  # Reserved for v1; v2 never fills model nodata.


class HealedReasonCode:
    """Stable reason codes of the healed-terrain stage."""

    TILE_INVALID = "HEALED_REM_TILE_INVALID"
    TERRAIN_INVALID = "HEALED_TERRAIN_INVALID"
    DATUM_MISMATCH = "HEALED_DATUM_MISMATCH"
    PARAMETER_INVALID = "HEALED_PARAMETER_INVALID"
    NOT_STAGED = "HEALED_INPUT_NOT_STAGED"
    MODEL_TERRAIN_INVALID = "HEALED_MODEL_TERRAIN_INVALID"


class HealedError(ValueError):
    """The inputs cannot be healed; ``code`` is the stable reason."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code


def lattice_difference(a: Any, b: Any) -> float:
    """Largest difference between two affine transforms, in cells of ``a``.

    Compares the cell size, rotation terms and origin. Two grids are the same lattice when this
    is within :data:`LATTICE_TOLERANCE`; no value is ever resampled to make them so.
    """
    cell = abs(a.a)
    return (
        max(abs(x - y) for x, y in zip(tuple(a)[:6], tuple(b)[:6], strict=True)) / cell
    )


def same_lattice(a: Any, b: Any) -> bool:
    """True when ``a`` and ``b`` differ by at most :data:`LATTICE_TOLERANCE` cells."""
    return lattice_difference(a, b) <= LATTICE_TOLERANCE


def vertical_scale(unit: str) -> float:
    """Multiplier from REM metres to ``unit`` (``ftUS``, ``ft`` or ``m``)."""
    if unit not in _SCALE:
        raise HealedError(
            HealedReasonCode.TERRAIN_INVALID,
            f"vertical unit {unit!r} is not ftUS, ft or m",
        )
    return _SCALE[unit]


def configure_terrio() -> None:
    """Set bounded defaults before Terrio initializes; runner settings take precedence.

    Rayon defaults to one worker, matching the prior sampler. Each output spills to local
    temporary storage above 64 MiB rather than buffering a whole reach in RAM.
    """
    os.environ.setdefault("RAYON_NUM_THREADS", "1")
    os.environ.setdefault("TERRIO_RIO_IN_MEMORY_LIMIT", str(64 * 1024 * 1024))


def _check_base(dataset: Any) -> None:
    transform = dataset.transform
    if (
        dataset.crs is None
        or not dataset.crs.is_projected
        or dataset.count != 1
        or dataset.dtypes != ("float32",)
        or transform.b != 0
        or transform.d != 0
        or transform.a <= 0
        or transform.e >= 0
        or not np.isclose(transform.a, -transform.e, rtol=1e-9)
        or dataset.nodata is None
    ):
        raise HealedError(
            HealedReasonCode.TERRAIN_INVALID,
            "the model terrain must be float32, one band on a north-up square projected grid with nodata",
        )


def _sample_block(sources: Sequence[Any], window: Any, base: Any, scale: float):
    """REM height (base units, NaN outside the footprint) and footprint for one window."""
    configure_terrio()
    from terrio.rio.enums import Resampling
    from terrio.rio.warp import reproject, transform_bounds
    from terrio.rio.windows import Window, from_bounds
    from terrio.rio.windows import bounds as window_bounds

    shape = (int(window.height), int(window.width))
    transform = base.window_transform(window)
    values = np.full(shape, np.nan, "float32")
    footprint = np.zeros(shape, bool)
    for source in sources:
        # Terrio's band adapter reads the whole source. Supply a bounded array instead:
        # densify the destination extent and retain four source cells for interpolation,
        # transformer approximation and rounding at the crop boundary.
        box = transform_bounds(
            base.crs, source.crs, *window_bounds(window, base.transform), densify_pts=21
        )
        crop = from_bounds(*box, transform=source.transform)
        left = max(0, int(np.floor(crop.col_off)) - 4)
        top = max(0, int(np.floor(crop.row_off)) - 4)
        right = min(source.width, int(np.ceil(crop.col_off + crop.width)) + 4)
        bottom = min(source.height, int(np.ceil(crop.row_off + crop.height)) + 4)
        if right <= left or bottom <= top:
            continue
        source_window = Window(left, top, right - left, bottom - top)
        source_data = source.read(1, window=source_window)
        source_transform = source.window_transform(source_window)
        near = np.full(shape, np.nan, "float32")
        bilinear = np.full(shape, np.nan, "float32")
        for out, method in (
            (near, Resampling.nearest),
            (bilinear, Resampling.bilinear),
        ):
            reproject(
                source=source_data,
                src_transform=source_transform,
                src_crs=source.crs,
                destination=out,
                src_nodata=source.nodata,
                dst_nodata=np.nan,
                dst_transform=transform,
                dst_crs=base.crs,
                resampling=method,
                num_threads=1,
            )
        here = np.isfinite(near) & ~footprint
        values[here] = np.where(np.isfinite(bilinear), bilinear, near)[here] * scale
        footprint |= here
    return values, footprint


def _has_rem(
    sources: Sequence[Any], window: Any, base: Any, valid: Sequence[Any]
) -> bool:
    """Cheap test that any REM footprint meets ``window`` (pre-read validity at REM resolution)."""
    if not sources:
        return False
    from terrio.rio.warp import transform_bounds
    from terrio.rio.windows import bounds as window_bounds

    west, south, east, north = transform_bounds(
        base.crs, sources[0].crs, *window_bounds(window, base.transform), densify_pts=5
    )
    for mask, mask_transform in valid:
        col0 = int(np.floor((west - mask_transform.c) / mask_transform.a)) - 2
        col1 = int(np.ceil((east - mask_transform.c) / mask_transform.a)) + 2
        row0 = int(np.floor((north - mask_transform.f) / mask_transform.e)) - 2
        row1 = int(np.ceil((south - mask_transform.f) / mask_transform.e)) + 2
        if mask[max(row0, 0) : max(row1, 0), max(col0, 0) : max(col1, 0)].any():
            return True
    return False


def _rem_validity(
    source: Any, base_bounds: tuple[float, float, float, float], base_crs: Any
):
    """Boolean REM footprint of ``source`` under the base extent, and its transform."""
    from terrio.rio.warp import transform_bounds
    from terrio.rio.windows import Window, from_bounds

    box = transform_bounds(base_crs, source.crs, *base_bounds, densify_pts=21)
    window = (
        from_bounds(*box, transform=source.transform).round_offsets().round_lengths()
    )
    window = window.intersection(Window(0, 0, source.width, source.height))
    if window.width <= 0 or window.height <= 0:
        return np.zeros((0, 0), bool), source.transform
    data = source.read(1, window=window)
    return data != source.nodata, source.window_transform(window)


def _percentiles(values: np.ndarray) -> dict[str, float | None]:
    if values.size == 0:
        return {"p01": None, "p05": None, "median": None, "p95": None, "p99": None}
    p = np.percentile(values, [1, 5, 50, 95, 99])
    return {
        k: round(float(v), 3)
        for k, v in zip(("p01", "p05", "median", "p95", "p99"), p, strict=True)
    }


def build_healed(
    base_path: Path,
    tile_paths: Sequence[Path],
    output: Path,
    *,
    vertical_unit: str = "ftUS",
    datum_tolerance_ft: float = 0.5,
) -> dict[str, Any]:
    """Write ``healed.tif``, ``validity.tif`` and ``difference.tif`` into ``output``; return QA.

    ``base_path`` is the exported model raster with terrain modifications applied (float32, square
    north-up projected grid, nodata set); ``tile_paths`` are REM section DSMs (EPSG:3081, m).
    ``output`` is a new scratch folder. The caller validates QA before using the result.
    Outputs are on the base grid, float32 (heights in ``vertical_unit``) with nodata -9999:

    - ``healed.tif``: max(base, sampled REM) where both are valid; else the base;
    - ``validity.tif``: uint8, 0 unchanged, 1 raised by REM; model nodata stays unchanged;
    - ``difference.tif``: nonnegative healed minus base where both inputs are valid, else nodata.

    Raises :class:`HealedError` for an unusable base. Never writes next to the inputs. The
    ``datum_check`` (REM minus base on corridor core, at least two REM cells inside the
    footprint edge) fails when the median is farther from zero than ``datum_tolerance_ft``.
    """
    configure_terrio()
    from scipy import ndimage
    from terrio import rio
    from terrio.rio.windows import Window

    started = time.perf_counter()
    scale = vertical_scale(vertical_unit)
    # `scale` converts REM metres to the base unit; QA statistics are reported in feet.
    ft_per_unit = {"ftUS": 1.0, "ft": 1.0, "m": FT_US_PER_M}[vertical_unit]
    output.mkdir(parents=True, exist_ok=True)
    profile: dict[str, Any] = {}
    sources = [rio.open(p) for p in tile_paths]
    try:
        with rio.open(base_path) as base:
            _check_base(base)
            for source in sources:
                if source.crs is None or source.crs.to_string() != "EPSG:3081":
                    raise HealedError(
                        HealedReasonCode.TILE_INVALID, "a REM tile is not EPSG:3081"
                    )
            profile = {
                "driver": "GTiff",
                "count": 1,
                "height": base.height,
                "width": base.width,
                "crs": base.crs,
                "transform": base.transform,
                "tiled": True,
                "blockxsize": 256,
                "blockysize": 256,
                "compress": "DEFLATE",
                "zlevel": 6,
            }
            cell = base.transform.a
            metres_per_unit = float(base.crs.linear_units_factor[1])
            core_distance = CORE_REM_CELLS * REM_CELL_M / metres_per_unit  # CRS units
            halo = int(np.ceil(core_distance / cell)) + 2
            validity = [_rem_validity(s, tuple(base.bounds), base.crs) for s in sources]
            t_prepare = time.perf_counter()
            float_profile = {
                **profile,
                "dtype": "float32",
                "nodata": NODATA,
                "predictor": 3,
            }
            counts = {
                "model_valid_cells": 0,
                "model_nodata_cells": 0,
                "rem_cells": 0,
                "rem_over_model_cells": 0,
                "rem_over_nodata_cells": 0,
                "raised_cells": 0,
                "ignored_lower_rem_cells": 0,
                "raised_cells_over_0p1_ft": 0,
                "lowered_cells_over_0p1_ft": 0,
                "blocks": 0,
                "blocks_with_rem": 0,
            }
            diffs_all: list[np.ndarray] = []
            diffs_core: list[np.ndarray] = []
            applied_all: list[np.ndarray] = []
            applied_core: list[np.ndarray] = []
            with (
                rio.open(output / "healed.tif", "w", **float_profile) as healed_out,
                rio.open(output / "difference.tif", "w", **float_profile) as diff_out,
                rio.open(
                    output / "validity.tif",
                    "w",
                    **{**profile, "dtype": "uint8", "predictor": 2},
                ) as mask_out,
            ):
                for row in range(0, base.height, BLOCK):
                    for col in range(0, base.width, BLOCK):
                        window = Window(
                            col,
                            row,
                            min(BLOCK, base.width - col),
                            min(BLOCK, base.height - row),
                        )
                        counts["blocks"] += 1
                        data = base.read(1, window=window)
                        base_ok = np.isfinite(data) & (data != base.nodata)
                        counts["model_valid_cells"] += int(base_ok.sum())
                        counts["model_nodata_cells"] += int((~base_ok).sum())
                        healed = np.where(base_ok, data, NODATA).astype("float32")
                        difference = np.full(healed.shape, NODATA, "float32")
                        mask = np.zeros(healed.shape, "uint8")
                        wide = Window(
                            max(col - halo, 0),
                            max(row - halo, 0),
                            min(col + window.width + halo, base.width)
                            - max(col - halo, 0),
                            min(row + window.height + halo, base.height)
                            - max(row - halo, 0),
                        )
                        if _has_rem(sources, window, base, validity):
                            counts["blocks_with_rem"] += 1
                            values, footprint = _sample_block(
                                sources, wide, base, scale
                            )
                            top, left = (
                                int(window.row_off - wide.row_off),
                                int(window.col_off - wide.col_off),
                            )
                            sel = (
                                slice(top, top + int(window.height)),
                                slice(left, left + int(window.width)),
                            )
                            inside = footprint[sel]
                            rem = values[sel]
                            over_data = inside & base_ok
                            fills = inside & ~base_ok
                            raises = over_data & (rem > data)
                            healed[raises] = rem[raises]
                            mask[raises] = MASK_REM
                            delta = rem[over_data] - data[over_data]
                            applied = healed[over_data] - data[over_data]
                            difference[over_data] = applied
                            counts["raised_cells"] += int(raises.sum())
                            counts["ignored_lower_rem_cells"] += int((delta < 0).sum())
                            counts["rem_cells"] += int(inside.sum())
                            counts["rem_over_model_cells"] += int(over_data.sum())
                            counts["rem_over_nodata_cells"] += int(fills.sum())
                            counts["raised_cells_over_0p1_ft"] += int(
                                (applied * ft_per_unit > 0.1).sum()
                            )
                            counts["lowered_cells_over_0p1_ft"] += int(
                                (applied * ft_per_unit < -0.1).sum()
                            )
                            if delta.size:
                                diffs_all.append(delta.astype("float32"))
                                applied_all.append(applied.astype("float32"))
                                # distance to the footprint edge, measured on the haloed window
                                edge = ndimage.distance_transform_edt(footprint) * cell
                                core = (edge[sel] >= core_distance) & over_data
                                diffs_core.append(
                                    (rem[core] - data[core]).astype("float32")
                                )
                                applied_core.append(
                                    (healed[core] - data[core]).astype("float32")
                                )
                        healed_out.write(healed, 1, window=window)
                        diff_out.write(difference, 1, window=window)
                        mask_out.write(mask, 1, window=window)
            grid = {
                "crs": base.crs.to_string(),
                "transform": list(base.transform)[:6],
                "shape": [base.height, base.width],
                "cell_size_crs_units": cell,
            }
    finally:
        for source in sources:
            source.close()
    all_ft = (
        np.concatenate(diffs_all) * ft_per_unit if diffs_all else np.empty(0, "float32")
    )
    core_ft = (
        np.concatenate(diffs_core) * ft_per_unit
        if diffs_core
        else np.empty(0, "float32")
    )
    applied_ft = (
        np.concatenate(applied_all) * ft_per_unit
        if applied_all
        else np.empty(0, "float32")
    )
    applied_core_ft = (
        np.concatenate(applied_core) * ft_per_unit
        if applied_core
        else np.empty(0, "float32")
    )
    median = float(np.median(core_ft)) if core_ft.size else None
    if median is None:
        datum = {"status": "not_evaluable", "reason": "no road core in the footprint"}
    else:
        datum = {
            "status": "pass" if abs(median) <= datum_tolerance_ft else "fail",
            "median_rem_minus_model_ft": round(median, 3),
            "tolerance_ft": datum_tolerance_ft,
            "core_cells": int(core_ft.size),
            "core_definition": (
                f"at least {CORE_REM_CELLS} REM cells ({CORE_REM_CELLS * REM_CELL_M:g} m) "
                "inside the REM footprint edge"
            ),
            **{k: v for k, v in _percentiles(core_ft).items() if k in ("p05", "p95")},
            "share_abs_over_1_ft": round(float((np.abs(core_ft) > 1.0).mean()), 4),
            "note": "median over all core cells, crossings included (robust to them)",
        }
    cell_area_ac = (cell * metres_per_unit) ** 2 / 4046.8564224
    qa = {
        "schema": SCHEMA,
        "grid": grid,
        "vertical": {
            "base_unit": vertical_unit,
            "rem_unit": "m",
            "rem_datum": "NAVD88 (assumed; the REM rasters do not declare it)",
            "metres_to_base_unit": scale,
        },
        "sampling": {
            "backend": "terrio",
            "backend_version": rio.__version__,
            "footprint": "nearest neighbour validity of the 3 m section DSM",
            "values": "bilinear over valid REM pixels, nearest where bilinear is undefined",
            "why": (
                "REM is 3 m on a ~1 m grid: bilinear avoids 3 x 3 stair-steps and invents "
                "nothing outside the footprint"
            ),
            "overlay": (
                "Raise only: max(model, sampled REM) where both are valid; "
                "lower/equal REM and model nodata stay unchanged; no ramp, no blend; "
                "the first tile in index order wins on overlap"
            ),
        },
        "counts": {
            **counts,
            "rem_area_acres": round(counts["rem_cells"] * cell_area_ac, 3),
        },
        "difference_ft": {
            "all_rem_cells": _percentiles(applied_ft),
            "core_cells": _percentiles(applied_core_ft),
            "max_raise_ft": round(float(applied_ft.max()), 3)
            if applied_ft.size
            else None,
            "max_lowering_ft": round(min(0.0, float(applied_ft.min())), 3)
            if applied_ft.size
            else None,
        },
        "datum_check": datum,
        "rem_difference_ft": {
            "all_rem_cells": _percentiles(all_ft),
            "core_cells": _percentiles(core_ft),
        },
        "healed_differs_from_model": counts["raised_cells"] > 0,
        "seconds": {
            "prepare": round(t_prepare - started, 2),
            "sample_and_write": round(time.perf_counter() - t_prepare, 2),
        },
        "mapping_only": True,
        "approval_state": "not_approved",
    }
    return qa
