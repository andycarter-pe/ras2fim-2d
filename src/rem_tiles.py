"""Select local 3 m TxDOT REM section DSMs for a model footprint.

Copyright (c) 2026, CLB Engineering Corporation. BSD-3-Clause; see NOTICE.
Adapted with permission from FIM Commander c0aefd016a66de22e19258cf654bfd65e7d0890d.
"""

from pathlib import Path
from typing import Any

import numpy as np

REM_CRS = "EPSG:3081"


def _dataset_has_data(dataset: Any, bounds: tuple[float, float, float, float]) -> bool:
    """Probe an already-open local/remote DSM without changing its backend."""
    from terrio.rio.windows import Window, from_bounds

    window = (
        from_bounds(*bounds, transform=dataset.transform)
        .round_offsets()
        .round_lengths()
    )
    window = window.intersection(Window(0, 0, dataset.width, dataset.height))
    if window.width <= 0 or window.height <= 0:
        return False
    step = 1024
    for row in range(0, int(window.height), step):
        for column in range(0, int(window.width), step):
            block = Window(
                window.col_off + column,
                window.row_off + row,
                min(step, window.width - column),
                min(step, window.height - row),
            )
            data = dataset.read(
                1,
                window=(
                    (int(block.row_off), int(block.row_off + block.height)),
                    (int(block.col_off), int(block.col_off + block.width)),
                ),
            )
            if (np.isfinite(data) & (data != dataset.nodata)).any():
                return True
    return False


def select_tiles(
    tile_root: Path, bounds: tuple[float, float, float, float], crs: Any
) -> tuple[list[Path], list[dict[str, Any]]]:
    """Choose the tiles of a ``rem-tiles`` unit that hold road data in the footprint.

    ``tile_root`` is the unit folder (tiles under ``<district>/dsm``). A tile is used when its
    box meets ``bounds`` (in ``crs``) and it has a valid pixel there; both tests read the
    local copy. Returns the paths in name order and a record of every tile with the reason
    it was skipped.
    """
    from terrio import rio
    from terrio.rio.warp import transform_bounds

    box = transform_bounds(crs, REM_CRS, *bounds, densify_pts=21)
    used: list[Path] = []
    considered: list[dict[str, Any]] = []
    for path in sorted(tile_root.glob("*/dsm/*.tif")):
        relative = path.relative_to(tile_root).as_posix()
        with rio.open(path) as dataset:
            if dataset.crs is None or dataset.crs.to_string() != REM_CRS:
                raise ValueError(f"REM tile must be {REM_CRS}: {relative}")
            west, south, east, north = dataset.bounds
            intersects = (
                west < box[2] and east > box[0] and south < box[3] and north > box[1]
            )
            has_data = intersects and _dataset_has_data(dataset, box)
        if not intersects:
            considered.append(
                {"section": relative, "used": False, "reason": "outside_footprint"}
            )
        elif has_data:
            used.append(path)
            considered.append({"section": relative, "used": True})
        else:
            considered.append(
                {
                    "section": relative,
                    "used": False,
                    "reason": "no_road_data_in_footprint",
                }
            )
    return used, considered
