"""Standalone healer regression tests, adapted from CLB Engineering's FIM Commander."""

from pathlib import Path

import numpy as np
import pytest
import rasterio
from rasterio.transform import Affine

import terrain_healer as healed

FT_US_PER_M = healed.FT_US_PER_M
CELL = 1.0 / 0.3048006096012192
SIDE = 300


def _origin() -> tuple[float, float]:
    """A point near Kerrville, in EPSG:2278 (ftUS) and EPSG:3081 (m), on whole-metre values."""
    from pyproj import Transformer

    x, y = Transformer.from_crs("EPSG:4326", "EPSG:2278", always_xy=True).transform(
        -99.14, 30.05
    )
    return round(x), round(y)


@pytest.fixture
def base_tif(tmp_path: Path) -> Path:
    """A 1,000 ft flat terrain on a 1 m (in ftUS) grid with a deep cut channel down the middle."""
    x0, y0 = _origin()
    data = np.full((SIDE, SIDE), 1000.0, "float32")
    data[:, 140:160] = (
        990.0  # the hydro-enforcement-like channel (kept as the un-cut base here)
    )
    path = tmp_path / "source.tif"
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=SIDE,
        width=SIDE,
        count=1,
        dtype="float32",
        crs="EPSG:2278",
        transform=Affine(CELL, 0, x0, 0, -CELL, y0),
        nodata=-9999.0,
        tiled=True,
        blockxsize=256,
        blockysize=256,
    ) as out:
        out.write(data, 1)
    return path


def _rem_tile(
    source: Path, base: Path, *, offset_ft: float, district: str = "15_SAT"
) -> Path:
    """A 3 m REM tile with a road strip crossing the base, ``offset_ft`` above 1,000 ft."""
    from rasterio.warp import transform_bounds

    with rasterio.open(base) as ds:
        west, south, east, north = transform_bounds(ds.crs, "EPSG:3081", *ds.bounds)
    x0, y1 = np.floor(west) - 30, np.ceil(north) + 30
    width = int((east - west + 60) // 3) + 1
    height = int((north - south + 60) // 3) + 1
    data = np.full((height, width), -9999.0, "float32")
    rows = slice(height // 3, height // 3 + 6)  # a 18 m strip west to east
    data[rows, :] = (1000.0 + offset_ft) / FT_US_PER_M
    path = source / district / "dsm" / "Section-1-Test.tif"
    path.parent.mkdir(parents=True)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype="float32",
        crs="EPSG:3081",
        transform=Affine(3.0, 0, x0, 0, -3.0, y1),
        nodata=-9999.0,
    ) as out:
        out.write(data, 1)
    return path


def test_overlay_raises_base_only_where_the_rem_is_higher(tmp_path, base_tif):
    tile = _rem_tile(tmp_path / "rem", base_tif, offset_ft=0.2)
    out = tmp_path / "out"
    qa = healed.build_healed(base_tif, [tile], out)
    with (
        rasterio.open(base_tif) as base,
        rasterio.open(out / "healed.tif") as heal,
        rasterio.open(out / "validity.tif") as mask,
        rasterio.open(out / "difference.tif") as diff,
    ):
        b, h, m, d = (r.read(1) for r in (base, heal, mask, diff))
        assert (
            heal.crs == base.crs
            and heal.transform == base.transform
            and h.shape == b.shape
        )
    on_rem = m == healed.MASK_REM
    assert on_rem.any() and not on_rem.all()
    np.testing.assert_array_equal(
        h[~on_rem], b[~on_rem]
    )  # the base is kept elsewhere, un-cut
    np.testing.assert_allclose(
        h[on_rem], 1000.2, atol=1e-3
    )  # REM replaces, no blend or ramp
    np.testing.assert_allclose(d[on_rem], h[on_rem] - b[on_rem], atol=1e-4)
    assert (d[~on_rem] == healed.NODATA).all()
    # the channel stays as delivered (the un-cut crop is the base) outside the strip
    assert h[:, 145][~on_rem[:, 145]].max() == 990.0
    assert qa["counts"]["rem_cells"] == int(on_rem.sum())
    assert qa["datum_check"]["status"] == "pass"
    assert qa["datum_check"]["median_rem_minus_model_ft"] == pytest.approx(
        0.2, abs=0.01
    )
    assert qa["healed_differs_from_model"] is True and qa["mapping_only"] is True


@pytest.mark.parametrize("offset_ft", [-0.25, -3.0])
def test_lower_rem_never_lowers_model_and_still_checks_candidate_datum(
    tmp_path, base_tif, offset_ft
):
    tile = _rem_tile(tmp_path / "rem", base_tif, offset_ft=offset_ft)
    with rasterio.open(base_tif, "r+") as base:
        data = np.full(base.shape, 1000.0, "float32")
        data[:10] = base.nodata
        base.write(data, 1)
    qa = healed.build_healed(base_tif, [tile], tmp_path / "out")
    with rasterio.open(tmp_path / "out/healed.tif") as result:
        np.testing.assert_array_equal(result.read(1), data)
    with rasterio.open(tmp_path / "out/validity.tif") as mask:
        assert not mask.read(1).any()
    with rasterio.open(tmp_path / "out/difference.tif") as difference:
        values = difference.read(1, masked=True)
        assert values.count() > 0 and np.all(values.compressed() == 0)
        assert values.mask[:10].all()
    assert qa["healed_differs_from_model"] is False
    assert qa["counts"]["raised_cells"] == 0
    assert qa["counts"]["ignored_lower_rem_cells"] > 0
    assert qa["counts"]["lowered_cells_over_0p1_ft"] == 0
    assert qa["difference_ft"]["max_lowering_ft"] == 0
    assert qa["datum_check"]["median_rem_minus_model_ft"] == pytest.approx(
        offset_ft, abs=0.05
    )
    assert qa["datum_check"]["status"] == ("pass" if offset_ft == -0.25 else "fail")


def test_mixed_rem_raises_only_and_preserves_model_nodata(tmp_path, base_tif):
    tile = _rem_tile(tmp_path / "rem", base_tif, offset_ft=0.2)
    with rasterio.open(base_tif, "r+") as base:
        data = base.read(1)
        data[:20] = base.nodata
        data[100:180] = 1002.0
        base.write(data, 1)
    qa = healed.build_healed(base_tif, [tile], tmp_path / "out")
    with (
        rasterio.open(tmp_path / "out/healed.tif") as result,
        rasterio.open(tmp_path / "out/validity.tif") as mask,
    ):
        h, m = result.read(1), mask.read(1)
    valid = data != healed.NODATA
    assert np.all(h[valid] >= data[valid])
    np.testing.assert_array_equal(h[~valid], data[~valid])
    np.testing.assert_array_equal(m == healed.MASK_REM, valid & (h > data))
    assert qa["counts"]["raised_cells"] == int((m == healed.MASK_REM).sum())
    assert qa["counts"]["ignored_lower_rem_cells"] > 0
    assert qa["schema"] == "ras2fim-2d-healed-terrain-qa/v1"


def test_datum_check_fails_when_the_rem_is_far_from_the_base(tmp_path, base_tif):
    tile = _rem_tile(tmp_path / "rem", base_tif, offset_ft=3.0)
    qa = healed.build_healed(base_tif, [tile], tmp_path / "out")
    assert qa["datum_check"]["status"] == "fail"
    assert qa["datum_check"]["median_rem_minus_model_ft"] == pytest.approx(
        3.0, abs=0.05
    )


def test_metre_base_unit_scales_the_rem_by_one(tmp_path):
    """A base in metres is compared in metres; the QA converts to feet."""
    from rasterio.warp import transform_bounds

    x0, y0 = _origin()
    from pyproj import Transformer

    mx, my = Transformer.from_crs("EPSG:2278", "EPSG:3081", always_xy=True).transform(
        x0, y0
    )
    base = tmp_path / "m.tif"
    with rasterio.open(
        base,
        "w",
        driver="GTiff",
        height=90,
        width=90,
        count=1,
        dtype="float32",
        crs="EPSG:3081",
        transform=Affine(1.0, 0, round(mx), 0, -1.0, round(my)),
        nodata=-9999.0,
    ) as out:
        out.write(np.full((90, 90), 300.0, "float32"), 1)
    with rasterio.open(base) as ds:
        west, _south, _east, north = transform_bounds(ds.crs, "EPSG:3081", *ds.bounds)
    tile = tmp_path / "rem" / "15_SAT" / "dsm" / "t.tif"
    tile.parent.mkdir(parents=True)
    with rasterio.open(
        tile,
        "w",
        driver="GTiff",
        height=40,
        width=40,
        count=1,
        dtype="float32",
        crs="EPSG:3081",
        transform=Affine(3.0, 0, west - 15, 0, -3.0, north + 15),
        nodata=-9999.0,
    ) as out:
        out.write(np.full((40, 40), 300.0, "float32"), 1)
    qa = healed.build_healed(base, [tile], tmp_path / "out", vertical_unit="m")
    assert qa["datum_check"]["median_rem_minus_model_ft"] == pytest.approx(
        0.0, abs=1e-3
    )
    assert qa["counts"]["rem_cells"] > 0


def test_a_base_that_is_not_a_square_north_up_grid_is_refused(tmp_path):
    bad = tmp_path / "bad.tif"
    with rasterio.open(
        bad,
        "w",
        driver="GTiff",
        height=4,
        width=4,
        count=1,
        dtype="float32",
        crs="EPSG:2278",
        transform=Affine(2.0, 0, 0, 0, -1.0, 0),
        nodata=-9999.0,
    ) as out:
        out.write(np.zeros((4, 4), "float32"), 1)
    with pytest.raises(
        healed.HealedError, match=healed.HealedReasonCode.TERRAIN_INVALID
    ):
        healed.build_healed(bad, [], tmp_path / "out")


@pytest.mark.parametrize("offset_ft", [-0.3, 0.2])
def test_terrio_sampling_matches_rasterio_with_overlap_holes_and_spilled_outputs(
    tmp_path, base_tif, monkeypatch, offset_ft
):
    """Independently check every haloed block using GDAL; validate Terrio writes with GDAL."""
    from rasterio.warp import reproject
    from terrio import rio

    first = _rem_tile(tmp_path / "first", base_tif, offset_ft=offset_ft)
    second = _rem_tile(tmp_path / "second", base_tif, offset_ft=5)
    with rasterio.open(first, "r+") as ds:
        data = ds.read(1)
        valid = data != ds.nodata
        rows, cols = np.indices(data.shape)
        data[valid] += (cols[valid] % 7) * 0.01
        data[(rows % 3 == 0) & (cols % 5 == 0)] = ds.nodata
        ds.write(data, 1)
    monkeypatch.setenv("TERRIO_RIO_IN_MEMORY_LIMIT", "0")
    monkeypatch.setattr(healed, "BLOCK", 128)
    actual_sample = healed._sample_block
    checked = []

    class WindowedSource:
        def __init__(self, dataset):
            self.dataset = dataset

        def __getattr__(self, name):
            return getattr(self.dataset, name)

        def read(self, *args, **kwargs):
            assert kwargs.get("window") is not None, "whole-section reads are forbidden"
            crop = kwargs["window"]
            # Blocks near the road must not decode the full section for every warp.
            assert crop.width < self.width or crop.height < self.height
            return self.dataset.read(*args, **kwargs)

    def checked_sample(sources, window, base, scale):
        actual, footprint = actual_sample(
            [WindowedSource(s) for s in sources], window, base, scale
        )
        expected = np.full(actual.shape, np.nan, "float32")
        expected_footprint = np.zeros(actual.shape, bool)
        for source in sources:
            with rasterio.open(source.name) as src:
                near = np.full(actual.shape, np.nan, "float32")
                bilinear = near.copy()
                for dest, method in (
                    (near, rasterio.enums.Resampling.nearest),
                    (bilinear, rasterio.enums.Resampling.bilinear),
                ):
                    reproject(
                        rasterio.band(src, 1),
                        dest,
                        src_nodata=src.nodata,
                        dst_nodata=np.nan,
                        dst_transform=base.window_transform(window),
                        dst_crs=base.crs.to_string(),
                        resampling=method,
                        num_threads=1,
                    )
                here = np.isfinite(near) & ~expected_footprint
                expected[here] = (
                    np.where(np.isfinite(bilinear), bilinear, near)[here] * scale
                )
                expected_footprint |= here
        np.testing.assert_array_equal(actual, expected)
        np.testing.assert_array_equal(footprint, expected_footprint)
        checked.append(True)
        return actual, footprint

    monkeypatch.setattr(healed, "_sample_block", checked_sample)
    out = tmp_path / "healed"
    qa = healed.build_healed(base_tif, [first, second], out)
    assert len(checked) > 1
    assert qa["sampling"]["backend"] == "terrio"
    assert qa["sampling"]["backend_version"] == "0.5.0"
    assert qa["counts"]["lowered_cells_over_0p1_ft"] == 0
    with rasterio.open(base_tif) as base:
        for name, dtype in (
            ("healed.tif", "float32"),
            ("difference.tif", "float32"),
            ("validity.tif", "uint8"),
        ):
            with rasterio.open(out / name) as written, rio.open(out / name) as other:
                assert written.crs == base.crs and written.transform == base.transform
                assert written.shape == base.shape and written.dtypes == (dtype,)
                np.testing.assert_array_equal(written.read(1), other.read(1))


def test_public_builder_without_rem_preserves_model_and_writes_all_outputs(
    tmp_path, base_tif
):
    out = tmp_path / "empty-rem"
    qa = healed.build_healed(base_tif, [], out)
    assert qa["counts"]["rem_cells"] == 0
    assert qa["datum_check"]["status"] == "not_evaluable"
    with rasterio.open(base_tif) as base, rasterio.open(out / "healed.tif") as result:
        np.testing.assert_array_equal(base.read(1), result.read(1))
    assert (out / "difference.tif").is_file() and (out / "validity.tif").is_file()
