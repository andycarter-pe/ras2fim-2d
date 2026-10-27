# Optional mapping terrain and engine

## What changes

Steps 5–7 can use a restored ground surface for depth, wet/dry mapping, and the
terrain embedded in the NetCDF catalog. Hydraulics, source HDFs, and original
terrain remain unchanged. The bundled healer is BSD-3-Clause code contributed by
CLB Engineering, adapted from FIM Commander (see [`NOTICE`](../NOTICE)); no private
package is required.

The active Dockerfile and conda environment use Python 3.14 (conda-forge). The Dockerfile copies
this checkout, rather than cloning main, so feature builds contain their changes.
`Dockerfile_old` is historical and does not provide this feature.

## Steps 5–7 version 2

The changes are in new files beside Andy's originals, which are unchanged:

- `src/generate_depth_wsel_rasters_mp_05_v2.py`
- `src/generate_and_simplify_rasters_mp_06_v2.py`
- `src/prep_for_netcdf_mp_07_v2.py`

Each is a copy of the original with a few added lines, the same entry function
name and arguments, and a `Revised - 2026.10.27 -- v2` header line.
`src/ras2fim-2d.py` selects them with three imports; to return to the originals,
change those three imports back (drop `_v2`). With the shipped `legacy` / `model`
settings, the v2 steps write the same files as the originals and add no NetCDF
attributes.

## Configuration and catalog defaults

The existing function arguments, CLI flags, product names, and default mapping
method remain available. Missing new settings mean `legacy` and `model`.

```ini
[05_raster_engine]
str_mapping_engine = legacy
int_mapping_threads = 8
int_mapping_max_working_mib = 1024
int_mapping_inflight_profiles = 1
str_mapping_terrain = model
str_rem_tiles =
str_terrain_vertical_unit = ftUS
flt_datum_tolerance_ft = 0.5
```

- `str_mapping_terrain`: `model` uses the original DEM; `healed` builds/reuses the
  raise-only mapping DEM. This choice is independent of mapping engine.
- `str_mapping_engine`: `legacy` runs the existing Python Step05 implementation;
  `fast_step05` calls `RasFim.map_ras2fim_results` using the same catalog WSEL columns;
  `rasmapper_sloping` uses native cell-corner reconstruction with caller-supplied
  catalog WSEL through the merged RasCommander API.
- `int_mapping_threads`: 1–64 Rust workers, default 8, for accelerated methods.
- `int_mapping_max_working_mib`: positive planned native mapping working budget,
  default 1024 MiB; not a hard RSS limit.
- `int_mapping_inflight_profiles`: 1–64 native profiles in flight, default 1.
- `str_rem_tiles`: local REM root, absolute or relative to this INI file, containing
  `<district>/dsm/<section>.tif`. Required for healing; there is no implicit download.
- `str_terrain_vertical_unit`: `ftUS` or `ft`, matching the source DEM and WSEL.
  The existing catalog/NetCDF workflow stores elevations in feet.
- `flt_datum_tolerance_ft`: finite, nonnegative median candidate-difference tolerance,
  default 0.5 feet. A failed check stops the healed mapping request.

To make healing the maintainer's default for new catalogs, set
`str_mapping_terrain = healed` and a mounted `str_rem_tiles` path in the distributed
configuration template (`src/config_global.ini`, and the Windows template if used).
An individual catalog can override with `model`. Existing INI files keep their
old behavior. To independently make accelerated Step05 the default, set
`str_mapping_engine = fast_step05` or `rasmapper_sloping`. To make both native
sloping and healed terrain the default, change both independent settings.
Use separate output directories for variants.

## Pull the published container

The image for this PR is published under the contributor fork:

```sh
docker pull ghcr.io/gpt-cmdr/ras2fim-2d:healed-native-20261027-r6
# or, pinned to the exact tested image:
docker pull ghcr.io/gpt-cmdr/ras2fim-2d@sha256:9ee963e954fea8c87fb9df18693038576be26d2cfc0ddbf50ce676244e4ab5c3
```

The `@sha256:...` registry digest differs from the Docker image ID and from the
SHA-256 of the compressed archive. The archive on the release is a transfer fallback.

How it is built and published (two workflows):

- `build-container.yml` runs on a pushed `healed-native-*` tag (or `workflow_dispatch`
  once on the default branch). On a GitHub runner it builds the Dockerfile, runs the
  tests, import and tool checks inside the image, and runs sample steps 5-8 from
  `s3://rasfim-2d-sample` with v1, v2 and `fast_step05`; v2 must match v1. It then
  attaches the Docker-save archive, `delivery.json`, `SHA256SUMS` and the logs to a
  DRAFT release. It pushes nothing to a registry.
- `publish-container.yml` runs when a person publishes that release. It checks the
  archive SHA-256 and the loaded image ID against `delivery.json`, then pushes that
  exact image to `ghcr.io/<owner>/<repository>:<tag>` with the repository-scoped
  Actions token, and writes the digest to the `registry-publication` artifact.
  Ordinary releases without `delivery.json` are unaffected.

The maintainer can run the same two workflows to publish under their own namespace.
Set the package to public once in GitHub package settings, then check an anonymous pull.

## Build and run the container

```sh
docker build -t ras2fim-2d:healed-terrio .
# Optional transfer to another Docker host:
docker save ras2fim-2d:healed-terrio | gzip > ras2fim-2d-healed-terrio.tar.gz
# On the receiving host:
gunzip -c ras2fim-2d-healed-terrio.tar.gz | docker load
```

The image includes the application source, healer, Terrio 0.5.0,
`ras-commander[fim]` pinned to merged feature commit
`777c7746ceb35fc3b1c29d3718935c9013ae009c`, fast-map 0.8.0 and fast-HDF 0.2.3.
The source commit matters: the PyPI 0.110.1 release does not include this API.
No new Rust build is required. The public source-archive URL is pinned in
`requirements-mapping-engine.txt` and the conda environment. Scientific and
native GIS libraries are solved for Python 3.14 from conda-forge in a separate
application environment, avoiding conflicts with the base image's conda Python. No Rust compiler,
private credentials, FIM Commander runtime, or HEC-RAS installation is required
for the healer or accelerated Step05 mapping. The image retains the existing
global-input downloads. Large REM inputs are mounted separately and are not baked in.
The Docker environment is not yet a fully locked scientific dependency set; retain
image digest, package inventory, and source snapshot with each delivered build.

A standalone terrain preparation lets an engineer inspect the surface first:

```sh
# config.ini sets str_mapping_terrain=healed and str_rem_tiles=/rem
docker run --rm --cpus=4 --memory=8g \
  -v /host/section-dsms:/rem:ro \
  -v /host/study:/study:ro -v /host/review:/review \
  ras2fim-2d:healed-terrio \
  python healed_mapping_terrain.py --config /study/config.ini \
  --terrain /study/model.tif --cache /review/terrain
```

Run Step05 with its existing CLI, substituting your actual inputs:

```sh
docker run --rm --cpus=8 --memory=16g \
  -v /host/section-dsms:/rem:ro \
  -v /host/study:/study:ro -v /host/output:/output \
  ras2fim-2d:healed-terrio \
  python generate_depth_wsel_rasters_mp_05_v2.py \
  -c /study/config.ini -i /study/plan.p01.hdf -t /study/model.tif \
  -g /study/hydraulic_results.gpkg -o /output/05_large_terrain/run
```

Continue the existing Steps 6–7 against the same output tree and configuration.
The original terrain argument stays the original DEM; Step7 resolves and verifies
the same cached healed selection, rather than building another variant.

For an existing Python 3.14 GIS environment instead of Docker:

```sh
python -m pip install -r requirements-healed-terrain.txt \
  -r requirements-mapping-engine.txt 'numpy<3' 'pandas<3'
```

## Inputs, resolution, and sampling

Use the native 3 m maintenance-section DSMs from the
[TxDOT Road Elevation Model](https://web.corral.tacc.utexas.edu/nfiedata/road3d/),
EPSG:3081 with heights in metres. These are road-corridor surfaces, including
bridge decks, not a continuous statewide bare-earth DEM. The district EPSG:3857
COG mosaics are resampled display products; do not substitute those for section DSMs.

The model DEM must be float32, single-band, projected, square, north-up, with NoData and
terrain modifications already rasterized. The output keeps its CRS, affine transform,
shape and cell size. A fine output grid does not create finer source-road detail.

The nearest-neighbor valid footprint of the REM admits candidate pixels. Inside
that footprint, elevations are bilinear, with nearest fallback where bilinear is
undefined. First section in sorted relative-path order wins on overlap. REM metres
are converted to configured model units, then only samples higher than valid model
ground are applied. Lower samples and model NoData remain unchanged. There is no
edge ramp, smoothing, channel carving or filling of model NoData. Bridge decks
remain decks; potential terrain obstructions require engineering review.

Processing uses 1024-cell blocks, a halo for the 6 m road-core diagnostic, bounded
source reads and a 64 MiB Terrio output spill threshold. Caller-provided Rayon/spill
settings take precedence. QA percentile samples can still consume memory on large
footprints; this is not a hard total-memory bound.

The REM rasters do not declare a vertical datum. NAVD88 is assumed, and a raw signed
REM-minus-model median is reported at least 6 m inside the corridor. This is a
consistency diagnostic, not a datum transformation or proof of datum compatibility.
No evaluable core produces `not_evaluable`, explicitly recorded for manual review.

## Artifacts, QA, and resuming

Step5 caches by source/selected-tile hashes, section names, units, tolerance and
method/backend version under `_healed_terrain/<SHA-256>/` beside per-run directories.

- `healed.tif`: float32 restored elevations on the original grid.
- `difference.tif`: float32 applied raise where both inputs are valid, zero where
  an eligible sample did not raise ground, NoData outside that comparison domain.
- `validity.tif`: uint8 mask, 1 where raised and 0 otherwise; no model NoData fill.
- `qa.json`: `ras2fim-2d-healed-terrain-qa/v1`, grid, units, sampling, counts,
  signed candidate and applied-change statistics, datum status, and timings.
- `review.html`: standalone original/healed/difference pixel comparison of the
  first raised 256-cell block, with a common elevation scale and QA text. Open it
  directly in a browser. It is a locator, not a full-extent review map.
- `mapping_terrain.json`: `ras2fim-2d-mapping-terrain/v1`, input/output SHA-256
  identities, settings, method, and explicit `mapping_only` / `not_approved` state.

Open original, healed and difference TIFFs together in QGIS or another GIS to
review full-extent road crossings, bridges, channels and corridor edges. Verify
source resolution and vertical compatibility before accepting a catalog. Process
success and raise-only arithmetic do not establish hydraulic acceptance.

Step5 writes `mapping_terrain.json` when terrain is healed and `mapping_engine.json`
when the engine is not `legacy`; the default writes neither, and an unmarked folder
reads as `legacy` / `model`. Step6 carries both selections to its shared output
folder. Step7 checks that identity before embedding the selected terrain and checks
its configured engine against Step6. When either option is not the default, Step7
adds NetCDF attributes numbered after `10_still_needed_values`: `11_terrain_kind`,
`12_mapping_engine`, and for healed terrain `13_mapping_terrain_sha256` and
`14_mapping_terrain_approval`. Source and cache artifacts are
hashed; changed cache files fail validation. Failed datum attempts retain
`failed-*` rasters, HTML and QA without a reusable success marker.

Changing terrain/engine/settings or source data requires a new output directory
and rerunning Steps 5–7. The fast adapter rejects existing mapped products and
propagates failures. It stages a complete mapping batch before transferring
products; transfer itself is not a multi-file atomic transaction. Remove a failed
partial output directory before retrying. No automatic legacy fallback occurs.
Fast outputs retain `09_*_depth_nudge.tif`, `10_*_wsel_nudge.tif`,
`11_*_floodplain_ar.gpkg`, and the public engine's `fim_manifest.json`.

## Native sloping catalog mapping

Select `str_mapping_engine = rasmapper_sloping` with either terrain mode. The
adapter calls `RasFim.map_sloping_profiles(wsel_profiles=...)` from merged
[RasCommander PR 710](https://github.com/gpt-cmdr/ras-commander/pull/710).
`fast_step05` remains the separate accelerated Delaunay method, preserving the
existing `grid_anchor=ras2fim` and coplanar-sample handling.

The native adapter accepts one mesh because existing Step4 analyzes only the
first mesh. It rejects multi-mesh catalogs rather than associating ambiguous IDs.
It reuses Step4's HDF wet-cell calculation, two-decimal WSEL rounding, and
rounded-depth-zero dry rule. It recomputes each existing discharge window before
spatial clipping, producing every physical and ghost cell ID. Cells proved never
wet are explicitly NaN/dry. Ghost values come from those same window calculations;
no whole-plan maximum or representative timestep is substituted. The retained
GeoPackage WSEL columns must match these results exactly at float32 precision,
otherwise mapping fails and asks for Step4 to be rerun with the matching HDF.

The legacy window slice ends at `(window_index + 1) * hour_count - 1`, exclusive.
This integration preserves that existing endpoint policy, including its omitted
last sample, and records exact sample ranges in `mapping_catalog.json`. It does
not change discharge selection or rounding. Profiles built from different cells'
peak times are aggregates; no simultaneous hydraulic state is implied.

Complete mesh values and topology remain available outside the output mask. The
GeoPackage's retained cell geometry limits output pixels. One selected original
or healed TIFF supplies mapping ground; original HDF geometry and hydraulic
properties remain unchanged. Plain cell-corner sloping is used; face weighting
and shallow enhancements are not enabled by this adapter.

Native rasters and the canonical profile-hash manifest remain in `native/` on the
selected terrain grid. Top-level 09/10 products are separate nearest-neighbor
warps to the configured catalog CRS and resolution on an aligned target grid.
Their 11 floodplain GeoPackages are regenerated from delivered depth, with the
existing `0 < depth <= 10000` threshold. This keeps Step6's expected CRS and names.
Native parity comparisons must use `native/`, before delivery reprojection.
`mapping_catalog.json` records HDF/GeoPackage/terrain hashes, window policy,
Commander commit, supplied-profile hashes, and delivery settings. Inputs are
rechecked before moving staged outputs. There is no method fallback.

This adds a full HDF WSEL read for complete profile preparation, using the same
existing Step4 helper. Allow memory for its depth/WSEL arrays in addition to the
Rust working budget. The engine timing alone is not an end-to-end speed guarantee.

See the [RasFim API](https://rascommander.info/ras/api/fim/) and
[qualification evidence](https://rascommander.info/ras/development/fim-mapping-performance/).
Selector equivalence and native qualification have distinct scopes. Arbitrary
window aggregates, changed terrain, and delivery reprojection require their own
review; this integration does not claim universal RASMapper pixel parity or speedup.

## Tests and attribution

```sh
python -m pip install pytest
python -m pytest -q
```

Public CI runs the standalone suite without private credentials. Synthetic tests
cover raise-only/NoData rules, datum rejection, cache integrity, Step5–7 hooks,
exact Terrio-vs-Rasterio sampling arrays, spilled outputs, and real Rust Step05
mapping against selected ground. The Rust smoke test replaces only HDF/GPKG
input ingestion. Native integration tests compare complete window profiles with
Step4, reject stale/multi-mesh inputs, and check native/delivery CRS and packaging
provenance. Full real-model catalog regression remains a separate qualification.

Code provenance: CLB Engineering's FIM Commander revision
`c0aefd016a66de22e19258cf654bfd65e7d0890d`, adapted with permission and contributed
under the BSD-3-Clause license; copyright and the list of contributed files are in
[`NOTICE`](../NOTICE). `LICENSE.md` is unchanged. The public mapping dependency carries
its own notices, including attribution to Andy's Step05 implementation. REM data:
Liu et al. (2025), [DOI 10.13139/ORNLNCCS/2574440](https://doi.org/10.13139/ORNLNCCS/2574440),
CC BY 4.0. Mounted data are not included in the code/container license.
