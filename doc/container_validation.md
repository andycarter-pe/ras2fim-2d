# Container validation

Image: `ghcr.io/gpt-cmdr/ras2fim-2d:healed-native-20261027-r6`
Digest: `ghcr.io/gpt-cmdr/ras2fim-2d@sha256:9ee963e954fea8c87fb9df18693038576be26d2cfc0ddbf50ce676244e4ab5c3`
Image ID: `sha256:a3ead9678d16945b2edf9f707f303922128da21a49d4cbc94752079687d8546e`
Source commit: `ab6310044aac288782693e7945456686d3cccb84`

Built and checked on a GitHub runner by `build-container.yml`
([run 40557944861](https://github.com/gpt-cmdr/ras2fim-2d/actions/runs/40557944861)),
published by `publish-container.yml`
([run 40559125768](https://github.com/gpt-cmdr/ras2fim-2d/actions/runs/40559125768)).
Logs, `delivery.json`, `SHA256SUMS`, package lists and the Docker-save archive are on the
[release](https://github.com/gpt-cmdr/ras2fim-2d/releases/tag/healed-native-20261027-r6).

Checked inside the image (Python 3.14.10, GDAL 3.13.3, rasterio 1.5.2, Terrio 0.5.0,
ras-commander at commit `777c7746ceb35fc3b1c29d3718935c9013ae009c`, fast-map 0.8.0,
fast-HDF 0.2.3):

- `pytest`: 28 passed.
- Steps 05-08 v1 and v2 files import; `RasFim.map_sloping_profiles` takes
  `wsel_profiles`; ras-commander is installed from the pinned commit;
  `ras2fim-2d.py --help`, `healed_mapping_terrain.py --help`, `nccopy` and
  `ras2fim-batch` are present.
- Andy's sample (`s3://rasfim-2d-sample`, steps 00-04 as inputs, `config_global_windows.ini`):
  steps 5-8 with v1 (the three imports switched back) and with v2 (shipped
  `legacy` / `model`) give the same 181 output files: TIFFs byte-identical,
  GeoPackage layers equal, NetCDF data and attributes equal except the
  `06_history` / `08_filter_history` timestamps. Result: MATCH.
- Same sample with `str_mapping_engine = fast_step05`: NetCDF WSEL arrays equal to legacy.
- Anonymous pull: GHCR token, manifest (tag and digest) and config blob fetched
  without credentials; the config digest equals the image ID above.

Not re-validated for this image:

- Windows test run, and the Windows vs Linux native-mapping bit comparison.
- Native `rasmapper_sloping` on the Davis DS Channel fixture (not public). On the
  ras2fim-2d sample it ran to completion outside the image, Python 3.14; its grid
  differs from legacy (one more row), as expected for a different method.
- Healed terrain on real REM tiles (needs the mounted REM library); covered here only
  by the synthetic tests.
- Steps 0-4, 2b, 2c and 3 (Windows / HEC-RAS).

Full real-model catalog regression and H&H acceptance remain separate work.
`fast_step05` and `rasmapper_sloping` are different methods; neither claims universal
RASMapper pixel parity or speedup.
