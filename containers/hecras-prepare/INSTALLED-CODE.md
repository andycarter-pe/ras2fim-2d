# Installed code and GitHub source

Source verified September 12, 2026.


The controller links identify source commit
`dc60b219091e85bcb4564eca45313475c39ce58a`, used for the CPU-aligned build.
Its installed scripts were checked against that source. The two-core setting passed Linux and Windows sample qualification. See the [image status](README.md#images).

| Installed file | Source and purpose |
|---|---|
| `/opt/hecras-prepare/prepare.py` | [prepare.py][build-controller]: Linux job controller. |
| `/opt/hecras-prepare/model_checks.py` | [model_checks.py](https://github.com/gpt-cmdr/ras2fim-2d/blob/dc60b219091e85bcb4564eca45313475c39ce58a/containers/hecras-prepare/model_checks.py): dependency, line-ending and HDF validation. |
| `/opt/hecras-prepare/windows_worker.py` | [windows_worker.py][build-worker]: calls [ras-commander][rc]. |
| `C:\ras2fim-runtime\verify_windows_runtime.py` | [profile_provenance_check_legacy.py][build-helper]: unused preparation helper retained in the runtime. |

The library is installed under `C:\Python311\Lib\site-packages\ras_commander`.
Exact installed source: [RasPrj.py][ras-prj], [RasPlan.py][plan-path],
[geom/GeomPreprocessor.py][clear-geom], [RasPreprocess.py][preprocess],
[RasTcu.py][tcu-status], [RasBco.py][bco], and [ComputeResults.py][result].
All 239 package files in the retained 0.99.2 wheel match public commit
`9e4217713e954236b0c16023e1815c6f2b7a5309` after normalizing line endings: 4 files match raw Git bytes exactly, 235 differ only by CRLF versus LF, and none have other differences. The
[complete installed package source](https://github.com/gpt-cmdr/ras-commander/tree/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander) is inspectable on GitHub.

## Build source

[Dockerfile][build-dockerfile] and [bundle_profile.py][build-exporter] now link
to the current build recipe and exporter. The [source snapshot][build-inputs]
also includes the manifest example and both runtime verification scripts.
The controller source commit is `dc60b219091e85bcb4564eca45313475c39ce58a`;
the installed Windows library source is `9e4217713e954236b0c16023e1815c6f2b7a5309`.
The external wheel, installed HEC-RAS runtime and registry data are not in Git.

See [How to Re-Create this Container](PREPARATION.md) for the build procedure and
[the inventory](runtime-inventory-current.json) for installed dependencies.

[rc]: https://rascommander.info/ras/
[ras-prj]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/RasPrj.py
[plan-path]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/RasPlan.py
[clear-geom]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/geom/GeomPreprocessor.py
[preprocess]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/RasPreprocess.py
[tcu-status]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/RasTcu.py
[bco]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/RasBco.py
[result]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/ComputeResults.py

[build-dockerfile]: https://github.com/gpt-cmdr/ras2fim-2d/blob/dc60b219091e85bcb4564eca45313475c39ce58a/containers/hecras-prepare/Dockerfile
[build-controller]: https://github.com/gpt-cmdr/ras2fim-2d/blob/dc60b219091e85bcb4564eca45313475c39ce58a/containers/hecras-prepare/prepare.py
[build-worker]: https://github.com/gpt-cmdr/ras2fim-2d/blob/dc60b219091e85bcb4564eca45313475c39ce58a/containers/hecras-prepare/windows_worker.py
[build-exporter]: https://github.com/gpt-cmdr/ras2fim-2d/blob/dc60b219091e85bcb4564eca45313475c39ce58a/containers/hecras-prepare/bundle_profile.py
[build-manifest]: https://github.com/gpt-cmdr/ras2fim-2d/blob/dc60b219091e85bcb4564eca45313475c39ce58a/containers/hecras-prepare/runtime-manifest.example.json
[build-helper]: https://github.com/gpt-cmdr/ras2fim-2d/blob/dc60b219091e85bcb4564eca45313475c39ce58a/containers/hecras-prepare/profile_provenance_check_legacy.py
[build-inputs]: https://github.com/gpt-cmdr/ras2fim-2d/blob/codex/phase2-linux-wine-preprocessing/containers/hecras-prepare/BUILD-INPUTS.md

## Input preparation and output validation

[model_checks.py](https://github.com/gpt-cmdr/ras2fim-2d/blob/dc60b219091e85bcb4564eca45313475c39ce58a/containers/hecras-prepare/model_checks.py), installed at `/opt/hecras-prepare/model_checks.py`, verifies
model dependencies, normalizes selected HEC-RAS input text to Windows CRLF,
and validates the 2D mesh and hydraulic tables in both output HDF files. It uses
the h5py and NumPy packages already installed in the Windows Python runtime.
The image Dockerfile copies this module with the controller and worker.

[test_preprocessing_inputs.py](https://github.com/gpt-cmdr/ras2fim-2d/blob/dc60b219091e85bcb4564eca45313475c39ce58a/containers/hecras-prepare/test_preprocessing_inputs.py) is a source-side Docker regression runner for
Windows and Linux callers. It retains separate disposable LF, CRLF, mixed and
missing-dependency test cases outside Git. It is not installed in the image.
