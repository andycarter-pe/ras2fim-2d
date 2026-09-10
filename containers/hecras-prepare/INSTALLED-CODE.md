# Installed code and GitHub source

Verified September 10, 2026.


The following links are pinned to the [current source snapshot][build-inputs]
in the public PR branch. The installed scripts match this source snapshot.

| Installed file | Source and purpose |
|---|---|
| `/opt/hecras-prepare/prepare.py` | [prepare.py][build-controller]: Linux job controller. |
| `/opt/hecras-prepare/windows_worker.py` | [windows_worker.py][build-worker]: calls [ras-commander][rc]. |
| `C:\ras2fim-runtime\verify_windows_runtime.py` | [profile_provenance_check_legacy.py][build-helper]: unused preparation helper retained in the runtime. |

The library is installed under `C:\Python311\Lib\site-packages\ras_commander`.
Public file references: [RasPrj.py][ras-prj], [RasPlan.py][plan-path],
[geom/GeomPreprocessor.py][clear-geom], [RasPreprocess.py][preprocess],
[RasTcu.py][tcu-status], [RasBco.py][bco], and [ComputeResults.py][result].
The preprocessing, TCU, and BCO files match installed source after line-ending
normalization. The other four are upstream references; the retained wheel
contains the exact installed files.

## Build source

[Dockerfile][build-dockerfile] and [bundle_profile.py][build-exporter] now link
to the current build recipe and exporter. The [source snapshot][build-inputs]
also includes the manifest example and both runtime verification scripts.
All files were read back from GitHub and compared with the working checkout.

The source snapshot commit is `3c5011dbe150d8311e45598401b16645c6e61932`. Public library references remain
upstream references where indicated above; publishing these container scripts
does not publish the external wheel, installed HEC-RAS runtime, or registry data.

See [How to Re-Create this Container](PREPARATION.md) for the build procedure and
[the inventory](runtime-inventory-20260910.json) for installed dependencies.

[rc]: https://rascommander.info/ras/
[ras-prj]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/RasPrj.py#L125
[plan-path]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/RasPlan.py#L787
[clear-geom]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/geom/GeomPreprocessor.py#L1152
[preprocess]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/RasPreprocess.py#L97
[tcu-status]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/RasTcu.py#L270
[bco]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/RasBco.py#L25
[result]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/ComputeResults.py#L185

[build-dockerfile]: https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/Dockerfile
[build-controller]: https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/prepare.py
[build-worker]: https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/windows_worker.py
[build-exporter]: https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/bundle_profile.py
[build-manifest]: https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/runtime-manifest.example.json
[build-helper]: https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/profile_provenance_check_legacy.py
[build-inputs]: https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/BUILD-INPUTS.md

## Input preparation and output validation

[model_checks.py](https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/model_checks.py), installed at `/opt/hecras-prepare/model_checks.py`, verifies
model dependencies, normalizes selected HEC-RAS input text to Windows CRLF,
and validates the 2D mesh and hydraulic tables in both output HDF files. It uses
the h5py and NumPy packages already installed in the Windows Python runtime.
The image Dockerfile copies this module with the controller and worker.

[test_preprocessing_inputs.py](https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/test_preprocessing_inputs.py) is a source-side Docker regression runner for
Windows and Linux callers. It retains separate disposable LF, CRLF, mixed and
missing-dependency test cases outside Git. It is not installed in the image.
