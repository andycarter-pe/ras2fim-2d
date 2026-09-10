# Inspectable container build inputs

This source snapshot contains the current container build scripts for review.
The Python controller, Windows worker, and retained profile helper match their
installed copies in the tested bundled images. Building still requires the
external prepared runtime described below.

| File | Purpose |
|---|---|
| [Dockerfile](Dockerfile) | Linux dependencies, installed runtime packaging, and entrypoint. |
| [prepare.py](prepare.py) | Linux job controller. |
| [windows_worker.py](windows_worker.py) | Calls the [ras-commander APIs](https://rascommander.info/ras/) under Wine. |
| [bundle_profile.py](bundle_profile.py) | Exports a prepared runtime into an external build context. |
| [runtime-manifest.example.json](runtime-manifest.example.json) | Manifest structure only; zero fingerprints are placeholders. |
| [verify_bundled_runtime.sh](verify_bundled_runtime.sh) | Fresh-copy runtime verification launcher. |
| [verify_windows_runtime.py](verify_windows_runtime.py) | Version-aware Windows installation and saved TCU check. |
| [profile_provenance_check_legacy.py](profile_provenance_check_legacy.py) | Exact source of the unused helper retained inside the installed profile. |

## External runtime profile

The actual `runtime.json` and `prefix/` are external build inputs, not files in
this source snapshot. The linked manifest example is not a usable release
manifest. A finalized profile supplies real file records and accepted TCU state.

```text
prepared-runtime/
  runtime.json
  prefix/
    drive_c/
      Program Files (x86)/HEC/HEC-RAS/<version>/Ras.exe
      Python311/python.exe
      ras2fim-runtime/ras_commander-0.99.2-py3-none-any.whl
    dosdevices/
      c: -> ../drive_c
      z: -> /
    system.reg
    user.reg
    userdef.reg
```

Vendor installers, installed runtimes, Wine registry contents, model data, and
solver outputs stay outside Git. This snapshot documents their layout without
including their contents. A complete empty-prefix installation recipe remains
unfinished; the current build starts with a retained prepared profile.

## Retained ras-commander wheel

The external runtime contains `ras_commander-0.99.2-py3-none-any.whl`, built from
retained commit `9e4217713e954236b0c16023e1815c6f2b7a5309`. The exact wheel is not
published in this source snapshot. The [public upstream source](https://github.com/gpt-cmdr/ras-commander)
provides library reference code; it does not replace the retained build input.

## Input preparation and output validation

`model_checks.py`, installed at `/opt/hecras-prepare/model_checks.py`, verifies
model dependencies, normalizes selected HEC-RAS input text to Windows CRLF,
and validates the 2D mesh and hydraulic tables in both output HDF files. It uses
the h5py and NumPy packages already installed in the Windows Python runtime.
The image Dockerfile copies this module with the controller and worker.

`test_preprocessing_inputs.py` is a source-side Docker regression runner for
Windows and Linux callers. It retains separate disposable LF, CRLF, mixed and
missing-dependency test cases outside Git. It is not installed in the image.
