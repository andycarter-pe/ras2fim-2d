# How to Re-Create this Container

[Installed code and GitHub source](INSTALLED-CODE.md) lists the packaged scripts,
links the current container source snapshot and public library references.

Companion to [the operating guide](OPERATION.md) and [release verification](RELEASE-VERIFICATION.md).


This section is for rebuilding the image. Routine users can start with the
installed image and the [run command](README.md#run).

**The repeatable build currently starts with a retained, prepared Wine profile.**
A complete installation recipe starting with an empty Wine prefix has not yet
been demonstrated. The steps below identify the inputs, build procedure, and
remaining work needed for that fully independent reconstruction.

### 1. Gather the build inputs

Use the [GitHub source snapshot][build-inputs] containing these files under
`containers/hecras-prepare/`. The source links require access to the private CLB
repository. External inputs link to their example or layout documentation.

| Input | Purpose |
|---|---|
| [Dockerfile][build-dockerfile] | Linux tools, user 1000, bundled Wine profile, and entrypoint. |
| [prepare.py][build-controller] | Linux job controller, private runtime setup, worker launch, and receipts. |
| [windows_worker.py][build-worker] | Calls the linked [ras-commander][rc] APIs under Wine. |
| [bundle_profile.py][build-exporter] | Exports a prepared runtime into the external build context. |
| External [runtime.json][build-manifest] and [prefix/][build-inputs] | Installed HEC-RAS 6.5 and dependencies. Links show the manifest example and layout; actual runtime contents stay outside Git. |
| Retained [ras_commander-0.99.2-py3-none-any.whl][build-inputs] | External wheel from build `9e4217713e954236b0c16023e1815c6f2b7a5309`; the link documents this input, not a wheel download. |

The base is `python:3.13.15-slim-trixie` pinned by its Dockerfile digest. The
Dockerfile also pins Wine, Xvfb, xauth, procps and tini. Windows package versions
are recorded in the inventory below. The public API links are reference
documentation; rebuilding the exact installed Python package requires the
retained wheel/build source, not an arbitrary current upstream checkout.

### 2. Prepare or recover the installed Wine profile

The retained profiles contain Windows Python 3.11.9, .NET 4.8, Visual C++
runtimes, GDI+, fonts, HEC-RAS, and accepted TCU state. Preserve the original;
work on a separate copy. Preserve `dosdevices/c: -> ../drive_c` and
`dosdevices/z: -> /` as links without following them during copies.

To create a replacement from an empty prefix, install those prerequisites,
the selected official HEC-RAS release, and the retained
[ras-commander][rc-github] wheel. Complete any installer and first-launch TCU
prompts, then close Wine cleanly and verify the saved state with
[RasTcu.status()][tcu-status]. The retained profile's dependency list is an
inventory; it does not yet supply all prerequisite download URLs or a tested
installer sequence.

Preparation details needed to match the current runtime:

- The retained installation wheel is
  `C:\ras2fim-runtime\ras_commander-0.99.2-py3-none-any.whl`. It was installed
  using Windows Python with `pip install --no-index --no-deps --no-cache-dir
  --disable-pip-version-check --force-reinstall --no-compile`. Its dependencies
  were already present. `--no-deps` does not install them in an empty prefix.
- Wine initialization used
  `WINEDLLOVERRIDES="mscoree,mshtml,winemenubuilder.exe,winedbg.exe="` for
  `wineboot`, followed by Winetricks prerequisites. This initialization setting
  must not disable `mscoree` during normal HEC-RAS 7.0.1 startup, which needs .NET.
- HEC-RAS 7.0.1's outer installer stalled on its .NET 4.8.1 web prerequisite.
  Its extracted official MSI was installed directly. The resulting .NET 4.8
  profile passed the retained sample. This is a recorded installation workaround.
- The 6.5 profile already held accepted TCU state. The 6.6 installer and first
  launch were accepted. The 7.0.1 application's first-launch terms were accepted.
  [RasTcu.accept()][tcu-accept] persisted acceptance where needed. Wine was
  stopped with `wineserver -k` followed by `wineserver -w` while the surrounding
  container process remained alive, allowing registry writes to finish.
- Saved TCU values are `650`, `660`, and `701` for the respective releases,
  scoped to the Windows user and full installed executable path. Fresh-copy
  checks verified them without invoking [RasTcu.accept()][tcu-accept] again.

Finalize `runtime.json` after installation: schema `ras-commander-runtime/v1`,
`kind: wine`, `hec_ras_version: "6.5"`, selected executable/Python/wheel
paths, build identity, and actual artifact records. A retained finalized profile
already supplies this manifest. If installed files change, regenerate the
manifest before export. A portable manifest-generation command remains part
of the reconstruction work; the example manifest's zero values are placeholders.

### 3. Export the profile and build the image

Run from the repository root on Linux. The source is the finalized profile
from step 2. The destination must be a new directory outside Git. Export removes
temporary installers, caches, files covered by the personal-file filters, and
unrelated HEC-RAS program directories while retaining vendor notices. Review the
exported contents before publishing.

```bash
python containers/hecras-prepare/bundle_profile.py \
  --source /controlled/hecras-6.5-wine \
  --destination /controlled/releases/hecras-6.5 \
  --version 6.5

docker buildx build --load --platform linux/amd64 \
  --build-context hecras_runtime=/controlled/releases/hecras-6.5 \
  --build-arg HEC_RAS_VERSION=6.5 \
  --build-arg RAS_COMMANDER_COMMIT=9e4217713e954236b0c16023e1815c6f2b7a5309 \
  --build-arg RAS_COMMANDER_WHEEL_SHA256=dc0c2f9baa9db66afee01e34eb00d7a04f53e7b82c3340c786b2e9aef1795932 \
  --tag local/hecras-6.5:review \
  --file containers/hecras-prepare/Dockerfile .
```

The external `hecras_runtime` context supplies the installed files; the Dockerfile
places them at `/runtime/wine-seed`, owned by root and readable by the worker.
The entrypoint is `tini -- python /opt/hecras-prepare/prepare.py`. A normal job
uses its own writable prefix and runs without downloads or installation.

### 4. Validate the reconstructed image

The repository's `verify_bundled_runtime.sh` and `verify_windows_runtime.py`
check the selected installation and saved TCU state from a fresh runtime copy.
Use a Linux Docker host with `/dev/ntsync`; run as host UID 1000 or root. The
audit directory must be new and outside Git:

```bash
bash containers/hecras-prepare/verify_bundled_runtime.sh \
  local/hecras-6.5:review 6.5 /scratch/hecras-audit-6.5-001
```

Then use the [run command](README.md#run) with `local/hecras-6.5:review` and a
disposable validation model. Require the expected three output files, a successful
receipt, no timeout, and no copied full-simulation result. Test with networking
disabled and no external runtime mount. Retain the build log, installed-software
inventory, TCU check, image identity, and sample-run evidence outside Git.

### Integrity checks recorded by the implementation

The runtime manifest contains SHA-256 file fingerprints. `prepare.py` reads the
declared runtime files at job startup and compares them with that inventory;
the Windows worker checks that its installed wheel matches the retained wheel.
The receipt also records fingerprints of the three generated output files,
which the host controller checks before staging them. These checks detect
changed runtime files or changed handoff files. They add file-reading work and
are an implementation choice for release traceability; HEC-RAS and Wine do not
need file hashing to perform preprocessing. They do not check hydraulic accuracy.

### Remaining work for a build starting from nothing

A fully independent rebuild still needs a tested empty-prefix installation
sequence, retrievable official installers and Python dependencies, locked
Winetricks/prerequisite versions, portable manifest creation, a distributable
validation fixture, and a clean-host rebuild test. The original HEC-RAS 6.5 setup
EXE identity is unavailable; a retained MSI-cache identity is recorded instead.
The prepared images are usable for the tested workflow while this reconstruction
work remains incomplete.

Use the corresponding HEC-RAS version and finalized profile to build 6.6 or 7.0.1.

### Installed software for reconstruction

The recorded inventory is common to the three tested bundled images except for the selected HEC-RAS installation. It contains 512 Linux packages, Linux Python and pip, 49 Windows Python distributions, and 24 Windows component registration records. This describes what is installed, not a minimal dependency list or an installer download lock file. Windows registration records can retain installation history.

| Component | Observed version or preparation |
|---|---|
| Linux base | Debian Trixie; Python base image pinned by digest in the Dockerfile. |
| Linux Python / pip | 3.13.15 / 26.2.1. |
| WineHQ stable | `11.0.0.0~trixie-1`, with 32-bit and 64-bit support. |
| Virtual display | Xvfb `2:21.1.16-1.3+deb13u3`, with xauth. |
| Process support | procps and tini, pinned in the Dockerfile. |
| Windows Python | 3.11.9, 64-bit, installed at `C:\Python311`. |
| Microsoft .NET | Framework 4.8, registered version `4.8.03761`; history includes the .NET 4.0 prerequisite. |
| Microsoft Visual C++ | 2010 x86/x64 `10.0.40219`; 2013 redistributables `12.0.30501.0` with `12.0.21005` components; 2015-2019 redistributables `14.24.28127.4` with `14.24.28127` components. |
| Windows drawing and fonts | Native GDI+ and core fonts installed through Winetricks; exact bootstrap installer revisions were not retained in this recipe. |
| Windows compatibility setting | Initialization history switches through Windows XP for prerequisite installation and ends at Windows 10. |
| [ras-commander][rc] | 0.99.2, installed from the retained installation wheel. |
| Windows Python dependencies | 49 distributions, including numpy, pandas, h5py, psutil, pywin32, pythonnet, and geospatial libraries. Every observed distribution version is listed below. |
| HEC-RAS | Only the selected 6.5, 6.6, or 7.0.1 program directory is retained in each release profile. |

Installed auxiliary packages such as `hms-commander` and `pathlib` were inherited
from the source environment. Their presence is recorded; this work does not
claim they are necessary for preprocessing or remove them without qualification.

### Windows Python distributions

These 49 distribution versions were the same in all three tested bundled runtimes.

| Package | Version |
|---|---|
| affine | `3.0.0` |
| attrs | `26.1.0` |
| certifi | `2026.7.22` |
| cffi | `2.1.1` |
| charset-normalizer | `3.5.1` |
| click | `8.4.2` |
| click-plugins | `1.1.1.2` |
| cligj | `0.7.2` |
| clr_loader | `0.2.10` |
| colorama | `0.4.6` |
| contourpy | `1.3.3` |
| cycler | `0.12.1` |
| fonttools | `4.63.0` |
| fsspec | `2026.7.0` |
| geopandas | `1.1.4` |
| h5py | `3.16.0` |
| hms-commander | `0.3.1` |
| idna | `3.19` |
| kiwisolver | `1.5.0` |
| matplotlib | `3.11.1` |
| numpy | `2.4.6` |
| packaging | `26.3` |
| pandas | `3.0.5` |
| pathlib | `1.0.1` |
| pefile | `2024.8.26` |
| pillow | `12.3.0` |
| pip | `24.0` |
| psutil | `7.2.2` |
| pycparser | `3.0` |
| pyogrio | `0.13.0` |
| pyparsing | `3.3.2` |
| pyproj | `3.7.2` |
| python-dateutil | `2.9.0.post0` |
| pythonnet | `3.0.5` |
| pywin32 | `312` |
| [ras-commander][rc] | `0.99.2` |
| rasterio | `1.4.4` |
| rasterstats | `0.21.0` |
| requests | `2.34.2` |
| rtree | `1.4.1` |
| scipy | `1.17.1` |
| setuptools | `65.5.0` |
| shapely | `2.1.2` |
| simplejson | `4.1.1` |
| six | `1.17.0` |
| tqdm | `4.70.0` |
| tzdata | `2026.3` |
| urllib3 | `2.7.0` |
| xarray | `2026.7.0` |

### Windows component registrations

These registration records describe the prepared profiles and can include installation history.

| Registered component | Version |
|---|---|
| Microsoft .NET Framework 4.8 | `4.8.03761` |
| Microsoft Visual C++ 2010  x64 Redistributable - 10.0.40219 | `10.0.40219` |
| Microsoft Visual C++ 2010  x86 Redistributable - 10.0.40219 | `10.0.40219` |
| Microsoft Visual C++ 2013 Redistributable (x64) - 12.0.30501 | `12.0.30501.0` |
| Microsoft Visual C++ 2013 Redistributable (x86) - 12.0.30501 | `12.0.30501.0` |
| Microsoft Visual C++ 2013 x64 Additional Runtime - 12.0.21005 | `12.0.21005` |
| Microsoft Visual C++ 2013 x64 Minimum Runtime - 12.0.21005 | `12.0.21005` |
| Microsoft Visual C++ 2013 x86 Additional Runtime - 12.0.21005 | `12.0.21005` |
| Microsoft Visual C++ 2013 x86 Minimum Runtime - 12.0.21005 | `12.0.21005` |
| Microsoft Visual C++ 2015-2019 Redistributable (x64) - 14.24.28127 | `14.24.28127.4` |
| Microsoft Visual C++ 2015-2019 Redistributable (x86) - 14.24.28127 | `14.24.28127.4` |
| Microsoft Visual C++ 2019 X64 Additional Runtime - 14.24.28127 | `14.24.28127` |
| Microsoft Visual C++ 2019 X64 Minimum Runtime - 14.24.28127 | `14.24.28127` |
| Microsoft Visual C++ 2019 X86 Additional Runtime - 14.24.28127 | `14.24.28127` |
| Microsoft Visual C++ 2019 X86 Minimum Runtime - 14.24.28127 | `14.24.28127` |
| Python 3.11.9 (64-bit) | `3.11.9150.0` |
| Python 3.11.9 Add to Path (64-bit) | `3.11.9150.0` |
| Python 3.11.9 Core Interpreter (64-bit) | `3.11.9150.0` |
| Python 3.11.9 Development Libraries (64-bit) | `3.11.9150.0` |
| Python 3.11.9 Executables (64-bit) | `3.11.9150.0` |
| Python 3.11.9 Standard Library (64-bit) | `3.11.9150.0` |
| Python 3.11.9 Tcl/Tk Support (64-bit) | `3.11.9150.0` |
| Python 3.11.9 Utility Scripts (64-bit) | `3.11.9150.0` |
| Python 3.11.9 pip Bootstrap (64-bit) | `3.11.9150.0` |

The [complete recorded inventory](runtime-inventory-20260910.json) includes all
Linux packages and Windows dependencies, custom code identities, and runtime settings.

### Custom code retained in the runtime

`prepare.py` and `windows_worker.py` are installed at `/opt/hecras-prepare/`.
An unused preparation helper remains at
`C:\ras2fim-runtime\verify_windows_runtime.py`; its exact source is
[profile_provenance_check_legacy.py](profile_provenance_check_legacy.py).
It assumes HEC-RAS 6.5 and is not called by normal jobs. Use the repository's
[verify_windows_runtime.py](verify_windows_runtime.py), launched by
[verify_bundled_runtime.sh](verify_bundled_runtime.sh), for current three-version
qualification. Installer and run evidence locations are in
[release verification](RELEASE-VERIFICATION.md).

[rc]: https://rascommander.info/ras/
[rc-github]: https://github.com/gpt-cmdr/ras-commander
[ras-prj]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/RasPrj.py#L125
[init]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/RasPrj.py#L2462
[plan-path]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/RasPlan.py#L787
[clear-geom]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/geom/GeomPreprocessor.py#L1152
[run-flags]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/RasPlan.py#L1394
[preprocess]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/RasPreprocess.py#L97
[tcu-status]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/RasTcu.py#L270
[tcu-accept]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/RasTcu.py#L449
[bco]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/RasBco.py#L25
[logging]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/RasBco.py#L95
[monitor]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/RasBco.py#L144
[terminate]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/RasPreprocess.py#L964
[result]: https://github.com/gpt-cmdr/ras-commander/blob/bab6179027fadfda2b143beccde42ce12e457a0f/ras_commander/ComputeResults.py#L185

[build-dockerfile]: https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/Dockerfile
[build-controller]: https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/prepare.py
[build-worker]: https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/windows_worker.py
[build-exporter]: https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/bundle_profile.py
[build-manifest]: https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/runtime-manifest.example.json
[build-helper]: https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/profile_provenance_check_legacy.py
[build-inputs]: https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/BUILD-INPUTS.md

## Reproduce the input-format regression checks

Run [test_preprocessing_inputs.py](https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/test_preprocessing_inputs.py) from a Windows or Linux Docker host with a
fresh external evidence directory. It makes separate disposable copies for LF,
CRLF, mixed line endings and missing dependencies. The first three must retain
the mesh and hydraulic tables; the missing-dependency case must fail without
changing the model files. The source sample is never modified.

```text
python containers/hecras-prepare/test_preprocessing_inputs.py --image rascommander/hec-ras-wine-precompute_6.5:v4 --source-project <absolute-spawned-sample.prj> --work-dir <new-external-evidence-directory>
```

The source project must have sibling `source_terrain` and `projection` folders
in its parent directory, as in `sample_data/sample_output/02_model_copies`.
The image Dockerfile includes [model_checks.py](https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/model_checks.py) beside `windows_worker.py`;
both are required build inputs. Run the same checks for each released runtime.
