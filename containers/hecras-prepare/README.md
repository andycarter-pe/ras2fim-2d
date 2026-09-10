# Linux/Wine HEC-RAS preprocessing

[Installed code and GitHub source](INSTALLED-CODE.md) lists the packaged scripts,
links the current container source snapshot and public library references.

The container provides Wine so [ras-commander][rc] can call installed Windows
HEC-RAS on Linux. It reads the model from a host folder mounted into the container
and writes preprocessing files through that same mount.

For Andy's review, start with [Wine, installation paths, and host data](OPERATION.md#wine-the-hec-ras-installation-and-host-data)
and the [API call sequence](OPERATION.md#exact-ras-commander-call-sequence).
[How to Re-Create this Container](PREPARATION.md) is a separate build guide with
installation notes, dependencies, and verification instructions. The
[release verification record](RELEASE-VERIFICATION.md) contains test evidence.

## Models created on Windows or Linux

Pull the current image before testing (`docker pull rascommander/hec-ras-wine-precompute_6.5:v4`, or the matching 6.6 image). Use a disposable model copy and mount its referenced terrain and projection
folders. Relative paths are resolved inside the container: `..\source_terrain`
and `..\projection` must lead to mounted folders. Mounting a model folder alone
does not expose its siblings. The terrain HDF must also have access to its
referenced TIFF files.

Use `--user root` for a command that works with both Docker Desktop and
Ubuntu WSL Windows-drive mounts. Docker Desktop 4.90 also passed the default
UID-1000 CRLF checks for HEC-RAS 6.5 and 6.6. Ubuntu's direct `/mnt/c` mount
produced a hidden HEC-RAS `Run-time error '75': Path/File access error` under
UID 1000; root passed there. Use a disposable model folder, since a failure
after HEC-RAS starts can leave partial outputs. See the
[host qualification results](RELEASE-VERIFICATION.md).

Before HEC-RAS starts, [model_checks.py][build-checks] reads the selected plan, geometry and
unsteady flow inputs, checks the existing 2D geometry HDF, and verifies the
projection, terrain HDF and referenced raster files. Missing dependencies fail
before any model file is modified.

The worker then normalizes the selected `.prj`, `.p##`, `.g##` and `.u##` text
files to Windows CRLF line endings, preserving all other bytes. LF, CRLF and
mixed inputs are accepted. This is necessary even though the container runs on
Linux: HEC-RAS is a Windows application under Wine. The generated `.x##` output
is still converted to Linux LF for the separate Linux calculation stage.

Before reporting success, the worker opens both the geometry HDF and temporary
plan HDF. They must retain the input's 2D area names and cell counts and contain
valid cell-volume and face-area elevation tables. Face counts may change when
HEC-RAS versions rebuild a mesh. The receipt records the normalized input names,
checked dependencies and HDF validation. A file's size, `File Type` attribute,
or presence of a `/Results` group is not proof of valid preprocessing.

The required starting model includes a populated geometry HDF, a `.rasmap`, a
projection file and a terrain HDF with its rasters. Qualification covers the
repository sample; it does not certify a complete hydraulic simulation.


## Images

| HEC-RAS | Image | Status |
|---|---|---|
| 6.5 | `rascommander/hec-ras-wine-precompute_6.5:v4` | Published, installed runtime and saved TCU state. |
| 6.6 | `rascommander/hec-ras-wine-precompute_6.6:v4` | Published, installed runtime and saved TCU state. |
| 7.0.1 | `rascommander/hec-ras-wine-precompute_7.0.1:v4` | Built and tested locally; publication pending. |

The 6.5 and 6.6 `latest` tags select their bundled releases. Normal jobs require
no separate HEC-RAS installation or external runtime mount.

## Run on a Windows drive

The following PowerShell command mounts the model and both dependency folders.
Set `$models` to the parent folder containing all three and `$name` to the
project filename without `.prj`. For HEC-RAS 6.6, use the matching 6.6 image.

```powershell
$models = 'E:\mac_test_output_Sept12\02_model_copies'
$name = '1919912_wb-2427466_wb-2427467_14-hr_100-cfs_to_11609-cfs'
docker pull rascommander/hec-ras-wine-precompute_6.5:v4
docker run --rm --user root --mount "type=bind,src=$models\$name,dst=/job" --mount "type=bind,src=$models\source_terrain,dst=/source_terrain,readonly" --mount "type=bind,src=$models\projection,dst=/projection,readonly" rascommander/hec-ras-wine-precompute_6.5:v4 prepare --project "/job/$name.prj" --plan 01 --timeout 900 --replace-generated
```

The native Windows Docker CLI and these mounts were qualified with Docker
Desktop 4.90 on CLB-08 for HEC-RAS 6.5 and 6.6. The PowerShell command also
passed with spaces in the host folder path. From a WSL shell, use Linux source
paths such as `/mnt/c/...`; `wsl -d Ubuntu -- docker ...` from PowerShell also
needs those WSL paths. See the [test record](RELEASE-VERIFICATION.md).

## Run on a Linux filesystem

```bash
docker pull rascommander/hec-ras-wine-precompute_6.5:v4

docker run --rm --network none --read-only \
  --user 1000:1000 --cpus 2 --memory 6g --memory-swap 6g \
  --cap-drop ALL --security-opt no-new-privileges \
  --tmpfs /tmp:rw,nosuid,nodev,size=512m,mode=1777 \
  --env RAS2FIM_JOB_ROOT=/job \
  --mount type=bind,src=/shared/model,dst=/job \
  --mount type=volume,dst=/run/ras-job \
  --device /dev/ntsync \
  rascommander/hec-ras-wine-precompute_6.5:v4 \
  prepare --project /job/model.prj --plan 01 \
          --timeout 900 --run-id example-001 --replace-generated
```

`/shared/model` is a disposable folder on the Linux Docker host, writable by
UID 1000. The container sees it as `/job`; generated files remain in that same
host folder after the container exits. Mount terrain/projection dependencies
where the model expects them. Use a unique run ID. `--replace-generated` can
remove an existing final plan HDF, so preserve wanted results elsewhere.

The Linux recipe uses `/dev/ntsync` when the host supplies it. CLB-08 WSL2
testing uses the simpler root-user command without passing this device. Read and agree to the [HEC-RAS terms](https://www.hec.usace.army.mil/confluence/rasdocs/rasum/6.6/terms-and-conditions-of-use)
before use. [init_ras_project()][init] uses `accept_tcu=True` with the prepared
profile's saved acceptance state.

For project `model`, plan 01 and geometry 01, the required outputs are
`model.p01.tmp.hdf`, `model.b01`, and `model.x01`. The receipt is
`.ras-commander/runs/<run-id>/prepare.json` under the model folder. Exit codes
are 0 for success, 1 for a recorded failure, and 2 for configuration/I/O errors.
Step 2b waits for the full batch before staging the next calculation stage.

## Build and inspection

Follow [How to Re-Create this Container](PREPARATION.md) to export a finalized
Wine profile and build with the external `hecras_runtime` context. The build
currently depends on retained installed profiles; a complete empty-prefix
installation procedure remains unfinished. The [inventory](runtime-inventory-20260910.json)
records installed software, and [release verification](RELEASE-VERIFICATION.md)
records exact image identities and test results.

Published Docker Hub documentation copies: [6.5](dockerhub/6.5.md),
[6.6](dockerhub/6.6.md), [7.0.1](dockerhub/7.0.1.md).

## Optional runtime override

Leave `str_wine_profile` blank to use the bundled installation. Advanced users
can supply a matching controlled profile at `/runtime/wine-seed` read-only.
It must declare the selected HEC-RAS version and matching [ras-commander][rc]
build in its [runtime manifest](runtime-manifest.example.json). The default
job uses the installation already included in the image.

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

[build-checks]: https://github.com/gpt-cmdr/ras2fim-2d/blob/3c5011dbe150d8311e45598401b16645c6e61932/containers/hecras-prepare/model_checks.py
