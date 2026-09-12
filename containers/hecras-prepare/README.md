# Linux/Wine HEC-RAS preprocessing

[Installed code and GitHub source](INSTALLED-CODE.md) lists the packaged scripts,
links the controller source and exact installed library source.

The container provides Wine so [ras-commander][rc] can call installed Windows
HEC-RAS on Linux. It reads the model from a host folder mounted into the container
and writes preprocessing files through that same mount.

For Andy's review, start with [Wine, installation paths, and host data](OPERATION.md#wine-the-hec-ras-installation-and-host-data)
and the [API call sequence](OPERATION.md#exact-ras-commander-call-sequence).
[How to Re-Create this Container](PREPARATION.md) is a separate build guide with
installation notes, dependencies, and verification instructions. The
[release verification record](https://github.com/gpt-cmdr/ras-commander/blob/codex/container-precompute-linux/containers/hecras-unsteady/RELEASE-CURRENT.md) contains test evidence.

The matching [native Linux image](https://github.com/gpt-cmdr/ras-commander/blob/codex/container-precompute-linux/containers/hecras-unsteady/README.md) completes the unsteady
calculation after preprocessing. The [Python operating guide](https://github.com/gpt-cmdr/ras-commander/blob/codex/container-precompute-linux/docs/user-guide/container-execution.md) and
[notebook](https://github.com/gpt-cmdr/ras-commander/blob/codex/container-precompute-linux/examples/512_docker_precompute_and_linux_compute.ipynb) show both stages through [ras-commander][rc].
See the [current release record](https://github.com/gpt-cmdr/ras-commander/blob/codex/container-precompute-linux/containers/hecras-unsteady/RELEASE-CURRENT.md) for full 266-hour Linux and one-hour Windows Docker Desktop qualification.

## Models created on Windows or Linux

Pull the matching current image before running. Use a disposable model copy and mount its referenced terrain and projection
folders. Relative paths are resolved inside the container: `..\source_terrain`
and `..\projection` must lead to mounted folders. Mounting a model folder alone
does not expose its siblings. The terrain HDF must also have access to its
referenced TIFF files.

Use `--user root` for the qualified Windows Docker Desktop route described below.

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

**Qualification status:** All three matching versions passed the full 266-hour Linux sample and the one-hour Windows Docker Desktop notebook, using two CPUs per container. Linux results contained 267 output times; Windows results contained two, with 6,548 finite water-surface values at every time. Live progress, resume and sequential batch checks passed on both hosts; the six Linux Wine LF/CRLF cases also passed. The images are published on Docker Hub, and anonymous pulls verified all six matching payloads.

| HEC-RAS | Image | Status |
|---|---|---|
| 6.5 | `rascommander/hec-ras-wine-precompute_6.5:v4` | Published; Linux and Windows sample passed. |
| 6.6 | `rascommander/hec-ras-wine-precompute_6.6:v4` | Published; Linux and Windows sample passed. |
| 7.0.1 | `rascommander/hec-ras-wine-precompute_7.0.1:v4` | Published; Linux and Windows sample passed. |

The `v4` and `latest` tags identify each published Wine image. Each image includes its matching installed runtime and saved TCU acceptance. Normal jobs require no separate HEC-RAS installation or external runtime mount.

## Run on a Windows drive

The following PowerShell command mounts the model and both dependency folders.
Set `$models` to the parent folder containing all three and `$name` to the
project filename without `.prj`. For HEC-RAS 6.6 or 7.0.1, use its matching image.

```powershell
$models = 'E:\mac_test_output_Sept12\02_model_copies'
$name = '1919912_wb-2427466_wb-2427467_14-hr_100-cfs_to_11609-cfs'
docker pull rascommander/hec-ras-wine-precompute_6.5:v4
docker run --rm --pull always --user root --cpus 2 --mount "type=bind,src=$models\$name,dst=/job" --mount "type=bind,src=$models\source_terrain,dst=/source_terrain,readonly" --mount "type=bind,src=$models\projection,dst=/projection,readonly" rascommander/hec-ras-wine-precompute_6.5:v4 prepare --project "/job/$name.prj" --plan 01 --num-cores 2 --timeout 900 --replace-generated
```

The published 6.5, 6.6 and 7.0.1 images passed the one-hour notebook on Windows Docker Desktop using root and a host folder containing spaces. Both stages used two CPUs, retained the prepared temporary HDF, and produced two output times with 6,548 finite water-surface values each. See the [current release record](https://github.com/gpt-cmdr/ras-commander/blob/codex/container-precompute-linux/containers/hecras-unsteady/RELEASE-CURRENT.md) for the full 266-hour Linux results and the one-hour Windows scope.
From a WSL shell, use Linux source paths such as `/mnt/c/...`; source paths must be visible to the selected Docker engine.

## Run on a Linux filesystem

```bash
docker pull rascommander/hec-ras-wine-precompute_6.5:v4

docker run --rm --pull always --network none --read-only \
  --user 1000:1000 --cpus 2 --memory 6g --memory-swap 6g \
  --cap-drop ALL --security-opt no-new-privileges \
  --tmpfs /tmp:rw,nosuid,nodev,size=512m,mode=1777 \
  --env RAS2FIM_JOB_ROOT=/job \
  --mount type=bind,src=/shared/model,dst=/job \
  --mount type=bind,src=/shared/source_terrain,dst=/source_terrain,readonly \
  --mount type=bind,src=/shared/projection,dst=/projection,readonly \
  --mount type=volume,dst=/run/ras-job \
  --device /dev/ntsync \
  rascommander/hec-ras-wine-precompute_6.5:v4 \
  prepare --project /job/model.prj --plan 01 \
          --num-cores 2 --timeout 900 --run-id example-001 --replace-generated
```

`/shared/model` is a disposable folder on the Linux Docker host, writable by
UID 1000. The container sees it as `/job`; generated files remain in that same
host folder after the container exits. Mount terrain/projection dependencies
where the model expects them. Use a unique run ID. `--replace-generated` can
remove an existing final plan HDF, so preserve wanted results elsewhere.

The Linux recipe uses `/dev/ntsync` when the host supplies it. Windows Docker Desktop qualification used the simpler root-user command without passing this device. Read and agree to the [HEC-RAS terms](https://www.hec.usace.army.mil/confluence/rasdocs/rasum/6.6/terms-and-conditions-of-use)
before use. [init_ras_project()][init] uses `accept_tcu=True` with the prepared
profile's saved acceptance state.

For project `model`, plan 01 and geometry 01, the required outputs are
`model.p01.tmp.hdf`, `model.b01`, and `model.x01`. The receipt is
`.ras-commander/runs/<run-id>/prepare.json` under the model folder. Exit codes
are 0 for success, 1 for a recorded failure, and 2 for configuration/I/O errors.
Step 2b waits for the full batch before staging the next calculation stage.

## CPU settings

The `prepare --num-cores N` option defaults to 2 and accepts integers from 1 to
8. Set Docker's `--cpus N` to the same count: Docker limits the container's CPU
time, while [RasPlan.set_num_cores()][plan-cores] and
[RasPlan.set_2d_flow_options()][plan-2d] write the selected count into the working
plan before HEC-RAS starts. The 2D call uses `include_default=True` to cover the default block and every
named mesh, inserting missing 2D settings; the worker reads the value back before
preprocessing. Each container still prepares one plan; batch scheduling happens on the
host. The receipt records the count as `arguments.num_cores`.

The current two-core setting passed Linux and Windows qualification. The commands above pass the same number to Docker and the worker.

Step 2b reads `int_prepare_cores_per_job` from `[03_run_hec_ras]`, defaults to 2,
and accepts 1 through 8. It passes the count to both Docker and the worker and
checks the receipt before staging the generated files.

## Build and inspection

Follow [How to Re-Create this Container](PREPARATION.md) to export a finalized
Wine profile and build with the external `hecras_runtime` context. The build
currently depends on retained installed profiles; a complete empty-prefix
installation procedure remains unfinished. The [inventory](runtime-inventory-current.json)
records installed software, and [release verification](https://github.com/gpt-cmdr/ras-commander/blob/codex/container-precompute-linux/containers/hecras-unsteady/RELEASE-CURRENT.md)
records exact image identities and test results.

Docker Hub documentation: [6.5](dockerhub/6.5.md),
[6.6](dockerhub/6.6.md), [7.0.1](dockerhub/7.0.1.md).

## Optional runtime override

Leave `str_wine_profile` blank to use the bundled installation. Advanced users
can supply a matching controlled profile at `/runtime/wine-seed` read-only.
It must declare the selected HEC-RAS version and matching [ras-commander][rc]
build in its [runtime manifest](https://github.com/gpt-cmdr/ras2fim-2d/blob/dc60b219091e85bcb4564eca45313475c39ce58a/containers/hecras-prepare/runtime-manifest.example.json). The default
job uses the installation already included in the image.

[rc]: https://rascommander.info/ras/
[rc-github]: https://github.com/gpt-cmdr/ras-commander
[ras-prj]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/RasPrj.py
[init]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/RasPrj.py
[plan-path]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/RasPlan.py
[clear-geom]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/geom/GeomPreprocessor.py
[run-flags]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/RasPlan.py
[plan-cores]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/RasPlan.py
[plan-2d]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/RasPlan.py
[preprocess]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/RasPreprocess.py
[tcu-status]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/RasTcu.py
[tcu-accept]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/RasTcu.py
[bco]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/RasBco.py
[logging]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/RasBco.py
[monitor]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/RasBco.py
[terminate]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/RasPreprocess.py
[result]: https://github.com/gpt-cmdr/ras-commander/blob/9e4217713e954236b0c16023e1815c6f2b7a5309/ras_commander/ComputeResults.py

[build-checks]: https://github.com/gpt-cmdr/ras2fim-2d/blob/dc60b219091e85bcb4564eca45313475c39ce58a/containers/hecras-prepare/model_checks.py
