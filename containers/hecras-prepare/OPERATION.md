# What happens inside the HEC-RAS preprocessing container

[Installed code and GitHub source](INSTALLED-CODE.md) lists the packaged scripts,
links the current container source snapshot and public library references.

Prepared for Andy's review. CLB Engineering Corporation, September 10, 2026.

The container provides Wine so [ras-commander][rc] can call installed Windows
HEC-RAS on Linux. It takes an existing model, prepares inputs for a separate
calculation stage, and leaves those files in the host's mounted working folder.

The bundled 6.5, 6.6 and 7.0.1 images are published with installed runtimes
and saved TCU acceptance. The operating sequence is the same
for each selected runtime. The path example below uses 6.5; other HEC-RAS
releases use their matching installation directories.

The matching [7.0.1 native Linux image](https://hub.docker.com/r/rascommander/hec-ras-linux-unsteady_7.0.1) completes the unsteady
calculation after preprocessing. The [Python operating guide](https://github.com/gpt-cmdr/ras-commander/blob/codex/container-precompute-linux/docs/user-guide/container-execution.md) and
[notebook](https://github.com/gpt-cmdr/ras-commander/blob/codex/container-precompute-linux/examples/512_docker_precompute_and_linux_compute.ipynb) show both stages through [ras-commander][rc].
See the [7.0.1 release record](https://github.com/gpt-cmdr/ras-commander/blob/codex/container-precompute-linux/containers/hecras-unsteady/RELEASE-7.0.1-20260911.md) for full Linux and Windows Docker
Desktop qualification.

## Models created on Windows or Linux

Pull the current image before testing (`docker pull rascommander/hec-ras-wine-precompute_6.5:v4`, or the matching 6.6 or 7.0.1 image). Use a disposable model copy and mount its referenced terrain and projection
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


## Wine, the HEC-RAS installation, and host data

[Wine](https://www.winehq.org/) supplies Windows application interfaces on Linux. It lets Windows Python
and the installed Windows HEC-RAS programs run inside this Linux container.
[ras-commander][rc] runs in Windows Python and calls HEC-RAS. HEC-RAS performs
the preprocessing; the container supplies the operating environment.

### Where HEC-RAS is installed

A Wine *prefix* is a directory holding a Windows-style C: drive, installed
programs, supporting DLLs, and Windows registry files. The image contains a
prepared prefix at `/runtime/wine-seed/prefix` and its settings in
`/runtime/wine-seed/runtime.json`.

| View | Installed executable path |
|---|---|
| Windows path seen by HEC-RAS and Python | `C:\Program Files (x86)\HEC\HEC-RAS\6.5\Ras.exe` |
| Linux path inside the image template | `/runtime/wine-seed/prefix/drive_c/Program Files (x86)/HEC/HEC-RAS/6.5/Ras.exe` |
| Linux path used by one running job | `/run/ras-job/<run-id>/wineprefix/drive_c/Program Files (x86)/HEC/HEC-RAS/6.5/Ras.exe` |

At startup the wrapper copies the prepared prefix into the job's private
scratch directory and sets `WINEPREFIX` to that copy. Wine can update its
registry and temporary state in that private copy. The worker reads the
installed template when making the copy; normal job updates use the copy.
Windows Python is `C:\Python311\python.exe`; Xvfb supplies a virtual screen
for the Windows application. No interactive desktop is needed for a prepared job.

### How inputs and outputs cross the host mount

See the [Docker bind-mount reference](https://docs.docker.com/engine/storage/bind-mounts/) for host-path behavior.

For example, `--mount type=bind,src=/shared/model,dst=/job` makes the Docker
host's `/shared/model` directory visible inside the container as `/job`.
These are **the same underlying files**. There is no upload at startup or
separate download at completion: opening `/job/model.prj` reads the host file,
and writing `/job/model.p01.tmp.hdf` writes into the host directory.

Wine's Z: drive maps the container's Linux filesystem, so `winepath -w` can
translate `/job/model.prj` to `Z:\job\model.prj`. That drive sees the container's
filesystem and mounted folders. It does not expose unmounted host directories.

```mermaid
flowchart LR
    H["Host: /shared/model"] <-->|"Docker bind mount: same files"| J["Container: /job"]
    J <-->|"Wine Z: drive mapping"| R["HEC-RAS: Z:/job/model.prj"]
    R --> O["Writes .p01.tmp.hdf, .b01 and .x01"]
    O --> J
    J --> P["Files remain in /shared/model on the host"]
```

| Location | What happens to its files |
|---|---|
| Host model folder mounted at `/job` | Read/write. Inputs, edited plans, generated files, logs, and receipts remain on the host after the container exits. |
| Referenced terrain/projection mounts | Usually read-only; their mounted paths must match the model's references. |
| `/runtime/wine-seed` | Installed template included in the image; the worker copies it for each job. |
| `/run/ras-job/<run-id>` | Private Wine prefix and intermediate worker files in the container's writable layer. |
| `/tmp` | Temporary storage in the container's writable layer for the command below. |

The bind-mount `src` path belongs to the machine running the Docker daemon.
When Docker is started over SSH, use that Linux host's path. A Windows
controller can access the same data through a network share. Its path and the
Linux host path must refer to the same shared folder.

On Linux filesystems, use a disposable model copy writable by UID 1000.
For Windows drive mounts, use `--user root`; see the [Windows command](README.md#run-on-a-windows-drive). The job edits preprocessing
settings and can replace generated files, including an existing final plan HDF
when `--replace-generated` is set. The documented `--rm` run removes the
container and its writable layer; the bind-mounted model folder remains.

## Exact ras-commander call sequence

The Windows worker uses these [ras-commander][rc] APIs. Every library function
shown in the diagram has a source link in the tables below. Links point to
public upstream references; the arguments and behavior described here were
checked against the installed 0.99.2 build. Its retained build identity is
recorded in the reconstruction section.

```mermaid
sequenceDiagram
    participant Host as Host model folder mounted at /job
    box Linux controller in container
        participant L as prepare.py
    end
    box Wine environment in same container
        participant W as Windows Python worker
        participant C as ras-commander
        participant H as HEC-RAS
    end
    L->>L: Select runtime and copy private WINEPREFIX
    L->>W: xvfb-run + wine + Windows Python
    W->>W: Check model dependencies and normalize inputs to CRLF
    W->>C: RasPrj()
    W->>C: init_ras_project(..., accept_tcu=True)
    C->>C: RasTcu.status() and RasTcu.accept() if needed
    W->>C: RasPlan.get_plan_path(...)
    W->>C: GeomPreprocessor.clear_geompre_files(...)
    W->>C: RasPlan.update_run_flags(..., geometry_preprocessor=True)
    W->>C: RasPreprocess.preprocess_plan(...)
    C->>C: BcoMonitor.enable_detailed_logging(plan_file)
    C->>H: Ras.exe -c with quoted project and plan paths
    H->>Host: Read model and write preprocessing files through mount
    C->>C: BcoMonitor(...).monitor_until_signal(process)
    Note over C,H: BCO signal or owned unsteady process + fresh files
    C->>H: RasPreprocess._terminate_process_tree(process) if running
    C-->>W: PreprocessResult
    W->>W: Validate 2D mesh and hydraulic tables in both HDF files
    W-->>L: worker-result.json
    L->>Host: Confirm output files and write prepare.json through mount
```

If Mermaid is unavailable in your viewer, the sequence is: Linux wrapper ->
Windows Python under Wine -> [ras-commander][rc] -> HEC-RAS -> files in the host
mount. The worker returns its result to the wrapper, which writes the job receipt.

| Direct worker API, in order | Actual arguments and purpose |
|---|---|
| [RasPrj()][ras-prj] | Creates the project object passed to subsequent calls as `ras_object`. |
| [init_ras_project()][init] | `project`, `ras_version=ras_executable`, `ras_object=ras_object`, `load_results_summary=False`, `load_hdf_metadata=False`, `hide_intro=True`, `accept_tcu=True`. Selects the installed executable and initializes project data. |
| [RasPlan.get_plan_path()][plan-path] | `plan, ras_object=ras_object`. Resolves the selected plan file; a missing plan stops the worker. |
| [GeomPreprocessor.clear_geompre_files()][clear-geom] | `plan_path, ras_object=ras_object`. Clears the matching geometry preprocessing cache and refreshes geometry metadata. |
| [RasPlan.update_run_flags()][run-flags] | `plan_path, geometry_preprocessor=True, ras_object=ras_object`. Enables geometry preprocessing in the working plan. |
| [RasPreprocess.preprocess_plan()][preprocess] | `plan, ras_object=ras_object, max_wait=timeout, clear_existing=replace_generated, fix_line_endings=True`. Starts and supervises HEC-RAS and returns its preprocessing result. |

| Calls and result used inside the library | Role |
|---|---|
| [RasTcu.status()][tcu-status] | Checks saved Terms and Conditions for Use acceptance for the selected executable. |
| [RasTcu.accept()][tcu-accept] | Called during initialization only when status is explicitly unaccepted and `accept_tcu=True`. Prepared images already have verified accepted state. |
| [BcoMonitor()][bco] and [BcoMonitor.enable_detailed_logging()][logging] | Configure calculation logging and create the readiness monitor. |
| [BcoMonitor.monitor_until_signal()][monitor] | Waits for the calculation-log signal or the owned-process/file readiness condition. |
| [RasPreprocess._terminate_process_tree()][terminate] | Internal cleanup helper; terminates the launcher and its descendants when still running. |
| [PreprocessResult][result] | Returns plan/geometry numbers, output paths, elapsed time, signal source, timeout state, and error information. |

The wrapper launches Windows Python with
`xvfb-run -a -s "-screen 0 1024x768x24" wine`, followed by the manifest's Python
path and translated worker arguments. Inside [RasPreprocess.preprocess_plan()][preprocess],
a Windows Python subprocess launches the full executable path as
`"Ras.exe" -c "<project.prj>" "<project.p01>"` with `shell=False`.

[BcoMonitor.monitor_until_signal()][monitor] watches for
`Starting Unsteady Flow Computations`. The alternate signal requires a
`RasUnsteady.exe` descendant of this job's launcher and three nonempty output
files changed from their prelaunch size/mtime baseline. The library uses
`psutil` for process ownership and cleanup. The readiness limit is `max_wait`;
the outer Windows worker subprocess has a separate `timeout + 120` second limit.

**The unsteady solver may start briefly before it is stopped.** This workflow
obtains preprocessing inputs for a separate calculation stage. A detected TCU
dialog blocks the job. The wrapper requires the expected plan, geometry and
output paths, an accepted `bco` or `owned_process_artifacts` readiness signal,
`timed_out=False`, `full_result_copied=False`, and successful HDF-content validation before reporting success.

## What remains on the host

For plan 01 and geometry 01, HEC-RAS produces:

| File in the mounted model folder | Purpose |
|---|---|
| `model.p01.tmp.hdf` | Temporary plan HDF produced during preprocessing. |
| `model.b01` | Plan binary input for the calculation stage. |
| `model.x01` | Geometry control input, converted to Linux LF line endings. |
| `.ras-commander/runs/<run-id>/prepare.json` | Job receipt recording success/failure, selected runtime, completion signal, and output details. |

Project names and plan/geometry numbers determine the filenames. Other working
files and calculation logs may also remain. Worker stdout/stderr are saved
beside the receipt when the worker subprocess returns normally; early validation
failures or forced timeouts can have less logging.

Exit code 0 means verified success, 1 a recorded job failure, and 2 a configuration
or I/O error. The Windows Step 2b controller waits for the full batch before it
stages files for the separate Linux calculation. Qualification here covers the
retained preprocessing sample; model suitability and compatibility with the
later native HEC-RAS 6.5 engine require their own checks.

## Run configuration and review sources

The supplied launcher runs as UID 1000 with two CPUs, 6 GiB RAM, networking
disabled, a read-only root, no added Linux capabilities, and no privilege
escalation. The qualified host supplies `/dev/ntsync` for Wine synchronization.
Hosts without it need separate qualification.

[Run instructions](README.md) provide the complete Docker command.
[How to Re-Create this Container](PREPARATION.md) covers installation, packaging,
software inventory, and build verification separately from this operating guide.
The [release verification record](RELEASE-VERIFICATION.md) contains the test evidence.

| Review file | Responsibility |
|---|---|
| [Dockerfile](Dockerfile) | Linux software, runtime packaging, user, and entrypoint. |
| [prepare.py](prepare.py) | Model request, runtime copy, worker launch, output checks, receipt. |
| [windows_worker.py](windows_worker.py) | The direct [ras-commander][rc] calls listed above. |
| [Step 2b controller](../../src/stage_hecras_for_linux_wine_02b.py) | Host mounts, Docker launch, batch completion, staging. |
| [Mermaid source](ras-commander-call-sequence.mmd) | Editable API sequence diagram. |

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
