# Preprocessing release verification — September 10, 2026

The 6.5 and 6.6 `v4` and `latest` images contain input normalization and HDF
validation. The 7.0.1 candidate passed the same checks locally; its binary
publication still awaits the separate approval recorded in AGENTS.md.

## Verified behavior

27 focused Python tests passed on Windows. Twelve real Docker/Wine regression
cases passed on the Linux worker: LF, CRLF, mixed line endings, and missing
terrain/projection mounts for each of 6.5, 6.6 and 7.0.1. The CRLF case ran as
the image's default UID 1000; the other cases ran as root to match Andy's command.

Each successful case retained all 6,548 cells in `Perimeter 1`, with cell-volume
and face-area elevation tables in both geometry and temporary plan HDF files.
Missing mounts produced a failed receipt and left all existing model files
unchanged. Dependency validation also checks raster files referenced by the
terrain HDF. The regression suite includes a missing-TIFF unit test.

| Runtime | LF temporary HDF | CRLF temporary HDF | Mixed temporary HDF | Missing mounts |
|---|---:|---:|---:|---|
| 6.5 | 2,680,668 bytes | 2,681,002 bytes | 2,680,668 bytes | Rejected without model changes |
| 6.6 | 2,681,108 bytes | 2,681,108 bytes | 2,681,108 bytes | Rejected without model changes |
| 7.0.1 | 4,676,101 bytes | 4,676,123 bytes | 4,676,101 bytes | Rejected without model changes |

The fixture was the retained repository 50-cfs sample, copied separately for
every run. This validates preprocessing, not a full hydraulic simulation or
Andy's exact 100-cfs model.

## Windows drive qualification on CLB-08

Eight additional cases passed against the published 6.5 and 6.6 image digests:
LF, CRLF, mixed inputs and missing dependencies for each version, all using
`--user root`. Model inputs and outputs were on Windows C:, bind-mounted from
`/mnt/c/Users/bill/ras2fim-reliability-20260910`. Docker Engine 29.1.3 ran inside
Ubuntu 26.04 LTS on WSL2, kernel `6.6.87.2-microsoft-standard-WSL2`, on Windows
11 Pro build 26200. Docker Engine was installed in the existing WSL distribution
for this test. Docker Desktop was not installed on CLB-08.

| Runtime | LF temporary HDF | CRLF temporary HDF | Mixed temporary HDF | Missing mounts |
|---|---:|---:|---:|---|
| 6.5 | 2,680,668 bytes | 2,681,002 bytes | 2,680,668 bytes | Rejected without model changes |
| 6.6 | 2,681,108 bytes | 2,681,108 bytes | 2,681,108 bytes | Rejected without model changes |

Independent inspection of all twelve successful geometry/temporary-plan HDF
files confirmed 6,548 cells and populated cell-volume and face-area elevation
tables. The actual spawning geometry writer was also executed under Linux
against a Windows-backed file and preserved CRLF. These runs did not pass a
`/dev/ntsync` device into the container.

The initial 6.5 matrix exposed a separate Windows-mount issue: UID 1000 timed
out with a hidden HEC-RAS `Run-time error '75': Path/File access error`, leaving
partial HDF files. A second UID-1000 run with LF inputs reproduced it. A root
run with CRLF passed on the same Windows mount, and UID 1000 passed when only
the writable project folder was moved to WSL's ext4 filesystem; terrain and
projection remained mounted from Windows C:. The failed runs are retained in
the evidence and are not counted as passes.

Use `--user root` for Windows drive mounts. UID 1000 remains qualified for
writable Linux model folders. A failure after HEC-RAS starts can leave partial
model outputs, so use disposable copies and require a successful, validated
receipt before staging computation. Windows Docker Desktop itself remains
unverified; WSL2 Docker with Windows-hosted data is the tested configuration.
HEC-RAS 7.0.1 has not been qualified on this Windows host.

### Repeat the Windows-drive cases

From WSL, run the checked-in runner with a fresh external work directory and a
source project whose parent has sibling `source_terrain` and `projection` folders:

```bash
python3 containers/hecras-prepare/test_preprocessing_inputs.py \
  --image rascommander/hec-ras-wine-precompute_6.5:v4 \
  --source-project /mnt/c/models/02_model_copies/model/model.prj \
  --work-dir /mnt/c/test-evidence/hecras-65-new \
  --container-user root
```

Repeat with the 6.6 image and another new work directory. Omitting
`--container-user` retains the Linux matrix's default-user CRLF case.
The runner records the resolved image ID, exact commands, receipts and each
case's outcome. Full evidence, including the initial failures and dialog
capture, remains outside Git under
`C:\Users\bill\ras2fim-reliability-20260910` on CLB-08.


## Implementation and build inputs

The final geometry and unsteady-flow writers explicitly use CRLF on either
operating system. Before HEC-RAS starts, the container checks model dependencies
and normalizes selected `.prj`, `.p##`, `.g##` and `.u##` inputs byte-for-byte apart
from their line endings. A success receipt requires populated 2D areas, preserved
cell counts and bounded, finite hydraulic-table data in both output HDF files.
Zero-length table rows for inactive cells are allowed; empty whole tables are not.
The generated `.x##` still uses LF for the separate Linux computation stage.

The source Dockerfile includes all three installed scripts. For this deployment,
the verified bundled runtime layers were reused and a source-only update layer
added after a full rebuild exceeded the worker's Docker storage capacity. The
exact update Dockerfile, base image identities, build commands and logs are
retained with the external release evidence. The full Dockerfile supports a
fresh build with the external runtime context on a host with sufficient space.

The release source is [commit 3c5011d](https://github.com/gpt-cmdr/ras2fim-2d/commit/3c5011dbe150d8311e45598401b16645c6e61932).
See [model_checks.py](model_checks.py), [test_preprocessing_inputs.py](test_preprocessing_inputs.py)
and [the reconstruction guide](PREPARATION.md). Runtime layers and vendor notices
were preserved. No installer, Wine profile, model or solver output was added to Git.

External evidence is retained under `/scratch/ras2fim-hecras/reliability-20260910`.
The current software inventory is [runtime-inventory-20260910.json](runtime-inventory-20260910.json).

## Published image identities

- 6.5: `sha256:50622b26c2c210c9034a38500aa8048e84c129dc0fbf2aa348fbc3c1e368c597` (`v4` and `latest`).
- 6.6: `sha256:b84b1cce0c73b8422c432a11167587a79a2bc2588b7a5d54904275dcf7af7e4a` (`v4` and `latest`).
