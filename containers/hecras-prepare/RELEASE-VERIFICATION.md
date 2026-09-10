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
Andy's exact 100-cfs model. Docker Desktop on the Windows workstation failed
during its own startup, so a Windows-host Docker Desktop end-to-end run is not
claimed. The same source-side regression runner supports Windows and Linux.

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
