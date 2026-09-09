# Linux/Wine HEC-RAS preprocessing

The `hecras-prepare` image runs HEC-RAS geometry preprocessing through
ras-commander and WineHQ Stable 11.0 on Linux. Each invocation prepares one
existing HEC-RAS plan and writes the three files consumed by the repository's
compute scripts:

- `<project>.p01.tmp.hdf`
- `<project>.b01`
- `<project>.x01`

The later unsteady simulation continues to use the existing
`civileng127/ras_v65:v0` compute image.

This Linux/Wine preprocessing integration was contributed by CLB Engineering
Corporation.

## Published images

Each Linux/amd64 image requires a runtime profile containing the same HEC-RAS
version.

| HEC-RAS | Docker Hub | Immutable image |
|---|---|---|
| 6.5 | [hec-ras-wine-precompute_6.5](https://hub.docker.com/r/rascommander/hec-ras-wine-precompute_6.5) | `rascommander/hec-ras-wine-precompute_6.5@sha256:0eb3aa3dbd1cb174651bf74196ac66b448acd0fb4dbdcb5852d664da521ea7ca` |
| 6.6 | [hec-ras-wine-precompute_6.6](https://hub.docker.com/r/rascommander/hec-ras-wine-precompute_6.6) | `rascommander/hec-ras-wine-precompute_6.6@sha256:4aa00301d6587f6fc073418655c585517e2a61a363c4898c144a21099f093c00` |

## Runtime profile

HEC-RAS, Windows Python, and the ras-commander wheel are supplied at run time in
a read-only Wine profile. The
[runtime manifest example](runtime-manifest.example.json) records their paths,
the HEC-RAS version, the ras-commander source commit and wheel hash, and hashes
for the executable artifacts. The image verifies these values before starting
HEC-RAS and clones the profile into private writable scratch for each job.

Create the profile with the official HEC-RAS installer and accept its displayed
Terms and Conditions for Use. The retained profile must contain that accepted
user state. The worker calls `init_ras_project(..., accept_tcu=True)` so
ras-commander transfers the state to the installed version and verifies it
before preprocessing.

Keep the runtime profile, installers, model data, and solver results outside the
repository and container image.

## Command-line interface

The container entry point exposes one command:

```text
hecras-prepare prepare --project PATH [--plan 01] [--timeout 300]
                       [--run-id ID] [--replace-generated]
```

| Option | Meaning |
|---|---|
| `--project PATH` | HEC-RAS `.prj` file inside the writable `/job` mount. |
| `--plan NUMBER` | Plan number to preprocess. Default: `01`. |
| `--timeout SECONDS` | Maximum HEC-RAS preprocessing time. Default: `300`. |
| `--run-id ID` | Optional receipt identifier. An identifier is generated when omitted. |
| `--replace-generated` | Clear and replace preprocessing files in a disposable project copy. |

A successful job writes its receipt to
`.ras-commander/runs/<run-id>/prepare.json` beneath the mounted project folder.
The receipt identifies the runtime and records the size and SHA-256 hash of each
output file.

```bash
docker run --rm --network none --read-only \
  --user 1000:1000 --cpus 2 --memory 6g \
  --tmpfs /tmp:rw,nosuid,nodev,size=512m,mode=1777 \
  --env RAS2FIM_JOB_ROOT=/job/project \
  --mount type=bind,src=/shared/model,dst=/job/project \
  --mount type=bind,src=/controlled/hecras-6.5-wine,dst=/runtime/wine-seed,readonly \
  --mount type=volume,dst=/run/ras-job \
  --device /dev/ntsync \
  rascommander/hec-ras-wine-precompute_6.5@sha256:0eb3aa3dbd1cb174651bf74196ac66b448acd0fb4dbdcb5852d664da521ea7ca \
  prepare --project /job/project/model.prj --plan 01 \
          --timeout 900 --run-id example-0001 --replace-generated
```

Omit `--device /dev/ntsync` when NT synchronization is unavailable. The Windows
Step 2b controller builds this command, starts the configured jobs in parallel,
waits for the complete batch, validates every receipt and hash, and then stages
the three output files for the existing compute scripts.

## Build

Build from the repository root. Licensed software is not part of the build
context.

```bash
docker buildx build --load --platform linux/amd64 \
  --build-arg HEC_RAS_VERSION=6.5 \
  --build-arg RAS_COMMANDER_COMMIT=9e4217713e954236b0c16023e1815c6f2b7a5309 \
  --build-arg RAS_COMMANDER_WHEEL_SHA256=dc0c2f9baa9db66afee01e34eb00d7a04f53e7b82c3340c786b2e9aef1795932 \
  --tag hecras-prepare:6.5 \
  --file containers/hecras-prepare/Dockerfile .
```

Use `HEC_RAS_VERSION=6.6` for the 6.6 image. At run time, mount the matching
profile and select an image pinned by its digest.
