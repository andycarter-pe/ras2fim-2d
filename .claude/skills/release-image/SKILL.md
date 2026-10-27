---
name: release-image
description: Build a new ras2fim-2d Docker image (civileng127/ras2fim2d with a dated YYYYMMDD tag, or the next vNN if Andy prefers), smoke-check it, propose the README docker pull/run and release-date changes as a diff, and stop to ask before any docker push. Andy starts this by typing /release-image; Claude never starts it on its own.
argument-hint: "[YYYYMMDD | vNN]"
disable-model-invocation: true
metadata:
  contributed-by: "CLB Engineering for Andy Carter, PE, 2026.10.27"
---

# /release-image for ras2fim-2d

Build and check a new image, update the README after Andy approves the diff, then stop before
the push. Work in this order and do not skip a stop.

Follow AGENTS.md if present. These rules apply even without it:
- **Never run `docker push`, `docker login` or `git push` unless Andy says "push" for this tag in
  this session.** A yes to the README diff is not a yes to push. Never push a `latest` tag.
- Never edit the Dockerfile, `Dockerfile_old`, code, configs or `environment_ras2fim2d.yml`. The
  only file this skill changes is README.md, and only after Andy approves the diff.
- Do not commit. Stage README.md and propose a short lowercase message; Andy commits.
- Report exactly what you ran. If a build or check fails, stop and show the error.

## 1. Learn the current scheme
- From README.md read every `civileng127/ras2fim2d:<tag>` (today `v01`): the setup bullet, the
  `docker pull` in Installation step 3 and its description bullet, the Phase 1 and Phase 4
  `docker run` blocks, and the "Docker (Reference)" sentence. Also read `**Release date:**`
  (written `YYYY-MM-DD`) and `**Version:**` under "Status".
- Read the published tags:
  `curl -s "https://hub.docker.com/v2/repositories/civileng127/ras2fim2d/tags?page_size=100"`.
  Also note the older repository `civileng127/ras2fim-2d` (one tag, `20250209`), which matches
  the Dockerfile header comment `civileng127/ras2fim-2d:v1`.
- New tag: the argument if given (`$ARGUMENTS`). If none, ask once, showing both choices:
  today's date as `YYYYMMDD` (Andy's scheme for lisflood2fim and the older ras2fim-2d image) or
  the next version after the README's (`v02` after `v01`). Suggest `YYYYMMDD`.
- If the tag already exists on Docker Hub, stop and ask: overwrite it, or pick another?
- Leave `civileng127/ras_v65` alone. It is the HEC-RAS worker image, built elsewhere.

## 2. Check that the build will contain what Andy expects
The Dockerfile does not copy this folder into the image. It runs
`git clone https://github.com/andycarter-pe/ras2fim-2d.git /ras2fim-2d` and downloads
`LimitingDischarge.csv` and `nextgen_12.gpkg` from the `ras2fim-2d-global-inputs` S3 bucket into
`/global_input`. So the image gets what is on GitHub `main`, not what is on disk.
- `git fetch origin`, then report `git status --short` and `git log --oneline origin/main..HEAD`.
- If `src/` has uncommitted or unpushed changes, stop and tell Andy they will not be in the
  image until he pushes them. Do not push for him.
- Docker caches the `git clone` layer, so a cached build can reuse an old clone. Build with
  `--no-cache` unless Andy says a cached build is fine.

## 3. Build
Show Andy the command and the tag, then run it from the repo root. Name the Dockerfile so
`Dockerfile_old` is never used by mistake:
```
docker build --no-cache -f Dockerfile -t civileng127/ras2fim2d:<tag> .
```
Run it in the background with output going to a log outside the repo (for example
`<temp>/release_image_ras2fim2d_<tag>.log`). Report the last 20 lines if it fails.

## 4. Smoke-check the image
Run each check and record the result:
```
docker run --rm civileng127/ras2fim2d:<tag> bash -c "python --version \
 && python -c 'import rasterio, geopandas, fiona, rioxarray, h5py, scipy, networkx, netCDF4, h5netcdf, tqdm; print(\"imports OK\")' \
 && which nccopy && ls -la /global_input && python ras2fim-2d.py -h | head -3 \
 && git -C /ras2fim-2d log -1 --format='%h %ad %s' --date=short"
```
- Python must be 3.8.12 (pinned in the Dockerfile). `/global_input` must hold
  `LimitingDischarge.csv` and `nextgen_fabric.gpkg`, both non-empty.
- `tqdm` is imported by `ras2fim-2d.py` but is not listed in the Dockerfile. If the import
  fails, stop and tell Andy; do not change the Dockerfile.
- The commit printed last must equal `git rev-parse --short origin/main`. If not, stop.
- Report `docker image ls civileng127/ras2fim2d:<tag>` (ID and size), and the size of the
  previous tag if it is on this machine.
- Offer a sample run on the new image before pushing (`/run-demo docker "(4,8)"` after the README
  points at the new tag). Do not start it without a yes.

## 5. Propose the README change, then apply it only after a yes
Change only these:
- Every `civileng127/ras2fim2d:<old>` listed in section 1, to the new tag.
- `**Release date:**` to today, in the same `YYYY-MM-DD` form. Ask before changing
  `**Version:**`; leave it if Andy does not say.

Write the proposed file to a temp path outside the repo (`<temp>/README_ras2fim2d_<tag>.md`) and show `diff -u README.md <temp>`.
README.md uses LF line endings; keep them. Ask: "Apply this README change?" Copy it over
README.md only after a yes. Then run `git diff --stat README.md` and `git add README.md`, and
propose a commit message in Andy's style, for example `docker container to <tag>`.

## 6. Stop before the push
Show Andy, and then stop and wait:
```
Built:   civileng127/ras2fim2d:<tag>  (<image id>, <size>)  from GitHub main <hash>
Checks:  Python <version> | imports OK | nccopy | /global_input OK | ras2fim-2d.py -h
README:  updated and staged | not changed
To publish, the command is:
  docker push civileng127/ras2fim2d:<tag>
Push now? (yes / no)
```
Run `docker push` only if Andy answers yes to this question in this session. After a push, check
the tag appears with the curl from section 1 and report it. Then remind Andy to commit and push
the README himself.
