# Agent notes for containers/hecras-prepare (contributed by CLB Engineering, 2026.10.27)

The repo-wide house rules live in the root AGENTS.md (added separately). These notes apply to this folder
and to src/stage_hecras_for_linux_wine_02b.py, which drives the container.

- containers/hecras-prepare/ builds the optional image that runs step 02b HEC-RAS preprocessing on a
  Linux Docker host. src/stage_hecras_for_linux_wine_02b.py drives it. Docs: containers/hecras-prepare/README.md.
- It is used only when [03_run_hec_ras] b_prepare_on_linux_wine = True. The default path,
  stage_hecras_for_linux_v66_02b.py on the Windows machine, must keep working unchanged.
- The container code is CLB-maintained and written in a different style. Do not restyle it, and do not
  copy that style into src/.
- tests/ covers only this container path (python -m pytest tests). Keep HEC-RAS installers, Wine
  profiles, model data and results out of Git.
