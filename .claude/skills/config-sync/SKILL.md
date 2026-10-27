---
name: config-sync
description: Add, retire, or rename one key in every ras2fim-2d config variant (src/config_global.ini and src/config_global_windows.ini) without disturbing comments, decorative rule lines, key order, blank lines or line endings. Reports the files touched and every line of code that reads the key. Use when a config key or [NN_section] value must change across the configs.
argument-hint: "add|retire|rename [NN_section] <key> [value or new_key]"
metadata:
  contributed-by: "CLB Engineering for Andy Carter, PE, 2026.10.27"
---

# /config-sync for ras2fim-2d

Follow AGENTS.md if present. These rules apply even without it:
- A key change goes into **every** config variant, or you say which variant was skipped and why.
- Never delete a line. Retired keys and old values stay in as commented-out lines.
- Never rename a key unless Andy asks for that rename in so many words. Keys are interfaces. This
  includes misspelled ones such as `flt_allowed_varianvce` in `[07_prep_netcdf]`.
- Don't reorder, re-indent, re-space or "tidy" any line you were not asked to change.
- Don't commit. Show the diff and propose a short lowercase commit message; Andy commits.

## The config variants

List them fresh each time (`git ls-files 'src/*.ini'`); today they are:

| File | Used for |
|---|---|
| `src/config_global.ini` | Linux Docker image (`-c config_global.ini`), paths under `/global_input/` |
| `src/config_global_windows.ini` | Windows run of step 02 / 02b, paths on `E:\ras2fim-2d-windows-resources\` |

Both are LF, ASCII, and end **without** a final newline. Some values differ on purpose
(`flt_min_q` 100.0 vs 50, `flt_max_allowed_gradient` 0.009 vs 0.027, the HEC-RAS controller).
When a new value is a path or a tuned number, ask for the value per variant.

Sections are numbered by step, each framed by its own rule line (`# >>>>>>`, `# ,,,,,,`,
`# ......`, `# ======`), with two blank lines between sections: `[global_input_files]`,
`[01_nextgen_hydrofabric]`, `[02_copy_geom]`, `[03_run_hec_ras]`, `[04_time_gradient]`,
`[05_raster_engine]` (empty), `[06_output_products]`, `[07_prep_netcdf]`, `[08_filter_netcdf]`.
Keys carry type prefixes (`str_`, `int_`, `flt_`, `b_`) because the step reads them with
`section.getfloat('flt_threshold', 0)`, `getint`, `getboolean` or `get`.

## Never write the configs back with configparser

`configparser.write()` drops every comment and rule line. Edit the files as text:
1. **Preferred:** the bundled script `${CLAUDE_SKILL_DIR}/config_sync.py` (standard library only,
   Python 3.8+, works in the `ras2fim2d` conda env or any Python). It reads bytes, keeps each
   file's line endings, encoding and missing final newline, refuses mixed-ending files, prints a
   diff, and checks that only the intended lines changed.
2. **Otherwise:** exact-text edits with the Edit tool, one file at a time, then the byte check
   below.

## Steps

1. **Pin down the request.** Action (add / retire / rename), section, key, value(s), and a short
   comment line in Andy's voice (lowercase is fine, e.g. `# minimum depth to keep in feet`).
   Match the section's key style and prefixes.
2. **Find every reader of the key** before changing anything:
   ```
   grep -rn -i "<key>" src/ README.md Dockerfile
   ```
   Also look for the section name (`'08_filter_netcdf'`) so you know which step reads it.
3. **Check the current state** of every variant:
   ```
   python "${CLAUDE_SKILL_DIR}/config_sync.py" -t show -f "src/config_global*.ini" -s <section> -k <key>
   ```
4. **Dry run, then write.** The script never writes without `-w True`.
   - Add (after the last key in the section, or after `-a <key>`):
     ```
     python "${CLAUDE_SKILL_DIR}/config_sync.py" -t add -f "src/config_global*.ini" -s 08_filter_netcdf -k flt_min_depth -v 0.1 -m "minimum depth to keep (ft)"
     ```
     For different values per variant, run once per file with `-f <that file>`.
   - Retire (comments the line out as `#key = value`; `-m` adds `# Retired - YYYY.MM.DD -- why`):
     `-t retire -s <section> -k <key> -m "<why>"`
   - Rename (only when asked): `-t rename -s <section> -k <old> -n <new> -m "<why>"`. The old line
     stays as a comment and the new key gets the same value on the next line.
   - A whole new section: give the block, in the file's style, on stdin with `-b -`
     (no temp file needed):
     ```
     python "${CLAUDE_SKILL_DIR}/config_sync.py" -t section -f "src/config_global*.ini" -s 09_check_netcdf -a 08_filter_netcdf -b - <<'EOF'
# .................
[09_check_netcdf]

# list of required global attributes
str_required_attrs = 00_stream_id
# .................
EOF
     ```
   Show Andy the dry-run diff. Write with `-w True` once it is what he asked for.
5. **Retire and rename need extra care.**
   - Retire: if any step still reads the key, stop and say so. ras2fim-2d steps read with a
     default (`section.getfloat('flt_threshold', 0)`), so a retired key does not crash; the step
     **silently** uses that default. Quote the default from the grep.
   - Rename: warn first that the old name lives on in configs outside this repo (TACC, other
     machines, older Docker images), then list every reader from step 2. Changing those readers is
     a Python edit: show each change, add a `# Revised - YYYY.MM.DD -- <what>` line under the
     file's header, and edit only after Andy agrees.
6. **Byte check** (required after Edit-tool edits; the script prints it for you):
   ```
   git diff --stat -- src/*.ini
   git diff -- src/*.ini
   python -c "import sys;[print(f, open(f,'rb').read().count(b'\r\n'), 'CRLF', open(f,'rb').read()[-1:]) for f in sys.argv[1:]]" src/config_global*.ini
   ```
   Expect 0 CRLF in these files and no final newline (last byte is not `\n`). When a key is added
   after the very last line, that old last line shows as changed in `git diff` only because it
   gains a line ending; its text is unchanged.

## Report

- Each variant: added / retired / renamed / skipped (with the reason, e.g. section missing).
- The grep of readers: file, line, and how it reads the key (`getfloat`, default value).
- Anything Andy must still decide (values per variant, readers to update).
- A proposed commit message, e.g. `config: add flt_min_depth to 08_filter_netcdf`.

## Supporting file

- `config_sync.py`: the byte-preserving editor used above. `python config_sync.py -h` lists the
  options.
