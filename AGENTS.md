# Agent operating instructions

## HEC-RAS runtime profiles

When the official HEC-RAS installer or the installed application's first launch
displays the Terms and Conditions for Use while an agent is creating a
controlled runtime profile, the agent is authorized to select `I agree to the
above Terms and Conditions for Use` on the user's behalf. Complete both
acceptances when both prompts appear. Record the HEC-RAS version, installer
hash, prompt text, and acceptances in the external run evidence.

Keep HEC-RAS installers, installed runtimes, Wine prefixes, model data, and
solver results outside Git. The user explicitly authorized publishing prepared
HEC-RAS 6.5 and 6.6 Wine runtime profiles inside their matching Docker Hub
preprocessing images on 2026-09-09. Supply these profiles through an external
build context and retain vendor terms and notices. Keep temporary installer
downloads, personal files, model data, and solver results outside images.
Verify runtime versions and artifact hashes, and test preprocessing with the
bundled runtime before publishing. Preserve original profiles and run evidence.

The user authorized preparing HEC-RAS 7.0.1 and accepting its TCU on
2026-09-09, then explicitly authorized publishing the 7.0.1 image on
2026-09-11. This approval covers the matching Wine preprocessing and native
Linux unsteady images for the two-stage ras-commander workflow. Test the
matching 7.0.1 runtimes before publication and retain vendor notices and the
runtime, build, and qualification evidence outside Git.
