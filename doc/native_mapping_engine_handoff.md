# Native mapping engine handoff: implemented and integrated

The requested caller-supplied `wsel_profiles=` API was implemented in
[RasCommander PR 710](https://github.com/gpt-cmdr/ras-commander/pull/710), merged
2026-10-27 as `d47ca7fb9188bdb810d95ab03245a0a676371965`.

This integration pins the immutable feature commit
`777c7746ceb35fc3b1c29d3718935c9013ae009c`; its development metadata is 0.110.1,
but the PyPI 0.110.1 release lacks this feature. Existing released fast-map 0.8.0
and fast-HDF 0.2.3 work without a Rust rebuild.

The consumer passes a wide table of exact `mesh_name`, mesh-local `cell_id`, and
named absolute WSEL columns. Every physical and ghost ID is present; explicit NaN
means dry. It binds source HDF and CRS attrs and passes the selected original or
healed TIFF through the existing terrain argument. Window aggregation stays in
ras2fim-2d. Existing calls/defaults and selector APIs remain available upstream.

See `healed_mapping_terrain.md` for the native adapter's complete-profile policy,
single-mesh limitation, existing rounding/window endpoint behavior, independent
terrain selection, output grid delivery, and qualification limits. The integration
uses plain sloping, without enhanced face weighting. No private package is needed.
