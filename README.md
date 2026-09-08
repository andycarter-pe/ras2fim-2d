# RAS2FIM-2D <img src="doc/Logo_CWE.png" align="right" alt="Center for Water and the Environment" height="80">

<br clear="right">

## *Flood Inundation Mapping using HEC-RAS 2D*

<img src="/doc/ras2fim2d-logo-20260907.png" align="right"
     alt="RAS2FIM-2D logo" width="160" height="160">

**RAS2FIM-2D** converts 2D HEC-RAS models into flood inundation mapping (FIM) libraries for National Water Model (NWM) NextGen stream segments represented within the HEC-RAS model's 2D computational area.

The workflow uses a "firehose" approach in which a series of constant-flow simulations is performed by applying flow to an internal HEC-RAS boundary condition named `Emitter1`. Flows are introduced near the upstream ends of the NWM NextGen hydrofabric catchments and the resulting inundation is clipped to the area associated with the corresponding NextGen `feature_id` (for example, `wb-2427467`).

The resulting FIM library contains water-surface elevations (WSEL) indexed by flow rate (cfs). The output NetCDF files also contain the terrain used by HEC-RAS to compute the water-surface elevations.

This project was developed in support of the National Weather Service under **Research Project NA22NWS4320003 / A25-0366-S018**.

<p align="center">
  <img src="/doc/ras2fim_animation.gif" alt="RAS2FIM-2D example" width="85%">
</p>

## Status

**Version:** 0.1 — Preliminary release  
**Release date:** 2026-09-07

## Technology

- Python 3.8.12
- HEC-RAS 2D
- Docker
- NetCDF
- GDAL / raster processing libraries

## Related Project

[NOAA-OWP/ras2fim](https://github.com/NOAA-OWP/ras2fim) — Flood Inundation Maps generated from HEC-RAS 1D models.

---

# HEC-RAS 2D Requirements

The base HEC-RAS model must meet the following requirements:

- **Internal Boundary Condition:** An internal boundary condition named `Emitter1` must be present within the 2D flow area.
- **Flow Hydrograph:** The flow hydrograph for `Emitter1` must be configured to **Use Simulation Time**.
- **Computation Interval:** Spawned simulations inherit the computation interval from the base HEC-RAS unsteady plan.
- **Projection:** The HEC-RAS projection file must be located in the same directory as the base HEC-RAS model.
- **Terrain:** The terrain file used by the HEC-RAS model must be accessible from the base HEC-RAS model directory.

---

# Sample Data

Sample HEC-RAS input data and RAS2FIM-2D output data are provided separately from the source repository.

### Browse the sample data

**[RAS2FIM-2D Sample Data](https://rasfim-2d-sample.s3.amazonaws.com/index.html)**

The sample data are hosted in an Amazon S3 bucket:

```text
s3://rasfim-2d-sample
```

The sample dataset includes representative HEC-RAS 2D input files and RAS2FIM-2D output files that can be used to test the workflow and examine the resulting FIM products.

> **Note:** The sample data are intentionally hosted separately from the GitHub repository because HEC-RAS model inputs and FIM outputs can be relatively large.

---

# Docker

RAS2FIM-2D uses Docker to provide a reproducible processing environment. The workflow also requires a containerized installation of **HEC-RAS 6.5**.

## HEC-RAS 6.5 Docker Image

A pre-built HEC-RAS 6.5 Docker image is available on Docker Hub:
**[Docker Hub HEC-RAS Linux v6.5](https://hub.docker.com/r/civileng127/ras_v65)**

Pull the image with:

```bash
docker pull civileng127/ras_v65:v0
```
This container provides the HEC-RAS 6.5 runtime to execute the 2D unsteady simulations in a Linux environment.

## Build the RAS2FIM-2D Image

A Dockerfile is included in this repository for building the RAS2FIM-2D processing environment.

Clone the repository and build the image:

```bash
docker build -t ras2fim2d .
```

## Run RAS2FIM-2D

After building the image, mount the directory containing the HEC-RAS model and output data into the container and run the RAS2FIM-2D workflow.

Refer to the example configuration and scripts in this repository for the required command-line arguments.

# Output

RAS2FIM-2D produces flood inundation mapping libraries containing water-surface elevations associated with multiple flow rates.

The primary output is a NetCDF file containing:

- Water-surface elevation (WSEL)
- Flow rate associated with each WSEL
- Terrain data used by the HEC-RAS simulation
- Spatial information required to interpret the results

The output can be used as a library for rapidly retrieving flood inundation information for NWM NextGen stream segments.

---

# License

See the repository license for terms of use.
