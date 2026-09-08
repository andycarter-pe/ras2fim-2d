# RAS2FIM-2D <img src="doc/Logo_CWE.png" align="right" alt="ras2fim2d agency" height="80"> <br> <br>
## <i>RAS2FIM-2D - Flood Inundation Mapping (FIM) using HEC-RAS 2D</i>

<img src="/doc/ras2fim2d-logo-20260907.png" align="right"
     alt="ras2fim2d logo" width="160" height="160">

**Description**:  RAS2FIM-2D convert 2D HEC-RAS models into a set of flood inundation mapping (fim) library of rasters (netCDF) for a corresponding National Water Model (NWM) NextGEN stream segments within the HEC-RAS model's 2D area.<br><br>
Simulations are run to a points of stability using a "firehose" method where a series of constant flows is set on an internal boundary condition ("Emitter1") <br><br>
Flows are emmitted at the upsteram ends of the NextGEN hydrofabric catchments and clipped to flooding around a given NextGEN id (Example: wb-2427467).  The ultimate output is a netCDF of 
multiple water surface elevations (WSEL) indexed to multiple flow rates in cfs.  The netCDF file also contains the source terrain on which the WSEL were determined.<br><br>
These scripts were developed in support of the National Weather Service (Research Project NA22NWS4320003 / A25-0366-S018).

<p align="center">
  <img src="/doc/ras2fim_animation.gif" alt="sample cross section" width="55%">
</p>

  - **Technology stack**: Scripts were all developed in Python 3.8.12<br>
  - **Status**:  Version 0.1- Preliminary release. (2026.09.07)<br>
  - **Related Projects**: Flood Inundation Maps from HEC-RAS 1D models  https://github.com/NOAA-OWP/ras2fim<br>
  
## HEC-RAS 2D Requirements
  - **Internal BC**: An internal boundary condition named "Emitter1" must be present in the 2D area<br>
  - **Flow Hydrograph**:  Flow hydrogrph for "Emitter1" must be set to 'Use Simulation Time'<br>
  - **Computation Interval**:  Spawned runs inherit computation interval from base HEC-RAS unsteady plan<br>
  - **Projection File**:  The projection file should be a folder where the base HEC-RAS exists<br>
  - **Terrain File**:  The terrain file should be a folder where the base HEC-RAS exists<br>
  
Note that a sample HEC-RAS input file is provided in this repository.
 
```
docker build -t ras2fim2d .
  
## Dockerfile
To build a container from this repository, clone to your local drive and build with the following command
```
docker build -t ras2fim2d .
```

## Docker Container
For convience, a container has been pre-built and pushed to DockerHub.  To pull this container to your machine...
```
docker pull civileng127/lisflood2fim:20260611
```
Run the containers demo: Note '/mnt/e/lisflood_dump' is the local directory where output will be saved
```
docker run -it \
-v /mnt/e/lisflood_dump:/mnt \
civileng127/lisflood2fim:20260611 \
bash -c "
