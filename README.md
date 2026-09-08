# RAS2FIM-2D <img src="doc/Logo_CWE.png" align="right" alt="ras2fim2d agency" height="80"> <br> <br>
## <i>RAS2FIM-2D - Flood Inundation Mapping (FIM) using HEC-RAS 2D</i>

<img src="/doc/ras2fim2d-logo-20260907.png" align="right"
     alt="ras2fim2d logo" width="160" height="160">

**Description**:  LISFLOOD2FIM creates Flood Inundation Maps (FIMs) indexed to various excess rainfall rates for a given catchment (watershed). The initial routine determines the stream network on which to apply the excess precipitation. This requires hydro-enforcement through dams and roadways. Input data needed for LISFLOOD-FP is then created, including parameters, terrain, boundary conditions, and stream-centric rainfall..<br><br>
Simulations are run using LISFLOOD-FP v8.1.0 to a point of “stability,” where inflow rainfall equals the outflow rate. This process is repeated for various excess rainfall intensities. The resulting “stable” flood depths for each intensity are aggregated into a single data cube (NetCDF), representing flood depth across multiple intensities over the watershed. <br><br>
These scripts were developed in support of the National Weather Service (Research Project NA22NWS4320003 / A25-0366-S018).

<p align="center">
  <img src="/doc/lisflood2fim.gif" alt="sample cross section" width="55%">
</p>

  - **Technology stack**: Scripts were all developed in Python 3.8.12<br>
  - **Status**:  Version 0.1- Preliminary release. (2026.09.07)<br>
  - **Related Projects**: Flood Inundation Maps from HEC-RAS 1D models  https://github.com/NOAA-OWP/ras2fim<br>
  
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