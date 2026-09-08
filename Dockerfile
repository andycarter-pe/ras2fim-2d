# ------------------------------------------------------------------------
# Dockerfile for running the RAS2FIM-2D workflow
# Dockhub Image name: civileng127/ras2fim-2d:v1
# Description: Create flood inundation products for NextGEN
# hydrofabric using HEC-RAS 2D floodplain models
#
# Base Image: continuumio/miniconda3:latest
#  -- Uses the Debian 11; 'bullseye'
# Version: 1.0.0
# Created by: Andy Carter, PE
# -- Center for Water and the Environment
# -- University of Texas at Austin
# Date: 2026-09-07
# License: BSD 3-Clause License
# ------------------------------------------------------------------------
FROM continuumio/miniconda3:latest

# --- System packages (single layer, clean up in same layer) ---
RUN apt-get -o Acquire::Check-Valid-Until=false update && \
    apt-get install -y --no-install-recommends git nano wget proj-bin && \
    rm -rf /var/lib/apt/lists/*

# --- Use the fast libmamba solver and conda-forge (strict) ---
# (libmamba is default on recent miniconda3, but set explicitly to be safe)
RUN conda install -n base -c conda-forge conda-libmamba-solver -y && \
    conda config --set solver libmamba && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

# --- ONE solve for the entire scientific/geospatial stack ---
# All C-extension geo libs come from conda-forge so PROJ/GEOS/GDAL
# stay consistent. No pip wheels for these = no PROJ conflicts.
RUN conda install -y \
        python=3.8.12 \
        gdal \
        rasterio \
        fiona=1.8.22 \
        geopandas=0.12.1 \
        rioxarray=0.13.1 \
        h5py=3.7.0 \
        scipy=1.9.3 \
        networkx=2.8.8 \
        matplotlib \
        netcdf4 \
        h5netcdf && \
    conda clean -a -y

ENV PATH /opt/conda/bin:$PATH

# --- Application code, directories, and global inputs ---
RUN git clone https://github.com/andycarter-pe/ras2fim-2d.git /ras2fim-2d && \
    mkdir /global_input /model_input /model_output && \
    wget https://ras2fim-2d-global-inputs.s3.us-east-1.amazonaws.com/LimitingDischarge.csv \
        -O /global_input/LimitingDischarge.csv && \
    wget https://ras2fim-2d-global-inputs.s3.us-east-1.amazonaws.com/nextgen_12.gpkg \
        -O /global_input/nextgen_fabric.gpkg

WORKDIR /ras2fim-2d/src
CMD [ "bash" ]