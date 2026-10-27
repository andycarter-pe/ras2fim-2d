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
# Revised: 2026-10-27 -- Python 3.14 env; optional healed terrain / mapping engine (CLB Engineering)
# License: BSD 3-Clause License
# ------------------------------------------------------------------------
FROM continuumio/miniconda3:latest

# --- System packages (single layer, clean up in same layer) ---
RUN apt-get -o Acquire::Check-Valid-Until=false update && \
    apt-get install -y --no-install-recommends git nano wget proj-bin && \
    rm -rf /var/lib/apt/lists/*

# --- Use the fast libmamba solver and conda-forge (strict) ---
# (libmamba is default on recent miniconda3, but set explicitly to be safe)
RUN conda install -n base --override-channels -c conda-forge conda-libmamba-solver -y && \
    conda config --set solver libmamba && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

# --- ONE solve for the entire scientific/geospatial stack ---
# All C-extension geo libs come from conda-forge so PROJ/GEOS/GDAL
# stay consistent. No pip wheels for these = no PROJ conflicts.
# Keep the base image Python/conda tools separate from the application runtime.
RUN conda create -p /opt/ras2fim-env --override-channels -c conda-forge -y \
        python=3.14 \
        gdal \
        rasterio \
        "fiona>=1.9" \
        "geopandas>=0.14,<2" \
        "rioxarray>=0.15,<1" \
        "h5py>=3.8" \
        "scipy>=1.11" \
        networkx=2.8.8 \
        matplotlib \
        "numpy<3" \
        "pandas<3" \
        pip \
        netcdf4 \
        h5netcdf && \
    conda clean -a -y

ENV PATH="/opt/ras2fim-env/bin:/opt/conda/bin:$PATH"
# The env is used by PATH, not activated, so point PROJ and GDAL at its data
ENV PROJ_DATA="/opt/ras2fim-env/share/proj" GDAL_DATA="/opt/ras2fim-env/share/gdal"

# Rust raster backend for the optional healed mapping stage.
RUN python -m pip install --no-deps terrio==0.5.0

# --- Application code, directories, and global inputs ---
COPY . /ras2fim-2d
RUN python -m pip install -r /ras2fim-2d/requirements-mapping-engine.txt "numpy<3" "pandas<3"
RUN mkdir -p /ras2fim-2d && \
    mkdir /global_input /model_input /model_output && \
    wget https://ras2fim-2d-global-inputs.s3.us-east-1.amazonaws.com/LimitingDischarge.csv \
        -O /global_input/LimitingDischarge.csv && \
    wget https://ras2fim-2d-global-inputs.s3.us-east-1.amazonaws.com/nextgen_12.gpkg \
        -O /global_input/nextgen_fabric.gpkg

WORKDIR /ras2fim-2d/src
CMD [ "bash" ]