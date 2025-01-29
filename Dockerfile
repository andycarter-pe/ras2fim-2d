# ------------------------------------------------------------------------
# Dockerfile for running the RAS2FIM-2D workflow
# Dockhub Image name: civileng127/ras2fim-2d:v0
# Description: Create flood inundation products for NextGEN
# hydrofabric using HEC-RAS 2D floodplain models
#
# Base Image: continuumio/miniconda3:latest
#  -- Uses the Debian 11; 'bullseye'
# Version: 1.0.0
# Created by: Andy Carter, PE 
# -- Center for Water and the Environment
# -- University of Texas at Austin
# Date: 2025-01-29
# License: BSD 3-Clause License
# ------------------------------------------------------------------------

# Example use: docker run -it -v E:\ras-docker-20240908:/ras/Linux_RAS_v65/mac-test civileng127/ras_v65:v0 /bin/bash -c "cd /ras/Linux_RAS_v65/mac-test && RasUnsteady sample_ras_name.p01,hdf.tmp x01"

# Use the Miniconda base image from Continuum
FROM continuumio/miniconda3:latest

# Install apt-get packages first (if any)
RUN apt-get update && apt-get install -y git nano wget

# Set the Python version you want (3.8.12 in this case)
RUN conda install python=3.8.12 -y

# Install the latest gdal via conda
RUN conda install gdal -y

# Install python libraries via pip
RUN pip install geopandas==0.12.1 rioxarray==0.13.1 h5py==3.7.0 networkx==2.8.8 scipy==1.9.3 matplotlib

# Clean up conda cache to reduce image size
RUN conda clean -a

# Remove unnecessary apt-get lists to reduce image size
RUN rm -rf /var/lib/apt/lists/*

# Set environment variables to ensure conda is available in PATH
ENV PATH /opt/conda/bin:$PATH

# Clone the desired GitHub repository
RUN git clone https://github.com/andycarter-pe/ras2fim-2d.git /ras2fim-2d

# Create necessary directories
RUN mkdir /global_input /model_input /model_output

# Download files from a public S3 bucket
RUN wget https://ras2fim-2d-global-inputs.s3.us-east-1.amazonaws.com/LimitingDischarge.csv -O /global_input/LimitingDischarge.csv && \
    wget https://ras2fim-2d-global-inputs.s3.us-east-1.amazonaws.com/nextgen_12.gpkg -O /global_input/nextgen_fabric.gpkg

# Set working directory
WORKDIR /ras2fim-2d/src

# Default command to run when container starts
CMD [ "bash" ]