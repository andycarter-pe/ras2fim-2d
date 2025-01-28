# Setting up ras2fim-2d code base to run in a miniconda environment

# Use the Miniconda base image from Continuum
FROM continuumio/miniconda3:latest

# Set the Python version you want (3.8.12 in this case)
RUN conda install python=3.8.12 -y

# Set environment variables to ensure conda is available in PATH
ENV PATH /opt/conda/bin:$PATH

# Install base packages using both pip and conda with -y to automatically confirm
RUN conda install -y gdal && \
    pip install geopandas rioxarray h5py networkx scipy matplotlib
	
# Install Git to clone repository and nano (text editor)
RUN apt-get update && apt-get install -y git nano s3fs

# Clone the desired GitHub repository
RUN git clone https://github.com/andycarter-pe/ras2fim-2d.git /ras2fim-2d

# Set working directory
WORKDIR /ras2fim-2d/src

# Clean up any unnecessary files to reduce image size
RUN conda clean -a

# Expose a port (optional for running a service later)
EXPOSE 8888

# Default command to run when container starts
CMD [ "bash" ]