# Official Ubuntu 22.04 base image
FROM ubuntu:22.04

# Set the environment variables to avoid user prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/London

# Set working directory
WORKDIR /app

# Copy codebase into image
COPY ./core /app/core
COPY ./README.md /app/README.md
COPY ./requirements.txt /app/requirements.txt

# Basic dependencies
RUN apt-get update && apt-get -y install \
    python3 \
    python3-pip \
    build-essential \
	unzip \
    curl \
    git \
    wget \
	libssl-dev \
    && apt-get clean

# Download utility functions
RUN git clone https://github.com/ellis-langford/ImageUtils.git utils

# Upgrade pip to the latest version
RUN python3 -m pip install --upgrade pip

# Install Python packages
RUN python3.10 -m pip install \
--trusted-host pypi.python.org -r /app/requirements.txt

# Install Jupyter
RUN pip install jupyterlab

# Expose JupyterLab's default port (8888)
EXPOSE 8888

# Revert working directory
WORKDIR /mnt/ellis

# Set the command to run JupyterLab when the container starts
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]
