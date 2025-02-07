# Use Ubuntu 20.04 as the base image with amd64 support
FROM --platform=linux/amd64 ubuntu:20.04

# Set environment variables
ENV DEBIAN_FRONTEND="noninteractive"

# Define build arguments
ARG LIBFABRIC_VERSION=1.18.1
ARG SPHAE_VERSION=1.4.5
ARG THREADS=8

# Install required packages and dependencies
RUN   apt -y update \
      && apt -y install build-essential wget doxygen gnupg gnupg2 curl apt-transport-https software-properties-common libgl1  \
      git vim gfortran libtool python3-venv ninja-build python3-pip \
      libnuma-dev python3-dev \
      && apt -y remove --purge --auto-remove cmake \
      && wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null \
      && apt-add-repository -y "deb https://apt.kitware.com/ubuntu/ focal-rc main" \
      && apt -y update

# Build and install libfabric
RUN (if [ -e /tmp/build ]; then rm -rf /tmp/build; fi;) \
      && mkdir -p /tmp/build \
      && cd /tmp/build \
      && wget https://github.com/ofiwg/libfabric/archive/refs/tags/v${LIBFABRIC_VERSION}.tar.gz \
      && tar xf v${LIBFABRIC_VERSION}.tar.gz \
      && cd libfabric-${LIBFABRIC_VERSION} \
      && ./autogen.sh \
      && ./configure \
      && make -j 16 \
      && make install

# Install Miniforge
RUN set -eux ; \
  curl -LO https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh ; \
  bash ./Miniforge3-* -b -p /opt/miniforge3 -s ; \
  rm -rf ./Miniforge3-*
ENV PATH /opt/miniforge3/bin:$PATH

# Install conda environment
RUN set -eux ; \
  mamba install -y -c conda-forge -c bioconda -c defaults sphae=${SPHAE_VERSION}=pyhdfd78af_0 python==3.11 
ENV PATH /opt/miniforge3/bin:$PATH
RUN conda clean -af -y

# Download test data
RUN git clone "https://github.com/linsalrob/sphae.git"

# Install Sphae databases (with dynamic threads)
RUN sphae install --threads ${THREADS} --conda-frontend mamba

# Environment settings for filtlong bug
ENV LC_ALL=C
ENV LANGUAGE=

# Create the directory if it doesn't exist, then list its contents
RUN mkdir -p sphae/tests/db && ls sphae/tests/db

# Create required conda environments without running
RUN sphae run --threads ${THREADS} --input sphae/tests/data/illumina-subset -k --use-conda --db_dir sphae/tests/db --conda-create-envs-only --use-singularity --sdm apptainer --use-conda
RUN sphae run --threads ${THREADS} --input sphae/tests/data/nanopore-subset --sequencing longread -k --conda-create-envs-only --db_dir sphae/tests/db --use-singularity --sdm apptainer --use-conda
RUN sphae annotate --threads ${THREADS} --input sphae/tests/data/genome -k --conda-create-envs-only --db_dir sphae/tests/db --use-singularity --sdm apptainer --use-conda

# Cleanup
RUN rm -rf sphae.out /tmp/* /var/tmp/* /var/lib/apt/lists/*

# Set default working directory
WORKDIR /sphae

# Set container entrypoint
CMD ["/bin/bash"]