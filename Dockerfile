FROM --platform=linux/amd64 ubuntu:20.04

ENV DEBIAN_FRONTEND="noninteractive"

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
  mamba install -y -c conda-forge -c bioconda -c defaults \
  sphae=${SPHAE_VERSION}=pyhdfd78af_0 python==3.11 
ENV PATH /opt/miniforge3/bin:$PATH
RUN conda clean -af -y

# Install Sphae databases (with dynamic threads)
RUN sphae install --threads ${THREADS} --conda-frontend mamba

# Environment settings for filtlong bug
ENV LC_ALL=C
ENV LANGUAGE=

# Download test data
RUN git clone "https://github.com/linsalrob/sphae.git"

#remove one of the test datasets
RUN rm -rf sphae/tests/data/illumina-subset/SRR16219309*

# Create required conda environments without running
RUN sphae run --threads ${THREADS} --input sphae/tests/data/illumina-subset --output example -k --conda-frontend mamba --conda-create-envs-only
RUN sphae run --threads ${THREADS} --input sphae/tests/data/nanopore-subset --sequencing longread --output examplelr -k --conda-frontend mamba --conda-create-envs-only

# Cleanup
RUN rm -rf example examplelr /tmp/* /var/tmp/* /var/lib/apt/lists/*
