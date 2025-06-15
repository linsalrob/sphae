FROM --platform=linux/amd64 ubuntu:20.04

ENV DEBIAN_FRONTEND="noninteractive"

ARG LIBFABRIC_VERSION=1.21.0

# Install required packages and dependencies
RUN   apt -y update \
      && apt -y install build-essential wget doxygen gnupg gnupg2 curl apt-transport-https software-properties-common libgl1  \
 git vim gfortran libtool python3-venv ninja-build python3-pip \
      libnuma-dev python3-dev \
      && apt -y remove --purge --auto-remove cmake \
      && wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null\
 | gpg --dearmor - | tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null \
      && apt-add-repository -y "deb https://apt.kitware.com/ubuntu/ jammy-rc main" \
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

ARG SPHAE_VERSION=1.4.6

#
# Install miniforge
#
RUN set -eux ; \
  curl -LO https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh ; \
  bash ./Miniforge3-* -b -p /opt/miniforge3 -s ; \
  rm -rf ./Miniforge3-*
ENV PATH /opt/miniforge3/bin:$PATH
#
# Install conda environment
# 
RUN set -eux ; \
  mamba install -y -c conda-forge -c bioconda -c defaults \
  sphae=${SPHAE_VERSION}=pyhdfd78af_0 python==3.11 
ENV PATH /opt/miniforge3/bin:$PATH
RUN conda clean -af -y

# install databases - needed or else sphae won't run - makes the container huge but nothing we can really do
RUN sphae install --threads 8 --conda-frontend mamba 

# for weird filtlong bug https://github.com/potree/PotreeConverter/issues/281
ENV LC_ALL=C
ENV LANGUAGE=

# download test data
RUN git clone "https://github.com/linsalrob/sphae.git"

# run with --conda-create-envs-only to create all required conda envs (without running itself for speedup)
RUN sphae run --threads 8 --input sphae/tests/data/illumina-subset --output example -k --conda-frontend mamba --conda-create-envs-only
RUN sphae run --threads 8 --input sphae/tests/data/nanopore-subset --sequencing longread --output examplelr -k --conda-frontend mamba --conda-create-envs-only
RUN sphae annotate --threads 8 --genome sphae/tests/data/genome --output exampleg --conda-frontend mamba --conda-create-envs-only
# cleanup
RUN rm -rf example
RUN rm -rf examplelr
RUN rm -rf exampleg
