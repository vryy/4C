FROM ubuntu:focal
LABEL project=4C

# Prevents tzdata asking for user feedback
ENV DEBIAN_FRONTEND=noninteractive

USER root

# Set locale information:. region and timezone
RUN apt-get update && apt-get install -y locales && \
    localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8
RUN locale-gen en_US.UTF-8  
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8

RUN apt-get update && apt-get install -y \
      build-essential \
      ccache \
      ffmpeg \
      git \
      libglu1-mesa \ 
      python3 \
      sudo \
      unzip \
      vim \
      wget

RUN apt-get update && apt-get install -y \
      cxxtest \
      doxygen \
      graphviz \
      texinfo \
      lcov \
      libblas-dev \
      libboost-all-dev \
      libhdf5-dev \
      libhdf5-openmpi-dev \
      libnetcdf-dev \
      libfftw3-dev \
      lld \
      python3 \
      python3.8-venv \
      python-is-python3 \
      liblapack-dev \
      libopenmpi-dev \
      libparmetis-dev \
      metis \
      ninja-build \
      libyaml-dev

# Create directory for dependencies
ARG NPROCS=12
ENV NPROCS=$NPROCS \
    INSTALL_DIR="/opt/4C-dependencies"
# Trilinos version: if unset the version defined in the install script is used
ARG TRILINOS_VERSION
RUN mkdir -p ${INSTALL_DIR}

COPY dependencies /dependencies
# Install cmake
RUN /dependencies/current/cmake/install.sh /usr/local

# Install SuiteSparse
RUN /dependencies/current/suitesparse/install.sh ${INSTALL_DIR}

# Install superLU_dist
RUN /dependencies/current/superlu_dist/install.sh ${INSTALL_DIR}

# Install qhull
RUN /dependencies/current/qhull/install.sh ${INSTALL_DIR}

 # Install CLN
RUN /dependencies/current/cln/install.sh ${INSTALL_DIR}

# Install Trilinos
RUN /dependencies/trilinos_develop/trilinos/install.sh ${INSTALL_DIR} ${TRILINOS_VERSION}

# Install (optional) backtrace library
RUN /dependencies/current/backtrace/install.sh ${INSTALL_DIR}

# Packages for testing
# Installation directory for dependencies concerning testing
ENV FOUR_C_TESTING_DEPENDENCIES_DIR="/opt/4C-dependencies-testing/"
RUN mkdir ${FOUR_C_TESTING_DEPENDENCIES_DIR}

# Install LLVM (use standard install directory)
RUN /dependencies/testing/llvm/install.sh /usr/local

# Install Mathjax
RUN /dependencies/testing/mathjax/install.sh ${FOUR_C_TESTING_DEPENDENCIES_DIR}

# Install Paraview
RUN /dependencies/testing/paraview/install.sh ${FOUR_C_TESTING_DEPENDENCIES_DIR}

# add and enable the default user
ARG USER=user
RUN adduser --disabled-password --shell '/usr/bin/bash' --gecos '' $USER
RUN adduser $USER sudo; echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

#make sure everything is in place
RUN chown -R $USER:$USER /home/$USER
USER $USER
ENV HOME=/home/$USER
ENV USER=$USER
# see https://github.com/open-mpi/ompi/issues/4948
ENV OMPI_MCA_btl_vader_single_copy_mechanism=none
WORKDIR $HOME
