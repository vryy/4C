FROM ubuntu:focal
LABEL project=baci

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
    INSTALL_DIR="/opt/baci-dependencies"
RUN mkdir -p ${INSTALL_DIR}

COPY dependencies /dependencies
# Install cmake
RUN /dependencies/cmake/install.sh /usr/local

# Install SuiteSparse 5.4.0
RUN /dependencies/suitesparse/install.sh ${INSTALL_DIR}

# Install superLU_dist 2.5
RUN /dependencies/superlu_dist/install.sh ${INSTALL_DIR}

# Install qhull 2012.1
RUN /dependencies/qhull/install.sh ${INSTALL_DIR}

 # Install CLN 1.3.4
RUN /dependencies/cln/install.sh ${INSTALL_DIR}

# Install Trilinos 2022_Q1
RUN /dependencies/trilinos/install.sh ${INSTALL_DIR}

# Packages for testing
# Installation directory for dependencies concerning testing
ENV BACI_TESTING_DEPENDENCIES_DIR="/opt/baci-dependencies-testing/"
RUN mkdir ${BACI_TESTING_DEPENDENCIES_DIR}

# Install clang-tidy
RUN /dependencies/testing/clang-tidy/install.sh ${BACI_TESTING_DEPENDENCIES_DIR}

# Install Mathjax
RUN /dependencies/testing/mathjax/install.sh ${BACI_TESTING_DEPENDENCIES_DIR}

# Install Paraview
RUN /dependencies/testing/paraview/install.sh ${BACI_TESTING_DEPENDENCIES_DIR}

# add and enable the default user
ARG USER=baci
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
