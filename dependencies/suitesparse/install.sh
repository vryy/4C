#!/bin/bash
# Install SuiteSparse
# Call with
# ./install.sh /path/to/install/dir

# Exit the script at the first failure
set -e

INSTALL_DIR="$1"
# Number of procs for building (default 4)
NPROCS=${NPROCS=4}
VERSION="5.4.0"
CHECKSUM="374dd136696c653e34ef3212dc8ab5b61d9a67a6791d5ec4841efb838e94dbd1"

wget --no-verbose https://people.engr.tamu.edu/davis/SuiteSparse/SuiteSparse-${VERSION}.tar.gz
# Verify checksum
if [ $CHECKSUM = `sha256sum SuiteSparse-${VERSION}.tar.gz | awk '{print $1}'` ]
then
  echo "Checksum matches"
else
  echo "Checksum does not match"
  exit 1
fi

tar -xzf SuiteSparse-${VERSION}.tar.gz
cd SuiteSparse/
make -j${NPROCS} library BLAS=-lblas
make install INSTALL=${INSTALL_DIR} BLAS=-lblas
cd ../
rm -rf SuiteSparse*
# Need to specify metis location?
