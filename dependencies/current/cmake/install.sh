#!/bin/bash
# Install cmake
# Call with
# ./install.sh /path/to/install/dir

# Exit the script at the first failure
set -e

INSTALL_DIR="$1"
# Number of procs for building (default 4)
NPROCS=${NPROCS=4}
VERSION="3.25.1"
CHECKSUM="6598da34f0e3a0f763809e25cfdd646aa1d5e4d133c4277821e63ae5cfe09457"

wget --no-verbose https://github.com/Kitware/CMake/releases/download/v${VERSION}/cmake-${VERSION}-linux-x86_64.sh
# Verify checksum
if [ $CHECKSUM = `sha256sum cmake-${VERSION}-linux-x86_64.sh | awk '{print $1}'` ]
then
  echo "Checksum matches"
else
  echo "Checksum does not match"
  exit 1
fi

chmod +x cmake-${VERSION}-linux-x86_64.sh
./cmake-${VERSION}-linux-x86_64.sh --prefix=${INSTALL_DIR} --skip-license
rm cmake-${VERSION}-linux-x86_64.sh
