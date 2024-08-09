#!/bin/bash
# Install the GNU Compiler Collection (GCC) in a standard configuration suitable for this project
# Call with
# ./install.sh /path/to/install/dir

# Exit the script at the first failure
set -e

INSTALL_DIR="$1"
# Number of procs for building (default 4)
NPROCS=${NPROCS=4}
# Release version 14.2.0
VERSION="04696df09633baf97cdbbdd6e9929b9d472161d3"

git clone git://gcc.gnu.org/git/gcc.git gcc-source
cd gcc-source
git checkout $VERSION
./contrib/download_prerequisites
cd ..
mkdir gcc-build && cd gcc-build
../gcc-source/configure --prefix=${INSTALL_DIR} --enable-languages=c,c++,fortran --disable-multilib
make -j${NPROCS}
make install
cd ../ && rm -rf gcc*


