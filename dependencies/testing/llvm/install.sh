#!/bin/bash
# Install relevant parts of the LLVM project (clang compiler and supporting tools)
# Call with
# ./install.sh /path/to/install/dir

# Exit the script at the first failure
set -e

INSTALL_DIR="$1"
# Number of procs for building (default 4)
NPROCS=${NPROCS=4}
VERSION="18.1.1"
CHECKSUM="62439f733311869dbbaf704ce2e02141d2a07092d952fc87ef52d1d636a9b1e4"

wget --no-verbose https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-${VERSION}.tar.gz
# Verify checksum
if [ $CHECKSUM = `sha256sum llvmorg-${VERSION}.tar.gz | awk '{print $1}'` ]
then
  echo "Checksum matches"
else
  echo "Checksum does not match"
  exit 1
fi

tar -xzf llvmorg-${VERSION}.tar.gz
cd llvm-project-llvmorg-${VERSION}
mkdir build && cd build
cmake -G "Unix Makefiles" \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DLLVM_ENABLE_PROJECTS="clang;clang-tools-extra" \
    -DLLVM_ENABLE_RUNTIMES="compiler-rt;openmp" \
    ../llvm
make -j${NPROCS} install
cd ../../ && rm -rf llvm*
