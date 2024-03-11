#!/bin/bash
# Install clang-tidy
# Call with
# ./install.sh /path/to/install/dir

# Exit the script at the first failure
set -e

# TODO: specify install dir
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
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DLLVM_ENABLE_PROJECTS="clang;clang-tools-extra" ../llvm
make -j${NPROCS} install clang-tidy
cd ../../ && rm -rf llvm*
