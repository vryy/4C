#!/bin/bash
# Install qhull
# Call with
# ./install.sh /path/to/install/dir

# Exit the script at the first failure
set -e

INSTALL_DIR="$1"
# Number of procs for building (default 4)
NPROCS=${NPROCS=4}
VERSION="2012.1"
CHECKSUM="a35ecaa610550b7f05c3ce373d89c30cf74b059a69880f03080c556daebcff88"

# Install qhull 2012.1
wget --no-verbose http://www.qhull.org/download/qhull-${VERSION}-src.tgz
# Verify checksum
if [ $CHECKSUM = `sha256sum qhull-${VERSION}-src.tgz | awk '{print $1}'` ]
then
  echo "Checksum matches"
else
  echo "Checksum does not match"
  exit 1
fi

tar -xzf qhull-${VERSION}-src.tgz
cd qhull-${VERSION}/build
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} ..
make -j${NPROCS} && make install
cd ../../ && rm -rf qhull*
