#!/bin/bash
# Install cln
# Call with
# ./install.sh /path/to/install/dir

# Exit the script at the first failure
set -e

INSTALL_DIR="$1"
# Number of procs for building (default 4)
NPROCS=${NPROCS=4}
VERSION="1.3.4"
CHECKSUM="2d99d7c433fb60db1e28299298a98354339bdc120d31bb9a862cafc5210ab748"

wget --no-verbose https://ginac.de/CLN/cln-${VERSION}.tar.bz2
# Verify checksum
if [ $CHECKSUM = `sha256sum cln-${VERSION}.tar.bz2 | awk '{print $1}'` ]
then
  echo "Checksum matches"
else
  echo "Checksum does not match"
  exit 1
fi

tar -xjf cln-${VERSION}.tar.bz2
cd cln-${VERSION}
./configure --prefix=${INSTALL_DIR}
make -j${NPROCS} && make install
cd ../ && rm -rf cln*
