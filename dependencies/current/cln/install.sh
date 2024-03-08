#!/bin/bash
# Install cln
# Call with
# ./install.sh /path/to/install/dir

# Exit the script at the first failure
set -e

INSTALL_DIR="$1"
# Number of procs for building (default 4)
NPROCS=${NPROCS=4}
VERSION="cln_1-3-4"
COMMIT="9b86a7fc69feb1b288469982001af565f73057eb"

git clone git://www.ginac.de/cln.git
cd cln
git checkout $VERSION

# Check commit sha1
if [ $COMMIT = `git rev-parse HEAD` ]
then
  echo "The tag $VERSION matches the commit sha1"
else
  echo "The tag $VERSION does not match the commit sha1"
  exit 1
fi

autoreconf -iv
./configure --prefix=${INSTALL_DIR}
make -j${NPROCS} && make install
cd ../ && rm -rf cln*
