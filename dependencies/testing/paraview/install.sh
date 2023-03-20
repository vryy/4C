#!/bin/bash
# Install mathjax
# Call with
# ./install.sh /path/to/install/dir

# Exit the script at the first failure
set -e

INSTALL_DIR="$1"
VERSION="5.5.2" # Not used as specific download path changes from version to version
CHECKSUM="b2ba52093bd04331217f32805d810198bccfa582de78d46185e95d0f39857772"

wget --no-verbose -O paraview.tar.gz 'https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.5&type=binary&os=Linux&downloadFile=ParaView-5.5.2-Qt5-MPI-Linux-64bit.tar.gz'
# Verify checksum
if [ $CHECKSUM = `sha256sum paraview.tar.gz | awk '{print $1}'` ]
then
  echo "Checksum matches"
else
  echo "Checksum does not match"
  exit 1
fi

tar -xzf paraview.tar.gz -C ${INSTALL_DIR}
