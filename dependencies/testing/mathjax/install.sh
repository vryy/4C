#!/bin/bash
# Install mathjax
# Call with
# ./install.sh /path/to/install/dir

# Exit the script at the first failure
set -e

INSTALL_DIR="$1"
VERSION="2.7.9"
CHECKSUM="c9167279c87da044f2ff910ad573a02ce90354cb59440ae568eb86e1630f65df"

wget --no-verbose https://github.com/mathjax/MathJax/archive/refs/tags/${VERSION}.tar.gz
# Verify checksum
if [ $CHECKSUM = `sha256sum ${VERSION}.tar.gz | awk '{print $1}'` ]
then
  echo "Checksum matches"
else
  echo "Checksum does not match"
  exit 1
fi

tar -xzf ${VERSION}.tar.gz -C ${INSTALL_DIR}
