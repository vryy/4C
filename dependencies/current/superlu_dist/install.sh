#!/bin/bash
# Install superLU_dist
# Call with
# ./install.sh /path/to/install/dir

# Exit the script at the first failure
set -e

INSTALL_DIR="$1"
# Number of procs for building (default 4)
NPROCS=${NPROCS=4}
VERSION="2.5"
CHECKSUM="ddba456a9aea3d589ade0b90ca197402ce75c665d344f4d8b735e80f5a821d28"

# Install superLU_dist 2.5
# Alternative mirror from web.archive https://web.archive.org/web/20221213023414/https://portal.nersc.gov/project/sparse/superlu/superlu_dist_2.5.tar.gz \
# Manual copy of header files to install directory because Makefile does not do it
# Additionally, Martin fixed a bug. So, we add the bug-fix to the super-hacky setup as well.
# > September 29, 2015:
# > Fixed bug in SuperLU_DIST_2.5/SRC/superlu_grid.c line 99; grid->comm set to MPI_COMM_WORLD instead of Bcomm because one later frees grid->comm if the communicator was cloned
wget --no-verbose https://portal.nersc.gov/project/sparse/superlu/superlu_dist_${VERSION}.tar.gz
# Verify checksum
if [ $CHECKSUM = `sha256sum superlu_dist_${VERSION}.tar.gz | awk '{print $1}'` ]
then
  echo "Checksum matches"
else
  echo "Checksum does not match"
  exit 1
fi

tar --no-same-owner -xzf superlu_dist_${VERSION}.tar.gz
cd SuperLU_DIST_${VERSION}/
sed -i -e '99s/Bcomm/MPI_COMM_WORLD/' SRC/superlu_grid.c
cp MAKE_INC/make.i386_linux make.inc
sed -i -e "s,\(^DSUPERLULIB\s*=\).*,\1 ${INSTALL_DIR}/lib/libsuperlu_dist_${VERSION}.a," \
     -e "s,\(^BLASLIB\s*=\).*,\1 -L/usr/lib -lblas," \
     -e "s,\(^METISLIB\s*=\).*,\1 -L/usr/lib -lmetis," \
     -e "s,\(^PARMETISLIB\s*=\).*,\1 -L/usr/lib -lparmetis," \
     -e "s,\(^CFLAGS.*\),\1 -fPIC," \
     -e "s,\(^NOOPTS.*\),\1 -O0 -fPIC," \
     -e "s,\(^F90FLAGS.*\),\1 -fPIC," \
     -e "s,\(^LOADER\s*=\).*,\1 mpicc," make.inc
make && make install
cp SRC/*h ${INSTALL_DIR}/include
cd ../
