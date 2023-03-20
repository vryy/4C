#!/bin/bash
# Install trilinos
# Call with
# ./install.sh /path/to/install/dir

# Exit the script at the first failure
set -e

INSTALL_DIR="$1"
# Number of procs for building (default 4)
NPROCS=${NPROCS=4}
#VERSION=""
#CHECKSUM=""


# Location of script to apply patches later
SCRIPT_DIR="`dirname "$0"`"
CMAKE_COMMAND=${INSTALL_DIR}/bin/cmake

git clone https://github.com/trilinos/Trilinos.git
cd Trilinos
git checkout cd1f2de822758156e2fc24a6cefeeb01dc83436d
git apply ${SCRIPT_DIR}/trilinos_for_baci.patch && git apply ${SCRIPT_DIR}/trilinos_for_baci_output.patch
cd .. && mkdir trilinos_build && cd trilinos_build

MPI_DIR=/usr
MPI_BIN_DIR=$MPI_DIR/bin

$CMAKE_COMMAND \
  -D CMAKE_BUILD_TYPE:STRING="RELEASE" \
  -D CMAKE_CXX_COMPILER:FILEPATH="$MPI_BIN_DIR/mpic++" \
  -D CMAKE_C_COMPILER:FILEPATH="$MPI_BIN_DIR/mpicc" \
  -D CMAKE_Fortran_COMPILER:FILEPATH="$MPI_BIN_DIR/mpif90" \
  -D CMAKE_INSTALL_PREFIX:STRING=$INSTALL_DIR \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
  -D CMAKE_COLOR_MAKEFILE:BOOL=ON \
  -D BUILD_SHARED_LIBS:BOOL=OFF \
  -D Gtest_SKIP_INSTALL:BOOL=TRUE \
  \
  -D Trilinos_ASSERT_MISSING_PACKAGES=OFF \
  -D Trilinos_VERBOSE_CONFIGURE:BOOL=OFF \
  -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
  -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
  -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
  -D Trilinos_ENABLE_TESTS:BOOL=OFF \
  -D Trilinos_ENABLE_EXAMPLES:BOOL=OFF \
  \
  -D Trilinos_ENABLE_Amesos:BOOL=ON \
  -D Trilinos_ENABLE_Amesos2:BOOL=OFF \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_AztecOO:BOOL=ON \
  -D Trilinos_ENABLE_Belos:BOOL=ON \
  -D Trilinos_ENABLE_Epetra:BOOL=ON \
    -D Epetra_HIDE_DEPRECATED_CODE:BOOL=ON \
  -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
    -D EpetraExt_BUILD_GRAPH_REORDERINGS:BOOL=ON \
    -D EpetraExt_ENABLE_HDF5:BOOL=OFF \
  -D Trilinos_ENABLE_Intrepid:BOOL=ON \
  -D Trilinos_ENABLE_Ifpack:BOOL=ON \
  -D Trilinos_ENABLE_Isorropia:BOOL=ON \
  -D Trilinos_ENABLE_Kokkos:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D Trilinos_ENABLE_Galeri:BOOL=OFF \
  -D Trilinos_ENABLE_MueLu:BOOL=ON \
    -D MueLu_ENABLE_EXAMPLES:BOOL=OFF \
    -D MueLu_ENABLE_TESTS:BOOL=OFF \
  -D Trilinos_ENABLE_NOX:BOOL=ON \
  -D Trilinos_ENABLE_Sacado:BOOL=ON \
  -D Trilinos_ENABLE_SEACASExodus:BOOL=ON \
  -D Trilinos_ENABLE_SEACASNemesis:BOOL=OFF \
  -D Trilinos_ENABLE_Shards:BOOL=ON \
  -D Trilinos_ENABLE_Teuchos:BOOL=ON \
    -D Teuchos_GLOBALLY_REDUCE_UNITTEST_RESULTS:BOOL=ON \
    -D Teuchos_HIDE_DEPRECATED_CODE:BOOL=ON \
  -D Trilinos_ENABLE_Tpetra:BOOL=ON \
    -D Tpetra_ENABLE_DEPRECATED_CODE:BOOL=OFF \
    -D Tpetra_INST_INT_INT:BOOL=ON \
  -D Trilinos_ENABLE_Xpetra:BOOL=ON \
    -D Xpetra_ENABLE_DEPRECATED_CODE:BOOL=OFF \
    -D Xpetra_ENABLE_Epetra:BOOL=ON \
    -D Xpetra_ENABLE_EpetraExt:BOOL=ON \
    -D Xpetra_ENABLE_TESTS:BOOL=OFF \
  -D Trilinos_ENABLE_Zoltan:BOOL=ON \
  -D Trilinos_ENABLE_Zoltan2:BOOL=OFF \
  \
  -D TPL_ENABLE_Boost:BOOL=ON \
  -D TPL_ENABLE_MPI:BOOL=ON \
  -D TPL_ENABLE_ParMETIS:BOOL=ON \
    -D ParMETIS_INCLUDE_DIRS:PATH="/usr/include" \
    -D ParMETIS_LIBRARY_DIRS:PATH="/usr/lib/libparmetis.so.4.0.3" \
  -D TPL_ENABLE_UMFPACK:BOOL=ON \
    -D UMFPACK_INCLUDE_DIRS:FILEPATH="$INSTALL_DIR/include" \
    -D UMFPACK_LIBRARY_DIRS:FILEPATH="$INSTALL_DIR/lib" \
  -D TPL_ENABLE_Netcdf:BOOL=ON \
  -D TPL_ENABLE_SuperLUDist:BOOL=ON \
    -D SuperLUDist_INCLUDE_DIRS:PATH="$INSTALL_DIR/include" \
    -D SuperLUDist_LIBRARY_DIRS:PATH="$INSTALL_DIR/lib" \
    -D SuperLUDist_LIBRARY_NAMES:STRING="superlu_dist_2.5" \
  -D TPL_ENABLE_DLlib:BOOL=OFF \
  ../Trilinos

make -j${NPROCS} install
cd ..
