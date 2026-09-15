@echo off
rem This file is part of 4C multiphysics licensed under the
rem GNU Lesser General Public License v3.0 or later.
rem
rem See the LICENSE.md file in the top-level for license information.
rem
rem SPDX-License-Identifier: LGPL-3.0-or-later

setlocal enabledelayedexpansion

set "INSTALL_DIR=%~1"

rem Number of procs for building (default 4 if not already set)
if "%NPROCS%"=="" set "NPROCS=4"
rem git sha from Trilinos repository:
set "VERSION=13519bfdf5504923a2b0160af7c83d837b5944a7"
set "USER_DIR=%USERPROFILE:\=\\%"
set "LIB_DIR=%USER_DIR%\\opt"
echo %LIB_DIR%

rem clone the repository
cd /d %USER_DIR%
git clone https://github.com/trilinos/Trilinos.git
cd Trilinos
git checkout %VERSION%
git apply %~dp0fixes.patch

rem compiling

cd ..

set "TRILINOS_HOME=%USER_DIR%\Trilinos"

mkdir "trilinos-build"

cd "trilinos-build"

cmake ^
  -G "Ninja" ^
  -D CMAKE_BUILD_TYPE="Release" ^
  -D CMAKE_CXX_STANDARD:STRING="17" ^
  -D CMAKE_C_COMPILER="cl.exe" ^
  -D CMAKE_CXX_COMPILER="cl.exe" ^
  -D CMAKE_CXX_FLAGS:STRING="%CMAKE_CXX_FLAGS% -D_WIN32 /EHsc /MP /bigobj" ^
  -D CMAKE_INSTALL_PREFIX="%INSTALL_DIR%" ^
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF ^
  -D CMAKE_COLOR_MAKEFILE:BOOL=ON ^
  -D BUILD_SHARED_LIBS:BOOL=OFF ^
  ^
  -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF ^
  -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON ^
  -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF ^
  -D Trilinos_ENABLE_TESTS:BOOL=OFF ^
  -D Trilinos_ENABLE_EXAMPLES:BOOL=OFF ^
  ^
  -D Trilinos_ASSERT_MISSING_PACKAGES=OFF ^
  -D Trilinos_ENABLE_Gtest:BOOL=OFF ^
  -D Trilinos_VERBOSE_CONFIGURE:BOOL=OFF ^
  ^
  -D Trilinos_ENABLE_Amesos:BOOL=ON ^
    -D Amesos_SHOW_DEPRECATED_WARNINGS:BOOL=OFF ^
  -D Trilinos_ENABLE_Amesos2:BOOL=ON ^
  -D Trilinos_ENABLE_AztecOO:BOOL=ON ^
    -D AztecOO_SHOW_DEPRECATED_WARNINGS:BOOL=OFF ^
  -D Trilinos_ENABLE_Anasazi:BOOL=ON ^
  -D Trilinos_ENABLE_Belos:BOOL=ON ^
  -D Trilinos_ENABLE_Epetra:BOOL=ON ^
    -D Epetra_SHOW_DEPRECATED_WARNINGS:BOOL=OFF ^
    -D Epetra_HIDE_DEPRECATED_CODE:BOOL=ON ^
  -D Trilinos_ENABLE_EpetraExt:BOOL=ON ^
    -D EpetraExt_BUILD_GRAPH_REORDERINGS:BOOL=ON ^
    -D EpetraExt_ENABLE_HDF5:BOOL=OFF ^
    -D EpetraExt_SHOW_DEPRECATED_WARNINGS:BOOL=OFF ^
  -D Trilinos_ENABLE_Intrepid:BOOL=ON ^
  -D Trilinos_ENABLE_Intrepid2:BOOL=ON ^
  -D Trilinos_ENABLE_Ifpack:BOOL=ON ^
    -D Ifpack_SHOW_DEPRECATED_WARNINGS:BOOL=OFF ^
  -D Trilinos_ENABLE_Ifpack2:BOOL=ON ^
  -D Trilinos_ENABLE_Isorropia:BOOL=ON ^
  -D Trilinos_ENABLE_Kokkos:BOOL=ON ^
  -D Trilinos_ENABLE_ML:BOOL=ON ^
    -D ML_SHOW_DEPRECATED_WARNINGS:BOOL=OFF ^
  -D Trilinos_ENABLE_Galeri:BOOL=OFF ^
  -D Trilinos_ENABLE_MueLu:BOOL=ON ^
    -D MueLu_ENABLE_EXAMPLES:BOOL=OFF ^
    -D MueLu_ENABLE_TESTS:BOOL=OFF ^
  -D Trilinos_ENABLE_NOX:BOOL=ON ^
    -D NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA:BOOL=OFF ^
    -D NOX_ENABLE_STRATIMIKOS_EPETRA_STACK:BOOL=OFF ^
  -D Trilinos_ENABLE_Sacado:BOOL=ON ^
  -D Trilinos_ENABLE_SEACASExodus:BOOL=ON ^
  -D Trilinos_ENABLE_SEACASNemesis:BOOL=OFF ^
  -D Trilinos_ENABLE_Shards:BOOL=ON ^
  -D Trilinos_ENABLE_Stratimikos:BOOL=ON ^
  -D Trilinos_ENABLE_Teko:BOOL=ON ^
  -D Trilinos_ENABLE_Teuchos:BOOL=ON ^
    -D Teuchos_GLOBALLY_REDUCE_UNITTEST_RESULTS:BOOL=ON ^
    -D Teuchos_HIDE_DEPRECATED_CODE:BOOL=ON ^
  -D Trilinos_ENABLE_Thyra:BOOL=ON ^
    -D Thyra_SHOW_DEPRECATED_WARNINGS:BOOL=OFF ^
    -D Trilinos_ENABLE_ThyraEpetraAdapters:BOOL=ON ^
    -D Trilinos_ENABLE_ThyraEpetraExtAdapters:BOOL=ON ^
  -D Trilinos_ENABLE_Tpetra:BOOL=ON ^
    -D Tpetra_ENABLE_DEPRECATED_CODE:BOOL=OFF ^
    -D Tpetra_INST_INT_INT:BOOL=ON ^
  -D Trilinos_ENABLE_Xpetra:BOOL=ON ^
    -D Xpetra_ENABLE_DEPRECATED_CODE:BOOL=ON ^
    -D Xpetra_ENABLE_Epetra:BOOL=ON ^
    -D Xpetra_ENABLE_EpetraExt:BOOL=ON ^
    -D Xpetra_ENABLE_TESTS:BOOL=OFF ^
    -D Xpetra_SHOW_DEPRECATED_WARNINGS:BOOL=OFF ^
  -D Trilinos_ENABLE_Zoltan:BOOL=ON ^
  -D Trilinos_ENABLE_Zoltan2:BOOL=ON ^
  ^
  -D Trilinos_MUST_FIND_ALL_TPL_LIBS=TRUE ^
  -D TPL_ENABLE_DLlib:BOOL=OFF ^
  ^
  -D TPL_ENABLE_BLAS:BOOL=ON ^
  -D TPL_BLAS_LIBRARIES="%LIB_DIR%\\lapack\\lib\\libblas.lib" ^
  -D TPL_ENABLE_LAPACK:BOOL=ON ^
  -D TPL_LAPACK_LIBRARIES="%LIB_DIR%\\lapack\\lib\\liblapack.lib" ^
  ^
  -D TPL_ENABLE_Boost:BOOL=OFF ^
  -D TPL_ENABLE_Netcdf:BOOL=ON ^
    -D Netcdf_INCLUDE_DIRS:PATH="%LIB_DIR%\\netcdf\\include" ^
    -D Netcdf_LIBRARY_DIRS:PATH="%LIB_DIR%\\netcdf\\lib" ^
    -D Netcdf_LIBRARY_NAMES="netcdf" ^
  -D TPL_ENABLE_MPI:BOOL=ON ^
  -D TPL_ENABLE_ParMETIS:BOOL=ON ^
    -D ParMETIS_INCLUDE_DIRS:PATH="%LIB_DIR%\\parmetis\\include" ^
    -D ParMETIS_LIBRARY_DIRS:PATH="%LIB_DIR%\\parmetis\\lib" ^
  -D TPL_ENABLE_UMFPACK:BOOL=ON ^
    -D UMFPACK_INCLUDE_DIRS:FILEPATH="%LIB_DIR%\\suitesparse\\include\\suitesparse" ^
    -D UMFPACK_LIBRARY_DIRS:FILEPATH="%LIB_DIR%\\suitesparse\\lib" ^
  -D TPL_ENABLE_AMD:BOOL=ON ^
    -D AMD_INCLUDE_DIRS:FILEPATH="%LIB_DIR%\\suitesparse\\include\\suitesparse" ^
    -D AMD_LIBRARY_DIRS:FILEPATH="%LIB_DIR%\\suitesparse\\lib" ^
  -D TPL_ENABLE_Cholmod:BOOL=ON ^
    -D Cholmod_INCLUDE_DIRS:FILEPATH="%LIB_DIR%\\suitesparse\\include\\suitesparse" ^
    -D Cholmod_LIBRARY_DIRS:FILEPATH="%LIB_DIR%\\suitesparse\\lib" ^
  -D TPL_ENABLE_SuperLUDist:BOOL=ON ^
    -D SuperLUDist_INCLUDE_DIRS:PATH="%LIB_DIR%\\superlu_dist\\include" ^
    -D SuperLUDist_LIBRARY_DIRS:PATH="%LIB_DIR%\\superlu_dist\\lib" ^
    -D SuperLUDist_LIBRARY_NAMES:STRING="superlu_dist" ^
  %TRILINOS_HOME%

if !ERRORLEVEL! neq 0 (
    echo ERROR: installing trilinos failed at configuration
    exit /b !ERRORLEVEL!
)

ninja install -j%NPROCS%

if !ERRORLEVEL! neq 0 (
    echo ERROR: installing trilinos failed at compilation
    exit /b !ERRORLEVEL!
)

rem

cd ..

rem Clean up downloaded and extracted artifacts
for /d %%D in (trilinos*) do rmdir /s /q "%%D"

if !ERRORLEVEL! neq 0 (
    echo ERROR: install.bat failed with exit code !ERRORLEVEL!
    exit /b !ERRORLEVEL!
)

endlocal
