# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# Exit the script at the first failure
$ErrorActionPreference = 'Stop'
$PSNativeCommandUseErrorActionPreference = $true

$env:DEP_DIR = "$env:USERPROFILE\opt"

$curDir = $PWD.Path

cd $curDir/dependencies/windows/parmetis

.\install.bat "$env:DEP_DIR\parmetis"

cd $curDir/dependencies/windows/lapack

.\install.bat "$env:DEP_DIR\lapack"

# cd $curDir/dependencies/windows/boost

# .\install.bat "$env:DEP_DIR\boost"

# cd $curDir/dependencies/windows/cln

# .\install.bat "$env:DEP_DIR\cln"

cd $curDir/dependencies/windows/zlib

.\install.bat "$env:DEP_DIR\zlib"

cd $curDir/dependencies/windows/hdf5

.\install.bat "$env:DEP_DIR\hdf5"

cd $curDir/dependencies/windows/netcdf

.\install.bat "$env:DEP_DIR\netcdf"

cd $curDir/dependencies/windows/qhull

.\install.bat "$env:DEP_DIR\qhull"

cd $curDir/dependencies/windows/fftw

.\install.bat "$env:DEP_DIR\fftw"

cd $curDir/dependencies/windows/suitesparse

.\install.bat "$env:DEP_DIR\suitesparse"

cd $curDir/dependencies/windows/superlu_dist

.\install.bat "$env:DEP_DIR\superlu_dist"

cd $curDir/dependencies/windows/trilinos

.\install.bat "$env:DEP_DIR\trilinos"
