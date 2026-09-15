# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# helper function to run and return on error
function Run($cmd, $args) {
    & $cmd $args
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
}

# Exit the script at the first failure
$ErrorActionPreference = 'Stop'

$env:DEP_DIR = "$env:USERPROFILE\opt"

$curDir = $PWD.Path

cd $curDir/dependencies/windows/parmetis

Run .\install.bat "$env:DEP_DIR\parmetis"

cd $curDir/dependencies/windows/lapack

Run .\install.bat "$env:DEP_DIR\lapack"

# cd $curDir/dependencies/windows/boost

# Run .\install.bat "$env:DEP_DIR\boost"

# cd $curDir/dependencies/windows/cln

# Run .\install.bat "$env:DEP_DIR\cln"

# cd $curDir/dependencies/windows/zlib

# Run .\install.bat "$env:DEP_DIR\zlib"

# cd $curDir/dependencies/windows/hdf5

# Run .\install.bat "$env:DEP_DIR\hdf5"

# cd $curDir/dependencies/windows/netcdf

# Run .\install.bat "$env:DEP_DIR\netcdf"

# cd $curDir/dependencies/windows/suitesparse

# Run .\install.bat "$env:DEP_DIR\suitesparse"

# cd $curDir/dependencies/windows/superlu_dist

# Run .\install.bat "$env:DEP_DIR\superlu_dist"

# cd $curDir/dependencies/windows/trilinos

# Run .\install.bat "$env:DEP_DIR\trilinos"
