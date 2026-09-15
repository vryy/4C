# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# Exit the script at the first failure
$ErrorActionPreference = 'Stop'

$env:DEP_DIR = "$env:USERPROFILE\opt"

$curDir = $PWD.Path

cd $curDir/dependencies/windows/parmetis

.\install.bat "$env:DEP_DIR\parmetis"

cd $curDir/dependencies/windows/lapack

.\install.bat "$env:DEP_DIR\lapack"

cd $curDir/dependencies/windows/boost

.\install.bat "$env:DEP_DIR\boost"
