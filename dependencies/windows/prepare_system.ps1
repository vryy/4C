# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# Exit the script at the first failure
$ErrorActionPreference = 'Stop'

cmake --version
ninja --version
cl.exe
$PSVersionTable.OS
$env:PROCESSOR_ARCHITECTURE

./dependencies/windows/msmpi/install.ps1
./dependencies/windows/intel/install.ps1
