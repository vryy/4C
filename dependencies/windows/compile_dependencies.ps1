# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# Exit the script at the first failure
$ErrorActionPreference = 'Stop'

set "DEP_DIR=%USERPROFILE%\opt"

cd dependencies/windows/parmetis

.\install.bat %DEP_DIR%\parmetis
