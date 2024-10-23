# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

#
# Set up the infrastructure that we need to manage compiler settings
#

# Define a target which pulls in all the compiler settings and definitions that we want to use to compile our own files.
# This target is filled with all compiler and linker features that are detected in the next steps.
add_library(four_c_private_compile_interface INTERFACE)
