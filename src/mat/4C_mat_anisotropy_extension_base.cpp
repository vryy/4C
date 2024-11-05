// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_anisotropy_extension_base.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_mat_anisotropy.hpp"
#include "4C_utils_shared_ptr_from_ref.hpp"

FOUR_C_NAMESPACE_OPEN

void Mat::BaseAnisotropyExtension::set_anisotropy(Mat::Anisotropy& anisotropy)
{
  anisotropy_ = Core::Utils::shared_ptr_from_ref(anisotropy);
}

FOUR_C_NAMESPACE_CLOSE
