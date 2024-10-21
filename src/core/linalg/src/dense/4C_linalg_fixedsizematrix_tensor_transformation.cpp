// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_fixedsizematrix_tensor_transformation.hpp"

FOUR_C_NAMESPACE_OPEN

template void Core::LinAlg::Tensor::tensor_rotation<3>(const Core::LinAlg::Matrix<3, 3>&,
    const Core::LinAlg::Matrix<3, 3>&, Core::LinAlg::Matrix<3, 3>&);
template void Core::LinAlg::Tensor::inverse_tensor_rotation<3>(const Core::LinAlg::Matrix<3, 3>&,
    const Core::LinAlg::Matrix<3, 3>&, Core::LinAlg::Matrix<3, 3>&);

FOUR_C_NAMESPACE_CLOSE
