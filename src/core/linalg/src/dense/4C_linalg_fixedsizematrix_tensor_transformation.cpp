/*! \file
\level 1
\brief Util functions for tensor transformations
*/

#include "4C_linalg_fixedsizematrix_tensor_transformation.hpp"

FOUR_C_NAMESPACE_OPEN

template void Core::LinAlg::Tensor::tensor_rotation<3>(const Core::LinAlg::Matrix<3, 3>&,
    const Core::LinAlg::Matrix<3, 3>&, Core::LinAlg::Matrix<3, 3>&);
template void Core::LinAlg::Tensor::inverse_tensor_rotation<3>(const Core::LinAlg::Matrix<3, 3>&,
    const Core::LinAlg::Matrix<3, 3>&, Core::LinAlg::Matrix<3, 3>&);

FOUR_C_NAMESPACE_CLOSE
