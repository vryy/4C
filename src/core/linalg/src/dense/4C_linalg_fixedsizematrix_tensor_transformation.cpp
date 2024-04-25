/*! \file
\level 1
\brief Util functions for tensor transformations
*/

#include "4C_linalg_fixedsizematrix_tensor_transformation.hpp"

FOUR_C_NAMESPACE_OPEN

template void CORE::LINALG::TENSOR::TensorRotation<3>(const CORE::LINALG::Matrix<3, 3>&,
    const CORE::LINALG::Matrix<3, 3>&, CORE::LINALG::Matrix<3, 3>&);
template void CORE::LINALG::TENSOR::InverseTensorRotation<3>(const CORE::LINALG::Matrix<3, 3>&,
    const CORE::LINALG::Matrix<3, 3>&, CORE::LINALG::Matrix<3, 3>&);

FOUR_C_NAMESPACE_CLOSE
