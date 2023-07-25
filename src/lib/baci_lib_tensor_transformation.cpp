/*! \file
\level 1
\brief Util functions for tensor transformations
*/

#include "baci_lib_tensor_transformation.H"

template void UTILS::TENSOR::TensorRotation<3>(const CORE::LINALG::Matrix<3, 3>&,
    const CORE::LINALG::Matrix<3, 3>&, CORE::LINALG::Matrix<3, 3>&);
template void UTILS::TENSOR::InverseTensorRotation<3>(const CORE::LINALG::Matrix<3, 3>&,
    const CORE::LINALG::Matrix<3, 3>&, CORE::LINALG::Matrix<3, 3>&);