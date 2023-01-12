/*! \file
\level 1
\brief Util functions for tensor transformations
*/

#include "lib_tensor_transformation.H"

template void UTILS::TENSOR::TensorRotation<3>(
    const LINALG::Matrix<3, 3>&, const LINALG::Matrix<3, 3>&, LINALG::Matrix<3, 3>&);
template void UTILS::TENSOR::InverseTensorRotation<3>(
    const LINALG::Matrix<3, 3>&, const LINALG::Matrix<3, 3>&, LINALG::Matrix<3, 3>&);