/*----------------------------------------------------------------------*/
/*! \file

\brief Provide 1D, 2D and 3D shape functions

The corresponding graphics and a detailed description can be found in
the Baci guide in the Convention chapter.

\level 0


*----------------------------------------------------------------------*/

#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <::DRT::Element::DiscretizationType distype, int probdim>
void CORE::DRT::UTILS::EvaluateShapeFunctionSpatialDerivativeInProbDim(
    CORE::LINALG::Matrix<probdim,
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& deriv_xyz,
    const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim,
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& deriv,
    const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        probdim>& xyze,
    const CORE::LINALG::Matrix<probdim, 1>& normal)
{
  const int nen = CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
  const int nsd_ele = CORE::DRT::UTILS::DisTypeToDim<distype>::dim;

  dsassert(nsd_ele != probdim,
      "This method is designed to be used if the dimension of the element is smaller than the "
      "problem dimension");

  CORE::LINALG::Matrix<probdim, nen> deriv_full;

  // transform the derivatives and Jacobians to the higher dimensional coordinates (problem
  // dimension)
  CORE::LINALG::Matrix<nsd_ele, probdim> dx_dr_red;
  CORE::LINALG::Matrix<probdim, probdim> dx_dr, dr_dx;
  dx_dr_red.MultiplyNN(deriv, xyze);

  for (unsigned i = 0; i < probdim; ++i)
  {
    for (unsigned j = 0; j < nsd_ele; ++j) dx_dr(j, i) = dx_dr_red(j, i);
    dx_dr(nsd_ele, i) = normal(i);
  }

  for (unsigned i = 0; i < nen; ++i)
  {
    for (unsigned j = 0; j < nsd_ele; ++j) deriv_full(j, i) = deriv(j, i);
    deriv_full(nsd_ele, i) = 0.0;
  }

  // special case: 1D element embedded in 3D problem
  if (nsd_ele == 1 and probdim == 3)
  {
    // compute second unit normal
    const double normalvec2_0 = dx_dr_red(0, 1) * normal(2) - normal(1) * dx_dr_red(0, 2);
    const double normalvec2_1 = dx_dr_red(0, 2) * normal(0) - normal(2) * dx_dr_red(0, 0);
    const double normalvec2_2 = dx_dr_red(0, 0) * normal(1) - normal(0) * dx_dr_red(0, 1);

    // norm
    const double norm2 = std::sqrt(
        normalvec2_0 * normalvec2_0 + normalvec2_1 * normalvec2_1 + normalvec2_2 * normalvec2_2);

    dx_dr(2, 0) = normalvec2_0 / norm2;
    dx_dr(2, 1) = normalvec2_1 / norm2;
    dx_dr(2, 2) = normalvec2_2 / norm2;

    for (unsigned i = 0; i < nen; i++) deriv_full(2, i) = 0.0;
  }

  dr_dx.Invert(dx_dr);

  // compute global spatial derivatives
  deriv_xyz.Multiply(dr_dx, deriv_full);
}

template void CORE::DRT::UTILS::EvaluateShapeFunctionSpatialDerivativeInProbDim<
    DRT::Element::DiscretizationType::line2, 2>(CORE::LINALG::Matrix<2, 2>& deriv_xyz,
    const CORE::LINALG::Matrix<1, 2>& deriv, const CORE::LINALG::Matrix<2, 2>& xyze,
    const CORE::LINALG::Matrix<2, 1>& normal);
template void CORE::DRT::UTILS::EvaluateShapeFunctionSpatialDerivativeInProbDim<
    DRT::Element::DiscretizationType::line2, 3>(CORE::LINALG::Matrix<3, 2>& deriv_xyz,
    const CORE::LINALG::Matrix<1, 2>& deriv, const CORE::LINALG::Matrix<2, 3>& xyze,
    const CORE::LINALG::Matrix<3, 1>& normal);
template void CORE::DRT::UTILS::EvaluateShapeFunctionSpatialDerivativeInProbDim<
    DRT::Element::DiscretizationType::line3, 2>(CORE::LINALG::Matrix<2, 3>& deriv_xyz,
    const CORE::LINALG::Matrix<1, 3>& deriv, const CORE::LINALG::Matrix<3, 2>& xyze,
    const CORE::LINALG::Matrix<2, 1>& normal);
template void CORE::DRT::UTILS::EvaluateShapeFunctionSpatialDerivativeInProbDim<
    DRT::Element::DiscretizationType::nurbs3, 2>(CORE::LINALG::Matrix<2, 3>& deriv_xyz,
    const CORE::LINALG::Matrix<1, 3>& deriv, const CORE::LINALG::Matrix<3, 2>& xyze,
    const CORE::LINALG::Matrix<2, 1>& normal);
template void CORE::DRT::UTILS::EvaluateShapeFunctionSpatialDerivativeInProbDim<
    DRT::Element::DiscretizationType::nurbs9, 3>(CORE::LINALG::Matrix<3, 9>& deriv_xyz,
    const CORE::LINALG::Matrix<2, 9>& deriv, const CORE::LINALG::Matrix<9, 3>& xyze,
    const CORE::LINALG::Matrix<3, 1>& normal);
template void CORE::DRT::UTILS::EvaluateShapeFunctionSpatialDerivativeInProbDim<
    DRT::Element::DiscretizationType::quad4, 3>(CORE::LINALG::Matrix<3, 4>& deriv_xyz,
    const CORE::LINALG::Matrix<2, 4>& deriv, const CORE::LINALG::Matrix<4, 3>& xyze,
    const CORE::LINALG::Matrix<3, 1>& normal);
template void CORE::DRT::UTILS::EvaluateShapeFunctionSpatialDerivativeInProbDim<
    DRT::Element::DiscretizationType::quad8, 3>(CORE::LINALG::Matrix<3, 8>& deriv_xyz,
    const CORE::LINALG::Matrix<2, 8>& deriv, const CORE::LINALG::Matrix<8, 3>& xyze,
    const CORE::LINALG::Matrix<3, 1>& normal);
template void CORE::DRT::UTILS::EvaluateShapeFunctionSpatialDerivativeInProbDim<
    DRT::Element::DiscretizationType::quad9, 3>(CORE::LINALG::Matrix<3, 9>& deriv_xyz,
    const CORE::LINALG::Matrix<2, 9>& deriv, const CORE::LINALG::Matrix<9, 3>& xyze,
    const CORE::LINALG::Matrix<3, 1>& normal);
template void CORE::DRT::UTILS::EvaluateShapeFunctionSpatialDerivativeInProbDim<
    DRT::Element::DiscretizationType::tri3, 3>(CORE::LINALG::Matrix<3, 3>& deriv_xyz,
    const CORE::LINALG::Matrix<2, 3>& deriv, const CORE::LINALG::Matrix<3, 3>& xyze,
    const CORE::LINALG::Matrix<3, 1>& normal);
template void CORE::DRT::UTILS::EvaluateShapeFunctionSpatialDerivativeInProbDim<
    DRT::Element::DiscretizationType::tri6, 3>(CORE::LINALG::Matrix<3, 6>& deriv_xyz,
    const CORE::LINALG::Matrix<2, 6>& deriv, const CORE::LINALG::Matrix<6, 3>& xyze,
    const CORE::LINALG::Matrix<3, 1>& normal);