/*----------------------------------------------------------------------*/
/*! \file

\brief little helpers for solid elements

\level 1
*----------------------------------------------------------------------*/

#include "solid_utils.H"
#include "../drt_lib/drt_element_integration_select.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/voigt_notation.H"

int STR::UTILS::DisTypeToNgpOptGaussRule(DRT::Element::DiscretizationType distype)
{
  switch (distype)
  {
    case DRT::Element::hex8:
      return DRT::UTILS::IntPointsAndWeights<3>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::hex8>::rule)
          .IP()
          .nquad;
    case DRT::Element::hex18:
      return DRT::UTILS::IntPointsAndWeights<3>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::hex18>::rule)
          .IP()
          .nquad;
    case DRT::Element::hex20:
      return DRT::UTILS::IntPointsAndWeights<3>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::hex20>::rule)
          .IP()
          .nquad;
    case DRT::Element::hex27:
      return DRT::UTILS::IntPointsAndWeights<3>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::hex27>::rule)
          .IP()
          .nquad;
    case DRT::Element::nurbs27:
      return DRT::UTILS::IntPointsAndWeights<3>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::nurbs27>::rule)
          .IP()
          .nquad;
    case DRT::Element::tet4:
      return DRT::UTILS::IntPointsAndWeights<3>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::tet4>::rule)
          .IP()
          .nquad;
    case DRT::Element::tet10:
      return DRT::UTILS::IntPointsAndWeights<3>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::tet10>::rule)
          .IP()
          .nquad;
    case DRT::Element::wedge6:
      return DRT::UTILS::IntPointsAndWeights<3>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::wedge6>::rule)
          .IP()
          .nquad;
    case DRT::Element::pyramid5:
      return DRT::UTILS::IntPointsAndWeights<3>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::pyramid5>::rule)
          .IP()
          .nquad;
    case DRT::Element::quad4:
      return DRT::UTILS::IntPointsAndWeights<2>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::quad4>::rule)
          .IP()
          .nquad;
    case DRT::Element::quad8:
      return DRT::UTILS::IntPointsAndWeights<2>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::quad8>::rule)
          .IP()
          .nquad;
    case DRT::Element::quad9:
      return DRT::UTILS::IntPointsAndWeights<2>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::quad9>::rule)
          .IP()
          .nquad;
    case DRT::Element::nurbs9:
      return DRT::UTILS::IntPointsAndWeights<2>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::nurbs9>::rule)
          .IP()
          .nquad;
    case DRT::Element::tri3:
      return DRT::UTILS::IntPointsAndWeights<2>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::tri3>::rule)
          .IP()
          .nquad;
    case DRT::Element::tri6:
      return DRT::UTILS::IntPointsAndWeights<2>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::tri6>::rule)
          .IP()
          .nquad;
    case DRT::Element::line2:
      return DRT::UTILS::IntPointsAndWeights<1>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::line2>::rule)
          .IP()
          .nquad;
    case DRT::Element::line3:
      return DRT::UTILS::IntPointsAndWeights<1>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::line3>::rule)
          .IP()
          .nquad;

    default:
      dserror("not implemented");
  }
  return -1;
}

void STR::UTILS::Pk2ToCauchy(const LINALG::Matrix<6, 1>& pk2, const LINALG::Matrix<3, 3>& defgrd,
    LINALG::Matrix<6, 1>& cauchy)
{
  static LINALG::Matrix<3, 3> pk2_matrix;
  ::UTILS::VOIGT::Stresses::VectorToMatrix(pk2, pk2_matrix);

  static LINALG::Matrix<3, 3> cauchy_matrix;
  cauchy_matrix.Multiply(defgrd, pk2_matrix);
  pk2_matrix.MultiplyNT(defgrd.Determinant(), cauchy_matrix, defgrd, 0.);

  ::UTILS::VOIGT::Stresses::MatrixToVector(pk2_matrix, cauchy);
}
