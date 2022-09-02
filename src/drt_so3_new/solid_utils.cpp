/*----------------------------------------------------------------------*/
/*! \file

\brief little helpers for solid elements

\level 1
*----------------------------------------------------------------------*/

#include "solid_utils.H"
#include "../drt_lib/drt_element_integration_select.H"
#include "../drt_lib/drt_element.H"

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
  VectorToMatrix(pk2, false, pk2_matrix);

  static LINALG::Matrix<3, 3> cauchy_matrix;
  cauchy_matrix.Multiply(defgrd, pk2_matrix);
  pk2_matrix.MultiplyNT(defgrd.Determinant(), cauchy_matrix, defgrd, 0.);

  MatrixToVector(pk2_matrix, false, cauchy);
}

void STR::UTILS::MatrixToVector(
    const LINALG::Matrix<3, 3>& matrix, const bool strain_like, LINALG::Matrix<6, 1>& vector)
{
  for (int i = 0; i < 3; ++i) vector(i) = matrix(i, i);
  vector(3) = matrix(0, 1) * (1. + strain_like);
  vector(4) = matrix(2, 1) * (1. + strain_like);
  vector(5) = matrix(0, 2) * (1. + strain_like);
}

void STR::UTILS::VectorToMatrix(
    const LINALG::Matrix<6, 1>& vector, const bool strain_like, LINALG::Matrix<3, 3>& matrix)
{
  for (int i = 0; i < 3; ++i) matrix(i, i) = vector(i);
  matrix(0, 1) = matrix(1, 0) = vector(3) / (1. + strain_like);
  matrix(2, 1) = matrix(1, 2) = vector(4) / (1. + strain_like);
  matrix(0, 2) = matrix(2, 0) = vector(5) / (1. + strain_like);
}
