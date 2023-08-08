/*----------------------------------------------------------------------*/
/*! \file

\brief little helpers for solid elements

\level 1
*----------------------------------------------------------------------*/

#include "solid_ele_utils.H"
#include "lib_element_integration_select.H"
#include "lib_element.H"
#include "lib_voigt_notation.H"
#include "lib_linedefinition.H"

int STR::UTILS::DisTypeToNgpOptGaussRule(DRT::Element::DiscretizationType distype)
{
  switch (distype)
  {
    case DRT::Element::hex8:
      return CORE::DRT::UTILS::IntPointsAndWeights<3>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::hex8>::rule)
          .IP()
          .nquad;
    case DRT::Element::hex18:
      return CORE::DRT::UTILS::IntPointsAndWeights<3>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::hex18>::rule)
          .IP()
          .nquad;
    case DRT::Element::hex20:
      return CORE::DRT::UTILS::IntPointsAndWeights<3>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::hex20>::rule)
          .IP()
          .nquad;
    case DRT::Element::hex27:
      return CORE::DRT::UTILS::IntPointsAndWeights<3>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::hex27>::rule)
          .IP()
          .nquad;
    case DRT::Element::nurbs27:
      return CORE::DRT::UTILS::IntPointsAndWeights<3>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::nurbs27>::rule)
          .IP()
          .nquad;
    case DRT::Element::tet4:
      return CORE::DRT::UTILS::IntPointsAndWeights<3>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::tet4>::rule)
          .IP()
          .nquad;
    case DRT::Element::tet10:
      return CORE::DRT::UTILS::IntPointsAndWeights<3>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::tet10>::rule)
          .IP()
          .nquad;
    case DRT::Element::wedge6:
      return CORE::DRT::UTILS::IntPointsAndWeights<3>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::wedge6>::rule)
          .IP()
          .nquad;
    case DRT::Element::pyramid5:
      return CORE::DRT::UTILS::IntPointsAndWeights<3>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::pyramid5>::rule)
          .IP()
          .nquad;
    case DRT::Element::quad4:
      return CORE::DRT::UTILS::IntPointsAndWeights<2>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::quad4>::rule)
          .IP()
          .nquad;
    case DRT::Element::quad8:
      return CORE::DRT::UTILS::IntPointsAndWeights<2>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::quad8>::rule)
          .IP()
          .nquad;
    case DRT::Element::quad9:
      return CORE::DRT::UTILS::IntPointsAndWeights<2>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::quad9>::rule)
          .IP()
          .nquad;
    case DRT::Element::nurbs9:
      return CORE::DRT::UTILS::IntPointsAndWeights<2>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::nurbs9>::rule)
          .IP()
          .nquad;
    case DRT::Element::tri3:
      return CORE::DRT::UTILS::IntPointsAndWeights<2>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::tri3>::rule)
          .IP()
          .nquad;
    case DRT::Element::tri6:
      return CORE::DRT::UTILS::IntPointsAndWeights<2>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::tri6>::rule)
          .IP()
          .nquad;
    case DRT::Element::line2:
      return CORE::DRT::UTILS::IntPointsAndWeights<1>(
          DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::line2>::rule)
          .IP()
          .nquad;
    case DRT::Element::line3:
      return CORE::DRT::UTILS::IntPointsAndWeights<1>(
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
  LINALG::Matrix<3, 3> S_matrix;
  ::UTILS::VOIGT::Stresses::VectorToMatrix(pk2, S_matrix);

  LINALG::Matrix<3, 3> FS;
  FS.MultiplyNN(defgrd, S_matrix);

  LINALG::Matrix<3, 3> cauchy_matrix;
  cauchy_matrix.MultiplyNT(1.0 / defgrd.Determinant(), FS, defgrd, 0.0);

  ::UTILS::VOIGT::Stresses::MatrixToVector(cauchy_matrix, cauchy);
}

LINALG::Matrix<6, 1> STR::UTILS::GreenLagrangeToEulerAlmansi(
    const LINALG::Matrix<6, 1>& gl, const LINALG::Matrix<3, 3>& defgrd)
{
  LINALG::Matrix<3, 3> invdefgrd(defgrd);
  invdefgrd.Invert();

  LINALG::Matrix<3, 3> E_matrix;
  ::UTILS::VOIGT::Strains::VectorToMatrix(gl, E_matrix);

  LINALG::Matrix<3, 3> iFTE;
  iFTE.MultiplyTN(invdefgrd, E_matrix);

  LINALG::Matrix<3, 3> ea_matrix;
  ea_matrix.MultiplyNN(iFTE, invdefgrd);

  LINALG::Matrix<6, 1> ea;
  ::UTILS::VOIGT::Strains::MatrixToVector(ea_matrix, ea);
  return ea;
}

int STR::UTILS::READELEMENT::ReadElementMaterial(DRT::INPUT::LineDefinition* linedef)
{
  int material = 0;
  linedef->ExtractInt("MAT", material);
  return material;
}

INPAR::STR::KinemType STR::UTILS::READELEMENT::ReadElementKinematicType(
    DRT::INPUT::LineDefinition* linedef)
{
  std::string kinem;
  linedef->ExtractString("KINEM", kinem);
  if (kinem == "nonlinear")
    return INPAR::STR::kinem_nonlinearTotLag;
  else if (kinem == "linear")
    return INPAR::STR::kinem_linear;
  else
  {
    dserror("unknown kinematic type %s", kinem.c_str());
    return INPAR::STR::kinem_vague;
  }
}

void STR::UTILS::READELEMENT::ReadAndSetEAS(DRT::INPUT::LineDefinition* linedef,
    ::STR::ELEMENTS::EASType& eastype, std::set<INPAR::STR::EleTech>& eletech)
{
  std::string type;
  linedef->ExtractString("EAS", type);
  if (type == "mild")
  {
    eastype = ::STR::ELEMENTS::EASType::eastype_h8_9;
    eletech.insert(INPAR::STR::EleTech::eas);
  }
  else if (type == "full")
  {
    eastype = ::STR::ELEMENTS::EASType::eastype_h8_21;
    eletech.insert(INPAR::STR::EleTech::eas);
  }
  else if (type == "none")
  {
    eastype = ::STR::ELEMENTS::EASType::soh8_easnone;
  }
  else
    dserror("unrecognized eas type for hex8: %s", type.c_str());
}

void STR::UTILS::NodalBlockInformationSolid(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  // todo: to this combined for 2D and 3D
  numdf = 3;
  dimns = 6;

  nv = 3;
}
