/*! \file

\brief Helpers for solid elements

\level 1
*/

#include "baci_solid_ele_utils.H"

#include "baci_io_linedefinition.H"
#include "baci_lib_element.H"
#include "baci_linalg_fixedsizematrix_voigt_notation.H"
#include "baci_solid_ele_properties.H"

BACI_NAMESPACE_OPEN


void STR::UTILS::Pk2ToCauchy(const CORE::LINALG::Matrix<6, 1>& pk2,
    const CORE::LINALG::Matrix<3, 3>& defgrd, CORE::LINALG::Matrix<6, 1>& cauchy)
{
  CORE::LINALG::Matrix<3, 3> S_matrix;
  CORE::LINALG::VOIGT::Stresses::VectorToMatrix(pk2, S_matrix);

  CORE::LINALG::Matrix<3, 3> FS;
  FS.MultiplyNN(defgrd, S_matrix);

  CORE::LINALG::Matrix<3, 3> cauchy_matrix;
  cauchy_matrix.MultiplyNT(1.0 / defgrd.Determinant(), FS, defgrd, 0.0);

  CORE::LINALG::VOIGT::Stresses::MatrixToVector(cauchy_matrix, cauchy);
}

CORE::LINALG::Matrix<6, 1> STR::UTILS::GreenLagrangeToEulerAlmansi(
    const CORE::LINALG::Matrix<6, 1>& gl, const CORE::LINALG::Matrix<3, 3>& defgrd)
{
  CORE::LINALG::Matrix<3, 3> invdefgrd(defgrd);
  invdefgrd.Invert();

  CORE::LINALG::Matrix<3, 3> E_matrix;
  CORE::LINALG::VOIGT::Strains::VectorToMatrix(gl, E_matrix);

  CORE::LINALG::Matrix<3, 3> iFTE;
  iFTE.MultiplyTN(invdefgrd, E_matrix);

  CORE::LINALG::Matrix<3, 3> ea_matrix;
  ea_matrix.MultiplyNN(iFTE, invdefgrd);

  CORE::LINALG::Matrix<6, 1> ea;
  CORE::LINALG::VOIGT::Strains::MatrixToVector(ea_matrix, ea);
  return ea;
}

int STR::UTILS::READELEMENT::ReadElementMaterial(INPUT::LineDefinition* linedef)
{
  int material = 0;
  linedef->ExtractInt("MAT", material);
  return material;
}

INPAR::STR::KinemType STR::UTILS::READELEMENT::ReadElementKinematicType(
    INPUT::LineDefinition* linedef)
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

DRT::ELEMENTS::ElementTechnology STR::UTILS::READELEMENT::ReadElementTechnology(
    INPUT::LineDefinition* linedef)
{
  std::string type;
  linedef->ExtractString("TECH", type);
  if (type == "fbar")
  {
    return DRT::ELEMENTS::ElementTechnology::fbar;
  }
  else if (type == "eas_full")
  {
    return DRT::ELEMENTS::ElementTechnology::eas_full;
  }
  else if (type == "eas_mild")
  {
    return DRT::ELEMENTS::ElementTechnology::eas_mild;
  }
  else if (type == "none")
  {
    return DRT::ELEMENTS::ElementTechnology::none;
  }
  else
    dserror("unrecognized element technology type %s", type.c_str());
}

DRT::ELEMENTS::PrestressTechnology STR::UTILS::READELEMENT::ReadPrestressTechnology(
    INPUT::LineDefinition* linedef)
{
  std::string type;
  linedef->ExtractString("PRESTRESS_TECH", type);
  if (type == "none")
  {
    return DRT::ELEMENTS::PrestressTechnology::none;
  }
  else if (type == "mulf")
  {
    return DRT::ELEMENTS::PrestressTechnology::mulf;
  }

  dserror("unrecognized prestress technology type %s", type.c_str());
}

void STR::UTILS::NodalBlockInformationSolid(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;

  nv = 3;
}

BACI_NAMESPACE_CLOSE
