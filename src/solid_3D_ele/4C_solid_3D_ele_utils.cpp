/*! \file

\brief Helpers for solid elements

\level 1
*/

#include "4C_solid_3D_ele_utils.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_lib_element.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_linalg_utils_densematrix_eigen.hpp"
#include "4C_solid_3D_ele_properties.hpp"

FOUR_C_NAMESPACE_OPEN


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

CORE::LINALG::Matrix<6, 1> STR::UTILS::GreenLagrangeToLogStrain(
    const CORE::LINALG::Matrix<6, 1>& gl)
{
  CORE::LINALG::Matrix<3, 3> E_matrix;
  CORE::LINALG::VOIGT::Strains::VectorToMatrix(gl, E_matrix);

  CORE::LINALG::Matrix<3, 3> pr_strain(true);  // squared principal strains
  CORE::LINALG::Matrix<3, 3> pr_dir(true);     // principal directions
  CORE::LINALG::SYEV(E_matrix, pr_strain, pr_dir);

  // compute principal logarithmic strains
  CORE::LINALG::Matrix<3, 3> pr_log_strain(true);
  for (int i = 0; i < 3; ++i) pr_log_strain(i, i) = std::log(std::sqrt(2 * pr_strain(i, i) + 1.0));

  // create logarithmic strain matrix
  CORE::LINALG::Matrix<3, 3> log_strain_matrix(true);
  CORE::LINALG::Matrix<3, 3> VH(false);
  VH.MultiplyNN(pr_dir, pr_log_strain);
  log_strain_matrix.MultiplyNT(VH, pr_dir);

  // convert to strain-like voigt notation
  CORE::LINALG::Matrix<6, 1> log_strain_voigt(true);
  CORE::LINALG::VOIGT::Strains::MatrixToVector(log_strain_matrix, log_strain_voigt);
  return log_strain_voigt;
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
    return INPAR::STR::KinemType::nonlinearTotLag;
  else if (kinem == "linear")
    return INPAR::STR::KinemType::linear;
  else
  {
    FOUR_C_THROW("unknown kinematic type %s", kinem.c_str());
    return INPAR::STR::KinemType::vague;
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
    FOUR_C_THROW("unrecognized element technology type %s", type.c_str());
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

  FOUR_C_THROW("unrecognized prestress technology type %s", type.c_str());
}

void STR::UTILS::NodalBlockInformationSolid(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;

  nv = 3;
}

FOUR_C_NAMESPACE_CLOSE
