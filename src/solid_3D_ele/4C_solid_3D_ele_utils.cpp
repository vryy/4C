/*! \file

\brief Helpers for solid elements

\level 1
*/

#include "4C_solid_3D_ele_utils.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_linalg_utils_densematrix_eigen.hpp"
#include "4C_solid_3D_ele_properties.hpp"

FOUR_C_NAMESPACE_OPEN


void STR::UTILS::pk2_to_cauchy(const Core::LinAlg::Matrix<6, 1>& pk2,
    const Core::LinAlg::Matrix<3, 3>& defgrd, Core::LinAlg::Matrix<6, 1>& cauchy)
{
  Core::LinAlg::Matrix<3, 3> S_matrix;
  Core::LinAlg::Voigt::Stresses::VectorToMatrix(pk2, S_matrix);

  Core::LinAlg::Matrix<3, 3> FS;
  FS.MultiplyNN(defgrd, S_matrix);

  Core::LinAlg::Matrix<3, 3> cauchy_matrix;
  cauchy_matrix.MultiplyNT(1.0 / defgrd.Determinant(), FS, defgrd, 0.0);

  Core::LinAlg::Voigt::Stresses::MatrixToVector(cauchy_matrix, cauchy);
}

Core::LinAlg::Matrix<6, 1> STR::UTILS::green_lagrange_to_euler_almansi(
    const Core::LinAlg::Matrix<6, 1>& gl, const Core::LinAlg::Matrix<3, 3>& defgrd)
{
  Core::LinAlg::Matrix<3, 3> invdefgrd(defgrd);
  invdefgrd.Invert();

  Core::LinAlg::Matrix<3, 3> E_matrix;
  Core::LinAlg::Voigt::Strains::VectorToMatrix(gl, E_matrix);

  Core::LinAlg::Matrix<3, 3> iFTE;
  iFTE.MultiplyTN(invdefgrd, E_matrix);

  Core::LinAlg::Matrix<3, 3> ea_matrix;
  ea_matrix.MultiplyNN(iFTE, invdefgrd);

  Core::LinAlg::Matrix<6, 1> ea;
  Core::LinAlg::Voigt::Strains::MatrixToVector(ea_matrix, ea);
  return ea;
}

Core::LinAlg::Matrix<6, 1> STR::UTILS::green_lagrange_to_log_strain(
    const Core::LinAlg::Matrix<6, 1>& gl)
{
  Core::LinAlg::Matrix<3, 3> E_matrix;
  Core::LinAlg::Voigt::Strains::VectorToMatrix(gl, E_matrix);

  Core::LinAlg::Matrix<3, 3> pr_strain(true);  // squared principal strains
  Core::LinAlg::Matrix<3, 3> pr_dir(true);     // principal directions
  Core::LinAlg::SYEV(E_matrix, pr_strain, pr_dir);

  // compute principal logarithmic strains
  Core::LinAlg::Matrix<3, 3> pr_log_strain(true);
  for (int i = 0; i < 3; ++i) pr_log_strain(i, i) = std::log(std::sqrt(2 * pr_strain(i, i) + 1.0));

  // create logarithmic strain matrix
  Core::LinAlg::Matrix<3, 3> log_strain_matrix(true);
  Core::LinAlg::Matrix<3, 3> VH(false);
  VH.MultiplyNN(pr_dir, pr_log_strain);
  log_strain_matrix.MultiplyNT(VH, pr_dir);

  // convert to strain-like voigt notation
  Core::LinAlg::Matrix<6, 1> log_strain_voigt(true);
  Core::LinAlg::Voigt::Strains::MatrixToVector(log_strain_matrix, log_strain_voigt);
  return log_strain_voigt;
}

int STR::UTILS::ReadElement::read_element_material(Input::LineDefinition* linedef)
{
  int material = 0;
  linedef->ExtractInt("MAT", material);
  return material;
}

Inpar::STR::KinemType STR::UTILS::ReadElement::read_element_kinematic_type(
    Input::LineDefinition* linedef)
{
  std::string kinem;
  linedef->ExtractString("KINEM", kinem);
  if (kinem == "nonlinear")
    return Inpar::STR::KinemType::nonlinearTotLag;
  else if (kinem == "linear")
    return Inpar::STR::KinemType::linear;
  else
  {
    FOUR_C_THROW("unknown kinematic type %s", kinem.c_str());
    return Inpar::STR::KinemType::vague;
  }
}

Discret::ELEMENTS::ElementTechnology STR::UTILS::ReadElement::read_element_technology(
    Input::LineDefinition* linedef)
{
  std::string type;
  linedef->ExtractString("TECH", type);
  if (type == "fbar")
  {
    return Discret::ELEMENTS::ElementTechnology::fbar;
  }
  else if (type == "eas_full")
  {
    return Discret::ELEMENTS::ElementTechnology::eas_full;
  }
  else if (type == "eas_mild")
  {
    return Discret::ELEMENTS::ElementTechnology::eas_mild;
  }
  else if (type == "none")
  {
    return Discret::ELEMENTS::ElementTechnology::none;
  }
  else
    FOUR_C_THROW("unrecognized element technology type %s", type.c_str());
}

Discret::ELEMENTS::PrestressTechnology STR::UTILS::ReadElement::read_prestress_technology(
    Input::LineDefinition* linedef)
{
  std::string type;
  linedef->ExtractString("PRESTRESS_TECH", type);
  if (type == "none")
  {
    return Discret::ELEMENTS::PrestressTechnology::none;
  }
  else if (type == "mulf")
  {
    return Discret::ELEMENTS::PrestressTechnology::mulf;
  }

  FOUR_C_THROW("unrecognized prestress technology type %s", type.c_str());
}

Discret::ELEMENTS::SolidElementProperties STR::UTILS::ReadElement::read_solid_element_properties(
    Input::LineDefinition* linedef)
{
  Discret::ELEMENTS::SolidElementProperties solid_properties{};
  // element technology
  if (linedef->HaveNamed("TECH"))
  {
    solid_properties.element_technology = STR::UTILS::ReadElement::read_element_technology(linedef);
  }

  // prestress technology
  if (linedef->HaveNamed("PRESTRESS_TECH"))
  {
    solid_properties.prestress_technology =
        STR::UTILS::ReadElement::read_prestress_technology(linedef);
  }
  // kinematic type
  solid_properties.kintype = STR::UTILS::ReadElement::read_element_kinematic_type(linedef);

  return solid_properties;
}

void STR::UTILS::nodal_block_information_solid(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;

  nv = 3;
}

FOUR_C_NAMESPACE_CLOSE
