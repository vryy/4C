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


void Solid::UTILS::pk2_to_cauchy(const Core::LinAlg::Matrix<6, 1>& pk2,
    const Core::LinAlg::Matrix<3, 3>& defgrd, Core::LinAlg::Matrix<6, 1>& cauchy)
{
  Core::LinAlg::Matrix<3, 3> S_matrix;
  Core::LinAlg::Voigt::Stresses::vector_to_matrix(pk2, S_matrix);

  Core::LinAlg::Matrix<3, 3> FS;
  FS.multiply_nn(defgrd, S_matrix);

  Core::LinAlg::Matrix<3, 3> cauchy_matrix;
  cauchy_matrix.multiply_nt(1.0 / defgrd.determinant(), FS, defgrd, 0.0);

  Core::LinAlg::Voigt::Stresses::matrix_to_vector(cauchy_matrix, cauchy);
}

Core::LinAlg::Matrix<6, 1> Solid::UTILS::green_lagrange_to_euler_almansi(
    const Core::LinAlg::Matrix<6, 1>& gl, const Core::LinAlg::Matrix<3, 3>& defgrd)
{
  Core::LinAlg::Matrix<3, 3> invdefgrd(defgrd);
  invdefgrd.invert();

  Core::LinAlg::Matrix<3, 3> E_matrix;
  Core::LinAlg::Voigt::Strains::vector_to_matrix(gl, E_matrix);

  Core::LinAlg::Matrix<3, 3> iFTE;
  iFTE.multiply_tn(invdefgrd, E_matrix);

  Core::LinAlg::Matrix<3, 3> ea_matrix;
  ea_matrix.multiply_nn(iFTE, invdefgrd);

  Core::LinAlg::Matrix<6, 1> ea;
  Core::LinAlg::Voigt::Strains::matrix_to_vector(ea_matrix, ea);
  return ea;
}

Core::LinAlg::Matrix<6, 1> Solid::UTILS::green_lagrange_to_log_strain(
    const Core::LinAlg::Matrix<6, 1>& gl)
{
  Core::LinAlg::Matrix<3, 3> E_matrix;
  Core::LinAlg::Voigt::Strains::vector_to_matrix(gl, E_matrix);

  Core::LinAlg::Matrix<3, 3> pr_strain(true);  // squared principal strains
  Core::LinAlg::Matrix<3, 3> pr_dir(true);     // principal directions
  Core::LinAlg::SYEV(E_matrix, pr_strain, pr_dir);

  // compute principal logarithmic strains
  Core::LinAlg::Matrix<3, 3> pr_log_strain(true);
  for (int i = 0; i < 3; ++i) pr_log_strain(i, i) = std::log(std::sqrt(2 * pr_strain(i, i) + 1.0));

  // create logarithmic strain matrix
  Core::LinAlg::Matrix<3, 3> log_strain_matrix(true);
  Core::LinAlg::Matrix<3, 3> VH(false);
  VH.multiply_nn(pr_dir, pr_log_strain);
  log_strain_matrix.multiply_nt(VH, pr_dir);

  // convert to strain-like voigt notation
  Core::LinAlg::Matrix<6, 1> log_strain_voigt(true);
  Core::LinAlg::Voigt::Strains::matrix_to_vector(log_strain_matrix, log_strain_voigt);
  return log_strain_voigt;
}

int Solid::UTILS::ReadElement::read_element_material(Input::LineDefinition* linedef)
{
  int material = 0;
  linedef->extract_int("MAT", material);
  return material;
}

Inpar::Solid::KinemType Solid::UTILS::ReadElement::read_element_kinematic_type(
    Input::LineDefinition* linedef)
{
  std::string kinem;
  linedef->extract_string("KINEM", kinem);
  if (kinem == "nonlinear")
    return Inpar::Solid::KinemType::nonlinearTotLag;
  else if (kinem == "linear")
    return Inpar::Solid::KinemType::linear;
  else
  {
    FOUR_C_THROW("unknown kinematic type %s", kinem.c_str());
    return Inpar::Solid::KinemType::vague;
  }
}

Discret::ELEMENTS::ElementTechnology Solid::UTILS::ReadElement::read_element_technology(
    Input::LineDefinition* linedef)
{
  std::string type;
  linedef->extract_string("TECH", type);
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

Discret::ELEMENTS::PrestressTechnology Solid::UTILS::ReadElement::read_prestress_technology(
    Input::LineDefinition* linedef)
{
  std::string type;
  linedef->extract_string("PRESTRESS_TECH", type);
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

Discret::ELEMENTS::SolidElementProperties Solid::UTILS::ReadElement::read_solid_element_properties(
    Input::LineDefinition* linedef)
{
  Discret::ELEMENTS::SolidElementProperties solid_properties{};
  // element technology
  if (linedef->has_named("TECH"))
  {
    solid_properties.element_technology =
        Solid::UTILS::ReadElement::read_element_technology(linedef);
  }

  // prestress technology
  if (linedef->has_named("PRESTRESS_TECH"))
  {
    solid_properties.prestress_technology =
        Solid::UTILS::ReadElement::read_prestress_technology(linedef);
  }
  // kinematic type
  solid_properties.kintype = Solid::UTILS::ReadElement::read_element_kinematic_type(linedef);

  return solid_properties;
}

void Solid::UTILS::nodal_block_information_solid(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;

  nv = 3;
}

FOUR_C_NAMESPACE_CLOSE
