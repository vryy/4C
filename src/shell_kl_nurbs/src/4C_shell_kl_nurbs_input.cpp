/*----------------------------------------------------------------------------*/
/*! \file
\brief Input for NURBS Kirchhoff-Love shell implementation

\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_shell_kl_nurbs.hpp"


FOUR_C_NAMESPACE_OPEN


/**
 *
 */
bool Discret::ELEMENTS::KirchhoffLoveShellNurbs::read_element(const std::string& eletype,
    const std::string& distype, const Core::IO::InputParameterContainer& container)
{
  // set discretization type
  set_dis_type(Core::FE::string_to_cell_type(distype));

  auto ngp = container.get<std::vector<int>>("GP");
  if (ngp[0] < 1 or ngp[1] < 1) FOUR_C_THROW("Number of Gauss points has to be a positive integer");

  for (unsigned int dim = 0; dim < 2; dim++)
    gaussrule_[dim] = Core::FE::num_gauss_points_to_gauss_rule<Core::FE::CellType::line2>(ngp[dim]);

  // read number of material model
  auto material = container.get<int>("MAT");
  set_material(0, Mat::factory(material));
  return true;
}

FOUR_C_NAMESPACE_CLOSE
