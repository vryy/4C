#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_thermo_element.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | read element and set required information                  gjb 01/08 |
 *----------------------------------------------------------------------*/
bool Thermo::Element::read_element(const std::string& eletype, const std::string& distype,
    const Core::IO::InputParameterContainer& container)
{
  // read number of material model
  int material_id = container.get<int>("MAT");
  set_material(0, Mat::factory(material_id));

  set_dis_type(Core::FE::string_to_cell_type(distype));

  return true;
}

FOUR_C_NAMESPACE_CLOSE
