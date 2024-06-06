/*----------------------------------------------------------------------*/
/*! \file
\brief element input routines
\level 1
*/

/*----------------------------------------------------------------------*
 | headers                                                    gjb 01/08 |
 *----------------------------------------------------------------------*/
#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_thermo_element.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | read element and set required information                  gjb 01/08 |
 *----------------------------------------------------------------------*/
bool Discret::ELEMENTS::Thermo::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(0, Mat::Factory(material));

  SetDisType(Core::FE::StringToCellType(distype));

  if (Shape() == Core::FE::CellType::nurbs27) SetNurbsElement() = true;

  return true;
}

FOUR_C_NAMESPACE_CLOSE
