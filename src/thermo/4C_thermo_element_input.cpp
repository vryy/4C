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
bool DRT::ELEMENTS::Thermo::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(0, MAT::Factory(material));

  SetDisType(CORE::FE::StringToCellType(distype));

  if (Shape() == CORE::FE::CellType::nurbs27) SetNurbsElement() = true;

  return true;
}

FOUR_C_NAMESPACE_CLOSE
