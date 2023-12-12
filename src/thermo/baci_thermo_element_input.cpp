/*----------------------------------------------------------------------*/
/*! \file
\brief element input routines
\level 1
*/

/*----------------------------------------------------------------------*
 | headers                                                    gjb 01/08 |
 *----------------------------------------------------------------------*/
#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_io_linedefinition.H"
#include "baci_thermo_element.H"

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | read element and set required information                  gjb 01/08 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Thermo::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  SetDisType(CORE::FE::StringToCellType(distype));

  if (Shape() == CORE::FE::CellType::nurbs27) SetNurbsElement() = true;

  return true;
}

BACI_NAMESPACE_CLOSE
