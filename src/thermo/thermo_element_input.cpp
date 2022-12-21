/*----------------------------------------------------------------------*/
/*! \file
\brief element input routines
\level 1
*/

/*----------------------------------------------------------------------*
 | headers                                                    gjb 01/08 |
 *----------------------------------------------------------------------*/
#include "thermo_element.H"
#include "linedefinition.H"
#include "utils_local_connectivity_matrices.H"


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

  SetDisType(DRT::StringToDistype(distype));

  if (Shape() == DRT::Element::nurbs27) SetNurbsElement() = true;

  return true;
}
