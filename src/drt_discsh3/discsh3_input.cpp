/*----------------------------------------------------------------------*/
/*! \file
\brief input routines for discsh3 element

\level 3

*/
/*----------------------------------------------------------------------*/
#include "discsh3.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::DiscSh3::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  //  int thickness = 0;
  linedef->ExtractDouble("THICK", thickness_);

  return true;
}
