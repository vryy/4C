/*----------------------------------------------------------------------*/
/*! \file

\brief ReadElement method of the fluid element implementation


\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_fluid_ele.hpp"
#include "4C_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Fluid::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  // set discretization type (setOptimalgaussrule is pushed into element
  // routine)
  SetDisType(CORE::FE::StringToCellType(distype));

  std::string na;
  linedef->ExtractString("NA", na);
  if (na == "ale" or na == "ALE" or na == "Ale")
  {
    is_ale_ = true;
  }
  else if (na == "euler" or na == "EULER" or na == "Euler")
    is_ale_ = false;
  else
    FOUR_C_THROW("Reading of fluid element failed: Euler/Ale");

  return true;
}

FOUR_C_NAMESPACE_CLOSE
