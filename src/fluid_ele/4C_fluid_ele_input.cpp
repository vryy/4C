/*----------------------------------------------------------------------*/
/*! \file

\brief ReadElement method of the fluid element implementation


\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_fluid_ele.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_material_factory.hpp"

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::Fluid::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->extract_int("MAT", material);
  SetMaterial(0, Mat::Factory(material));

  // set discretization type (setOptimalgaussrule is pushed into element
  // routine)
  SetDisType(Core::FE::StringToCellType(distype));

  std::string na;
  linedef->extract_string("NA", na);
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
