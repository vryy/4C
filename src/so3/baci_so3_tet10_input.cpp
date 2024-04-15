/*----------------------------------------------------------------------*/
/*! \file
\brief Solid Tet10 Element
\level 2
*----------------------------------------------------------------------*/

#include "baci_io_linedefinition.hpp"
#include "baci_mat_so3_material.hpp"
#include "baci_so3_tet10.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::SoTet10::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  SolidMaterial()->Setup(NUMGPT_SOTET10, linedef);

  std::string buffer;
  linedef->ExtractString("KINEM", buffer);

  // geometrically linear
  if (buffer == "linear")
  {
    kintype_ = INPAR::STR::KinemType::linear;
    dserror("Reading of SO_TET10 element failed only nonlinear kinematics implemented");
  }
  // geometrically non-linear with Total Lagrangean approach
  else if (buffer == "nonlinear")
    kintype_ = INPAR::STR::KinemType::nonlinearTotLag;
  // geometrically non-linear with Updated Lagrangean approach
  else
    dserror("Reading of SO_TET10 element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  return true;
}

FOUR_C_NAMESPACE_CLOSE
