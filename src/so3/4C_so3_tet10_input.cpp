/*----------------------------------------------------------------------*/
/*! \file
\brief Solid Tet10 Element
\level 2
*----------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_tet10.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::SoTet10::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->extract_int("MAT", material);
  SetMaterial(0, Mat::Factory(material));

  SolidMaterial()->setup(NUMGPT_SOTET10, linedef);

  std::string buffer;
  linedef->extract_string("KINEM", buffer);

  // geometrically linear
  if (buffer == "linear")
  {
    kintype_ = Inpar::STR::KinemType::linear;
    FOUR_C_THROW("Reading of SO_TET10 element failed only nonlinear kinematics implemented");
  }
  // geometrically non-linear with Total Lagrangean approach
  else if (buffer == "nonlinear")
    kintype_ = Inpar::STR::KinemType::nonlinearTotLag;
  // geometrically non-linear with Updated Lagrangean approach
  else
    FOUR_C_THROW("Reading of SO_TET10 element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  return true;
}

FOUR_C_NAMESPACE_CLOSE
