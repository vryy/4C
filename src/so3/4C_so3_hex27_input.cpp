/*----------------------------------------------------------------------*/
/*! \file
\brief tri-quadratic displacement based solid element
\level 1

*----------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_hex27.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::SoHex27::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  Teuchos::RCP<CORE::MAT::Material> mat = Material();

  SolidMaterial()->Setup(NUMGPT_SOH27, linedef);

  std::string buffer;
  linedef->ExtractString("KINEM", buffer);

  if (buffer == "linear")
  {
    kintype_ = INPAR::STR::KinemType::linear;
  }
  else if (buffer == "nonlinear")
  {
    kintype_ = INPAR::STR::KinemType::nonlinearTotLag;
  }
  else
    FOUR_C_THROW("Reading SO_HEX27 element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (SolidMaterial()->UsesExtendedUpdate())
    FOUR_C_THROW("This element currently does not support the extended update call.");

  return true;
}

FOUR_C_NAMESPACE_CLOSE
