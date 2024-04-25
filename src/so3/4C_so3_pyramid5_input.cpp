/*----------------------------------------------------------------------*/
/*! \file
\brief pyramid shaped solid element
\level 2


*----------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_pyramid5.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::SoPyramid5::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  Teuchos::RCP<MAT::Material> mat = Material();

  SolidMaterial()->Setup(NUMGPT_SOP5, linedef);

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
    FOUR_C_THROW("Reading SO_PYRAMID5 element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (SolidMaterial()->UsesExtendedUpdate())
    FOUR_C_THROW("This element currently does not support the extended update call.");

  // only for linear SVK materials and small strain plastic materials
  bool admissibl_mat = false;
  if ((mat->MaterialType() == INPAR::MAT::m_stvenant) or
      (mat->MaterialType() == INPAR::MAT::m_thermostvenant) or
      (mat->MaterialType() == INPAR::MAT::m_pllinelast) or
      (mat->MaterialType() == INPAR::MAT::m_thermopllinelast) or
      (mat->MaterialType() == INPAR::MAT::m_elpldamage))
    admissibl_mat = true;

  // check for SVK material if geometrically linear
  if ((kintype_ == INPAR::STR::KinemType::linear) and (admissibl_mat == false))
    FOUR_C_THROW("ERROR: Only linear elasticity (SVK) for geometrically linear pyramid5 element");

  return true;
}

FOUR_C_NAMESPACE_CLOSE
