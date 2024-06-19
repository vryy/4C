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
bool Discret::ELEMENTS::SoPyramid5::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->extract_int("MAT", material);
  SetMaterial(0, Mat::Factory(material));

  Teuchos::RCP<Core::Mat::Material> mat = Material();

  SolidMaterial()->setup(NUMGPT_SOP5, linedef);

  std::string buffer;
  linedef->extract_string("KINEM", buffer);

  if (buffer == "linear")
  {
    kintype_ = Inpar::STR::KinemType::linear;
  }
  else if (buffer == "nonlinear")
  {
    kintype_ = Inpar::STR::KinemType::nonlinearTotLag;
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
  if ((mat->MaterialType() == Core::Materials::m_stvenant) or
      (mat->MaterialType() == Core::Materials::m_thermostvenant) or
      (mat->MaterialType() == Core::Materials::m_pllinelast) or
      (mat->MaterialType() == Core::Materials::m_thermopllinelast) or
      (mat->MaterialType() == Core::Materials::m_elpldamage))
    admissibl_mat = true;

  // check for SVK material if geometrically linear
  if ((kintype_ == Inpar::STR::KinemType::linear) and (admissibl_mat == false))
    FOUR_C_THROW("ERROR: Only linear elasticity (SVK) for geometrically linear pyramid5 element");

  return true;
}

FOUR_C_NAMESPACE_CLOSE
