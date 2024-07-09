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
bool Discret::ELEMENTS::SoPyramid5::read_element(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material_id = 0;
  linedef->extract_int("MAT", material_id);
  set_material(0, Mat::Factory(material_id));

  Teuchos::RCP<Core::Mat::Material> mat = material();

  solid_material()->setup(NUMGPT_SOP5, linedef);

  std::string buffer;
  linedef->extract_string("KINEM", buffer);

  if (buffer == "linear")
  {
    kintype_ = Inpar::Solid::KinemType::linear;
  }
  else if (buffer == "nonlinear")
  {
    kintype_ = Inpar::Solid::KinemType::nonlinearTotLag;
  }
  else
    FOUR_C_THROW("Reading SO_PYRAMID5 element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  solid_material()->valid_kinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (solid_material()->uses_extended_update())
    FOUR_C_THROW("This element currently does not support the extended update call.");

  // only for linear SVK materials and small strain plastic materials
  bool admissibl_mat = false;
  if ((mat->material_type() == Core::Materials::m_stvenant) or
      (mat->material_type() == Core::Materials::m_thermostvenant) or
      (mat->material_type() == Core::Materials::m_pllinelast) or
      (mat->material_type() == Core::Materials::m_thermopllinelast) or
      (mat->material_type() == Core::Materials::m_elpldamage))
    admissibl_mat = true;

  // check for SVK material if geometrically linear
  if ((kintype_ == Inpar::Solid::KinemType::linear) and (admissibl_mat == false))
    FOUR_C_THROW("ERROR: Only linear elasticity (SVK) for geometrically linear pyramid5 element");

  return true;
}

FOUR_C_NAMESPACE_CLOSE
