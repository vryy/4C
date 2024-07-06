/*----------------------------------------------------------------------*/
/*! \file
\brief Solid Tet4 Element
\level 1
*----------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_tet4.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::SoTet4::read_element(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material_id = 0;
  linedef->extract_int("MAT", material_id);
  set_material(0, Mat::Factory(material_id));

  Teuchos::RCP<Core::Mat::Material> mat = material();

  solid_material()->setup(NUMGPT_SOTET4, linedef);

  std::string buffer;
  linedef->extract_string("KINEM", buffer);

  // geometrically linear
  if (buffer == "linear")
  {
    kintype_ = Inpar::Solid::KinemType::linear;
  }
  // geometrically non-linear with Total Lagrangean approach
  else if (buffer == "nonlinear")
  {
    kintype_ = Inpar::Solid::KinemType::nonlinearTotLag;
  }
  else
  {
    FOUR_C_THROW("Reading of SO_TET4 element failed KINEM unknown");
  }

  // check if material kinematics is compatible to element kinematics
  solid_material()->valid_kinematics(kintype_);

  return true;
}

FOUR_C_NAMESPACE_CLOSE
