/*----------------------------------------------------------------------*/
/*! \file
\brief
\level 1


*----------------------------------------------------------------------*/


#include "4C_so3_shw6.hpp"  //**
#include "4C_mat_so3_material.hpp"
#include "4C_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::SoShw6::read_element(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material_id = 0;
  linedef->extract_int("MAT", material_id);
  set_material(0, Mat::Factory(material_id));

  solid_material()->setup(NUMGPT_WEG6, linedef);

  std::string buffer;
  linedef->extract_string("KINEM", buffer);


  // geometrically non-linear with Total Lagrangean approach
  if (buffer == "nonlinear")
  {
    kintype_ = Inpar::Solid::KinemType::nonlinearTotLag;
  }
  // geometrically linear
  else if (buffer == "linear")
  {
    kintype_ = Inpar::Solid::KinemType::linear;
    FOUR_C_THROW("Reading of SOLIDSHW6 element failed onlz nonlinear kinetmatics implemented");
  }

  // geometrically non-linear with Updated Lagrangean approach
  else
    FOUR_C_THROW("Reading of SOLIDSHW6 element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  solid_material()->valid_kinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (solid_material()->uses_extended_update())
    FOUR_C_THROW("This element currently does not support the extended update call.");

  linedef->extract_string("EAS", buffer);

  // full sohw6 EAS technology
  if (buffer == "soshw6")
  {
    eastype_ = soshw6_easpoisthick;  // EAS to allow linear thickness strain
    neas_ = soshw6_easpoisthick;     // number of eas parameters
    soshw6_easinit();
  }
  // no EAS technology
  else if (buffer == "none")
  {
    eastype_ = soshw6_easnone;
    neas_ = 0;  // number of eas parameters
    std::cout << "Warning: Solid-Shell Wegde6 without EAS" << std::endl;
  }
  else
    FOUR_C_THROW("Reading of SOLIDSHW6 EAS technology failed");

  // check for automatically align material space optimally with parameter space
  optimal_parameterspace_map_ = false;
  nodes_rearranged_ = false;
  if (linedef->has_named("OPTORDER")) optimal_parameterspace_map_ = true;

  return true;
}

FOUR_C_NAMESPACE_CLOSE
