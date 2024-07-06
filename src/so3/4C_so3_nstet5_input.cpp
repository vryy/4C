/*----------------------------------------------------------------------*/
/*! \file

\brief NStet5 element

\level 2


*----------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_mat_elasthyper.hpp"
#include "4C_so3_nstet5.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::NStet5::read_element(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material_id = 0;
  linedef->extract_int("MAT", material_id);
  set_material(0, Mat::Factory(material_id));

  if (material()->material_type() == Core::Materials::m_elasthyper)
  {
    Mat::ElastHyper* elahy = dynamic_cast<Mat::ElastHyper*>(material().get());
    elahy->setup(0, linedef);
  }

  std::string buffer;
  linedef->extract_string("KINEM", buffer);
  if (buffer == "linear")
  {
    // kintype_ not yet implemented for nstet5
    // kintype_ = sonstet5_linear;
    FOUR_C_THROW("Reading of SO_NSTET5 element failed only nonlinear kinematics implemented");
  }
  else if (buffer == "nonlinear")
  {
    // kintype_ not yet implemented for nstet5
    // kintype_ = sonstet5_nonlinear;
  }
  else
    FOUR_C_THROW("Reading SO_NSTET5 element failed KINEM unknown");

  return true;
}

FOUR_C_NAMESPACE_CLOSE
