/*----------------------------------------------------------------------*/
/*! \file
\brief

\level 3


\brief Nonlinear Membrane Finite Element input

*----------------------------------------------------------------------*/
#include "4C_io_linedefinition.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_membrane.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  ReadElement                                            fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
bool Discret::ELEMENTS::Membrane<distype>::ReadElement(
    const std::string& eletype, const std::string& eledistype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->extract_int("MAT", material);
  SetMaterial(0, Mat::Factory(material));

  // set up of materials with GP data (e.g., history variables)
  SolidMaterial()->setup(intpoints_.nquad, linedef);


  // read element thickness
  linedef->extract_double("THICK", thickness_);
  if (thickness_ <= 0) FOUR_C_THROW("Membrane element thickness needs to be > 0");

  // initialize current thickness at all gp
  for (int i = 0; i < intpoints_.nquad; ++i) cur_thickness_[i] = thickness_;


  // temporary variable for read-in
  std::string buffer;


  // reduced dimension assumption
  linedef->extract_string("STRESS_STRAIN", buffer);
  if (buffer == "plane_stress")
  {
    planetype_ = plane_stress;
  }
  else if (buffer == "plane_strain")
  {
    FOUR_C_THROW("Membrane not intended for plane strain evaluation");
  }
  else
    FOUR_C_THROW("Reading STRESS_STRAIN state failed");

  return true;
}

template class Discret::ELEMENTS::Membrane<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::Membrane<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::Membrane<Core::FE::CellType::quad4>;
template class Discret::ELEMENTS::Membrane<Core::FE::CellType::quad9>;

FOUR_C_NAMESPACE_CLOSE
