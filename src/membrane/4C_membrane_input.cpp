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
 |  read_element                                            fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
bool Discret::ELEMENTS::Membrane<distype>::read_element(const std::string& eletype,
    const std::string& eledistype, const Core::IO::InputParameterContainer& container)
{
  // read number of material model
  int material_id = container.get<int>("MAT");
  set_material(0, Mat::Factory(material_id));

  // set up of materials with GP data (e.g., history variables)
  solid_material()->setup(intpoints_.nquad, container);

  // read element thickness
  thickness_ = container.get<double>("THICK");
  if (thickness_ <= 0) FOUR_C_THROW("Membrane element thickness needs to be > 0");

  // initialize current thickness at all gp
  for (int i = 0; i < intpoints_.nquad; ++i) cur_thickness_[i] = thickness_;

  // reduced dimension assumption
  std::string buffer = container.get<std::string>("STRESS_STRAIN");
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
