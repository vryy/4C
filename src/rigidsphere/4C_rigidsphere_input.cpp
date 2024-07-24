/*----------------------------------------------------------------------------*/
/*! \file

\brief spherical particle element for brownian dynamics

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_rigidsphere.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::Rigidsphere::read_element(const std::string& eletype,
    const std::string& distype, const Core::IO::InputParameterContainer& container)
{
  // currently only rotationally symmetric profiles for beam --> Iyy = Izz
  radius_ = container.get<double>("RADIUS");
  rho_ = container.get<double>("DENSITY");

  return (true);
}

FOUR_C_NAMESPACE_CLOSE
