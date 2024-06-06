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
bool Discret::ELEMENTS::Rigidsphere::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // currently only rotationally symmetric profiles for beam --> Iyy = Izz
  linedef->ExtractDouble("RADIUS", radius_);
  linedef->ExtractDouble("DENSITY", rho_);

  return (true);
}

FOUR_C_NAMESPACE_CLOSE
