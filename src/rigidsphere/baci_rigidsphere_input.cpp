/*----------------------------------------------------------------------------*/
/*! \file

\brief spherical particle element for brownian dynamics

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "baci_io_linedefinition.H"
#include "baci_rigidsphere.H"

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Rigidsphere::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // currently only rotationally symmetric profiles for beam --> Iyy = Izz
  linedef->ExtractDouble("RADIUS", radius_);
  linedef->ExtractDouble("DENSITY", rho_);

  return (true);
}

BACI_NAMESPACE_CLOSE
