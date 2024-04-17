/*----------------------------------------------------------------------------*/
/*! \file
\brief three dimensional total Lagrange truss element

\level 3


*/
/*---------------------------------------------------------------------------*/

#include "baci_io_linedefinition.hpp"
#include "baci_truss3.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Truss3::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  linedef->ExtractDouble("CROSS", crosssec_);

  std::string buffer;
  linedef->ExtractString("KINEM", buffer);

  if (buffer == "totlag")  // geometrically non-linear with Total Lagrangean approach
    kintype_ = KinematicType::tr3_totlag;
  else if (buffer == "engstr")  // geometrically non-linear approach with engineering strains
    kintype_ = KinematicType::tr3_engstrain;
  else
    dserror("Reading of Truss3 element failed because of unknown kinematic type!");

  return true;
}

/*------------------------------------------------------------------------*
 | Set cross section area                           (public) mueller 03/12|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::SetCrossSec(const double& crosssec) { crosssec_ = crosssec; }

FOUR_C_NAMESPACE_CLOSE
