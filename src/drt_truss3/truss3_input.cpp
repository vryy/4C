/*----------------------------------------------------------------------------*/
/*! \file
\brief three dimensional total Lagrange truss element

\level 3


*/
/*---------------------------------------------------------------------------*/

#include "truss3.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Truss3::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  linedef->ExtractDouble("CROSS", crosssec_);

  std::string buffer;
  linedef->ExtractString("KINEM", buffer);

  // geometrically non-linear with Total Lagrangean approach
  if (buffer == "totlag") kintype_ = tr3_totlag;

  // geometrically non-linear approach with engineering strains
  else if (buffer == "engstr")
    kintype_ = tr3_engstrain;

  else
    dserror("Reading of Truss3 element failed because of unknown kinematic type!");

  return true;
}

/*------------------------------------------------------------------------*
 | Set cross section area                           (public) mueller 03/12|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::SetCrossSec(const double& crosssec)
{
  crosssec_ = crosssec;
  return;
}
