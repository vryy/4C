/*!----------------------------------------------------------------------
\file truss3cl_input.cpp
 \brief three dimensional interpolated total Lagrange hybrid beam-truss element
 (can be connected to beam3 elements)

<pre>
Maintainer: Dhrubajyoti Mukherjee
            mukherjee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15270
</pre>

*----------------------------------------------------------------------*/

#include "truss3cl.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Truss3CL::ReadElement(
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
    dserror("Reading of Torsion2 element failed because of unknown kinematic type!");

  linedef->ExtractDoubleVector("BPOS", mybindingposition_);

  return true;
}

/*------------------------------------------------------------------------*
 | Set cross section area                         (public) mukherjee 01/14|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3CL::SetCrossSec(const double& crosssec)
{
  crosssec_ = crosssec;
  return;
}

/*------------------------------------------------------------------------*
 | Set internodal Binding Position                (public) mukherjee 01/14|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3CL::SetBindingPosition(const double& pos1, const double& pos2)
{
  mybindingposition_.resize(2);
  mybindingposition_[0] = pos1;
  mybindingposition_[1] = pos2;

  return;
}
