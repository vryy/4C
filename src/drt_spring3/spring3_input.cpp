/*!----------------------------------------------------------------------
\file spring3_input.cpp
\brief three dimensional spring element

<pre>
Maintainer: Dhrubajyoti Mukherjee
            mukherjee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15270
</pre>

*----------------------------------------------------------------------*/

#include "spring3.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Spring3::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  // Gruesse aus Barcelona, Martin und Dhruba
  linedef->ExtractDouble("CROSS", crosssec_);


  // set nodal tridas according to input file
  Qnew_.resize(NumNode());
  Qold_.resize(NumNode());

  Qold_ = Qnew_;

  return true;
}
/*------------------------------------------------------------------------*
 | Set cross section area                         (public) mukherjee 04/15|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Spring3::SetCrossSec(const double& crosssec)
{
  crosssec_ = crosssec;
  return;
}
