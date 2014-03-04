/*!----------------------------------------------------------------------
\file rigidsphere_input.cpp
\brief

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*----------------------------------------------------------------------*/

#include "rigidsphere.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Rigidsphere::ReadElement(const std::string& eletype,
                                                  const std::string& distype,
                                                  DRT::INPUT::LineDefinition* linedef)
{
 //currently only rotationally symmetric profiles for beam --> Iyy = Izz
 linedef->ExtractDouble("RADIUS",radius_);
 linedef->ExtractDouble("DENSITY",rho_);

  return (true);
}
