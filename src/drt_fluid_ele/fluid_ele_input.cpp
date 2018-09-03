/*!----------------------------------------------------------------------
\file fluid_ele_input.cpp
\brief

<pre>
Maintainer: Ursula Rasthofer & Volker Gravemeier
            {rasthofer,vgravem}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236/-245
</pre>

*----------------------------------------------------------------------*/

#include "fluid_ele.H"
#include "../drt_lib/drt_linedefinition.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Fluid::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  // set discretization type (setOptimalgaussrule is pushed into element
  // routine)
  SetDisType(DRT::StringToDistype(distype));

  std::string na;
  linedef->ExtractString("NA", na);
  if (na == "ale" or na == "ALE" or na == "Ale")
  {
    is_ale_ = true;
  }
  else if (na == "euler" or na == "EULER" or na == "Euler")
    is_ale_ = false;
  else
    dserror("Reading of fluid element failed: Euler/Ale");

  return true;
}
