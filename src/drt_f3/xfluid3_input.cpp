/*!----------------------------------------------------------------------
\file xfluid3_input.cpp
\brief input stuff

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "xfluid3.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::XFluid3::ReadElement(const std::string& eletype,
                                         const std::string& distype,
                                         DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  std::string na;
  linedef->ExtractString("NA",na);
  if (na=="ale" or na=="ALE" or na=="Ale")
  {
    is_ale_ = true;
  }
  else if (na=="euler" or na=="EULER" or na=="Euler")
    is_ale_ = false;
  else
    dserror("Reading of XFLUID3 element failed: Euler/Ale");

  return true;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
