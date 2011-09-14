/*!----------------------------------------------------------------------
\file fluid3_input.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "fluid3.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_linedefinition.H"

extern struct _GENPROB     genprob;

using namespace DRT::UTILS;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Fluid3::ReadElement(const std::string& eletype,
                                        const std::string& distype,
                                        DRT::INPUT::LineDefinition* linedef)
{
  // Only Fluid3 elements are allowed -> the support for Fluid2 elements is under dvelopment
  // if (genprob.ndim!=3)
  //   dserror("Problem defined as %dd, but Fluid3 does not support 2D elements yet",genprob.ndim);

  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  // Set Discretization Type
  // setOptimalgaussrule is pushed into the element routine
  SetDisType(DRT::StringToDistype(distype));

  std::string na;
  linedef->ExtractString("NA",na);
  if (na=="ale" or na=="ALE" or na=="Ale")
  {
    is_ale_ = true;
  }
  else if (na=="euler" or na=="EULER" or na=="Euler")
    is_ale_ = false;
  else
    dserror("Reading of FLUID3 element failed: Euler/Ale");

  return true;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
