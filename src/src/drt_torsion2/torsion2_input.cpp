/*!----------------------------------------------------------------------
\file torsion2_input.cpp
\brief two dimensional torsion spring element

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_TORSION2
#ifdef CCADISCRET

#include "torsion2.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/linalg_fixedsizematrix.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Torsion2::ReadElement(const std::string& eletype,
                                          const std::string& distype,
                                          DRT::INPUT::LineDefinition* linedef)
{
  //read constant of torsion spring
  linedef->ExtractDouble("SPRING",springconstant_);
  return true;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_TORSION2
