/*!----------------------------------------------------------------------
\file truss3_input.cpp
\brief three dimensional total Lagrange truss element

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_TORSION3
#ifdef CCADISCRET

#include "torsion3.H"
#include "../drt_lib/drt_linedefinition.H"



/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Torsion3::ReadElement(const std::string& eletype,
                                          const std::string& distype,
                                          DRT::INPUT::LineDefinition* linedef)
{
  //read constant of torsion spring
  linedef->ExtractDouble("SPRING",springconstant_);
  return true;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_TORSION3
