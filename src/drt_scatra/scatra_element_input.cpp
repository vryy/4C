/*----------------------------------------------------------------------*/
/*!
\file scatra_element_input.cpp
\brief

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#if defined(D_FLUID2) || defined(D_FLUID3)
#ifdef CCADISCRET

#include "scatra_element.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*
 | read element input (public)                                          |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Transport::ReadElement(const std::string& eletype,
                                           const std::string& distype,
                                           DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  SetDisType(DRT::StringToDistype(distype));

  return true;
}


#endif  // #ifdef CCADISCRET
#endif  // #if defined(D_FLUID2) || defined(D_FLUID3)
