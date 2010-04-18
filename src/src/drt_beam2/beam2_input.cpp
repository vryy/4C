/*!----------------------------------------------------------------------
\file beam2_input.cpp
\brief

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BEAM2
#ifdef CCADISCRET

#include "beam2.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |  read element input (public)                              cyron 01/08|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam2::ReadElement(const std::string& eletype,
                                       const std::string& distype,
                                       DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);
  
  linedef->ExtractDouble("CROSS",crosssec_);

  double shear_correction = 0.0;
  linedef->ExtractDouble("SHEARCORR",shear_correction);
  crosssecshear_ = crosssec_ * shear_correction;
  
  linedef->ExtractDouble("INERMOM",mominer_);
   
  return true;
} // Beam2::ReadElement()


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM2
