/*!----------------------------------------------------------------------
\file beam2_input.cpp
\brief

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*----------------------------------------------------------------------*/

#include "beam2r.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam2r::ReadElement(const std::string& eletype,
                                        const std::string& distype,
                                        DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  //Set the applied Gaussrule ( It can be proven that we need 1 GP less than nodes to integrate exact )
  //note: we use a static cast for the enumeration here cf. Practical C++ Programming p.185
  gaussrule_ = static_cast<enum DRT::UTILS::GaussRule1D>(NumNode()-1);

  linedef->ExtractDouble("CROSS",crosssec_); 
  
  double shear_correction = 0.0;
  linedef->ExtractDouble("SHEARCORR",shear_correction);
  crosssecshear_ = crosssec_ * shear_correction;

  linedef->ExtractDouble("INERMOM",mominer_);
   
  return true;
}

