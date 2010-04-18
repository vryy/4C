/*!----------------------------------------------------------------------
\file smoothrod_input.cpp
\brief

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SMOOTHROD
#ifdef CCADISCRET

#include "smoothrod.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Smoothrod::ReadElement(const std::string& eletype,
                                       const std::string& distype,
                                       DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  linedef->ExtractDouble("CROSS",crosssec_);

  /*read beam moments of inertia of area; currently the beam3 element works only with rotationally symmetric
   * crosssection so that the moment of inertia of area around both principal can be expressed by one input
   * number I_; however, the implementation itself is a general one and works also for other cases; the only
   * point which has to be made sure is that the nodal triad T_ is initialized in the registration process
   * (->beam3.cpp) in such a way that t1 is the unit vector along the beam axis and t2 and t3 are the principal
   * axes with moment of inertia of area Iyy_ and Izz_, respectively; so a modification to more general kinds of
   * cross sections can be done easily by allowing for more complex input right here and by calculating an approxipate
   * initial nodal triad in the frame of the registration; */

  linedef->ExtractDouble("MOMIN",Iyy_);
  linedef->ExtractDouble("MOMIN",Izz_);
  linedef->ExtractDouble("MOMINPOL",Irr_);

  return true;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SMOOTHROD
