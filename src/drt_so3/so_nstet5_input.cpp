/*!----------------------------------------------------------------------**###
\file so_nstet5_input.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_nstet5.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::NStet5::ReadElement(const std::string& eletype,
                                      const std::string& distype,
                                      DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);
  std::string buffer;
  linedef->ExtractString("KINEM",buffer);
  if (buffer=="linear")
  {
    // kintype_ not yet implemented for nstet5
    //kintype_ = sonstet5_linear;
    dserror("Reading of SO_NSTET5 element failed only nonlinear kinematics implemented");
  }
  else if (buffer=="nonlinear")
  {
    // kintype_ not yet implemented for nstet5
    //kintype_ = sonstet5_nonlinear;
  }
  else dserror ("Reading SO_NSTET5 element failed KINEM unknown");

  return true;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
