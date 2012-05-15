/*!----------------------------------------------------------------------
\file truss2_input.cpp
\brief three dimensional total Lagrange truss element

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/

#include "truss2.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Truss2::ReadElement(const std::string&          eletype,
                                        const std::string&          distype,
                                        DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  linedef->ExtractDouble("CROSS",crosssec_);

  std::string buffer;
  linedef->ExtractString("KINEM",buffer);

  // geometrically non-linear with Total Lagrangean approach
  if (buffer=="totlag")
    kintype_ = tr2_totlag;

  // geometrically non-linear approach with engineering strains
  else if (buffer=="engstr")
    kintype_ = tr2_engstrain;

  else
    dserror("Reading of Torsion2 element failed because of unknown kinematic type!");

  return true;
}
