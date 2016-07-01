/*!----------------------------------------------------------------------
\file membrane_input.cpp
\brief

\level 3

<pre>
\maintainer Fabian Bräu
            braeu@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

\brief Nonlinear Membrane Finite Element input

*----------------------------------------------------------------------*/
#include "membrane.H"
#include "../drt_mat/so3_material.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*
 |  ReadElement                                            fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::Membrane<distype>::ReadElement(const std::string& eletype,
                                                   const std::string& eledistype,
                                                   DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  // set up of materials with GP data (e.g., history variables)
  SolidMaterial()->Setup(intpoints_.nquad,linedef);


  // read element thickness
  linedef->ExtractDouble("THICK",thickness_);
  if  (thickness_<=0) dserror("Membrane element thickness needs to be > 0");

  // initialize current thickness at all gp
  for(int i=0;i<intpoints_.nquad;++i)
    curr_thickness_[i] = thickness_;


  // temporary variable for read-in
   std::string buffer;


  // reduced dimension assumption
  linedef->ExtractString("STRESS_STRAIN",buffer);
  if (buffer=="plane_stress")
  {
    planetype_ = plane_stress;
  }
  else if (buffer=="plane_strain")
  {
    dserror("Membrane not intended for plane strain evaluation");
  }
  else dserror("Reading STRESS_STRAIN state failed");

  return true;
}

template class DRT::ELEMENTS::Membrane<DRT::Element::tri3>;
template class DRT::ELEMENTS::Membrane<DRT::Element::tri6>;
template class DRT::ELEMENTS::Membrane<DRT::Element::quad4>;
template class DRT::ELEMENTS::Membrane<DRT::Element::quad9>;
