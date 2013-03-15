/*!----------------------------------------------------------------------
\file so_weg6_input.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/

#include "so_weg6.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_mat/so3_material.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_weg6::ReadElement(const std::string& eletype,
                                         const std::string& distype,
                                         DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());
  so3mat->Setup(NUMGPT_WEG6, linedef);

  std::string buffer;
  linedef->ExtractString("KINEM",buffer);
  if (buffer=="linear")
  {
    kintype_ = sow6_linear;
    dserror("Reading of SO_WEG6 element failed only nonlinear kinematics implemented");
  }
  else if (buffer=="nonlinear")
  {
    kintype_ = sow6_nonlinear;
  }
  else dserror ("Reading SO_WEG6 element failed KINEM unknwon");

  return true;
}
