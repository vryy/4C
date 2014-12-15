/*!----------------------------------------------------------------------**###
\file so_tet4_input.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
writen by : Alexander Volf
            alexander.volf@mytum.de
</pre>

*----------------------------------------------------------------------*/

#include "so_tet4.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_mat/so3_material.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_tet4::ReadElement(const std::string& eletype,
                                         const std::string& distype,
                                         DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  Teuchos::RCP<MAT::Material> mat = Material();

  Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());
  so3mat->Setup(NUMGPT_SOTET4, linedef);

  std::string buffer;
  linedef->ExtractString("KINEM",buffer);

  // geometrically linear
  if(buffer=="linear")
  {
    kintype_ = INPAR::STR::kinem_linear;
    dserror("Reading of SO_TET4 element failed only nonlinear kinematics implemented");
  }
  // geometrically non-linear with Total Lagrangean approach
  else if (buffer=="nonlinear")    kintype_ = INPAR::STR::kinem_nonlinearTotLag;
  else dserror("Reading of SO_TET4 element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  so3mat->ValidKinematics(kintype_);

  return true;
}


