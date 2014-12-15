/*!----------------------------------------------------------------------
\file so_hex8fbar_input.cpp
\brief

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/

#include "so_hex8fbar.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_mat/so3_material.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_hex8fbar::ReadElement(const std::string& eletype,
                                             const std::string& distype,
                                             DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  // set up of materials with GP data (e.g., history variables)
  Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());
  so3mat->Setup(NUMGPT_SOH8, linedef);

  // temporary variable for read-in
   std::string buffer;

  // read kinematic flag
  linedef->ExtractString("KINEM",buffer);
  if (buffer=="linear")
  {
    dserror("Only nonlinear kinematics for SO_HEX8FBAR implemented!");
  }
  else if (buffer=="nonlinear")
  {
    kintype_ = INPAR::STR::kinem_nonlinearTotLag;
  }
  else dserror ("Reading SO_HEX8FBAR element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  so3mat->ValidKinematics(kintype_);

  return true;
}

