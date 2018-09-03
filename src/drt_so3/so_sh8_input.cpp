/*!----------------------------------------------------------------------
\file so_sh8_input.cpp
\brief

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/


#include "so_sh8.H"
#include "../drt_mat/so3_material.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_sh8::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  SolidMaterial()->Setup(NUMGPT_SOH8, linedef);

  // temporary variable for read-in
  std::string buffer;


  // read kinematic flag
  linedef->ExtractString("KINEM", buffer);
  if (buffer == "linear")
  {
    // kintype_ = soh8_linear;
    dserror("Only nonlinear kinematics for SO_SH8 implemented!");
  }
  else if (buffer == "nonlinear")
  {
    kintype_ = INPAR::STR::kinem_nonlinearTotLag;
  }
  else
    dserror("Reading SO_HEX8p1j1 element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  // read EAS technology flag
  linedef->ExtractString("EAS", buffer);

  // full EAS technology
  if (buffer == "sosh8")
  {
    eastype_ = soh8_eassosh8;
    neas_ = 7;  // number of eas parameters for EAS_SOSH8
    soh8_easinit();
  }
  // no EAS technology
  else if (buffer == "none")
  {
    eastype_ = soh8_easnone;
    neas_ = 0;  // number of eas parameters for EAS_SOSH8
  }
  else
    dserror("Reading of SO_SH8 EAS technology failed");

  // read ANS technology flag
  linedef->ExtractString("ANS", buffer);
  if (buffer == "sosh8")
  {
    anstype_ = anssosh8;
  }
  // no ANS technology
  else if (buffer == "none")
  {
    anstype_ = ansnone;
  }
  else
    dserror("Reading of SO_SH8 ANS technology failed");

  linedef->ExtractString("THICKDIR", buffer);
  nodes_rearranged_ = false;

  // global X
  if (buffer == "xdir") thickdir_ = globx;
  // global Y
  else if (buffer == "ydir")
    thickdir_ = globy;
  // global Z
  else if (buffer == "zdir")
    thickdir_ = globz;
  // find automatically through Jacobian of Xrefe
  else if (buffer == "auto")
    thickdir_ = autoj;
  // local r
  else if (buffer == "rdir")
    thickdir_ = enfor;
  // local s
  else if (buffer == "sdir")
    thickdir_ = enfos;
  // local t
  else if (buffer == "tdir")
    thickdir_ = enfot;
  // no noderearrangement
  else if (buffer == "none")
  {
    thickdir_ = none;
    nodes_rearranged_ = true;
  }
  else
    dserror("Reading of SO_SH8 thickness direction failed");

  return true;
}
