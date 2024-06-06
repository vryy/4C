/*----------------------------------------------------------------------*/
/*! \file
\brief 8-node solid shell element
\level 2
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */

#include "4C_io_linedefinition.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_sh8p8.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::SoSh8p8::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(0, Mat::Factory(material));

  SolidMaterial()->Setup(NUMGPT_SOH8, linedef);

  // a temprorary variable for read-in
  std::string buffer;
  // read kinematic flag
  linedef->ExtractString("KINEM", buffer);
  if (buffer == "linear")
  {
    FOUR_C_THROW("Only nonlinear kinematics for SO_SH8P8 implemented!");
  }
  else if (buffer == "nonlinear")
  {
    kintype_ = Inpar::STR::KinemType::nonlinearTotLag;
  }
  else
    FOUR_C_THROW("Reading SO_SH8P8 element failed unknown KINEM Type");

  // we expect kintype to be total lagrangian
  kintype_ = Inpar::STR::KinemType::nonlinearTotLag;

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (SolidMaterial()->UsesExtendedUpdate())
    FOUR_C_THROW("This element currently does not support the extended update call.");


  // read EAS technology flag
  linedef->ExtractString("EAS", buffer);

  if (buffer == "sosh8")
  {
    eastype_ = soh8_eassosh8;
    neas_ = NUMEAS_SOSH8_;
  }
  else if (buffer == "atype")
  {
    eastype_ = soh8_easa;
    neas_ = NUMEAS_A_;
  }
  else if (buffer == "None")
  {
    eastype_ = soh8_easnone;
    neas_ = 0;
  }
  else if (buffer == "none")
  {
    eastype_ = soh8_easnone;
    neas_ = 0;
  }
  else
    FOUR_C_THROW("Reading of SO_SH8P8 EAS type failed");

  if (eastype_ != soh8_easnone)
  {
    eas_init();
  }

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
    FOUR_C_THROW("Reading of SO_SH8P8 thickness direction failed");

  linedef->ExtractString("STAB", buffer);
  if (buffer == "Aff")
    stab_ = stab_affine;
  else if (buffer == "NonAff")
    stab_ = stab_nonaffine;
  else if (buffer == "SpatAff")
    stab_ = stab_spatialaffine;
  else if (buffer == "Spat")
    stab_ = stab_spatial;
  else if (buffer == "PureDisp")
    stab_ = stab_puredisp;
  else
    FOUR_C_THROW("Reading of SO_SH8P8 stabilisation failed");

  linedef->ExtractString("ANS", buffer);
  if (buffer == "Later")
    ans_ = ans_lateral;
  else if (buffer == "OnSpot")
    ans_ = ans_onspot;
  else if (buffer == "None")
    ans_ = ans_none;
  else
    FOUR_C_THROW("Reading of SO_SH8P8 ANS type failed");

  // Linearization
  linedef->ExtractString("LIN", buffer);
  if (buffer == "One")
    lin_ = lin_one;
  else if (buffer == "Half")
    lin_ = lin_half;
  else if (buffer == "Sixth")
    lin_ = lin_sixth;
  else
    FOUR_C_THROW("Reading of SO_SH8P8 LIN type failed");

  // Isochoric way
  linedef->ExtractString("ISO", buffer);
  if (buffer == "Mat")
    iso_ = iso_material;
  else if (buffer == "Enf")
    iso_ = iso_enforced;
  else
    FOUR_C_THROW("Reading of SO_SH8P8 ISO type failed");

  return true;
}

FOUR_C_NAMESPACE_CLOSE
