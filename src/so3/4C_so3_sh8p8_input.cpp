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
bool Discret::ELEMENTS::SoSh8p8::read_element(const std::string& eletype,
    const std::string& distype, const Core::IO::InputParameterContainer& container)
{
  // read number of material model
  int material_id = container.get<int>("MAT");
  set_material(0, Mat::Factory(material_id));

  solid_material()->setup(NUMGPT_SOH8, container);


  // read kinematic flag
  std::string kinem = container.get<std::string>("KINEM");

  if (kinem == "linear")
  {
    FOUR_C_THROW("Only nonlinear kinematics for SO_SH8P8 implemented!");
  }
  else if (kinem == "nonlinear")
  {
    kintype_ = Inpar::Solid::KinemType::nonlinearTotLag;
  }
  else
    FOUR_C_THROW("Reading SO_SH8P8 element failed unknown KINEM Type");

  // we expect kintype to be total lagrangian
  kintype_ = Inpar::Solid::KinemType::nonlinearTotLag;

  // check if material kinematics is compatible to element kinematics
  solid_material()->valid_kinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (solid_material()->uses_extended_update())
    FOUR_C_THROW("This element currently does not support the extended update call.");


  // read EAS technology flag
  std::string eas = container.get<std::string>("EAS");

  if (eas == "sosh8")
  {
    eastype_ = soh8_eassosh8;
    neas_ = NUMEAS_SOSH8_;
  }
  else if (eas == "atype")
  {
    eastype_ = soh8_easa;
    neas_ = NUMEAS_A_;
  }
  else if (eas == "None")
  {
    eastype_ = soh8_easnone;
    neas_ = 0;
  }
  else if (eas == "none")
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


  // read THICKDIR flag
  std::string thickdir = container.get<std::string>("THICKDIR");
  nodes_rearranged_ = false;

  // global X
  if (thickdir == "xdir") thickdir_ = globx;
  // global Y
  else if (thickdir == "ydir")
    thickdir_ = globy;
  // global Z
  else if (thickdir == "zdir")
    thickdir_ = globz;
  // find automatically through Jacobian of Xrefe
  else if (thickdir == "auto")
    thickdir_ = autoj;
  // local r
  else if (thickdir == "rdir")
    thickdir_ = enfor;
  // local s
  else if (thickdir == "sdir")
    thickdir_ = enfos;
  // local t
  else if (thickdir == "tdir")
    thickdir_ = enfot;
  // no noderearrangement
  else if (thickdir == "none")
  {
    thickdir_ = none;
    nodes_rearranged_ = true;
  }
  else
    FOUR_C_THROW("Reading of SO_SH8P8 thickness direction failed");

  // read STAB flag
  auto stab = container.get<std::string>("STAB");

  if (stab == "Aff")
    stab_ = stab_affine;
  else if (stab == "NonAff")
    stab_ = stab_nonaffine;
  else if (stab == "SpatAff")
    stab_ = stab_spatialaffine;
  else if (stab == "Spat")
    stab_ = stab_spatial;
  else if (stab == "PureDisp")
    stab_ = stab_puredisp;
  else
    FOUR_C_THROW("Reading of SO_SH8P8 stabilisation failed");


  // read ANS flag
  std::string ans = container.get<std::string>("ANS");

  if (ans == "Later")
    ans_ = ans_lateral;
  else if (ans == "OnSpot")
    ans_ = ans_onspot;
  else if (ans == "None")
    ans_ = ans_none;
  else
    FOUR_C_THROW("Reading of SO_SH8P8 ANS type failed");

  // read linearization flag
  std::string lin = container.get<std::string>("LIN");

  if (lin == "One")
    lin_ = lin_one;
  else if (lin == "Half")
    lin_ = lin_half;
  else if (lin == "Sixth")
    lin_ = lin_sixth;
  else
    FOUR_C_THROW("Reading of SO_SH8P8 LIN type failed");

  // read isochoric flag
  std::string iso = container.get<std::string>("ISO");

  if (iso == "Mat")
    iso_ = iso_material;
  else if (iso == "Enf")
    iso_ = iso_enforced;
  else
    FOUR_C_THROW("Reading of SO_SH8P8 ISO type failed");

  return true;
}

FOUR_C_NAMESPACE_CLOSE
