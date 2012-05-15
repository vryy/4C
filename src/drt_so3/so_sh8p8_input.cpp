/*----------------------------------------------------------------------*/
/*!
\file so_sh8p8_input.cpp
\brief

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */

#include "so_sh8p8.H"
#include "../drt_mat/artwallremod.H"
#include "../drt_mat/anisotropic_balzani.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/visconeohooke.H"
#include "../drt_mat/viscogenmax.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_mat/aaaraghavanvorp_damage.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_sh8p8::ReadElement(const std::string& eletype,
                                          const std::string& distype,
                                          DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  // special element-dependent input of material parameters
  switch (Material()->MaterialType())
  {
  case INPAR::MAT::m_artwallremod:
  {
    MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(Material().get());
    remo->Setup(NUMGPT_SOH8, this->Id(), linedef);
    break;
  }
  case INPAR::MAT::m_anisotropic_balzani:
  {
    MAT::AnisotropicBalzani* balz = static_cast <MAT::AnisotropicBalzani*>(Material().get());
    balz->Setup(linedef);
    break;
  }
  case INPAR::MAT::m_viscoanisotropic:
  {
    MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(Material().get());
    visco->Setup(NUMGPT_SOH8, linedef);
    break;
  }
  case INPAR::MAT::m_visconeohooke:
  {
    MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(Material().get());
    visco->Setup(NUMGPT_SOH8);
    break;
  }
  case INPAR::MAT::m_viscogenmax:
  {
    MAT::ViscoGenMax* viscogenmax = static_cast <MAT::ViscoGenMax*>(Material().get());
    viscogenmax->Setup(NUMGPT_SOH8,linedef);
    break;
  }
  case INPAR::MAT::m_elasthyper:
  {
    MAT::ElastHyper* elahy = static_cast <MAT::ElastHyper*>(Material().get());
    elahy->Setup(linedef);
    break;
  }
  case INPAR::MAT::m_aaaraghavanvorp_damage:
  {
    double strength = 0.0; // section for extracting the element strength
    linedef->ExtractDouble("STRENGTH",strength);
    MAT::AAAraghavanvorp_damage* aaadamage = static_cast <MAT::AAAraghavanvorp_damage*>(Material().get());
    aaadamage->Setup(NUMGPT_SOH8,strength);
    //aaadamage->Setup(NUMGPT_SOH8);
  }
  default:
    // Do nothing. Simple material.
    break;
  }

  // a temprorary variable for read-in
  std::string buffer;
  // read kinematic flag
  linedef->ExtractString("KINEM",buffer);
  if (buffer=="linear")
  {
    kintype_ = soh8_linear;
    dserror ("Only nonlinear kinematics for SO_SH8P8 implemented!");
  }
  else if (buffer=="nonlinear")
  {
    kintype_ = soh8_nonlinear;
  }
  else dserror ("Reading SO_SH8P8 element failed unknown KINEM Type");


  // we expect kintype to be total lagrangian
  kintype_ = soh8_nonlinear;
  

  // read EAS technology flag
  linedef->ExtractString("EAS",buffer);

  if (buffer=="sosh8")
  {
    eastype_ = soh8_eassosh8;
    neas_ = NUMEAS_SOSH8_;
  }
  else if (buffer=="atype")
  {
    eastype_ = soh8_easa;
    neas_ = NUMEAS_A_;
  }
  else if (buffer=="None")
  {
    eastype_ = soh8_easnone;
    neas_ = 0;
  }
  else if (buffer=="none")
  {
    eastype_ = soh8_easnone;
    neas_ = 0;
  }
  else
    dserror("Reading of SO_SH8P8 EAS type failed");

  if (eastype_ != soh8_easnone)
  {
    EasInit();
  }

  linedef->ExtractString("THICKDIR",buffer);
  nodes_rearranged_ = false;

  // global X
  if      (buffer=="xdir")    thickdir_ = globx;
  // global Y
  else if (buffer=="ydir")    thickdir_ = globy;
  // global Z
  else if (buffer=="zdir")    thickdir_ = globz;
  // find automatically through Jacobian of Xrefe
  else if (buffer=="auto")    thickdir_ = autoj;
  // local r
  else if (buffer=="rdir")    thickdir_ = enfor;
  // local s
  else if (buffer=="sdir")    thickdir_ = enfos;
  // local t
  else if (buffer=="tdir")    thickdir_ = enfot;
  // no noderearrangement
  else if (buffer=="none")
  {
    thickdir_ = none;
    nodes_rearranged_ = true;
  }
  else dserror("Reading of SO_SH8P8 thickness direction failed");

  linedef->ExtractString("STAB",buffer);
  if (buffer=="Aff")
    stab_ = stab_affine;
  else if (buffer=="NonAff")
    stab_ = stab_nonaffine;
  else if (buffer=="SpatAff")
    stab_ = stab_spatialaffine;
  else if (buffer=="Spat")
    stab_ = stab_spatial;
  else if (buffer=="PureDisp")
    stab_ = stab_puredisp;
  else
    dserror("Reading of SO_SH8P8 stabilisation failed");

  linedef->ExtractString("ANS",buffer);
  if (buffer=="Later")
    ans_ = ans_lateral;
  else if (buffer=="OnSpot")
    ans_ = ans_onspot;
  else if (buffer=="None")
    ans_ = ans_none;
  else
    dserror("Reading of SO_SH8P8 ANS type failed");

  // Linearization
  linedef->ExtractString("LIN",buffer);
  if (buffer=="One")
    lin_ = lin_one;
  else if (buffer=="Half")
    lin_ = lin_half;
  else if (buffer=="Sixth")
    lin_ = lin_sixth;
  else
    dserror("Reading of SO_SH8P8 LIN type failed");

  // Isochoric way
  linedef->ExtractString("ISO",buffer);
  if (buffer=="Mat")
    iso_ = iso_material;
  else if (buffer=="Enf")
    iso_ = iso_enforced;
  else
    dserror("Reading of SO_SH8P8 ISO type failed");

  return true;
}
