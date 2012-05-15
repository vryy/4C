/*!----------------------------------------------------------------------
\file so_sh8_input.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/


#include "so_sh8.H"
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
bool DRT::ELEMENTS::So_sh8::ReadElement(const std::string& eletype,
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


  // temporary variable for read-in
  std::string buffer;


  // read kinematic flag
  linedef->ExtractString("KINEM",buffer);
  if (buffer=="linear")
  {
   //kintype_ = soh8_linear;
   dserror ("Only nonlinear kinematics for SO_SH8 implemented!");
  }
  else if (buffer=="nonlinear")
  {
   kintype_ = soh8_nonlinear;
  }
  else dserror ("Reading SO_HEX8p1j1 element failed KINEM unknown");

  // we expect kintype to be total lagrangian
  kintype_ = soh8_nonlinear;

  // read EAS technology flag
  linedef->ExtractString("EAS",buffer);

  // full EAS technology
  if      (buffer=="sosh8")
  {
    eastype_ = soh8_eassosh8;
    neas_ = 7;               // number of eas parameters for EAS_SOSH8
    soh8_easinit();
  }
  // no EAS technology
  else if (buffer=="none")
  {
    eastype_ = soh8_easnone;
    neas_ = 0;               // number of eas parameters for EAS_SOSH8
  }
  else
    dserror("Reading of SO_SH8 EAS technology failed");

  // read ANS technology flag
  linedef->ExtractString("ANS",buffer);
  if      (buffer=="sosh8")
  {
    anstype_ = anssosh8;
  }
  // no ANS technology
  else if (buffer=="none")
  {
    anstype_ = ansnone;
  }
  else
    dserror("Reading of SO_SH8 ANS technology failed");
  
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
  else dserror("Reading of SO_SH8 thickness direction failed");

  return true;
}


