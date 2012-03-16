/*!----------------------------------------------------------------------
\file so_hex8_input.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "so_hex8.H"
#include "../drt_mat/artwallremod.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/visconeohooke.H"
#include "../drt_mat/viscogenmax.H"
#include "../drt_mat/charmm.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_mat/aaaraghavanvorp_damage.H"
#include "../drt_mat/holzapfelcardiovascular.H"
#include "../drt_mat/humphreycardiovascular.H"
#include "../drt_mat/plasticneohooke.H"
#include "../drt_mat/plasticlinelast.H"
#include "../drt_mat/thermoplasticlinelast.H"
#include "../drt_mat/robinson.H"
#include "../drt_mat/growth_ip.H"
#include "../drt_mat/constraintmixture.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_hex8::ReadElement(const std::string& eletype,
                                         const std::string& distype,
                                         DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);

  SetMaterial(material);

  // special element-dependent input of material parameters
  if (Material()->MaterialType() == INPAR::MAT::m_artwallremod)
  {
    MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(Material().get());
    remo->Setup(NUMGPT_SOH8, this->Id(), linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_viscoanisotropic){
    MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(Material().get());
    visco->Setup(NUMGPT_SOH8, linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_visconeohooke){
    MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(Material().get());
    visco->Setup(NUMGPT_SOH8);
  } else if (Material()->MaterialType() == INPAR::MAT::m_viscogenmax){
    MAT::ViscoGenMax* viscogenmax = static_cast <MAT::ViscoGenMax*>(Material().get());
    viscogenmax->Setup(NUMGPT_SOH8,linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_charmm){
    MAT::CHARMM* charmm = static_cast <MAT::CHARMM*>(Material().get());
    charmm->Setup(data_);
  } else if (Material()->MaterialType() == INPAR::MAT::m_elasthyper){
    MAT::ElastHyper* elahy = static_cast <MAT::ElastHyper*>(Material().get());
    elahy->Setup(linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_aaaraghavanvorp_damage){
    double strength = 0.0; // section for extracting the element strength
    linedef->ExtractDouble("STRENGTH",strength);
    MAT::AAAraghavanvorp_damage* aaadamage = static_cast <MAT::AAAraghavanvorp_damage*>(Material().get());
    aaadamage->Setup(NUMGPT_SOH8,strength);
  } else if (Material()->MaterialType() == INPAR::MAT::m_holzapfelcardiovascular){
    MAT::HolzapfelCardio* holzcard = static_cast <MAT::HolzapfelCardio*>(Material().get());
    holzcard->Setup(NUMGPT_SOH8, linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_plneohooke){
    MAT::PlasticNeoHooke* plastic = static_cast <MAT::PlasticNeoHooke*>(Material().get());
    plastic->Setup(NUMGPT_SOH8);
  } else if (Material()->MaterialType() == INPAR::MAT::m_pllinelast){
    MAT::PlasticLinElast* plastic = static_cast <MAT::PlasticLinElast*>(Material().get());
    plastic->Setup(NUMGPT_SOH8);
  } else if (Material()->MaterialType() == INPAR::MAT::m_thermopllinelast){
    MAT::ThermoPlasticLinElast* plastic = static_cast <MAT::ThermoPlasticLinElast*>(Material().get());
    plastic->Setup(NUMGPT_SOH8);
  } else if (Material()->MaterialType() == INPAR::MAT::m_vp_robinson){
    MAT::Robinson* robinson = static_cast <MAT::Robinson*>(Material().get());
    robinson->Setup(NUMGPT_SOH8, linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_humphreycardiovascular){
    MAT::HumphreyCardio* humcard = static_cast <MAT::HumphreyCardio*>(Material().get());
    humcard->Setup(NUMGPT_SOH8, linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_growth){
    MAT::Growth* grow = static_cast <MAT::Growth*>(Material().get());
    grow->Setup(NUMGPT_SOH8, linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_constraintmixture){
    MAT::ConstraintMixture* comix = static_cast <MAT::ConstraintMixture*>(Material().get());
    comix->Setup(NUMGPT_SOH8, linedef);
  }

  // temporary variable for read-in
   std::string buffer;

  // read kinematic flag
  linedef->ExtractString("KINEM",buffer);
  if (buffer=="linear")
  {
   kintype_ = soh8_linear;
  }
  else if (buffer=="nonlinear")
  {
   kintype_ = soh8_nonlinear;
  }
  else dserror ("Reading SO_HEX8 element failed");

  // read EAS technology flag
  linedef->ExtractString("EAS",buffer);

  // full EAS technology
  if      (buffer=="full")
  {
    eastype_ = soh8_easfull;
    neas_ = 21;               // number of eas parameters for full EAS
    soh8_easinit();
  }
  // mild EAS technology
  else if (buffer=="mild")
  {
    eastype_ = soh8_easmild;
    neas_ = 9;               // number of eas parameters for mild EAS
    soh8_easinit();
  }
  // no EAS technology
  else if (buffer=="none")
  {
    eastype_ = soh8_easnone;
    neas_ = 0;
  }
  else
    dserror("Reading of SO_HEX8 EAS technology failed");

  return true;
}

#endif  // #ifdef CCADISCRET
