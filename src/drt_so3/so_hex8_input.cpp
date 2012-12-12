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
#include "../drt_mat/damage.H"
#include "../drt_mat/growth_ip.H"
#include "../drt_mat/constraintmixture.H"
#include "../drt_mat/structporo.H"
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

  Teuchos::RCP<MAT::Material> mat = Material();

  if(mat->MaterialType() == INPAR::MAT::m_structporo)
  {
    MAT::StructPoro* actmat = static_cast<MAT::StructPoro*>(mat.get());
    //setup is done in so3_poro
    //actmat->Setup(NUMGPT_SOH8);
    mat = actmat->GetMaterial();
  }

  switch (mat->MaterialType()){
  // special element-dependent input of material parameters
    case  INPAR::MAT::m_artwallremod:
    {
      MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(mat.get());
      remo->Setup(NUMGPT_SOH8, this->Id(), linedef);
    }
      break;
    case INPAR::MAT::m_viscoanisotropic:
    {
      MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(mat.get());
      visco->Setup(NUMGPT_SOH8, linedef);
    }
      break;
    case INPAR::MAT::m_visconeohooke:
    {
      MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(mat.get());
      visco->Setup(NUMGPT_SOH8);
    }
      break;
    case  INPAR::MAT::m_viscogenmax:
    {
      MAT::ViscoGenMax* viscogenmax = static_cast <MAT::ViscoGenMax*>(mat.get());
      viscogenmax->Setup(NUMGPT_SOH8,linedef);
    }
      break;
    case INPAR::MAT::m_charmm:
    {
      MAT::CHARMM* charmm = static_cast <MAT::CHARMM*>(mat.get());
      charmm->Setup(data_);
    }
      break;
    case INPAR::MAT::m_elasthyper:
    {
      MAT::ElastHyper* elahy = static_cast <MAT::ElastHyper*>(mat.get());
      elahy->Setup(linedef);
    }
      break;
    case  INPAR::MAT::m_aaaraghavanvorp_damage:
    {
      double strength = 0.0; // section for extracting the element strength
      linedef->ExtractDouble("STRENGTH",strength);
      MAT::AAAraghavanvorp_damage* aaadamage = static_cast <MAT::AAAraghavanvorp_damage*>(mat.get());
      aaadamage->Setup(NUMGPT_SOH8,strength);
    }
      break;
    case INPAR::MAT::m_holzapfelcardiovascular:
    {
      MAT::HolzapfelCardio* holzcard = static_cast <MAT::HolzapfelCardio*>(mat.get());
      holzcard->Setup(NUMGPT_SOH8, linedef);
    }
      break;
    case INPAR::MAT::m_plneohooke:
    {
      MAT::PlasticNeoHooke* plastic = static_cast <MAT::PlasticNeoHooke*>(mat.get());
      plastic->Setup(NUMGPT_SOH8);
    }
      break;
    case  INPAR::MAT::m_pllinelast:
    {
      MAT::PlasticLinElast* plastic = static_cast <MAT::PlasticLinElast*>(mat.get());
      plastic->Setup(NUMGPT_SOH8);
    }
      break;
    case  INPAR::MAT::m_thermopllinelast:
    {
      MAT::ThermoPlasticLinElast* plastic = static_cast <MAT::ThermoPlasticLinElast*>(mat.get());
      plastic->Setup(NUMGPT_SOH8);
    }
      break;
    case  INPAR::MAT::m_vp_robinson:
    {
      MAT::Robinson* robinson = static_cast <MAT::Robinson*>(mat.get());
      robinson->Setup(NUMGPT_SOH8, linedef);
    }
      break;
    case  INPAR::MAT::m_elpldamage:
    {
      MAT::Damage* damage = static_cast <MAT::Damage*>(mat.get());
      damage->Setup(NUMGPT_SOH8);
    }
      break;
    case  INPAR::MAT::m_humphreycardiovascular:
    {
      MAT::HumphreyCardio* humcard = static_cast <MAT::HumphreyCardio*>(mat.get());
      humcard->Setup(NUMGPT_SOH8, linedef);
    }
      break;
    case INPAR::MAT::m_growth:
    {
      MAT::Growth* grow = static_cast <MAT::Growth*>(mat.get());
      grow->Setup(NUMGPT_SOH8, linedef);
    }
      break;
    case INPAR::MAT::m_constraintmixture:
    {
      MAT::ConstraintMixture* comix = static_cast <MAT::ConstraintMixture*>(mat.get());
      comix->Setup(NUMGPT_SOH8, linedef);
    }
      break;
    default:
      break;
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

