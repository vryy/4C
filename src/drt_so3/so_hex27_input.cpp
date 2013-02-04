/*!----------------------------------------------------------------------
\file so_hex27_input.cpp
\brief

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#include "so_hex27.H"
#include "../drt_mat/artwallremod.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/visconeohooke.H"
#include "../drt_mat/viscogenmax.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_mat/charmm.H"
#include "../drt_mat/aaaraghavanvorp_damage.H"
#include "../drt_mat/holzapfelcardiovascular.H"
#include "../drt_mat/humphreycardiovascular.H"
#include "../drt_mat/structporo.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_hex27::ReadElement(const std::string& eletype,
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

  // special element-dependent input of material parameters
  switch (mat->MaterialType())
  {
  case INPAR::MAT::m_artwallremod:
  {
    MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(mat.get());
    remo->Setup(NUMGPT_SOH27, this->Id(), linedef);
    break;
  }
  case INPAR::MAT::m_viscoanisotropic:
  {
    MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(mat.get());
    visco->Setup(NUMGPT_SOH27, linedef);
    break;
  }
  case INPAR::MAT::m_visconeohooke:
  {
    MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(mat.get());
    visco->Setup(NUMGPT_SOH27);
    break;
  }
  case INPAR::MAT::m_viscogenmax:
  {
    MAT::ViscoGenMax* viscogenmax = static_cast <MAT::ViscoGenMax*>(mat.get());
    viscogenmax->Setup(NUMGPT_SOH27, linedef);
    break;
  }
  case INPAR::MAT::m_charmm:
  {
    MAT::CHARMM* charmm = static_cast <MAT::CHARMM*>(mat.get());
    charmm->Setup(data_);
    break;
  }
  case INPAR::MAT::m_aaaraghavanvorp_damage:
  {
    double strength = 0.0; // section for extracting the element strength
    linedef->ExtractDouble("STRENGTH",strength);
    MAT::AAAraghavanvorp_damage* aaadamage = static_cast <MAT::AAAraghavanvorp_damage*>(mat.get());
    aaadamage->Setup(NUMGPT_SOH27,strength);
    //aaadamage->Setup(NUMGPT_SOH27);
    break;
  }
  case INPAR::MAT::m_holzapfelcardiovascular:
  {
  	MAT::HolzapfelCardio* holzcard = static_cast <MAT::HolzapfelCardio*>(mat.get());
  	holzcard->Setup(NUMGPT_SOH27, linedef);
  	break;
  }
  case INPAR::MAT::m_humphreycardiovascular:
  {
  	MAT::HumphreyCardio* humcard = static_cast <MAT::HumphreyCardio*>(mat.get());
  	humcard->Setup(NUMGPT_SOH27, linedef);
  	break;
  }
  case INPAR::MAT::m_elasthyper
  :{
    MAT::ElastHyper* elahy = static_cast <MAT::ElastHyper*>(mat.get());
    elahy->Setup(linedef);
    break;
  }
  default:
    // Do nothing. Simple material.
    break;
  }

  std::string buffer;
  linedef->ExtractString("KINEM",buffer);

  if (buffer=="linear")
  {
    kintype_ = soh27_linear;
  }
  else if (buffer=="nonlinear")
  {
    kintype_ = soh27_nonlinear;
  }
  else dserror ("Reading SO_HEX27 element failed KINEM unknown");

  // check for SVK material if geometrically linear
  if (kintype_==soh27_linear && mat->MaterialType()!=INPAR::MAT::m_stvenant)
    dserror("ERROR: Only linear elasticity (SVK) for geometrically linear hex27 element");

  return true;
}


