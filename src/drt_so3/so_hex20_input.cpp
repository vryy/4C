/*!----------------------------------------------------------------------
\file so_hex20_input.cpp
\brief

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "so_hex20.H"
#include "../drt_mat/artwallremod.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/visconeohooke.H"
#include "../drt_mat/viscogenmax.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_mat/charmm.H"
#include "../drt_mat/aaaraghavanvorp_damage.H"
#include "../drt_mat/holzapfelcardiovascular.H"
#include "../drt_mat/humphreycardiovascular.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_hex20::ReadElement(const std::string& eletype,
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
    remo->Setup(NUMGPT_SOH20, this->Id(), linedef);
    break;
  }
  case INPAR::MAT::m_viscoanisotropic:
  {
    MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(Material().get());
    visco->Setup(NUMGPT_SOH20, linedef);
    break;
  }
  case INPAR::MAT::m_visconeohooke:
  {
    MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(Material().get());
    visco->Setup(NUMGPT_SOH20);
    break;
  }
  case INPAR::MAT::m_viscogenmax:
  {
    MAT::ViscoGenMax* viscogenmax = static_cast <MAT::ViscoGenMax*>(Material().get());
    viscogenmax->Setup(NUMGPT_SOH20,linedef);
    break;
  }
  case INPAR::MAT::m_charmm:
  {
    MAT::CHARMM* charmm = static_cast <MAT::CHARMM*>(Material().get());
    charmm->Setup(data_);
    break;
  }
  case INPAR::MAT::m_aaaraghavanvorp_damage:
  {
    double strength = 0.0; // section for extracting the element strength
    linedef->ExtractDouble("STRENGTH",strength);
    MAT::AAAraghavanvorp_damage* aaadamage = static_cast <MAT::AAAraghavanvorp_damage*>(Material().get());
    aaadamage->Setup(NUMGPT_SOH20,strength);
    //aaadamage->Setup(NUMGPT_SOH20);
    break;
  }
  case INPAR::MAT::m_holzapfelcardiovascular:
  {
  	MAT::HolzapfelCardio* holzcard = static_cast <MAT::HolzapfelCardio*>(Material().get());
   	holzcard->Setup(NUMGPT_SOH20, linedef);
   	break;
  }
  case INPAR::MAT::m_humphreycardiovascular:
  {
  	MAT::HumphreyCardio* humcard = static_cast <MAT::HumphreyCardio*>(Material().get());
  	humcard->Setup(NUMGPT_SOH20, linedef);
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
    kintype_ = soh20_linear;
  }
  else if (buffer=="nonlinear")
  {
    kintype_ = soh20_nonlinear;
  }
  else dserror ("Reading SO_HEX20 element failed KINEM unknown");

  // check for SVK material if geometrically linear
  if (kintype_==soh20_linear && Material()->MaterialType()!=INPAR::MAT::m_stvenant)
    dserror("ERROR: Only linear elasticity (SVK) for geometrically linear hex20 element");

  return true;
}


#endif  // #ifdef CCADISCRET
