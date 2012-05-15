/*!----------------------------------------------------------------------
\file so_weg6_input.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/

#include "so_weg6.H"
#include "../drt_mat/artwallremod.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/holzapfelcardiovascular.H"
#include "../drt_mat/humphreycardiovascular.H"
#include "../drt_mat/growth_ip.H"
#include "../drt_mat/constraintmixture.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_mat/elasthyper.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_weg6::ReadElement(const std::string& eletype,
                                         const std::string& distype,
                                         DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  // special element-dependent input of material parameters
  if (Material()->MaterialType() == INPAR::MAT::m_artwallremod){
    MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(Material().get());
    remo->Setup(NUMGPT_WEG6, this->Id(), linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_viscoanisotropic){
    MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(Material().get());
    visco->Setup(NUMGPT_WEG6, linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_holzapfelcardiovascular){
    MAT::HolzapfelCardio* holzcard = static_cast <MAT::HolzapfelCardio*>(Material().get());
    holzcard->Setup(NUMGPT_WEG6, linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_humphreycardiovascular){
    MAT::HumphreyCardio* humcard = static_cast <MAT::HumphreyCardio*>(Material().get());
    humcard->Setup(NUMGPT_WEG6, linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_growth){
    MAT::Growth* grow = static_cast <MAT::Growth*>(Material().get());
    grow->Setup(NUMGPT_WEG6, linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_constraintmixture){
    MAT::ConstraintMixture* comix = static_cast <MAT::ConstraintMixture*>(Material().get());
    comix->Setup(NUMGPT_WEG6, linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_elasthyper){
    MAT::ElastHyper* elahy = static_cast <MAT::ElastHyper*>(Material().get());
    elahy->Setup(linedef);
  }

  std::string buffer;
  linedef->ExtractString("KINEM",buffer);
  if (buffer=="linear")
  {
    kintype_ = sow6_linear;
    dserror("Reading of SO_WEG6 element failed only nonlinear kinematics implemented");
  }
  else if (buffer=="nonlinear")
  {
    kintype_ = sow6_nonlinear;
  }
  else dserror ("Reading SO_WEG6 element failed KINEM unknwon");

  return true;
}
