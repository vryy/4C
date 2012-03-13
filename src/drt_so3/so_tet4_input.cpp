/*!----------------------------------------------------------------------**###
\file so_tet4_input.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
writen by : Alexander Volf
			alexander.volf@mytum.de
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "so_tet4.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_mat/holzapfelcardiovascular.H"
#include "../drt_mat/humphreycardiovascular.H"
#include "../drt_mat/constraintmixture.H"
#include "../drt_mat/elasthyper.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_tet4::ReadElement(const std::string& eletype,
                                         const std::string& distype,
                                         DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  if (Material()->MaterialType() == INPAR::MAT::m_holzapfelcardiovascular){
    MAT::HolzapfelCardio* holzcard = static_cast <MAT::HolzapfelCardio*>(Material().get());
    holzcard->Setup(NUMGPT_SOTET4, linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_humphreycardiovascular){
    MAT::HumphreyCardio* humcard = static_cast <MAT::HumphreyCardio*>(Material().get());
    humcard->Setup(NUMGPT_SOTET4, linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_constraintmixture){
    MAT::ConstraintMixture* comix = static_cast <MAT::ConstraintMixture*>(Material().get());
    comix->Setup(NUMGPT_SOTET4, linedef);
  } else if (Material()->MaterialType() == INPAR::MAT::m_elasthyper){
    MAT::ElastHyper* elahy = static_cast <MAT::ElastHyper*>(Material().get());
    elahy->Setup(linedef);
  }

  std::string buffer;
  linedef->ExtractString("KINEM",buffer);

  // geometrically linear
  if(buffer=="linear")
  {
    kintype_ = so_tet4_linear;
    dserror("Reading of SO_TET4 element failed only nonlinear kinematics implemented");
  }
  // geometrically non-linear with Total Lagrangean approach
  else if (buffer=="nonlinear")    kintype_ = so_tet4_nonlinear;
  else dserror("Reading of SO_TET4 element failed KINEM unknown");

  return true;
}


#endif  // #ifdef CCADISCRET
