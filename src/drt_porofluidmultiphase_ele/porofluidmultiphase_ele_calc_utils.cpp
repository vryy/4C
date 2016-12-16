/*----------------------------------------------------------------------*/
/*!
 \file porofluidmultiphase_ele_calc_utils.cpp

 \brief helpful methods for the porofluidmultiphase element

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/


#include "porofluidmultiphase_ele_calc_utils.H"

#include "../drt_mat/fluidporo_multiphase.H"
#include "../drt_mat/fluidporo_singlephase.H"
#include "../drt_mat/fluidporo_multiphase_reactions.H"
#include "../drt_mat/fluidporo_multiphase_singlereaction.H"

/*-----------------------------------------------------------------------------------------*
 * get the single phase material from the element multiphase reaction material   vuong 08/16 |
*-------------------------------------------------------------------------------------------*/
MAT::FluidPoroSingleReaction&
POROFLUIDMULTIPHASE::ELEUTILS::GetSingleReactionMatFromMultiReactionsMaterial(
    const MAT::FluidPoroMultiPhaseReactions& multiphasereacmat,
    int phasenum)
{
  // get the single phase material by its ID
  const int matid = multiphasereacmat.ReacID(phasenum);
  Teuchos::RCP< MAT::Material> singlemat = multiphasereacmat.MaterialById(matid);

  // safety check and cast
  if(singlemat->MaterialType() != INPAR::MAT::m_fluidporo_singlereaction)
    dserror("only poro singleraction material valid");

  return static_cast<MAT::FluidPoroSingleReaction&>(*singlemat);
}

/*----------------------------------------------------------------------------------*
 * get the single phase material from the element multiphase material    vuong 08/16 |
*-----------------------------------------------------------------------------------*/
const MAT::FluidPoroSinglePhase&
POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMultiMaterial(
    const MAT::FluidPoroMultiPhase& multiphasemat,
    int phasenum)
{
  // get the single phase material by its ID
  const int matid = multiphasemat.MatID(phasenum);
  Teuchos::RCP< MAT::Material> singlemat = multiphasemat.MaterialById(matid);

  // safety check and cast
  if(singlemat->MaterialType() != INPAR::MAT::m_fluidporo_singlephase)
    dserror("only poro singlephase material valid");

  return static_cast<const MAT::FluidPoroSinglePhase&>(*singlemat);
}

/*------------------------------------------------------------------------*
 *  get the single phase material from the element material   vuong 08/16 |
*-------------------------------------------------------------------------*/
const MAT::FluidPoroSinglePhase&
POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMaterial(
    const MAT::Material& material,
    int phasenum)
{
  //safety check
  if(material.MaterialType() != INPAR::MAT::m_fluidporo_multiphase and
     material.MaterialType() != INPAR::MAT::m_fluidporo_multiphase_reactions)
    dserror("only poro multiphase material valid");

  //cast
  const MAT::FluidPoroMultiPhase& multiphasemat =
      static_cast<const MAT::FluidPoroMultiPhase&>(material);

  return GetSinglePhaseMatFromMultiMaterial(multiphasemat,phasenum);
}
