/*----------------------------------------------------------------------*/
/*! \file
 \brief helpful methods for the porofluidmultiphase element

   \level 3

 *----------------------------------------------------------------------*/


#include "4C_porofluidmultiphase_ele_calc_utils.hpp"

#include "4C_mat_fluidporo_multiphase.hpp"
#include "4C_mat_fluidporo_multiphase_reactions.hpp"
#include "4C_mat_fluidporo_multiphase_singlereaction.hpp"
#include "4C_mat_fluidporo_singlephase.hpp"

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------*
 * get the single phase material from the element multiphase reaction material   vuong 08/16 |
 *-------------------------------------------------------------------------------------------*/
Mat::FluidPoroSingleReaction&
POROFLUIDMULTIPHASE::ElementUtils::GetSingleReactionMatFromMultiReactionsMaterial(
    const Mat::FluidPoroMultiPhaseReactions& multiphasereacmat, int phasenum)
{
  // get the single phase material by its ID
  const int matid = multiphasereacmat.ReacID(phasenum);
  Teuchos::RCP<Core::Mat::Material> singlemat = multiphasereacmat.MaterialById(matid);

  // safety check and cast
  if (singlemat->MaterialType() != Core::Materials::m_fluidporo_singlereaction)
    FOUR_C_THROW("only poro singleraction material valid");

  return static_cast<Mat::FluidPoroSingleReaction&>(*singlemat);
}

/*----------------------------------------------------------------------------------*
 * get the single phase material from the element multiphase material    vuong 08/16 |
 *-----------------------------------------------------------------------------------*/
const Mat::FluidPoroSinglePhase&
POROFLUIDMULTIPHASE::ElementUtils::GetSinglePhaseMatFromMultiMaterial(
    const Mat::FluidPoroMultiPhase& multiphasemat, int phasenum)
{
  // get the single phase material by its ID
  const int matid = multiphasemat.MatID(phasenum);
  Teuchos::RCP<Core::Mat::Material> singlemat = multiphasemat.MaterialById(matid);

  // safety check and cast
  if (singlemat->MaterialType() != Core::Materials::m_fluidporo_singlephase)
    FOUR_C_THROW("check at position %i/%i failed, only poro singlephase material valid",
        phasenum + 1, multiphasemat.NumMat());

  return static_cast<const Mat::FluidPoroSinglePhase&>(*singlemat);
}

/*------------------------------------------------------------------------*
 *  get the single phase material from the element material   vuong 08/16 |
 *-------------------------------------------------------------------------*/
const Mat::FluidPoroSinglePhase& POROFLUIDMULTIPHASE::ElementUtils::GetSinglePhaseMatFromMaterial(
    const Core::Mat::Material& material, int phasenum)
{
  // safety check
  if (material.MaterialType() != Core::Materials::m_fluidporo_multiphase and
      material.MaterialType() != Core::Materials::m_fluidporo_multiphase_reactions)
    FOUR_C_THROW("only poro multiphase material valid");

  // cast
  const Mat::FluidPoroMultiPhase& multiphasemat =
      static_cast<const Mat::FluidPoroMultiPhase&>(material);

  return GetSinglePhaseMatFromMultiMaterial(multiphasemat, phasenum);
}

/*---------------------------------------------------------------------------------------*
 * get the single volfrac material from the element multiphase material kremheller 08/17 |
 *----------------------------------------------------------------------------------------*/
const Mat::FluidPoroSingleVolFrac&
POROFLUIDMULTIPHASE::ElementUtils::GetSingleVolFracMatFromMultiMaterial(
    const Mat::FluidPoroMultiPhase& multiphasemat, int volfracnum)
{
  // get the single phase material by its ID
  const int matid = multiphasemat.MatID(volfracnum);
  Teuchos::RCP<Core::Mat::Material> singlemat = multiphasemat.MaterialById(matid);

  // safety check and cast
  if (singlemat->MaterialType() != Core::Materials::m_fluidporo_singlevolfrac)
    FOUR_C_THROW("check at position %i/%i failed, only poro single vol fraction material valid",
        volfracnum + 1, multiphasemat.NumMat());

  return static_cast<const Mat::FluidPoroSingleVolFrac&>(*singlemat);
}

/*-------------------------------------------------------------------------------*
 *  get the single volfrac material from the element material   kremheller 08/17 |
 *--------------------------------------------------------------------------------*/
const Mat::FluidPoroSingleVolFrac&
POROFLUIDMULTIPHASE::ElementUtils::GetSingleVolFracMatFromMaterial(
    const Core::Mat::Material& material, int volfracnum)
{
  // safety check
  if (material.MaterialType() != Core::Materials::m_fluidporo_multiphase and
      material.MaterialType() != Core::Materials::m_fluidporo_multiphase_reactions)
    FOUR_C_THROW("only poro multiphase material valid");

  // cast
  const Mat::FluidPoroMultiPhase& multiphasemat =
      static_cast<const Mat::FluidPoroMultiPhase&>(material);

  return GetSingleVolFracMatFromMultiMaterial(multiphasemat, volfracnum);
}

/*-------------------------------------------------------------------------------------------------*
 * get the volume fraction pressure material from the element multiphase material kremheller 02/18 |
 *--------------------------------------------------------------------------------------------------*/
const Mat::FluidPoroVolFracPressure&
POROFLUIDMULTIPHASE::ElementUtils::GetVolFracPressureMatFromMultiMaterial(
    const Mat::FluidPoroMultiPhase& multiphasemat, int volfracnum)
{
  // get the single phase material by its ID
  const int matid = multiphasemat.MatID(volfracnum);
  Teuchos::RCP<Core::Mat::Material> singlemat = multiphasemat.MaterialById(matid);

  // safety check and cast
  if (singlemat->MaterialType() != Core::Materials::m_fluidporo_volfracpressure)
    FOUR_C_THROW("check at position %i/%i failed, only poro single vol fraction material valid",
        volfracnum + 1, multiphasemat.NumMat());

  return static_cast<const Mat::FluidPoroVolFracPressure&>(*singlemat);
}

/*-----------------------------------------------------------------------------------------*
 *  get the volume fraction pressure material from the element material   kremheller 02/18 |
 *------------------------------------------------------------------------------------------*/
const Mat::FluidPoroVolFracPressure&
POROFLUIDMULTIPHASE::ElementUtils::GetVolFracPressureMatFromMaterial(
    const Core::Mat::Material& material, int volfracnum)
{
  // safety check
  if (material.MaterialType() != Core::Materials::m_fluidporo_multiphase and
      material.MaterialType() != Core::Materials::m_fluidporo_multiphase_reactions)
    FOUR_C_THROW("only poro multiphase material valid");

  // cast
  const Mat::FluidPoroMultiPhase& multiphasemat =
      static_cast<const Mat::FluidPoroMultiPhase&>(material);

  return GetVolFracPressureMatFromMultiMaterial(multiphasemat, volfracnum);
}

FOUR_C_NAMESPACE_CLOSE
