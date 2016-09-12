/*----------------------------------------------------------------------*/
/*!
 \file porofluid_reactionevaluator.cpp

 \brief evaluation class for reactive terms (mass sources and sinks)

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/



#include "porofluid_reactionevaluator.H"

#include "porofluidmultiphase_ele_calc_utils.H"

#include "../drt_mat/fluidporo_multiphase_reactions.H"
#include "../drt_mat/fluidporo_singlephase.H"
#include "../drt_mat/fluidporo_multiphase_singlereaction.H"
#include "porofluid_phasemanager.H"


/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::PoroFluidReactionEvaluator::PoroFluidReactionEvaluator()
: isevaluated_(false)
{
  return;
}


/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidReactionEvaluator::EvaluateGPState(
    const PoroFluidPhaseManager&       phasemanager,
    const MAT::Material&               material,
    const double&                      porosity,
    const std::vector<double>&         scalar
    )
{
  if(material.MaterialType() != INPAR::MAT::m_fluidporo_multiphase_reactions)
    dserror("Invalid material! Only MAT_FluidPoroMultiPhaseReactions material valid for reaction evaluation!");

  // cast the material to multiphase material
  const MAT::FluidPoroMultiPhaseReactions& multiphasemat =
      dynamic_cast<const MAT::FluidPoroMultiPhaseReactions&>(material);

  // get number of reactions
  const int numreactions = multiphasemat.NumReac();
  // get number of phases
  const int numphases = phasemanager.NumPhases();

  reacterms_.clear();
  reactermsderivs_.clear();
  reacterms_.resize(numphases,0.0);
  reactermsderivs_.resize(numphases,std::vector<double>(numphases,0.0));
  isreactive_.clear();
  isreactive_.resize(numphases,false);

  for(int ireac=0; ireac<numreactions; ireac++)
  {
    // get the single phase material
    MAT::FluidPoroSingleReaction& singlephasemat =
        POROFLUIDMULTIPHASE::ELEUTILS::GetSingleReactionMatFromMultiReactionsMaterial(multiphasemat,ireac);

    // evaluate the reaction
    singlephasemat.EvaluateReaction(
        reacterms_,
        reactermsderivs_,
        phasemanager.Pressure(),
        phasemanager.Saturation(),
        porosity,
        scalar
        );

    for(int iphase=0; iphase<numphases; iphase++)
    {
      isreactive_[iphase] = isreactive_[iphase] or singlephasemat.IsReactive(iphase);
    }
  }

  isevaluated_=true;

  return;
}

/*----------------------------------------------------------------------*
 | zero all values at GP                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidReactionEvaluator::ClearGPState()
{
  // this is just a safety call. All quantities should be recomputed
  // in the next EvaluateGPState call anyway, but you never know ...

  reacterms_.clear();
  reactermsderivs_.clear();
  isreactive_.clear();

  isevaluated_=false;
  return;
}

/*----------------------------------------------------------------------*
 | get the reaction term                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::PoroFluidReactionEvaluator::GetReacTerm(int phasenum) const
{
  if(not isevaluated_)
    dserror("EvaluateGPState was not called!");

  return reacterms_[phasenum];
}

/*----------------------------------------------------------------------*
 | get the derivative of the reaction term                   vuong 08/16 |
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::PoroFluidReactionEvaluator::GetReacDeriv(int phasenum, int doftoderive) const
{
  if(not isevaluated_)
    dserror("EvaluateGPState was not called!");

  return reactermsderivs_[phasenum][doftoderive];
}
