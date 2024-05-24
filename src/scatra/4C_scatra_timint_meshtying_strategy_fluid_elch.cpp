/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid-fluid meshtying strategy for electrochemistry problems

\level 2



*----------------------------------------------------------------------*/
#include "4C_scatra_timint_meshtying_strategy_fluid_elch.hpp"

#include "4C_fluid_meshtying.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyFluidElch::MeshtyingStrategyFluidElch(SCATRA::ScaTraTimIntElch* elchtimint)
    : MeshtyingStrategyFluid(elchtimint)
{
  return;
}  // SCATRA::MeshtyingStrategyFluidElch::MeshtyingStrategyFluidElch


/*----------------------------------------------------------------------------*
 | perform setup of fluid-fluid meshtying (electrochemistry)       fang 12/14 |
 *----------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyFluidElch::InitMeshtying()
{
  // Important: Meshtying for electrochemistry problems is not well tested!
  // safety check
  if (CORE::UTILS::IntegralValue<INPAR::FLUID::MeshTying>(*(scatratimint_->ScatraParameterList()),
          "MESHTYING") == INPAR::FLUID::condensed_bmat_merged and
      elch_tim_int()->EquPot() == INPAR::ELCH::equpot_enc)
    FOUR_C_THROW(
        "In the context of meshtying, the ion-transport system including the electroneutrality "
        "condition cannot be solved in a block matrix!");

  // call setup in base class
  SCATRA::MeshtyingStrategyFluid::InitMeshtying();

  return;
}

void SCATRA::MeshtyingStrategyFluidElch::SetupMeshtying()
{
  // safety check
  if (scatratimint_->NumScal() < 1)
    FOUR_C_THROW("Number of transported scalars not correctly set!");

  // define coupling
  // standard case: all dofs (transported scalars and electric potential) are coupled
  // =>   coupleddof = [1, 1, ..., 1, 1]
  // special case: only potential is coupled
  // =>   coupleddof = [0, 0, ..., 0, 1]
  std::vector<int> coupleddof(scatratimint_->NumScal() + 1, 1);
  if (CORE::UTILS::IntegralValue<int>(*(elch_tim_int()->ElchParameterList()), "ONLYPOTENTIAL"))
    for (int i = 0; i < scatratimint_->NumScal(); ++i) coupleddof[i] = 0;

  meshtying_->SetupMeshtying(coupleddof);
}

/*------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   fang 02/16 |
 *------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyFluidElch::init_conv_check_strategy()
{
  convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyStdElch(
      scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));

  return;
}  // SCATRA::MeshtyingStrategyFluidElch::init_conv_check_strategy

FOUR_C_NAMESPACE_CLOSE
