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
ScaTra::MeshtyingStrategyFluidElch::MeshtyingStrategyFluidElch(ScaTra::ScaTraTimIntElch* elchtimint)
    : MeshtyingStrategyFluid(elchtimint)
{
  return;
}  // ScaTra::MeshtyingStrategyFluidElch::MeshtyingStrategyFluidElch


/*----------------------------------------------------------------------------*
 | perform setup of fluid-fluid meshtying (electrochemistry)       fang 12/14 |
 *----------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyFluidElch::init_meshtying()
{
  // Important: Meshtying for electrochemistry problems is not well tested!
  // safety check
  if (Core::UTILS::IntegralValue<Inpar::FLUID::MeshTying>(*(scatratimint_->scatra_parameter_list()),
          "MESHTYING") == Inpar::FLUID::condensed_bmat_merged and
      elch_tim_int()->equ_pot() == Inpar::ElCh::equpot_enc)
    FOUR_C_THROW(
        "In the context of meshtying, the ion-transport system including the electroneutrality "
        "condition cannot be solved in a block matrix!");

  // call setup in base class
  ScaTra::MeshtyingStrategyFluid::init_meshtying();

  return;
}

void ScaTra::MeshtyingStrategyFluidElch::setup_meshtying()
{
  // safety check
  if (scatratimint_->num_scal() < 1)
    FOUR_C_THROW("Number of transported scalars not correctly set!");

  // define coupling
  // standard case: all dofs (transported scalars and electric potential) are coupled
  // =>   coupleddof = [1, 1, ..., 1, 1]
  // special case: only potential is coupled
  // =>   coupleddof = [0, 0, ..., 0, 1]
  std::vector<int> coupleddof(scatratimint_->num_scal() + 1, 1);
  if (Core::UTILS::IntegralValue<int>(*(elch_tim_int()->elch_parameter_list()), "ONLYPOTENTIAL"))
    for (int i = 0; i < scatratimint_->num_scal(); ++i) coupleddof[i] = 0;

  meshtying_->setup_meshtying(coupleddof);
}

/*------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   fang 02/16 |
 *------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyFluidElch::init_conv_check_strategy()
{
  convcheckstrategy_ = Teuchos::rcp(new ScaTra::ConvCheckStrategyStdElch(
      scatratimint_->scatra_parameter_list()->sublist("NONLINEAR")));

  return;
}  // ScaTra::MeshtyingStrategyFluidElch::init_conv_check_strategy

FOUR_C_NAMESPACE_CLOSE
