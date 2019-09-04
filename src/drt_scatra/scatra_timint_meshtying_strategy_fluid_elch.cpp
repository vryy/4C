/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid-fluid meshtying strategy for electrochemistry problems

\level 2

\maintainer Christoph Schmidt


*----------------------------------------------------------------------*/
#include "scatra_timint_meshtying_strategy_fluid_elch.H"

#include "../drt_fluid/fluid_meshtying.H"

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
  if (DRT::INPUT::IntegralValue<INPAR::FLUID::MeshTying>(*(scatratimint_->ScatraParameterList()),
          "MESHTYING") == INPAR::FLUID::condensed_bmat_merged and
      ElchTimInt()->EquPot() == INPAR::ELCH::equpot_enc)
    dserror(
        "In the context of meshtying, the ion-transport system including the electroneutrality "
        "condition cannot be solved in a block matrix!");

  // call setup in base class
  SCATRA::MeshtyingStrategyFluid::InitMeshtying();

  return;
}  // SCATRA::MeshtyingStrategyFluidElch::SetupMeshtying


/*------------------------------------------------------------------------------------*
 | initialize system matrix for fluid-fluid meshtying (electrochemistry)   fang 12/14 |
 *------------------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseOperator> SCATRA::MeshtyingStrategyFluidElch::InitSystemMatrix() const
{
  // safety check
  if (scatratimint_->NumScal() < 1) dserror("Number of transported scalars not correctly set!");

  // define coupling
  // standard case: all dofs (transported scalars and electric potential) are coupled
  // =>   coupleddof = [1, 1, ..., 1, 1]
  // special case: only potential is coupled
  // =>   coupleddof = [0, 0, ..., 0, 1]
  std::vector<int> coupleddof(scatratimint_->NumScal() + 1, 1);
  if (DRT::INPUT::IntegralValue<int>(*(ElchTimInt()->ElchParameterList()), "ONLYPOTENTIAL"))
    for (int i = 0; i < scatratimint_->NumScal(); ++i) coupleddof[i] = 0;

  return meshtying_->Setup(coupleddof);
}  // SCATRA::MeshtyingStrategyFluidElch::InitSystemMatrix


/*------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   fang 02/16 |
 *------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyFluidElch::InitConvCheckStrategy()
{
  convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyStdElch(
      scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));

  return;
}  // SCATRA::MeshtyingStrategyFluidElch::InitConvCheckStrategy
