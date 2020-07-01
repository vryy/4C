/*----------------------------------------------------------------------*/
/*! \file

\brief Standard solution strategy for electrochemistry problems (without meshtying)

\level 2


*----------------------------------------------------------------------*/
#include "scatra_timint_meshtying_strategy_std_elch.H"

#include "../drt_scatra/scatra_utils_splitstrategy.H"

#include "../linalg/linalg_utils_sparse_algebra_create.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyStdElch::MeshtyingStrategyStdElch(SCATRA::ScaTraTimIntElch* elchtimint)
    : MeshtyingStrategyStd(elchtimint)
{
  return;
}  // SCATRA::MeshtyingStrategyStdElch::MeshtyingStrategyStdElch


/*----------------------------------------------------------------------*
 | initialize system matrix for electrochemistry problems    fang 12/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseOperator> SCATRA::MeshtyingStrategyStdElch::InitSystemMatrix() const
{
  Teuchos::RCP<LINALG::SparseOperator> systemmatrix;

  if (DRT::INPUT::IntegralValue<int>(*(ElchTimInt()->ElchParameterList()), "BLOCKPRECOND"))
  {
    // safety checks
    if (ElchTimInt()->EquPot() == INPAR::ELCH::equpot_undefined)
      dserror("Type of closing equation for electric potential not correctly set!");
    if (ElchTimInt()->EquPot() != INPAR::ELCH::equpot_enc)
      dserror("Special ELCH assemble strategy for block-matrix will not assemble A_11 block!");
    if (scatratimint_->NumScal() < 1) dserror("Number of transported scalars not correctly set!");

    // initial guess for non-zeros per row: 27 neighboring nodes for hex8
    // this is enough! A higher guess would require too much memory!
    // A more precise guess for every submatrix would read:
    // A_00: 27*1,  A_01: 27*1,  A_10: 27*numscal due to electroneutrality, A_11: EMPTY matrix !!!!!
    // usage of a split strategy that makes use of the ELCH-specific sparsity pattern
    LINALG::MapExtractor splitter;
    LINALG::CreateMapExtractorFromDiscretization(
        *(scatratimint_->Discretization()), scatratimint_->NumScal(), splitter);
    systemmatrix = Teuchos::rcp(
        new LINALG::BlockSparseMatrix<SCATRA::SplitStrategy>(splitter, splitter, 27, false, true));
    Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrix<SCATRA::SplitStrategy>>(systemmatrix)
        ->SetNumScal(scatratimint_->NumScal());
  }

  else
  {
    // initialize standard (stabilized) system matrix (and save its graph)
    systemmatrix = Teuchos::rcp(
        new LINALG::SparseMatrix(*(scatratimint_->Discretization()->DofRowMap()), 27, false, true));
  }

  return systemmatrix;
}  // SCATRA::MeshtyingStrategyStdElch::InitSystemMatrix


/*------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   fang 02/16 |
 *------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyStdElch::InitConvCheckStrategy()
{
  if (ElchTimInt()->MacroScale())
    convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyStdMacroScaleElch(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));
  else
    convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyStdElch(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));

  return;
}  // SCATRA::MeshtyingStrategyStdElch::InitConvCheckStrategy
