/*----------------------------------------------------------------------*/
/*! \file

\brief Standard solution strategy for electrochemistry problems (without meshtying)

\level 2


*----------------------------------------------------------------------*/
#include "baci_scatra_timint_meshtying_strategy_std_elch.H"

#include "baci_linalg_utils_sparse_algebra_create.H"
#include "baci_scatra_utils_splitstrategy.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyStdElch::MeshtyingStrategyStdElch(SCATRA::ScaTraTimIntElch* elchtimint)
    : MeshtyingStrategyStd(elchtimint)
{
}  // SCATRA::MeshtyingStrategyStdElch::MeshtyingStrategyStdElch


/*----------------------------------------------------------------------*
 | initialize system matrix for electrochemistry problems    fang 12/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseOperator> SCATRA::MeshtyingStrategyStdElch::InitSystemMatrix()
    const
{
  Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix;

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
    CORE::LINALG::MapExtractor splitter;
    CORE::LINALG::CreateMapExtractorFromDiscretization(
        *(scatratimint_->Discretization()), scatratimint_->NumScal(), splitter);
    systemmatrix = Teuchos::rcp(new CORE::LINALG::BlockSparseMatrix<SCATRA::SplitStrategy>(
        splitter, splitter, 27, false, true));
    Teuchos::rcp_dynamic_cast<CORE::LINALG::BlockSparseMatrix<SCATRA::SplitStrategy>>(systemmatrix)
        ->SetNumScal(scatratimint_->NumScal());
  }

  else
  {
    // initialize standard (stabilized) system matrix (and save its graph)
    switch (scatratimint_->MatrixType())
    {
      case CORE::LINALG::MatrixType::sparse:
      {
        systemmatrix = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
            *scatratimint_->Discretization()->DofRowMap(), 27, false, true));
        break;
      }

      case CORE::LINALG::MatrixType::block_condition:
      case CORE::LINALG::MatrixType::block_condition_dof:
      {
        systemmatrix = Teuchos::rcp(
            new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
                *scatratimint_->BlockMaps(), *scatratimint_->BlockMaps(), 81, false, true));

        break;
      }

      default:
      {
        dserror("Unknown matrix type of ScaTra field");
        break;
      }
    }
  }

  return systemmatrix;
}  // SCATRA::MeshtyingStrategyStdElch::InitSystemMatrix


/*------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   fang 02/16 |
 *------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyStdElch::InitConvCheckStrategy()
{
  if (ElchTimInt()->MacroScale())
  {
    convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyStdMacroScaleElch(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));
  }
  else
  {
    convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyStdElch(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));
  }
}  // SCATRA::MeshtyingStrategyStdElch::InitConvCheckStrategy
