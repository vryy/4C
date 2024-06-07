/*----------------------------------------------------------------------*/
/*! \file

\brief Standard solution strategy for electrochemistry problems (without meshtying)

\level 2


*----------------------------------------------------------------------*/
#include "4C_scatra_timint_meshtying_strategy_std_elch.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_scatra_utils_splitstrategy.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
ScaTra::MeshtyingStrategyStdElch::MeshtyingStrategyStdElch(ScaTra::ScaTraTimIntElch* elchtimint)
    : MeshtyingStrategyStd(elchtimint)
{
}  // ScaTra::MeshtyingStrategyStdElch::MeshtyingStrategyStdElch


/*----------------------------------------------------------------------*
 | initialize system matrix for electrochemistry problems    fang 12/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseOperator> ScaTra::MeshtyingStrategyStdElch::init_system_matrix()
    const
{
  Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix;

  if (Core::UTILS::IntegralValue<int>(*(elch_tim_int()->ElchParameterList()), "BLOCKPRECOND"))
  {
    // safety checks
    if (elch_tim_int()->EquPot() == Inpar::ElCh::equpot_undefined)
      FOUR_C_THROW("Type of closing equation for electric potential not correctly set!");
    if (elch_tim_int()->EquPot() != Inpar::ElCh::equpot_enc)
      FOUR_C_THROW("Special ELCH assemble strategy for block-matrix will not assemble A_11 block!");
    if (scatratimint_->NumScal() < 1)
      FOUR_C_THROW("Number of transported scalars not correctly set!");

    // initial guess for non-zeros per row: 27 neighboring nodes for hex8
    // this is enough! A higher guess would require too much memory!
    // A more precise guess for every submatrix would read:
    // A_00: 27*1,  A_01: 27*1,  A_10: 27*numscal due to electroneutrality, A_11: EMPTY matrix !!!!!
    // usage of a split strategy that makes use of the ELCH-specific sparsity pattern
    Core::LinAlg::MapExtractor splitter;
    Core::LinAlg::CreateMapExtractorFromDiscretization(
        *(scatratimint_->discretization()), scatratimint_->NumScal(), splitter);
    systemmatrix = Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<ScaTra::SplitStrategy>(
        splitter, splitter, 27, false, true));
    Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrix<ScaTra::SplitStrategy>>(systemmatrix)
        ->SetNumScal(scatratimint_->NumScal());
  }

  else
  {
    // initialize standard (stabilized) system matrix (and save its graph)
    switch (scatratimint_->MatrixType())
    {
      case Core::LinAlg::MatrixType::sparse:
      {
        systemmatrix = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
            *scatratimint_->discretization()->dof_row_map(), 27, false, true));
        break;
      }

      case Core::LinAlg::MatrixType::block_condition:
      case Core::LinAlg::MatrixType::block_condition_dof:
      {
        systemmatrix = Teuchos::rcp(
            new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
                *scatratimint_->BlockMaps(), *scatratimint_->BlockMaps(), 81, false, true));

        break;
      }

      default:
      {
        FOUR_C_THROW("Unknown matrix type of ScaTra field");
        break;
      }
    }
  }

  return systemmatrix;
}  // ScaTra::MeshtyingStrategyStdElch::init_system_matrix


/*------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   fang 02/16 |
 *------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyStdElch::init_conv_check_strategy()
{
  if (elch_tim_int()->MacroScale())
  {
    convcheckstrategy_ = Teuchos::rcp(new ScaTra::ConvCheckStrategyStdMacroScaleElch(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));
  }
  else
  {
    convcheckstrategy_ = Teuchos::rcp(new ScaTra::ConvCheckStrategyStdElch(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));
  }
}  // ScaTra::MeshtyingStrategyStdElch::init_conv_check_strategy

FOUR_C_NAMESPACE_CLOSE
