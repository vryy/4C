/*----------------------------------------------------------------------*/
/*! \file

\brief Standard solution strategy for electrochemistry problems (without meshtying)

\level 2


*----------------------------------------------------------------------*/
#include "4C_scatra_timint_meshtying_strategy_std_elch.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_scatra_utils_splitstrategy.hpp"

#include <Teuchos_ParameterList.hpp>

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

  // initialize standard (stabilized) system matrix (and save its graph)
  switch (scatratimint_->matrix_type())
  {
    case Core::LinAlg::MatrixType::sparse:
    {
      systemmatrix = Teuchos::make_rcp<Core::LinAlg::SparseMatrix>(
          *scatratimint_->discretization()->dof_row_map(), 27, false, true);
      break;
    }

    case Core::LinAlg::MatrixType::block_condition:
    case Core::LinAlg::MatrixType::block_condition_dof:
    {
      systemmatrix = Teuchos::make_rcp<
          Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>>(

          *scatratimint_->block_maps(), *scatratimint_->block_maps(), 81, false, true);

      break;
    }

    default:
    {
      FOUR_C_THROW("Unknown matrix type of ScaTra field");
      break;
    }
  }

  return systemmatrix;
}  // ScaTra::MeshtyingStrategyStdElch::init_system_matrix


/*------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   fang 02/16 |
 *------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyStdElch::init_conv_check_strategy()
{
  if (elch_tim_int()->macro_scale())
  {
    convcheckstrategy_ = Teuchos::make_rcp<ScaTra::ConvCheckStrategyStdMacroScaleElch>(
        scatratimint_->scatra_parameter_list()->sublist("NONLINEAR"));
  }
  else
  {
    convcheckstrategy_ = Teuchos::make_rcp<ScaTra::ConvCheckStrategyStdElch>(
        scatratimint_->scatra_parameter_list()->sublist("NONLINEAR"));
  }
}  // ScaTra::MeshtyingStrategyStdElch::init_conv_check_strategy

FOUR_C_NAMESPACE_CLOSE
