/*----------------------------------------------------------------------*/
/*! \file
\brief Application of contact contributions strategy for monolithic/partitioning SSI

\level 2

 */
/*----------------------------------------------------------------------*/

#include "4C_ssi_contact_strategy.hpp"

#include "4C_contact_nitsche_strategy_ssi.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_ssi_utils.hpp"

FOUR_C_NAMESPACE_OPEN

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::ContactStrategyBase::ContactStrategyBase(
    Teuchos::RCP<CONTACT::NitscheStrategySsi> contact_nitsche_strategy,
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps)
    : contact_strategy_nitsche_(std::move(contact_nitsche_strategy)), ssi_maps_(std::move(ssi_maps))
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::ContactStrategySparse::ContactStrategySparse(
    Teuchos::RCP<CONTACT::NitscheStrategySsi> contact_nitsche_strategy,
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps)
    : ContactStrategyBase(contact_nitsche_strategy, ssi_maps)
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::ContactStrategyBlock::ContactStrategyBlock(
    Teuchos::RCP<CONTACT::NitscheStrategySsi> contact_nitsche_strategy,
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps)
    : ContactStrategyBase(contact_nitsche_strategy, ssi_maps)
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategyBase::apply_contact_to_scatra_residual(
    Teuchos::RCP<Epetra_Vector> scatra_residual)
{
  scatra_residual->Update(
      1.0, *nitsche_strategy_ssi()->GetRhsBlockPtr(CONTACT::VecBlockType::scatra), 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategySparse::apply_contact_to_scatra_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatra_scatra_matrix)
{
  auto scatra_scatra_matrix_sparse =
      Core::LinAlg::CastToSparseMatrixAndCheckSuccess(scatra_scatra_matrix);

  const auto& scatra_scatra_sparsematrix =
      nitsche_strategy_ssi()->GetMatrixBlockPtr(CONTACT::MatBlockType::scatra_scatra);

  scatra_scatra_matrix_sparse->Add(*scatra_scatra_sparsematrix, false, 1.0, 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategyBlock::apply_contact_to_scatra_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatra_scatra_matrix)
{
  auto scatra_scatra_matrix_block =
      Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(scatra_scatra_matrix);

  // get scatra-scatra block matrix and complete split matrix
  const auto& scatra_scatra_blockmatrix =
      nitsche_strategy_ssi()
          ->GetMatrixBlockPtr(CONTACT::MatBlockType::scatra_scatra)
          ->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *ssi_maps()->BlockMapScaTra(), *ssi_maps()->BlockMapScaTra());
  scatra_scatra_blockmatrix->Complete();

  scatra_scatra_matrix_block->Add(*scatra_scatra_blockmatrix, false, 1.0, 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategySparse::apply_contact_to_scatra_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatra_structure_matrix)
{
  auto scatra_structure_matrix_sparse =
      Core::LinAlg::CastToSparseMatrixAndCheckSuccess(scatra_structure_matrix);
  scatra_structure_matrix_sparse->UnComplete();

  const auto& scatra_struct_matrix =
      nitsche_strategy_ssi()->GetMatrixBlockPtr(CONTACT::MatBlockType::scatra_displ);

  scatra_structure_matrix_sparse->Add(*scatra_struct_matrix, false, 1.0, 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategyBlock::apply_contact_to_scatra_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatra_structure_matrix)
{
  auto scatra_structure_matrix_block =
      Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(scatra_structure_matrix);

  // get scatra-structure block matrix and complete split matrix
  const auto& scatra_struct_blockmatrix =
      nitsche_strategy_ssi()
          ->GetMatrixBlockPtr(CONTACT::MatBlockType::scatra_displ)
          ->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *ssi_maps()->BlockMapStructure(), *ssi_maps()->BlockMapScaTra());
  scatra_struct_blockmatrix->Complete();

  scatra_structure_matrix_block->Add(*scatra_struct_blockmatrix, false, 1.0, 1.0);
}


/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategySparse::apply_contact_to_structure_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> structure_scatra_matrix)
{
  auto structure_scatra_matrix_sparse =
      Core::LinAlg::CastToSparseMatrixAndCheckSuccess(structure_scatra_matrix);
  structure_scatra_matrix_sparse->UnComplete();

  const auto& struct_scatra_matrix =
      nitsche_strategy_ssi()->GetMatrixBlockPtr(CONTACT::MatBlockType::displ_scatra);

  structure_scatra_matrix_sparse->Add(*struct_scatra_matrix, false, 1.0, 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategyBlock::apply_contact_to_structure_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> structure_scatra_matrix)
{
  auto structure_scatra_matrix_block =
      Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(structure_scatra_matrix);

  // get structure-scatra block matrix and complete split matrix
  const auto& struct_scatra_blockmatrix =
      nitsche_strategy_ssi()
          ->GetMatrixBlockPtr(CONTACT::MatBlockType::displ_scatra)
          ->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *ssi_maps()->BlockMapScaTra(), *ssi_maps()->BlockMapStructure());
  struct_scatra_blockmatrix->Complete();

  structure_scatra_matrix_block->Add(*struct_scatra_blockmatrix, false, 1.0, 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
Teuchos::RCP<SSI::ContactStrategyBase> SSI::BuildContactStrategy(
    Teuchos::RCP<CONTACT::NitscheStrategySsi> contact_nitsche_strategy,
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, Core::LinAlg::MatrixType matrixtype_scatra)
{
  Teuchos::RCP<SSI::ContactStrategyBase> contact_strategy(Teuchos::null);

  switch (matrixtype_scatra)
  {
    case Core::LinAlg::MatrixType::block_condition:
    case Core::LinAlg::MatrixType::block_condition_dof:
    {
      contact_strategy =
          Teuchos::rcp(new SSI::ContactStrategyBlock(contact_nitsche_strategy, ssi_maps));
      break;
    }
    case Core::LinAlg::MatrixType::sparse:
    {
      contact_strategy =
          Teuchos::rcp(new SSI::ContactStrategySparse(contact_nitsche_strategy, ssi_maps));
      break;
    }

    default:
    {
      FOUR_C_THROW("unknown matrix type of ScaTra field");
      break;
    }
  }

  return contact_strategy;
}

FOUR_C_NAMESPACE_CLOSE
