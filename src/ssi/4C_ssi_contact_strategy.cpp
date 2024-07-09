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
      1.0, *nitsche_strategy_ssi()->get_rhs_block_ptr(CONTACT::VecBlockType::scatra), 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategySparse::apply_contact_to_scatra_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatra_scatra_matrix)
{
  auto scatra_scatra_matrix_sparse =
      Core::LinAlg::CastToSparseMatrixAndCheckSuccess(scatra_scatra_matrix);

  const auto& scatra_scatra_sparsematrix =
      nitsche_strategy_ssi()->get_matrix_block_ptr(CONTACT::MatBlockType::scatra_scatra);

  scatra_scatra_matrix_sparse->add(*scatra_scatra_sparsematrix, false, 1.0, 1.0);
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
          ->get_matrix_block_ptr(CONTACT::MatBlockType::scatra_scatra)
          ->split<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *ssi_maps()->block_map_sca_tra(), *ssi_maps()->block_map_sca_tra());
  scatra_scatra_blockmatrix->complete();

  scatra_scatra_matrix_block->add(*scatra_scatra_blockmatrix, false, 1.0, 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategySparse::apply_contact_to_scatra_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatra_structure_matrix)
{
  auto scatra_structure_matrix_sparse =
      Core::LinAlg::CastToSparseMatrixAndCheckSuccess(scatra_structure_matrix);
  scatra_structure_matrix_sparse->un_complete();

  const auto& scatra_struct_matrix =
      nitsche_strategy_ssi()->get_matrix_block_ptr(CONTACT::MatBlockType::scatra_displ);

  scatra_structure_matrix_sparse->add(*scatra_struct_matrix, false, 1.0, 1.0);
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
          ->get_matrix_block_ptr(CONTACT::MatBlockType::scatra_displ)
          ->split<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *ssi_maps()->block_map_structure(), *ssi_maps()->block_map_sca_tra());
  scatra_struct_blockmatrix->complete();

  scatra_structure_matrix_block->add(*scatra_struct_blockmatrix, false, 1.0, 1.0);
}


/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategySparse::apply_contact_to_structure_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> structure_scatra_matrix)
{
  auto structure_scatra_matrix_sparse =
      Core::LinAlg::CastToSparseMatrixAndCheckSuccess(structure_scatra_matrix);
  structure_scatra_matrix_sparse->un_complete();

  const auto& struct_scatra_matrix =
      nitsche_strategy_ssi()->get_matrix_block_ptr(CONTACT::MatBlockType::displ_scatra);

  structure_scatra_matrix_sparse->add(*struct_scatra_matrix, false, 1.0, 1.0);
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
          ->get_matrix_block_ptr(CONTACT::MatBlockType::displ_scatra)
          ->split<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *ssi_maps()->block_map_sca_tra(), *ssi_maps()->block_map_structure());
  struct_scatra_blockmatrix->complete();

  structure_scatra_matrix_block->add(*struct_scatra_blockmatrix, false, 1.0, 1.0);
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
