/*----------------------------------------------------------------------*/
/*! \file
\brief Application of contact contributions strategy for monolithic SSI

\level 2

 */
/*----------------------------------------------------------------------*/

#include "baci_ssi_monolithic_contact_strategy.H"

#include "baci_contact_nitsche_strategy_ssi.H"
#include "baci_linalg_blocksparsematrix.H"
#include "baci_ssi_monolithic.H"
#include "baci_ssi_utils.H"


/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::ContactStrategyBase::ContactStrategyBase(
    Teuchos::RCP<CONTACT::CoNitscheStrategySsi> contact_nitsche_strategy,
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps)
    : contact_strategy_nitsche_(std::move(contact_nitsche_strategy)), ssi_maps_(std::move(ssi_maps))
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::ContactStrategySparse::ContactStrategySparse(
    Teuchos::RCP<CONTACT::CoNitscheStrategySsi> contact_nitsche_strategy,
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps)
    : ContactStrategyBase(contact_nitsche_strategy, ssi_maps)
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::ContactStrategyBlock::ContactStrategyBlock(
    Teuchos::RCP<CONTACT::CoNitscheStrategySsi> contact_nitsche_strategy,
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps)
    : ContactStrategyBase(contact_nitsche_strategy, ssi_maps)
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategyBase::ApplyContactToScatraResidual(
    Teuchos::RCP<Epetra_Vector> scatra_residual)
{
  scatra_residual->Update(
      1.0, *CoNitscheStrategySsi()->GetRhsBlockPtr(DRT::UTILS::VecBlockType::scatra), 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategySparse::ApplyContactToScatraScatra(
    Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_scatra_matrix)
{
  auto scatra_scatra_matrix_sparse =
      CORE::LINALG::CastToSparseMatrixAndCheckSuccess(scatra_scatra_matrix);

  const auto& scatra_scatra_sparsematrix =
      CoNitscheStrategySsi()->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::scatra_scatra);

  scatra_scatra_matrix_sparse->Add(*scatra_scatra_sparsematrix, false, 1.0, 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategyBlock::ApplyContactToScatraScatra(
    Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_scatra_matrix)
{
  auto scatra_scatra_matrix_block =
      CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(scatra_scatra_matrix);

  // get scatra-scatra block matrix and complete split matrix
  const auto& scatra_scatra_blockmatrix =
      CoNitscheStrategySsi()
          ->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::scatra_scatra)
          ->Split<CORE::LINALG::DefaultBlockMatrixStrategy>(
              *SSIMaps()->BlockMapScaTra(), *SSIMaps()->BlockMapScaTra());
  scatra_scatra_blockmatrix->Complete();

  scatra_scatra_matrix_block->Add(*scatra_scatra_blockmatrix, false, 1.0, 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategySparse::ApplyContactToScatraStructure(
    Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_structure_matrix)
{
  auto scatra_structure_matrix_sparse =
      CORE::LINALG::CastToSparseMatrixAndCheckSuccess(scatra_structure_matrix);
  scatra_structure_matrix_sparse->UnComplete();

  const auto& scatra_struct_matrix =
      CoNitscheStrategySsi()->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::scatra_displ);

  scatra_structure_matrix_sparse->Add(*scatra_struct_matrix, false, 1.0, 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategyBlock::ApplyContactToScatraStructure(
    Teuchos::RCP<CORE::LINALG::SparseOperator> scatra_structure_matrix)
{
  auto scatra_structure_matrix_block =
      CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(scatra_structure_matrix);

  // get scatra-structure block matrix and complete split matrix
  const auto& scatra_struct_blockmatrix =
      CoNitscheStrategySsi()
          ->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::scatra_displ)
          ->Split<CORE::LINALG::DefaultBlockMatrixStrategy>(
              *SSIMaps()->BlockMapStructure(), *SSIMaps()->BlockMapScaTra());
  scatra_struct_blockmatrix->Complete();

  scatra_structure_matrix_block->Add(*scatra_struct_blockmatrix, false, 1.0, 1.0);
}


/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategySparse::ApplyContactToStructureScatra(
    Teuchos::RCP<CORE::LINALG::SparseOperator> structure_scatra_matrix)
{
  auto structure_scatra_matrix_sparse =
      CORE::LINALG::CastToSparseMatrixAndCheckSuccess(structure_scatra_matrix);
  structure_scatra_matrix_sparse->UnComplete();

  const auto& struct_scatra_matrix =
      CoNitscheStrategySsi()->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::displ_scatra);

  structure_scatra_matrix_sparse->Add(*struct_scatra_matrix, false, 1.0, 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategyBlock::ApplyContactToStructureScatra(
    Teuchos::RCP<CORE::LINALG::SparseOperator> structure_scatra_matrix)
{
  auto structure_scatra_matrix_block =
      CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(structure_scatra_matrix);

  // get structure-scatra block matrix and complete split matrix
  const auto& struct_scatra_blockmatrix =
      CoNitscheStrategySsi()
          ->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::displ_scatra)
          ->Split<CORE::LINALG::DefaultBlockMatrixStrategy>(
              *SSIMaps()->BlockMapScaTra(), *SSIMaps()->BlockMapStructure());
  struct_scatra_blockmatrix->Complete();

  structure_scatra_matrix_block->Add(*struct_scatra_blockmatrix, false, 1.0, 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
Teuchos::RCP<SSI::ContactStrategyBase> SSI::BuildContactStrategy(
    Teuchos::RCP<CONTACT::CoNitscheStrategySsi> contact_nitsche_strategy,
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, CORE::LINALG::MatrixType matrixtype_scatra)
{
  Teuchos::RCP<SSI::ContactStrategyBase> contact_strategy(Teuchos::null);

  switch (matrixtype_scatra)
  {
    case CORE::LINALG::MatrixType::block_condition:
    case CORE::LINALG::MatrixType::block_condition_dof:
    {
      contact_strategy =
          Teuchos::rcp(new SSI::ContactStrategyBlock(contact_nitsche_strategy, ssi_maps));
      break;
    }
    case CORE::LINALG::MatrixType::sparse:
    {
      contact_strategy =
          Teuchos::rcp(new SSI::ContactStrategySparse(contact_nitsche_strategy, ssi_maps));
      break;
    }

    default:
    {
      dserror("unknown matrix type of ScaTra field");
      break;
    }
  }

  return contact_strategy;
}