/*----------------------------------------------------------------------*/
/*! \file
\brief Application of contact contributions strategy for monolithic SSI

\level 2

 */
/*----------------------------------------------------------------------*/

#include "ssi_monolithic_contact_strategy.H"
#include "ssi_monolithic.H"
#include "ssi_utils.H"

#include "../drt_contact/contact_nitsche_strategy_ssi.H"

#include "../linalg/linalg_blocksparsematrix.H"


/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::ContactStrategyBase::ContactStrategyBase(const SSI::SSIMono& ssi_mono) : ssi_mono_(ssi_mono) {}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::ContactStrategySparse::ContactStrategySparse(const SSI::SSIMono& ssi_mono)
    : ContactStrategyBase(ssi_mono)
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::ContactStrategyBlock::ContactStrategyBlock(const SSI::SSIMono& ssi_mono)
    : ContactStrategyBase(ssi_mono)
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategyBase::ApplyContactToScatraResidual(
    Teuchos::RCP<Epetra_Vector> scatra_residual)
{
  scatra_residual->Update(1.0,
      *SSIMono().CoNitscheStrategySsi()->GetRhsBlockPtr(DRT::UTILS::VecBlockType::scatra), 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategySparse::ApplyContactToScatraScatra(
    Teuchos::RCP<LINALG::SparseOperator> scatra_scatra_matrix)
{
  auto scatra_scatra_matrix_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(scatra_scatra_matrix);

  const auto& scatra_scatra_sparsematrix =
      SSIMono().CoNitscheStrategySsi()->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::scatra_scatra);

  scatra_scatra_matrix_sparse->Add(*scatra_scatra_sparsematrix, false, 1.0, 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategyBlock::ApplyContactToScatraScatra(
    Teuchos::RCP<LINALG::SparseOperator> scatra_scatra_matrix)
{
  auto scatra_scatra_matrix_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(scatra_scatra_matrix);

  // get scatra-scatra block matrix and complete split matrix
  const auto& scatra_scatra_blockmatrix =
      SSIMono()
          .CoNitscheStrategySsi()
          ->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::scatra_scatra)
          ->Split<LINALG::DefaultBlockMatrixStrategy>(
              *SSIMono().MapsScatra(), *SSIMono().MapsScatra());
  scatra_scatra_blockmatrix->Complete();

  scatra_scatra_matrix_block->Add(*scatra_scatra_blockmatrix, false, 1.0, 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategySparse::ApplyContactToScatraStructure(
    Teuchos::RCP<LINALG::SparseOperator> scatra_structure_matrix)
{
  auto scatra_structure_matrix_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(scatra_structure_matrix);
  scatra_structure_matrix_sparse->UnComplete();

  const auto& scatra_struct_matrix =
      SSIMono().CoNitscheStrategySsi()->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::scatra_displ);

  scatra_structure_matrix_sparse->Add(*scatra_struct_matrix, false, 1.0, 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategyBlock::ApplyContactToScatraStructure(
    Teuchos::RCP<LINALG::SparseOperator> scatra_structure_matrix)
{
  auto scatra_structure_matrix_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(scatra_structure_matrix);

  // get scatra-structure block matrix and complete split matrix
  const auto& scatra_struct_blockmatrix =
      SSIMono()
          .CoNitscheStrategySsi()
          ->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::scatra_displ)
          ->Split<LINALG::DefaultBlockMatrixStrategy>(
              *SSIMono().MapStructure(), *SSIMono().MapsScatra());
  scatra_struct_blockmatrix->Complete();

  scatra_structure_matrix_block->Add(*scatra_struct_blockmatrix, false, 1.0, 1.0);
}


/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategySparse::ApplyContactToStructureScatra(
    Teuchos::RCP<LINALG::SparseOperator> structure_scatra_matrix)
{
  auto structure_scatra_matrix_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(structure_scatra_matrix);
  structure_scatra_matrix_sparse->UnComplete();

  const auto& struct_scatra_matrix =
      SSIMono().CoNitscheStrategySsi()->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::displ_scatra);

  structure_scatra_matrix_sparse->Add(*struct_scatra_matrix, false, 1.0, 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::ContactStrategyBlock::ApplyContactToStructureScatra(
    Teuchos::RCP<LINALG::SparseOperator> structure_scatra_matrix)
{
  auto structure_scatra_matrix_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(structure_scatra_matrix);

  // get structure-scatra block matrix and complete split matrix
  const auto& struct_scatra_blockmatrix =
      SSIMono()
          .CoNitscheStrategySsi()
          ->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::displ_scatra)
          ->Split<LINALG::DefaultBlockMatrixStrategy>(
              *SSIMono().MapsScatra(), *SSIMono().MapStructure());
  struct_scatra_blockmatrix->Complete();

  structure_scatra_matrix_block->Add(*struct_scatra_blockmatrix, false, 1.0, 1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
Teuchos::RCP<SSI::ContactStrategyBase> SSI::BuildContactStrategy(
    const SSI::SSIMono& ssi_mono, LINALG::MatrixType matrixtype_scatra)
{
  Teuchos::RCP<SSI::ContactStrategyBase> contact_strategy(Teuchos::null);

  switch (matrixtype_scatra)
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      contact_strategy = Teuchos::rcp(new SSI::ContactStrategyBlock(ssi_mono));
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      contact_strategy = Teuchos::rcp(new SSI::ContactStrategySparse(ssi_mono));
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