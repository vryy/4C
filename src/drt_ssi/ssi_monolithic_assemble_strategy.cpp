/*----------------------------------------------------------------------*/
/*! \file
\brief Assemble strategy for monolithic SSI
\level 2

 */
/*----------------------------------------------------------------------*/
#include "ssi_monolithic_assemble_strategy.H"

#include "ssi_monolithic.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"

#include "../drt_contact/contact_nitsche_strategy_ssi.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_locsys.H"

#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_matrixtransform.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategyBase::AssembleStrategyBase(const SSI::SSIMono& ssi_mono) : ssi_mono_(ssi_mono)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategyBlock::AssembleStrategyBlock(const SSI::SSIMono& ssi_mono)
    : AssembleStrategyBase(ssi_mono),
      block_position_scatra_(Teuchos::null),
      block_position_scatra_manifold_(Teuchos::null),
      position_structure_(-1)
{
  block_position_scatra_ = SSIMono().GetBlockPositions(SSI::Subproblem::scalar_transport);
  position_structure_ = SSIMono().GetBlockPositions(SSI::Subproblem::structure)->at(0);
  if (SSIMono().IsScaTraManifold())
    block_position_scatra_manifold_ = SSIMono().GetBlockPositions(SSI::Subproblem::manifold);

  if (block_position_scatra_ == Teuchos::null) dserror("Cannot get position of scatra blocks");
  if (position_structure_ == -1) dserror("Cannot get position of structure block");
  if (SSIMono().IsScaTraManifold() and block_position_scatra_manifold_ == Teuchos::null)
    dserror("Cannot get position of scatra manifold blocks");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategyBlockBlock::AssembleStrategyBlockBlock(const SSI::SSIMono& ssi_mono)
    : AssembleStrategyBlock(ssi_mono)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategyBlockSparse::AssembleStrategyBlockSparse(const SSI::SSIMono& ssi_mono)
    : AssembleStrategyBlock(ssi_mono)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategySparse::AssembleStrategySparse(const SSI::SSIMono& ssi_mono)
    : AssembleStrategyBase(ssi_mono)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScatraScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatra_scatra_matrix)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatra_scatra_matrix_block =
      LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatra_scatra_matrix);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(BlockPositionScaTra()->size()); ++jblock)
    {
      auto& systemmatrix_block_iscatra_jscatra = systemmatrix_block->Matrix(
          BlockPositionScaTra()->at(iblock), BlockPositionScaTra()->at(jblock));

      systemmatrix_block_iscatra_jscatra.Add(
          scatra_scatra_matrix_block->Matrix(iblock, jblock), false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleScatraScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatra_scatra_matrix)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatra_scatra_matrix_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(scatra_scatra_matrix);

  auto& systemmatrix_block_scatra_scatra =
      systemmatrix_block->Matrix(BlockPositionScaTra()->at(0), BlockPositionScaTra()->at(0));

  systemmatrix_block_scatra_scatra.Add(*scatra_scatra_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleScatraScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatra_scatra_matrix)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatra_scatra_matrix_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(scatra_scatra_matrix);

  // add scalar transport system matrix to global system matrix
  systemmatrix_sparse->Add(*scatra_scatra_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleStructureStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseMatrix> structure_structure_matrix)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto& systemmatrix_block_struct_struct =
      systemmatrix_block->Matrix(PositionStructure(), PositionStructure());

  systemmatrix_block_struct_struct.Add(*structure_structure_matrix, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleStructureStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseMatrix> structure_structure_matrix)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto& systemmatrix_block_struct_struct =
      systemmatrix_block->Matrix(PositionStructure(), PositionStructure());

  systemmatrix_block_struct_struct.Add(*structure_structure_matrix, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleStructureStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseMatrix> structure_structure_matrix)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  systemmatrix_sparse->Add(*structure_structure_matrix, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScatraStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatra_structure_matrix)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatra_structure_matrix_block =
      LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatra_structure_matrix);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
  {
    auto& systemmatrix_block_iscatra_struct =
        systemmatrix_block->Matrix(BlockPositionScaTra()->at(iblock), PositionStructure());

    systemmatrix_block_iscatra_struct.Add(
        scatra_structure_matrix_block->Matrix(iblock, 0), false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleScatraStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatra_structure_matrix)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatra_structure_matrix_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(scatra_structure_matrix);

  auto& systemmatrix_block_scatra_struct =
      systemmatrix_block->Matrix(BlockPositionScaTra()->at(0), PositionStructure());

  systemmatrix_block_scatra_struct.Add(*scatra_structure_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleScatraStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatra_structure_matrix)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatra_structure_matrix_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(scatra_structure_matrix);

  systemmatrix_sparse->Add(*scatra_structure_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScatraScatramanifold(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatra_scatramanifold_matrix)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatra_scatramanifold_matrix_block =
      LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatra_scatramanifold_matrix);

  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(BlockPositionScaTraManifold()->size()); ++jblock)
    {
      systemmatrix_block
          ->Matrix(BlockPositionScaTra()->at(iblock), BlockPositionScaTraManifold()->at(jblock))
          .Add(scatra_scatramanifold_matrix_block->Matrix(iblock, jblock), false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleScatraScatramanifold(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatra_scatramanifold_matrix)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatra_scatramanifold_matrix_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(scatra_scatramanifold_matrix);

  systemmatrix_block->Matrix(BlockPositionScaTra()->at(0), BlockPositionScaTraManifold()->at(0))
      .Add(*scatra_scatramanifold_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleScatraScatramanifold(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatra_scatramanifold_matrix)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatra_scatramanifold_matrix_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(scatra_scatramanifold_matrix);

  systemmatrix_sparse->Add(*scatra_scatramanifold_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleStructureScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> structure_scatra_matrix)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto structure_scatra_matrix_block =
      LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(structure_scatra_matrix);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
  {
    auto& systemmatrix_block_struct_iscatra =
        systemmatrix_block->Matrix(PositionStructure(), BlockPositionScaTra()->at(iblock));
    systemmatrix_block_struct_iscatra.Add(
        structure_scatra_matrix_block->Matrix(0, iblock), false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleStructureScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> structure_scatra_matrix)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto structure_scatra_matrix_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(structure_scatra_matrix);

  auto& systemmatrix_block_struct_scatra =
      systemmatrix_block->Matrix(PositionStructure(), BlockPositionScaTra()->at(0));
  systemmatrix_block_struct_scatra.Add(*structure_scatra_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleStructureScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> structure_scatra_matrix)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto structure_scatra_matrix_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(structure_scatra_matrix);

  systemmatrix_sparse->Add(*structure_scatra_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScatramanifoldScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_scatra_matrix)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatramanifold_scatra_matrix_block =
      LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatramanifold_scatra_matrix);

  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTraManifold()->size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(BlockPositionScaTra()->size()); ++jblock)
    {
      systemmatrix_block
          ->Matrix(BlockPositionScaTraManifold()->at(iblock), BlockPositionScaTra()->at(jblock))
          .Add(scatramanifold_scatra_matrix_block->Matrix(iblock, jblock), false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleScatramanifoldScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_scatra_matrix)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatramanifold_scatra_matrix_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(scatramanifold_scatra_matrix);

  systemmatrix_block->Matrix(BlockPositionScaTraManifold()->at(0), BlockPositionScaTra()->at(0))
      .Add(*scatramanifold_scatra_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleScatramanifoldScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_scatra_matrix)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatramanifold_scatra_matrix_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(scatramanifold_scatra_matrix);

  systemmatrix_sparse->Add(*scatramanifold_scatra_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScatramanifoldScatramanifold(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_scatramanifold_matrix)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatramanifold_scatramanifold_matrix_block =
      LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatramanifold_scatramanifold_matrix);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTraManifold()->size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(BlockPositionScaTraManifold()->size()); ++jblock)
    {
      auto& systemmatrix_block_iscatramanifold_jscatramanifold = systemmatrix_block->Matrix(
          BlockPositionScaTraManifold()->at(iblock), BlockPositionScaTraManifold()->at(jblock));

      systemmatrix_block_iscatramanifold_jscatramanifold.Add(
          scatramanifold_scatramanifold_matrix_block->Matrix(iblock, jblock), false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleScatramanifoldScatramanifold(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_scatramanifold_matrix)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatramanifold_scatramanifold_matrix_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(scatramanifold_scatramanifold_matrix);

  auto& systemmatrix_block_scatramanifold_scatramanifold = systemmatrix_block->Matrix(
      BlockPositionScaTraManifold()->at(0), BlockPositionScaTraManifold()->at(0));

  systemmatrix_block_scatramanifold_scatramanifold.Add(
      *scatramanifold_scatramanifold_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleScatramanifoldScatramanifold(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_scatramanifold_matrix)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatramanifold_scatramanifold_matrix_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(scatramanifold_scatramanifold_matrix);

  systemmatrix_sparse->Add(*scatramanifold_scatramanifold_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScatramanifoldStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_structure_matrix)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatramanifold_structure_matrix_block =
      LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatramanifold_structure_matrix);

  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTraManifold()->size()); ++iblock)
  {
    auto& systemmatrix_block_iscatramanifold_struct =
        systemmatrix_block->Matrix(BlockPositionScaTraManifold()->at(iblock), PositionStructure());
    systemmatrix_block_iscatramanifold_struct.Add(
        scatramanifold_structure_matrix_block->Matrix(iblock, 0), false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleScatramanifoldStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_structure_matrix)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatramanifold_structure_matrix_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(scatramanifold_structure_matrix);

  auto& systemmatrix_block_scatramanifold_struct =
      systemmatrix_block->Matrix(BlockPositionScaTraManifold()->at(0), PositionStructure());
  systemmatrix_block_scatramanifold_struct.Add(
      *scatramanifold_structure_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleScatramanifoldStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatramanifold_structure_matrix)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatramanifold_structure_matrix_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(scatramanifold_structure_matrix);

  systemmatrix_sparse->Add(*scatramanifold_structure_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::AssembleRHS(Teuchos::RCP<Epetra_Vector> rhs,
    Teuchos::RCP<const Epetra_Vector> rhs_scatra, Teuchos::RCP<const Epetra_Vector> rhs_structure,
    Teuchos::RCP<const Epetra_Vector> rhs_manifold,
    Teuchos::RCP<const Epetra_Vector> rhs_manifold_scatra_coupling_manifold_side,
    Teuchos::RCP<const Epetra_Vector> rhs_manifold_scatra_coupling_scatra_side)
{
  SSIMono().MapsSubProblems()->InsertVector(
      rhs_scatra, SSIMono().GetProblemPosition(SSI::Subproblem::scalar_transport), rhs);

  if (SSIMono().IsScaTraManifold())
  {
    SSIMono().MapsSubProblems()->InsertVector(
        rhs_manifold, SSIMono().GetProblemPosition(SSI::Subproblem::manifold), rhs);

    SSIMono().MapsSubProblems()->InsertVector(
        rhs_manifold, SSIMono().GetProblemPosition(SSI::Subproblem::manifold), rhs);

    SSIMono().MapsSubProblems()->AddVector(rhs_manifold_scatra_coupling_manifold_side,
        SSIMono().GetProblemPosition(Subproblem::manifold), rhs);

    SSIMono().MapsSubProblems()->AddVector(rhs_manifold_scatra_coupling_scatra_side,
        SSIMono().GetProblemPosition(Subproblem::scalar_transport), rhs);
  }

  SSIMono().MapsSubProblems()->AddVector(
      rhs_structure, SSIMono().GetProblemPosition(SSI::Subproblem::structure), rhs, -1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
Teuchos::RCP<SSI::AssembleStrategyBase> SSI::BuildAssembleStrategy(const SSI::SSIMono& ssi_mono,
    LINALG::MatrixType matrixtype_ssi, LINALG::MatrixType matrixtype_scatra)
{
  Teuchos::RCP<SSI::AssembleStrategyBase> assemblestrategy = Teuchos::null;

  switch (matrixtype_ssi)
  {
    case LINALG::MatrixType::block_field:
    {
      switch (matrixtype_scatra)
      {
        case LINALG::MatrixType::block_condition:
        case LINALG::MatrixType::block_condition_dof:
        {
          assemblestrategy = Teuchos::rcp(new SSI::AssembleStrategyBlockBlock(ssi_mono));
          break;
        }
        case LINALG::MatrixType::sparse:
        {
          assemblestrategy = Teuchos::rcp(new SSI::AssembleStrategyBlockSparse(ssi_mono));
          break;
        }

        default:
        {
          dserror("unknown matrix type of ScaTra field");
          break;
        }
      }
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      assemblestrategy = Teuchos::rcp(new SSI::AssembleStrategySparse(ssi_mono));
      break;
    }
    default:
    {
      dserror("unknown matrix type of SSI problem");
      break;
    }
  }

  return assemblestrategy;
}