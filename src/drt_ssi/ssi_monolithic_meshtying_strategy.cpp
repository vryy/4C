/*----------------------------------------------------------------------*/
/*! \file
\brief Mesh tying strategy for monolithic SSI

\level 2

 */
/*----------------------------------------------------------------------*/

#include "ssi_monolithic_meshtying_strategy.H"

#include "ssi_monolithic.H"
#include "ssi_utils.H"
#include "Epetra_Map.h"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/adapter_coupling.H"

#include "../drt_lib/drt_locsys.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_matrixtransform.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::MeshtyingStrategyBase::MeshtyingStrategyBase(const SSI::SSIMono& ssi_mono)
    : temp_scatra_struct_domain_mat_(Teuchos::null),
      temp_scatra_struct_interface_mat_(Teuchos::null),
      temp_scatramanifold_struct_mat_(Teuchos::null),
      temp_struct_scatra_mat_(Teuchos::null),
      mapstructurecondensed_(
          ssi_mono.SSIInterfaceMeshtying() ? ssi_mono.MapStructureCondensed() : Teuchos::null),
      mapstructureslave_(
          ssi_mono.SSIInterfaceMeshtying() ? ssi_mono.MapsCoupStruct()->Map(1) : Teuchos::null),
      mapstructureslave3domainintersection_(
          (ssi_mono.SSIInterfaceMeshtying() and ssi_mono.Meshtying3DomainIntersection())
              ? ssi_mono.MapsCoupStruct3DomainIntersection()->Map(1)
              : Teuchos::null),
      meshtying_3_domain_intersection_(
          ssi_mono.SSIInterfaceMeshtying() and ssi_mono.Meshtying3DomainIntersection()),
      slave_side_converter_(
          ssi_mono.SSIInterfaceMeshtying() ? ssi_mono.SlaveSideConverter() : Teuchos::null),
      ssi_mono_(ssi_mono)
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::MeshtyingStrategySparse::MeshtyingStrategySparse(
    const SSI::SSIMono& ssi_mono, Teuchos::RCP<const Epetra_Map> interface_map_scatra)
    : MeshtyingStrategyBase(ssi_mono), interface_map_scatra_(interface_map_scatra)
{
  temp_scatra_struct_domain_mat_ =
      SSI::UTILS::SSIMatrices::SetupSparseMatrix(ssi_mono.ScaTraField()->DofRowMap());

  if (ssi_mono.SSIInterfaceMeshtying())
    temp_scatra_struct_interface_mat_ =
        SSI::UTILS::SSIMatrices::SetupSparseMatrix(interface_map_scatra);

  if (ssi_mono.IsScaTraManifold())
    temp_scatramanifold_struct_mat_ =
        SSI::UTILS::SSIMatrices::SetupSparseMatrix(ssi_mono.ScaTraManifold()->DofRowMap());

  temp_struct_scatra_mat_ =
      SSI::UTILS::SSIMatrices::SetupSparseMatrix(ssi_mono.StructureField()->DofRowMap());
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::MeshtyingStrategyBlock::MeshtyingStrategyBlock(const SSI::SSIMono& ssi_mono)
    : MeshtyingStrategyBase(ssi_mono),
      block_position_scatra_(Teuchos::null),
      block_position_scatra_manifold_(Teuchos::null),
      position_structure_(-1)
{
  block_position_scatra_ = SSIMono().GetBlockPositions(SSI::Subproblem::scalar_transport);
  position_structure_ = SSIMono().GetBlockPositions(SSI::Subproblem::structure)->at(0);
  if (ssi_mono.IsScaTraManifold())
    block_position_scatra_manifold_ = ssi_mono.GetBlockPositions(SSI::Subproblem::manifold);

  // safety checks
  if (block_position_scatra_ == Teuchos::null) dserror("Cannot get position of scatra blocks");
  if (position_structure_ == -1) dserror("Cannot get position of structure block");
  if (ssi_mono.IsScaTraManifold() and block_position_scatra_manifold_ == Teuchos::null)
    dserror("Cannot get position of scatra manifold blocks");
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::MeshtyingStrategyBlockSparse::MeshtyingStrategyBlockSparse(
    const SSI::SSIMono& ssi_mono, Teuchos::RCP<const Epetra_Map> interface_map_scatra)
    : MeshtyingStrategyBlock(ssi_mono), interface_map_scatra_(interface_map_scatra)
{
  temp_scatra_struct_domain_mat_ =
      SSI::UTILS::SSIMatrices::SetupSparseMatrix(ssi_mono.ScaTraField()->DofRowMap());

  if (ssi_mono.SSIInterfaceMeshtying())
    temp_scatra_struct_interface_mat_ =
        SSI::UTILS::SSIMatrices::SetupSparseMatrix(interface_map_scatra);

  if (ssi_mono.IsScaTraManifold())
    temp_scatramanifold_struct_mat_ =
        SSI::UTILS::SSIMatrices::SetupSparseMatrix(ssi_mono.ScaTraManifold()->DofRowMap());

  temp_struct_scatra_mat_ =
      SSI::UTILS::SSIMatrices::SetupSparseMatrix(ssi_mono.StructureField()->DofRowMap());
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::MeshtyingStrategyBlockBlock::MeshtyingStrategyBlockBlock(
    const SSI::SSIMono& ssi_mono, Teuchos::RCP<const Epetra_Map> interface_map_scatra)
    : MeshtyingStrategyBlock(ssi_mono)
{
  temp_scatra_struct_domain_mat_ = SSI::UTILS::SSIMatrices::SetupBlockMatrix(
      Teuchos::rcpFromRef(ssi_mono.ScaTraField()->BlockMaps()), ssi_mono.MapStructure());

  if (ssi_mono.SSIInterfaceMeshtying())
  {
    const auto block_map_scatra_interface =
        SSI::UTILS::SSIMatrices::GetScaTraInterfaceBlockMap(ssi_mono, interface_map_scatra);

    temp_scatra_struct_interface_mat_ = SSI::UTILS::SSIMatrices::SetupBlockMatrix(
        block_map_scatra_interface, ssi_mono.MapStructure());
  }

  if (ssi_mono.IsScaTraManifold())
  {
    // structure dofs on manifold discretization
    const auto map_structure_manifold =
        Teuchos::rcp(new LINALG::MultiMapExtractor(*ssi_mono.MapStructureOnScaTraManifold()->Map(0),
            std::vector<Teuchos::RCP<const Epetra_Map>>(
                1, ssi_mono.MapStructureOnScaTraManifold()->Map(0))));

    temp_scatramanifold_struct_mat_ = SSI::UTILS::SSIMatrices::SetupBlockMatrix(
        Teuchos::rcpFromRef(ssi_mono.ScaTraManifold()->BlockMaps()), map_structure_manifold);
  }

  temp_struct_scatra_mat_ =
      SSI::UTILS::SSIMatrices::SetupBlockMatrix(ssi_mono.MapStructure(), ssi_mono.MapsScatra());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::MeshtyingStrategyBase::ApplyMeshtyingToStructureMatrix(
    LINALG::SparseMatrix& ssi_structure_matrix,
    Teuchos::RCP<const LINALG::SparseMatrix> structure_matrix)
{
  /* Transform and assemble the structure matrix into the ssi structure matrix block by block:
   * S_m: structure interior and master side dofs
   * S_ss: structure slave surface dofs
   * S_sl: structure slave line dofs
   *
   *       S_m  S_ss  S_sl
   *       --------------
   * S_m  |  a |  b |  c |
   * S_ss |  d |  e |  f |
   * S_sl |  g |  h |  i |
   *       --------------
   */
  // assemble derivatives of interior & master dofs w.r.t. interior & master dofs (block a)
  LINALG::MatrixLogicalSplitAndTransform()(*structure_matrix, *MapStructureCondensed(),
      *MapStructureCondensed(), 1.0, nullptr, nullptr, ssi_structure_matrix, true, true);

  // assemble derivatives of surface slave dofs w.r.t. master & interior dofs (block d)
  LINALG::MatrixLogicalSplitAndTransform()(*structure_matrix, *MapStructureSlave(),
      *MapStructureCondensed(), 1.0, &StructureSlaveConverter(), nullptr, ssi_structure_matrix,
      true, true);

  // assemble derivatives of master & interior w.r.t. surface slave dofs (block b)
  LINALG::MatrixLogicalSplitAndTransform()(*structure_matrix, *MapStructureCondensed(),
      *MapStructureSlave(), 1.0, nullptr, &StructureSlaveConverter(), ssi_structure_matrix, true,
      true);

  // assemble derivatives of surface slave dofs w.r.t. surface slave dofs (block e)
  LINALG::MatrixLogicalSplitAndTransform()(*structure_matrix, *MapStructureSlave(),
      *MapStructureSlave(), 1.0, &StructureSlaveConverter(), &StructureSlaveConverter(),
      ssi_structure_matrix, true, true);

  if (Meshtying3DomainIntersection())
  {
    // assemble derivatives of line slave dofs w.r.t. master & interior (block g)
    LINALG::MatrixLogicalSplitAndTransform()(*structure_matrix,
        *MapStructureSlave3DomainIntersection(), *MapStructureCondensed(), 1.0,
        &StructureSlaveConverter3DomainIntersection(), nullptr, ssi_structure_matrix, true, true);

    // assemble derivatives of master & interior w.r.t. line slave dofs (block c)
    LINALG::MatrixLogicalSplitAndTransform()(*structure_matrix, *MapStructureCondensed(),
        *MapStructureSlave3DomainIntersection(), 1.0, nullptr,
        &StructureSlaveConverter3DomainIntersection(), ssi_structure_matrix, true, true);

    // assemble derivatives of line slave dof w.r.t. line slave dofs (block i)
    LINALG::MatrixLogicalSplitAndTransform()(*structure_matrix,
        *MapStructureSlave3DomainIntersection(), *MapStructureSlave3DomainIntersection(), 1.0,
        &StructureSlaveConverter3DomainIntersection(),
        &StructureSlaveConverter3DomainIntersection(), ssi_structure_matrix, true, true);

    // assemble derivatives of surface slave dofs w.r.t. line slave dofs (block f)
    LINALG::MatrixLogicalSplitAndTransform()(*structure_matrix, *MapStructureSlave(),
        *MapStructureSlave3DomainIntersection(), 1.0, &StructureSlaveConverter(),
        &StructureSlaveConverter3DomainIntersection(), ssi_structure_matrix, true, true);

    // assemble derivatives of line slave dofs w.r.t. surface slave dofs (block h)
    LINALG::MatrixLogicalSplitAndTransform()(*structure_matrix,
        *MapStructureSlave3DomainIntersection(), *MapStructureSlave(), 1.0,
        &StructureSlaveConverter3DomainIntersection(), &StructureSlaveConverter(),
        ssi_structure_matrix, true, true);
  }

  FinalizeMeshtyingStructureMatrix(ssi_structure_matrix);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::MeshtyingStrategyBase::ApplyMeshtyingToXXXStructure(
    LINALG::SparseMatrix& ssi_xxx_structure_matrix,
    const LINALG::SparseMatrix& xxx_structure_matrix)
{
  // assemble derivatives of x w.r.t. structure master & interior dofs
  LINALG::MatrixLogicalSplitAndTransform()(xxx_structure_matrix, xxx_structure_matrix.RangeMap(),
      *MapStructureCondensed(), 1.0, nullptr, nullptr, ssi_xxx_structure_matrix, true, true);

  // assemble derivatives of x w.r.t. structure surface slave dofs
  LINALG::MatrixLogicalSplitAndTransform()(xxx_structure_matrix, xxx_structure_matrix.RangeMap(),
      *MapStructureSlave(), 1.0, nullptr, &StructureSlaveConverter(), ssi_xxx_structure_matrix,
      true, true);

  if (Meshtying3DomainIntersection())
  {
    // assemble derivatives of x w.r.t. structure line slave dofs
    LINALG::MatrixLogicalSplitAndTransform()(xxx_structure_matrix, xxx_structure_matrix.RangeMap(),
        *MapStructureSlave3DomainIntersection(), 1.0, nullptr,
        &StructureSlaveConverter3DomainIntersection(), ssi_xxx_structure_matrix, true, true);
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
Epetra_Vector SSI::MeshtyingStrategyBase::ApplyMeshtyingToStructureRHS(
    Teuchos::RCP<const Epetra_Vector> structure_rhs)
{
  // make copy of structure right-hand side vector
  Epetra_Vector rhs_structure(*structure_rhs);

  // transform slave-side part of structure right-hand side vector to master side
  const auto rhs_structure_only_slave_dofs =
      SSIMono().MapsCoupStruct()->ExtractVector(rhs_structure, 1);

  const auto rhs_structure_only_master_dofs =
      SSIMono().InterfaceCouplingAdapterStructure()->SlaveToMaster(rhs_structure_only_slave_dofs);

  auto rhs_structure_master =
      SSIMono().MapsCoupStruct()->InsertVector(rhs_structure_only_master_dofs, 2);

  if (Meshtying3DomainIntersection())
  {
    const auto rhs_structure_3_domain_intersection_only_slave_dofs =
        SSIMono().MapsCoupStruct3DomainIntersection()->ExtractVector(rhs_structure, 1);

    const auto rhs_structure_3_domain_intersection_only_master_dofs =
        SSIMono().InterfaceCouplingAdapterStructure3DomainIntersection()->SlaveToMaster(
            rhs_structure_3_domain_intersection_only_slave_dofs);

    const auto rhs_structure_3_domain_intersection_master =
        SSIMono().MapsCoupStruct3DomainIntersection()->InsertVector(
            rhs_structure_3_domain_intersection_only_master_dofs, 2);

    rhs_structure_master->Update(1.0, *rhs_structure_3_domain_intersection_master, 1.0);
  }

  // assemble master-side part of structure right-hand side vector
  rhs_structure.Update(1.0, *rhs_structure_master, 1.0);

  // zero out slave-side part of structure right-hand side vector
  SSIMono().MapsCoupStruct()->PutScalar(rhs_structure, 1, 0.0);
  if (Meshtying3DomainIntersection())
    SSIMono().MapsCoupStruct3DomainIntersection()->PutScalar(rhs_structure, 1, 0.0);

  return rhs_structure;
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::MeshtyingStrategySparse::ApplyMeshtyingToScatraManifoldStructure(
    Teuchos::RCP<LINALG::SparseOperator> manifold_structure_matrix)
{
  temp_scatramanifold_struct_mat_->Zero();
  auto temp_scatramanifold_struct_sparse_matrix =
      LINALG::CastToSparseMatrixAndCheckSuccess(temp_scatramanifold_struct_mat_);
  auto manifold_structure_sparse_matrix =
      LINALG::CastToSparseMatrixAndCheckSuccess(manifold_structure_matrix);

  // apply mesh tying contributions to temp matrix and complete it
  ApplyMeshtyingToXXXStructure(
      *temp_scatramanifold_struct_sparse_matrix, *manifold_structure_sparse_matrix);
  temp_scatramanifold_struct_sparse_matrix->Complete(
      *SSIMono().StructureField()->DofRowMap(), *SSIMono().ScaTraManifold()->DofRowMap());

  // uncomplete matrix, add mesh tying entries stored in temp matrix to matrix and complete again
  manifold_structure_sparse_matrix->UnComplete();
  manifold_structure_sparse_matrix->Add(*temp_scatramanifold_struct_sparse_matrix, false, 1.0, 0.0);
  manifold_structure_sparse_matrix->Complete(
      *SSIMono().StructureField()->DofRowMap(), *SSIMono().ScaTraManifold()->DofRowMap());
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::MeshtyingStrategyBlockSparse::ApplyMeshtyingToScatraManifoldStructure(
    Teuchos::RCP<LINALG::SparseOperator> manifold_structure_matrix)
{
  temp_scatramanifold_struct_mat_->Zero();
  auto temp_scatramanifold_struct_sparse_matrix =
      LINALG::CastToSparseMatrixAndCheckSuccess(temp_scatramanifold_struct_mat_);
  auto manifold_structure_sparse_matrix =
      LINALG::CastToSparseMatrixAndCheckSuccess(manifold_structure_matrix);

  // apply mesh tying contributions to temp matrix and complete it
  ApplyMeshtyingToXXXStructure(
      *temp_scatramanifold_struct_sparse_matrix, *manifold_structure_sparse_matrix);
  temp_scatramanifold_struct_sparse_matrix->Complete(
      *SSIMono().StructureField()->DofRowMap(), *SSIMono().ScaTraManifold()->DofRowMap());

  // uncomplete matrix, add mesh tying entries stored in temp matrix to matrix and complete again
  manifold_structure_sparse_matrix->UnComplete();
  manifold_structure_sparse_matrix->Add(*temp_scatramanifold_struct_sparse_matrix, false, 1.0, 0.0);
  manifold_structure_sparse_matrix->Complete(
      *SSIMono().StructureField()->DofRowMap(), *SSIMono().ScaTraManifold()->DofRowMap());
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::MeshtyingStrategyBlockBlock::ApplyMeshtyingToScatraManifoldStructure(
    Teuchos::RCP<LINALG::SparseOperator> manifold_structure_matrix)
{
  temp_scatramanifold_struct_mat_->Zero();
  auto temp_scatramanifold_struct_block_matrix =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(temp_scatramanifold_struct_mat_);
  auto manifold_structure_matrix_block_matrix =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(manifold_structure_matrix);

  // apply mesh tying contributions to temp matrix blocks and complete the resulting matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTraManifold()->size()); ++iblock)
  {
    ApplyMeshtyingToXXXStructure(temp_scatramanifold_struct_block_matrix->Matrix(iblock, 0),
        manifold_structure_matrix_block_matrix->Matrix(iblock, 0));
  }
  temp_scatramanifold_struct_block_matrix->Complete();

  // uncomplete matrix, add mesh tying entries stored in temp matrix to matrix and complete again
  manifold_structure_matrix_block_matrix->UnComplete();
  manifold_structure_matrix_block_matrix->Add(
      *temp_scatramanifold_struct_block_matrix, false, 1.0, 0.0);
  manifold_structure_matrix_block_matrix->Complete();
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::MeshtyingStrategySparse::ApplyMeshtyingToScatraStructure(
    Teuchos::RCP<LINALG::SparseOperator> scatra_structure_domain,
    Teuchos::RCP<LINALG::SparseOperator> scatra_structure_interface)
{
  // apply mesh tying to the domain contributions
  if (scatra_structure_domain != Teuchos::null)
  {
    temp_scatra_struct_domain_mat_->Zero();
    auto temp_scatra_struct_domain_sparse_matrix =
        LINALG::CastToSparseMatrixAndCheckSuccess(temp_scatra_struct_domain_mat_);
    auto scatra_structure_domain_sparse =
        LINALG::CastToSparseMatrixAndCheckSuccess(scatra_structure_domain);

    ApplyMeshtyingToXXXStructure(
        *temp_scatra_struct_domain_sparse_matrix, *scatra_structure_domain_sparse);
    temp_scatra_struct_domain_sparse_matrix->Complete(
        *SSIMono().StructureField()->DofRowMap(), *SSIMono().ScaTraField()->DofRowMap());

    // uncomplete matrix, add mesh tying entries stored in temp matrix to matrix and complete again
    scatra_structure_domain_sparse->UnComplete();
    scatra_structure_domain_sparse->Add(*temp_scatra_struct_domain_sparse_matrix, false, 1.0, 0.0);
    scatra_structure_domain_sparse->Complete(
        *SSIMono().StructureField()->DofRowMap(), *SSIMono().ScaTraField()->DofRowMap());
  }

  // apply mesh tying to the interface contributions
  if (SSIMono().SSIInterfaceMeshtying())
  {
    temp_scatra_struct_interface_mat_->Zero();
    auto temp_scatra_struct_interface_sparse_matrix =
        LINALG::CastToSparseMatrixAndCheckSuccess(temp_scatra_struct_interface_mat_);
    auto scatra_structure_interface_sparse =
        LINALG::CastToSparseMatrixAndCheckSuccess(scatra_structure_interface);

    ApplyMeshtyingToXXXStructure(
        *temp_scatra_struct_interface_sparse_matrix, *scatra_structure_interface_sparse);
    temp_scatra_struct_interface_sparse_matrix->Complete(
        *SSIMono().StructureField()->DofRowMap(), *interface_map_scatra_);

    // uncomplete matrix, add mesh tying entries stored in temp matrix to matrix and complete again
    scatra_structure_interface_sparse->UnComplete();
    scatra_structure_interface_sparse->Add(
        *temp_scatra_struct_interface_sparse_matrix, false, 1.0, 0.0);
    scatra_structure_interface_sparse->Complete(
        *SSIMono().StructureField()->DofRowMap(), *interface_map_scatra_);
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::MeshtyingStrategyBlockSparse::ApplyMeshtyingToScatraStructure(
    Teuchos::RCP<LINALG::SparseOperator> scatra_structure_domain,
    Teuchos::RCP<LINALG::SparseOperator> scatra_structure_interface)
{
  // apply mesh tying to the domain contributions
  if (scatra_structure_domain != Teuchos::null)
  {
    temp_scatra_struct_domain_mat_->Zero();
    auto temp_scatra_struct_domain_sparse_matrix =
        LINALG::CastToSparseMatrixAndCheckSuccess(temp_scatra_struct_domain_mat_);
    auto scatra_structure_domain_sparse =
        LINALG::CastToSparseMatrixAndCheckSuccess(scatra_structure_domain);

    ApplyMeshtyingToXXXStructure(
        *temp_scatra_struct_domain_sparse_matrix, *scatra_structure_domain_sparse);
    temp_scatra_struct_domain_sparse_matrix->Complete(
        *SSIMono().StructureField()->DofRowMap(), *SSIMono().ScaTraField()->DofRowMap());

    // uncomplete matrix, add mesh tying entries stored in temp matrix to matrix and complete again
    scatra_structure_domain_sparse->UnComplete();
    scatra_structure_domain_sparse->Add(*temp_scatra_struct_domain_sparse_matrix, false, 1.0, 0.0);
    scatra_structure_domain_sparse->Complete(
        *SSIMono().StructureField()->DofRowMap(), *SSIMono().ScaTraField()->DofRowMap());
  }

  // apply mesh tying to the interface contributions
  if (SSIMono().SSIInterfaceMeshtying())
  {
    temp_scatra_struct_interface_mat_->Zero();
    auto temp_scatra_struct_interface_sparse_matrix =
        LINALG::CastToSparseMatrixAndCheckSuccess(temp_scatra_struct_interface_mat_);
    auto scatra_structure_interface_sparse =
        LINALG::CastToSparseMatrixAndCheckSuccess(scatra_structure_interface);

    ApplyMeshtyingToXXXStructure(
        *temp_scatra_struct_interface_sparse_matrix, *scatra_structure_interface_sparse);
    temp_scatra_struct_interface_sparse_matrix->Complete(
        *SSIMono().StructureField()->DofRowMap(), *interface_map_scatra_);

    // uncomplete matrix, add mesh tying entries stored in temp matrix to matrix and complete again
    scatra_structure_interface_sparse->UnComplete();
    scatra_structure_interface_sparse->Add(
        *temp_scatra_struct_interface_sparse_matrix, false, 1.0, 0.0);
    scatra_structure_interface_sparse->Complete(
        *SSIMono().StructureField()->DofRowMap(), *interface_map_scatra_);
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::MeshtyingStrategyBlockBlock::ApplyMeshtyingToScatraStructure(
    Teuchos::RCP<LINALG::SparseOperator> scatra_structure_domain,
    Teuchos::RCP<LINALG::SparseOperator> scatra_structure_interface)
{
  // apply mesh tying to the domain contributions
  if (scatra_structure_domain != Teuchos::null)
  {
    temp_scatra_struct_domain_mat_->Zero();
    auto temp_scatra_struct_domain_block_sparse_matrix =
        LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(temp_scatra_struct_domain_mat_);
    auto scatra_structure_domain_block_sparse =
        LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(scatra_structure_domain);

    // apply mesh tying for all blocks
    for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
    {
      ApplyMeshtyingToXXXStructure(temp_scatra_struct_domain_block_sparse_matrix->Matrix(iblock, 0),
          scatra_structure_domain_block_sparse->Matrix(iblock, 0));
    }
    temp_scatra_struct_domain_block_sparse_matrix->Complete();

    // uncomplete matrix, add mesh tying entries stored in temp matrix to matrix and complete again
    scatra_structure_domain_block_sparse->UnComplete();
    scatra_structure_domain_block_sparse->Add(
        *temp_scatra_struct_domain_block_sparse_matrix, false, 1.0, 0.0);
    scatra_structure_domain_block_sparse->Complete();
  }

  // apply mesh tying to the interface contributions
  if (SSIMono().SSIInterfaceMeshtying())
  {
    temp_scatra_struct_interface_mat_->Zero();
    auto temp_scatra_structure_interface_block_sparse_matrix =
        LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(temp_scatra_struct_interface_mat_);
    auto scatra_structure_interface_block_sparse =
        LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(scatra_structure_interface);

    // apply mesh tying for all blocks
    for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
    {
      ApplyMeshtyingToXXXStructure(
          temp_scatra_structure_interface_block_sparse_matrix->Matrix(iblock, 0),
          scatra_structure_interface_block_sparse->Matrix(iblock, 0));
    }
    temp_scatra_structure_interface_block_sparse_matrix->Complete();

    // uncomplete matrix, add mesh tying entries stored in temp matrix to matrix and complete again
    scatra_structure_interface_block_sparse->UnComplete();
    scatra_structure_interface_block_sparse->Add(
        *temp_scatra_structure_interface_block_sparse_matrix, false, 1.0, 0.0);
    scatra_structure_interface_block_sparse->Complete();
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::MeshtyingStrategySparse::ApplyMeshtyingToStructureScatra(
    Teuchos::RCP<LINALG::SparseOperator> structure_scatra_matrix)
{
  temp_struct_scatra_mat_->Zero();
  auto temp_struct_scatra_sparse_matrix =
      LINALG::CastToSparseMatrixAndCheckSuccess(temp_struct_scatra_mat_);
  auto struct_scatra_sparse_matrix =
      LINALG::CastToSparseMatrixAndCheckSuccess(structure_scatra_matrix);

  // apply mesh tying contributions to temp matrix and complete it
  ApplyMeshtyingToStructureXXX(*temp_struct_scatra_sparse_matrix, *struct_scatra_sparse_matrix);
  temp_struct_scatra_sparse_matrix->Complete(
      *SSIMono().ScaTraField()->DofRowMap(), *SSIMono().StructureField()->DofRowMap());

  // uncomplete matrix, add mesh tying entries stored in temp matrix to matrix and complete again
  struct_scatra_sparse_matrix->UnComplete();
  struct_scatra_sparse_matrix->Add(*temp_struct_scatra_sparse_matrix, false, 1.0, 0.0);
  struct_scatra_sparse_matrix->Complete(
      *SSIMono().ScaTraField()->DofRowMap(), *SSIMono().StructureField()->DofRowMap());
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::MeshtyingStrategyBlockSparse::ApplyMeshtyingToStructureScatra(
    Teuchos::RCP<LINALG::SparseOperator> structure_scatra_matrix)
{
  temp_struct_scatra_mat_->Zero();
  auto temp_struct_scatra_sparse_matrix =
      LINALG::CastToSparseMatrixAndCheckSuccess(temp_struct_scatra_mat_);
  auto struct_scatra_sparse_matrix =
      LINALG::CastToSparseMatrixAndCheckSuccess(structure_scatra_matrix);

  // apply mesh tying contributions to temp matrix and complete it
  ApplyMeshtyingToStructureXXX(*temp_struct_scatra_sparse_matrix, *struct_scatra_sparse_matrix);
  temp_struct_scatra_sparse_matrix->Complete(
      *SSIMono().ScaTraField()->DofRowMap(), *SSIMono().StructureField()->DofRowMap());

  // uncomplete matrix, add mesh tying entries stored in temp matrix to matrix and complete again
  struct_scatra_sparse_matrix->UnComplete();
  struct_scatra_sparse_matrix->Add(*temp_struct_scatra_sparse_matrix, false, 1.0, 0.0);
  struct_scatra_sparse_matrix->Complete(
      *SSIMono().ScaTraField()->DofRowMap(), *SSIMono().StructureField()->DofRowMap());
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::MeshtyingStrategyBlockBlock::ApplyMeshtyingToStructureScatra(
    Teuchos::RCP<LINALG::SparseOperator> structure_scatra_matrix)
{
  temp_struct_scatra_mat_->Zero();
  auto temp_struct_scatra_block_matrix =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(temp_struct_scatra_mat_);
  auto structure_scatra_block_matrix =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(structure_scatra_matrix);

  // apply mesh tying contributions to temp matrix blocks and complete the resulting matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
  {
    ApplyMeshtyingToStructureXXX(temp_struct_scatra_block_matrix->Matrix(0, iblock),
        structure_scatra_block_matrix->Matrix(0, iblock));
  }
  temp_struct_scatra_block_matrix->Complete();

  // uncomplete matrix, add mesh tying entries stored in temp matrix to matrix and complete again
  structure_scatra_block_matrix->UnComplete();
  structure_scatra_block_matrix->Add(*temp_struct_scatra_block_matrix, false, 1.0, 0.0);
  structure_scatra_block_matrix->Complete();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::MeshtyingStrategyBase::ApplyMeshtyingToStructureXXX(
    LINALG::SparseMatrix& ssi_structure_xxx_matrix,
    const LINALG::SparseMatrix& structure_xxx_matrix)
{
  // assemble derivatives of structure master & interior dofs w.r.t. scatra dofs
  LINALG::MatrixLogicalSplitAndTransform()(structure_xxx_matrix, *MapStructureCondensed(),
      structure_xxx_matrix.DomainMap(), 1.0, nullptr, nullptr, ssi_structure_xxx_matrix, true,
      true);

  // assemble derivatives of structure surface slave dofs & interior dofs w.r.t. scatra dofs
  LINALG::MatrixLogicalSplitAndTransform()(structure_xxx_matrix, *MapStructureSlave(),
      structure_xxx_matrix.DomainMap(), 1.0, &StructureSlaveConverter(), nullptr,
      ssi_structure_xxx_matrix, true, true);

  if (Meshtying3DomainIntersection())
  {
    // assemble derivatives of structure surface line dofs & interior dofs w.r.t. scatra dofs
    LINALG::MatrixLogicalSplitAndTransform()(structure_xxx_matrix,
        *MapStructureSlave3DomainIntersection(), structure_xxx_matrix.DomainMap(), 1.0,
        &StructureSlaveConverter3DomainIntersection(), nullptr, ssi_structure_xxx_matrix, true,
        true);
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::MeshtyingStrategyBase::FinalizeMeshtyingStructureMatrix(
    LINALG::SparseMatrix& ssi_structure_matrix)
{
  // map for slave side structure degrees of freedom
  Teuchos::RCP<const Epetra_Map> slavemaps = Teuchos::null;
  if (Meshtying3DomainIntersection())
  {
    slavemaps = LINALG::MultiMapExtractor::MergeMaps({SSIMono().MapsCoupStruct()->Map(1),
        SSIMono().MapsCoupStruct3DomainIntersection()->Map(1)});
  }
  else
    slavemaps = SSIMono().MapsCoupStruct()->Map(1);

  // subject slave-side rows of structure system matrix to pseudo Dirichlet conditions to finalize
  // structure mesh tying
  const double one(1.0);
  for (int doflid_slave = 0; doflid_slave < slavemaps->NumMyElements(); ++doflid_slave)
  {
    // extract global ID of current slave-side row
    const int dofgid_slave = slavemaps->GID(doflid_slave);
    if (dofgid_slave < 0) dserror("Local ID not found!");

    // apply pseudo Dirichlet conditions to filled matrix, i.e., to local row and column indices
    if (ssi_structure_matrix.Filled())
    {
      const int rowlid_slave = ssi_structure_matrix.RowMap().LID(dofgid_slave);
      if (rowlid_slave < 0) dserror("Global ID not found!");
      if (ssi_structure_matrix.EpetraMatrix()->ReplaceMyValues(
              rowlid_slave, 1, &one, &rowlid_slave))
        dserror("ReplaceMyValues failed!");
    }

    // apply pseudo Dirichlet conditions to unfilled matrix, i.e., to global row and column indices
    else if (ssi_structure_matrix.EpetraMatrix()->InsertGlobalValues(
                 dofgid_slave, 1, &one, &dofgid_slave))
      dserror("InsertGlobalValues failed!");
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
Teuchos::RCP<SSI::MeshtyingStrategyBase> SSI::BuildMeshtyingStrategy(const SSI::SSIMono& ssi_mono,
    LINALG::MatrixType matrixtype_ssi, LINALG::MatrixType matrixtype_scatra,
    Teuchos::RCP<const Epetra_Map> interface_map_scatra)
{
  Teuchos::RCP<SSI::MeshtyingStrategyBase> meshtying_strategy = Teuchos::null;

  switch (matrixtype_ssi)
  {
    case LINALG::MatrixType::block_field:
    {
      switch (matrixtype_scatra)
      {
        case LINALG::MatrixType::block_condition:
        case LINALG::MatrixType::block_condition_dof:
        {
          meshtying_strategy = Teuchos::rcp(
              new SSI::MeshtyingStrategyBlockBlock(ssi_mono, std::move(interface_map_scatra)));
          break;
        }
        case LINALG::MatrixType::sparse:
        {
          meshtying_strategy = Teuchos::rcp(
              new SSI::MeshtyingStrategyBlockSparse(ssi_mono, std::move(interface_map_scatra)));
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
      meshtying_strategy =
          Teuchos::rcp(new SSI::MeshtyingStrategySparse(ssi_mono, std::move(interface_map_scatra)));
      break;
    }
    default:
    {
      dserror("unknown matrix type of SSI problem");
      break;
    }
  }

  return meshtying_strategy;
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
ADAPTER::CouplingSlaveConverter& SSI::MeshtyingStrategyBase::StructureSlaveConverter() const
{
  return slave_side_converter_->InterfaceCouplingAdapterStructureSlaveConverter();
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
ADAPTER::CouplingSlaveConverter&
SSI::MeshtyingStrategyBase::StructureSlaveConverter3DomainIntersection() const
{
  return slave_side_converter_
      ->InterfaceCouplingAdapterStructureSlaveConverter3DomainIntersection();
}