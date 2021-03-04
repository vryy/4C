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
SSI::AssembleStrategyBase::AssembleStrategyBase(const Teuchos::RCP<const SSI::SSIMono> ssi_mono)
    : ssi_mono_(ssi_mono),
      mapstructurecondensed_(
          ssi_mono->SSIInterfaceMeshtying() ? ssi_mono->MapStructureCondensed() : Teuchos::null),
      mapstructureslave_(
          ssi_mono->SSIInterfaceMeshtying() ? ssi_mono->MapsCoupStruct()->Map(1) : Teuchos::null),
      mapstructureslave3domainintersection_(
          (ssi_mono->SSIInterfaceMeshtying() and ssi_mono->Meshtying3DomainIntersection())
              ? ssi_mono->MapsCoupStruct3DomainIntersection()->Map(1)
              : Teuchos::null),
      slave_side_converter_(
          ssi_mono->SSIInterfaceMeshtying() ? ssi_mono->SlaveSideConverter() : Teuchos::null),
      meshtying_3_domain_intersection_(
          ssi_mono->SSIInterfaceMeshtying() and ssi_mono->Meshtying3DomainIntersection())
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategyBlock::AssembleStrategyBlock(const Teuchos::RCP<const SSI::SSIMono> ssi_mono)
    : AssembleStrategyBase(ssi_mono),
      block_position_scatra_(Teuchos::null),
      block_position_scatra_manifold_(Teuchos::null),
      position_structure_(-1)
{
  block_position_scatra_ = ssi_mono_->GetBlockPositions(SSI::Subproblem::scalar_transport);
  position_structure_ = ssi_mono_->GetBlockPositions(SSI::Subproblem::structure)->at(0);
  if (ssi_mono_->IsScaTraManifold())
    block_position_scatra_manifold_ = ssi_mono_->GetBlockPositions(SSI::Subproblem::manifold);

  if (block_position_scatra_ == Teuchos::null) dserror("Cannot get position of scatra blocks");
  if (position_structure_ == -1) dserror("Cannot get position of structure block");
  if (ssi_mono_->IsScaTraManifold() and block_position_scatra_manifold_ == Teuchos::null)
    dserror("Cannot get position of scatra manifold blocks");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategyBlockBlock::AssembleStrategyBlockBlock(
    const Teuchos::RCP<const SSI::SSIMono> ssi_mono)
    : AssembleStrategyBlock(ssi_mono)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategyBlockSparse::AssembleStrategyBlockSparse(
    const Teuchos::RCP<const SSI::SSIMono> ssi_mono)
    : AssembleStrategyBlock(ssi_mono)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategySparse::AssembleStrategySparse(const Teuchos::RCP<const SSI::SSIMono> ssi_mono)
    : AssembleStrategyBase(ssi_mono)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatradomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatradomain_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(scatradomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(BlockPositionScaTra()->size()); ++jblock)
    {
      auto& systemmatrix_block_iscatra_jscatra = systemmatrix_block->Matrix(
          BlockPositionScaTra()->at(iblock), BlockPositionScaTra()->at(jblock));

      systemmatrix_block_iscatra_jscatra.Add(
          scatradomain_block->Matrix(iblock, jblock), false, 1.0, 1.0);
    }
  }

  if (ssi_mono_->SSIInterfaceContact())
  {
    // get scatra-scatra-interface block matrix
    const auto& scatra_scatra_interface_blockmatrix =
        ssi_mono_->CoNitscheStrategySsi()
            ->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::scatra_scatra)
            ->Split<LINALG::DefaultBlockMatrixStrategy>(
                *ssi_mono_->MapsScatra(), *ssi_mono_->MapsScatra());
    scatra_scatra_interface_blockmatrix->Complete();

    // assemble it into the system matrix
    for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
    {
      for (int jblock = 0; jblock < static_cast<int>(BlockPositionScaTra()->size()); ++jblock)
      {
        auto& systemmatrix_block_iscatra_jscatra = systemmatrix_block->Matrix(
            BlockPositionScaTra()->at(iblock), BlockPositionScaTra()->at(jblock));

        // get relevant block and complete it, such that it can be added to system matrix block
        auto& scatra_scatra_interface_block_i_j =
            scatra_scatra_interface_blockmatrix->Matrix(iblock, jblock);
        systemmatrix_block_iscatra_jscatra.Add(scatra_scatra_interface_block_i_j, false, 1.0, 1.0);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatradomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatradomain_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(scatradomain);

  auto& systemmatrix_block_scatra_scatra =
      systemmatrix_block->Matrix(BlockPositionScaTra()->at(0), BlockPositionScaTra()->at(0));

  systemmatrix_block_scatra_scatra.Add(*scatradomain_sparse, false, 1.0, 1.0);

  if (ssi_mono_->SSIInterfaceContact())
  {
    systemmatrix_block_scatra_scatra.Add(*(ssi_mono_->CoNitscheStrategySsi()->GetMatrixBlockPtr(
                                             DRT::UTILS::MatBlockType::scatra_scatra)),
        false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatradomain)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatradomain_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(scatradomain);

  // add scalar transport system matrix to global system matrix
  systemmatrix_sparse->Add(*scatradomain_sparse, false, 1.0, 1.0);

  if (ssi_mono_->SSIInterfaceContact())
  {
    systemmatrix_sparse->Add(*(ssi_mono_->CoNitscheStrategySsi()->GetMatrixBlockPtr(
                                 DRT::UTILS::MatBlockType::scatra_scatra)),
        false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseMatrix> structuredomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  // add entire block or assemble slave side to master side
  if (ssi_mono_->SSIInterfaceMeshtying())
  {
    AssembleStructureMeshtying(
        systemmatrix_block->Matrix(PositionStructure(), PositionStructure()), structuredomain);
  }
  else
  {
    auto& systemmatrix_block_struct_struct =
        systemmatrix_block->Matrix(PositionStructure(), PositionStructure());

    systemmatrix_block_struct_struct.Add(*structuredomain, false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseMatrix> structuredomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  // add entire block or assemble slave side to master side
  if (ssi_mono_->SSIInterfaceMeshtying())
  {
    AssembleStructureMeshtying(
        systemmatrix_block->Matrix(PositionStructure(), PositionStructure()), structuredomain);
  }
  else
  {
    auto& systemmatrix_block_struct_struct =
        systemmatrix_block->Matrix(PositionStructure(), PositionStructure());

    systemmatrix_block_struct_struct.Add(*structuredomain, false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseMatrix> structuredomain)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);

  // add entire block or assemble slave side to master side
  if (ssi_mono_->SSIInterfaceMeshtying())
    AssembleStructureMeshtying(*systemmatrix_sparse, structuredomain);
  else
    systemmatrix_sparse->Add(*structuredomain, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::AssembleStructureMeshtying(
    LINALG::SparseMatrix& systemmatrix_structure,
    Teuchos::RCP<LINALG::SparseMatrix> structurematrix)
{
  /* Transform and assemble the structural matrix in the global system matrix block by block:
   * S_m: structural interior and master side dofs
   * S_ss: structural slave surface dofs
   * S_sl: structural slave line dofs
   *
   *       S_m  S_ss  S_sl
   *       --------------
   * S_m  |  a |  b |  c |
   * S_ss |  d |  e |  f |
   * S_sl |  g |  h |  i |
   *       --------------
   */
  // assemble derivs. of interior & master dofs w.r.t. interior & master dofs (block a)
  LINALG::MatrixLogicalSplitAndTransform()(*structurematrix, *MapStructureCondensed(),
      *MapStructureCondensed(), 1.0, nullptr, nullptr, systemmatrix_structure, true, true);

  // assemble  derivs. of surface slave dofs w.r.t. master & interior dofs (block d)
  LINALG::MatrixLogicalSplitAndTransform()(*structurematrix, *MapStructureSlave(),
      *MapStructureCondensed(), 1.0, &StructureSlaveConverter(), nullptr, systemmatrix_structure,
      true, true);

  // assemble derivs. of master & interior w.r.t. surface slave dofs (block b)
  LINALG::MatrixLogicalSplitAndTransform()(*structurematrix, *MapStructureCondensed(),
      *MapStructureSlave(), 1.0, nullptr, &StructureSlaveConverter(), systemmatrix_structure, true,
      true);

  // assemble derivs. of surface slave dofs w.r.t. surface slave dofs (block e)
  LINALG::MatrixLogicalSplitAndTransform()(*structurematrix, *MapStructureSlave(),
      *MapStructureSlave(), 1.0, &StructureSlaveConverter(), &StructureSlaveConverter(),
      systemmatrix_structure, true, true);

  if (Meshtying3DomainIntersection())
  {
    // assemble derivs. of line slave dofs w.r.t. master & interior (block g)
    LINALG::MatrixLogicalSplitAndTransform()(*structurematrix,
        *MapStructureSlave3DomainIntersection(), *MapStructureCondensed(), 1.0,
        &StructureSlaveConverter3DomainIntersection(), nullptr, systemmatrix_structure, true, true);

    // assemble derivs. of master & interior w.r.t. line slave dofs (block c)
    LINALG::MatrixLogicalSplitAndTransform()(*structurematrix, *MapStructureCondensed(),
        *MapStructureSlave3DomainIntersection(), 1.0, nullptr,
        &StructureSlaveConverter3DomainIntersection(), systemmatrix_structure, true, true);

    // assemble derivs. of line slave dof w.r.t. line slave dofs (block i)
    LINALG::MatrixLogicalSplitAndTransform()(*structurematrix,
        *MapStructureSlave3DomainIntersection(), *MapStructureSlave3DomainIntersection(), 1.0,
        &StructureSlaveConverter3DomainIntersection(),
        &StructureSlaveConverter3DomainIntersection(), systemmatrix_structure, true, true);

    // assemble derivs. of surface slave dofs w.r.t. line slave dofs (block f)
    LINALG::MatrixLogicalSplitAndTransform()(*structurematrix, *MapStructureSlave(),
        *MapStructureSlave3DomainIntersection(), 1.0, &StructureSlaveConverter(),
        &StructureSlaveConverter3DomainIntersection(), systemmatrix_structure, true, true);

    // assemble derivs. of line slave dofs w.r.t. surface slave dofs (block h)
    LINALG::MatrixLogicalSplitAndTransform()(*structurematrix,
        *MapStructureSlave3DomainIntersection(), *MapStructureSlave(), 1.0,
        &StructureSlaveConverter3DomainIntersection(), &StructureSlaveConverter(),
        systemmatrix_structure, true, true);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScatraStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatrastructuredomain_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(scatrastructuredomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
  {
    // add entire block or assemble slave side to master side
    if (ssi_mono_->SSIInterfaceMeshtying())
    {
      AssembleXXXStructureMeshtying(
          systemmatrix_block->Matrix(BlockPositionScaTra()->at(iblock), PositionStructure()),
          scatrastructuredomain_block->Matrix(iblock, 0));

      auto scatrastructureinterface_block =
          LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(scatrastructureinterface);

      AssembleXXXStructureMeshtying(
          systemmatrix_block->Matrix(BlockPositionScaTra()->at(iblock), PositionStructure()),
          scatrastructureinterface_block->Matrix(iblock, 0));
    }
    else
    {
      auto& systemmatrix_block_iscatra_struct =
          systemmatrix_block->Matrix(BlockPositionScaTra()->at(iblock), PositionStructure());

      systemmatrix_block_iscatra_struct.Add(
          scatrastructuredomain_block->Matrix(iblock, 0), false, 1.0, 1.0);
    }
  }

  if (ssi_mono_->SSIInterfaceContact())
  {
    // get scatra-structure-interface block matrix
    const auto& scatra_struct_interface_blockmatrix =
        ssi_mono_->CoNitscheStrategySsi()
            ->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::scatra_displ)
            ->Split<LINALG::DefaultBlockMatrixStrategy>(
                *ssi_mono_->MapStructure(), *ssi_mono_->MapsScatra());
    scatra_struct_interface_blockmatrix->Complete();

    // assemble it into the system matrix
    for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
    {
      auto& systemmatrix_block_iscatra_struct =
          systemmatrix_block->Matrix(BlockPositionScaTra()->at(iblock), PositionStructure());

      // get relevant block and complete it, such that it can be added to system matrix block
      auto& iscatra_struct_interface_block = scatra_struct_interface_blockmatrix->Matrix(iblock, 0);
      systemmatrix_block_iscatra_struct.Add(iscatra_struct_interface_block, false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleScatraStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatrastructuredomain_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(scatrastructuredomain);

  // add entire block or assemble slave side to master side
  if (ssi_mono_->SSIInterfaceMeshtying())
  {
    AssembleXXXStructureMeshtying(
        systemmatrix_block->Matrix(BlockPositionScaTra()->at(0), PositionStructure()),
        *scatrastructuredomain_sparse);

    auto scatrastructureinterface_sparse =
        LINALG::CastToSparseMatrixAndCheckSuccess(scatrastructureinterface);

    AssembleXXXStructureMeshtying(
        systemmatrix_block->Matrix(BlockPositionScaTra()->at(0), PositionStructure()),
        *scatrastructureinterface_sparse);
  }
  else
  {
    auto& systemmatrix_block_scatra_struct =
        systemmatrix_block->Matrix(BlockPositionScaTra()->at(0), PositionStructure());

    systemmatrix_block_scatra_struct.Add(*scatrastructuredomain_sparse, false, 1.0, 1.0);
  }

  if (ssi_mono_->SSIInterfaceContact())
  {
    const auto& scatra_struct_interface_matrix =
        ssi_mono_->CoNitscheStrategySsi()->GetMatrixBlockPtr(
            DRT::UTILS::MatBlockType::scatra_displ);

    auto& systemmatrix_block_scatra_struct =
        systemmatrix_block->Matrix(BlockPositionScaTra()->at(0), PositionStructure());

    systemmatrix_block_scatra_struct.Add(*scatra_struct_interface_matrix, false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleScatraStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatrastructuredomain_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(scatrastructuredomain);

  // add entire block or assemble slave side to master side
  if (ssi_mono_->SSIInterfaceMeshtying())
  {
    AssembleXXXStructureMeshtying(*systemmatrix_sparse, *scatrastructuredomain_sparse);

    auto scatrastructureinterface_sparse =
        LINALG::CastToSparseMatrixAndCheckSuccess(scatrastructureinterface);

    AssembleXXXStructureMeshtying(*systemmatrix_sparse, *scatrastructureinterface_sparse);
  }
  else
    systemmatrix_sparse->Add(*scatrastructuredomain_sparse, false, 1.0, 1.0);

  if (ssi_mono_->SSIInterfaceContact())
  {
    systemmatrix_sparse->Add(*(ssi_mono_->CoNitscheStrategySsi()->GetMatrixBlockPtr(
                                 DRT::UTILS::MatBlockType::scatra_displ)),
        false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::AssembleXXXStructureMeshtying(
    LINALG::SparseMatrix& systemmatrix_xxx_structure,
    const LINALG::SparseMatrix& xxx_structurematrix)
{
  // assemble derivs. of x w.r.t. structural master & interior dofs
  LINALG::MatrixLogicalSplitAndTransform()(xxx_structurematrix, xxx_structurematrix.RangeMap(),
      *MapStructureCondensed(), 1.0, nullptr, nullptr, systemmatrix_xxx_structure, true, true);

  // assemble derivs. of x w.r.t. structural surface slave dofs
  LINALG::MatrixLogicalSplitAndTransform()(xxx_structurematrix, xxx_structurematrix.RangeMap(),
      *MapStructureSlave(), 1.0, nullptr, &StructureSlaveConverter(), systemmatrix_xxx_structure,
      true, true);

  if (Meshtying3DomainIntersection())
  {
    // assemble derivs. of x w.r.t. structural line slave dofs
    LINALG::MatrixLogicalSplitAndTransform()(xxx_structurematrix, xxx_structurematrix.RangeMap(),
        *MapStructureSlave3DomainIntersection(), 1.0, nullptr,
        &StructureSlaveConverter3DomainIntersection(), systemmatrix_xxx_structure, true, true);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleStructureScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> structurescatradomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto structurescatradomain_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(structurescatradomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
  {
    // add entire block or assemble slave side to master side
    if (ssi_mono_->SSIInterfaceMeshtying())
    {
      AssembleStructureXXXMeshtying(
          systemmatrix_block->Matrix(PositionStructure(), BlockPositionScaTra()->at(iblock)),
          structurescatradomain_block->Matrix(0, iblock));
    }
    else
    {
      auto& systemmatrix_block_struct_iscatra =
          systemmatrix_block->Matrix(PositionStructure(), BlockPositionScaTra()->at(iblock));

      systemmatrix_block_struct_iscatra.Add(
          structurescatradomain_block->Matrix(0, iblock), false, 1.0, 1.0);
    }
  }

  if (ssi_mono_->SSIInterfaceContact())
  {
    // get structure-scatra-interface block matrix
    const auto& struct_scatra_interface_blockmatrix =
        ssi_mono_->CoNitscheStrategySsi()
            ->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::displ_scatra)
            ->Split<LINALG::DefaultBlockMatrixStrategy>(
                *ssi_mono_->MapsScatra(), *ssi_mono_->MapStructure());
    struct_scatra_interface_blockmatrix->Complete();

    // assemble it into the system matrix
    for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
    {
      auto& systemmatrix_block_struct_iscatra =
          systemmatrix_block->Matrix(PositionStructure(), BlockPositionScaTra()->at(iblock));

      // get relevant block and complete it, such that it can be added to system matrix block
      auto& struct_iscatra_interface_block = struct_scatra_interface_blockmatrix->Matrix(0, iblock);
      systemmatrix_block_struct_iscatra.Add(struct_iscatra_interface_block, false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleStructureScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> structurescatradomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto structurescatradomain_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(structurescatradomain);

  // add entire block or assemble slave side to master side
  if (ssi_mono_->SSIInterfaceMeshtying())
  {
    AssembleStructureXXXMeshtying(
        systemmatrix_block->Matrix(PositionStructure(), BlockPositionScaTra()->at(0)),
        *structurescatradomain_sparse);
  }
  else
  {
    auto& systemmatrix_block_struct_scatra =
        systemmatrix_block->Matrix(PositionStructure(), BlockPositionScaTra()->at(0));

    systemmatrix_block_struct_scatra.Add(*structurescatradomain_sparse, false, 1.0, 1.0);
  }

  if (ssi_mono_->SSIInterfaceContact())
  {
    const auto& struct_scatra_interface_matrix =
        ssi_mono_->CoNitscheStrategySsi()->GetMatrixBlockPtr(
            DRT::UTILS::MatBlockType::displ_scatra);

    auto& systemmatrix_block_struct_scatra =
        systemmatrix_block->Matrix(PositionStructure(), BlockPositionScaTra()->at(0));

    systemmatrix_block_struct_scatra.Add(*struct_scatra_interface_matrix, false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleStructureScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> structurescatradomain)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto structurescatradomain_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(structurescatradomain);

  // add entire block or assemble slave side to master side
  if (ssi_mono_->SSIInterfaceMeshtying())
    AssembleStructureXXXMeshtying(*systemmatrix_sparse, *structurescatradomain_sparse);
  else
    systemmatrix_sparse->Add(*structurescatradomain_sparse, false, 1.0, 1.0);

  // add contact contributions
  if (ssi_mono_->SSIInterfaceContact())
  {
    systemmatrix_sparse->Add(*(ssi_mono_->CoNitscheStrategySsi()->GetMatrixBlockPtr(
                                 DRT::UTILS::MatBlockType::displ_scatra)),
        false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::AssembleStructureXXXMeshtying(
    LINALG::SparseMatrix& systemmatrix_structure_xxx,
    const LINALG::SparseMatrix& structure_xxx_matrix)
{
  // assemble derivs. of structural master & interior dofs w.r.t. scatra dofs
  LINALG::MatrixLogicalSplitAndTransform()(structure_xxx_matrix, *MapStructureCondensed(),
      structure_xxx_matrix.DomainMap(), 1.0, nullptr, nullptr, systemmatrix_structure_xxx, true,
      true);

  // assemble derivs. of structural surface slave dofs & interior dofs w.r.t. scatra dofs
  LINALG::MatrixLogicalSplitAndTransform()(structure_xxx_matrix, *MapStructureSlave(),
      structure_xxx_matrix.DomainMap(), 1.0, &StructureSlaveConverter(), nullptr,
      systemmatrix_structure_xxx, true, true);

  if (Meshtying3DomainIntersection())
  {
    // assemble derivs. of structural surface line dofs & interior dofs w.r.t. scatra dofs
    LINALG::MatrixLogicalSplitAndTransform()(structure_xxx_matrix,
        *MapStructureSlave3DomainIntersection(), structure_xxx_matrix.DomainMap(), 1.0,
        &StructureSlaveConverter3DomainIntersection(), nullptr, systemmatrix_structure_xxx, true,
        true);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScaTraManifold(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> manifolddomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto manifolddomain_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(manifolddomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTraManifold()->size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(BlockPositionScaTraManifold()->size()); ++jblock)
    {
      auto& systemmatrix_block_iscatramanifold_jscatramanifold = systemmatrix_block->Matrix(
          BlockPositionScaTraManifold()->at(iblock), BlockPositionScaTraManifold()->at(jblock));

      systemmatrix_block_iscatramanifold_jscatramanifold.Add(
          manifolddomain_block->Matrix(iblock, jblock), false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleScaTraManifold(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> manifolddomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto manifolddomain_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(manifolddomain);

  auto& systemmatrix_block_scatramanifold_scatramanifold = systemmatrix_block->Matrix(
      BlockPositionScaTraManifold()->at(0), BlockPositionScaTraManifold()->at(0));

  systemmatrix_block_scatramanifold_scatramanifold.Add(*manifolddomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleScaTraManifold(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> manifolddomain)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto manifolddomain_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(manifolddomain);

  systemmatrix_sparse->Add(*manifolddomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScaTraManifoldStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> manifoldstructuredomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto manifoldstructuredomain_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(manifoldstructuredomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTraManifold()->size()); ++iblock)
  {
    // add entire block or assemble slave side to master side
    if (ssi_mono_->SSIInterfaceMeshtying())
    {
      AssembleXXXStructureMeshtying(
          systemmatrix_block->Matrix(
              BlockPositionScaTraManifold()->at(iblock), PositionStructure()),
          manifoldstructuredomain_block->Matrix(iblock, 0));
    }
    else
    {
      auto& systemmatrix_block_iscatramanifold_struct = systemmatrix_block->Matrix(
          BlockPositionScaTraManifold()->at(iblock), PositionStructure());

      systemmatrix_block_iscatramanifold_struct.Add(
          manifoldstructuredomain_block->Matrix(iblock, 0), false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleScaTraManifoldStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> manifoldstructuredomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto manifoldstructuredomain_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(manifoldstructuredomain);

  // add entire block or assemble slave side to master side
  if (ssi_mono_->SSIInterfaceMeshtying())
  {
    AssembleXXXStructureMeshtying(
        systemmatrix_block->Matrix(BlockPositionScaTraManifold()->at(0), PositionStructure()),
        *manifoldstructuredomain_sparse);
  }
  else
  {
    auto& systemmatrix_block_scatramanifold_struct =
        systemmatrix_block->Matrix(BlockPositionScaTraManifold()->at(0), PositionStructure());

    systemmatrix_block_scatramanifold_struct.Add(*manifoldstructuredomain_sparse, false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleScaTraManifoldStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> manifoldstructuredomain)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto manifoldstructuredomain_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(manifoldstructuredomain);

  // add entire block or assemble slave side to master side
  if (ssi_mono_->SSIInterfaceMeshtying())
    AssembleXXXStructureMeshtying(*systemmatrix_sparse, *manifoldstructuredomain_sparse);
  else
    systemmatrix_sparse->Add(*manifoldstructuredomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::ApplyMeshtyingSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  ApplyMeshtyingSysMat(systemmatrix_block->Matrix(PositionStructure(), PositionStructure()));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::ApplyMeshtyingSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  ApplyMeshtyingSysMat(systemmatrix_block->Matrix(PositionStructure(), PositionStructure()));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::ApplyMeshtyingSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  ApplyMeshtyingSysMat(*systemmatrix_sparse);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::ApplyMeshtyingSysMat(LINALG::SparseMatrix& systemmatrix_structure)
{
  // map for slave side structural degrees of freedom
  Teuchos::RCP<const Epetra_Map> slavemaps;
  if (Meshtying3DomainIntersection())
  {
    slavemaps = LINALG::MultiMapExtractor::MergeMaps({ssi_mono_->MapsCoupStruct()->Map(1),
        ssi_mono_->MapsCoupStruct3DomainIntersection()->Map(1)});
  }
  else
    slavemaps = ssi_mono_->MapsCoupStruct()->Map(1);

  // subject slave-side rows of structural system matrix to pseudo Dirichlet conditions to
  // finalize structural meshtying
  const double one(1.0);
  for (int doflid_slave = 0; doflid_slave < slavemaps->NumMyElements(); ++doflid_slave)
  {
    // extract global ID of current slave-side row
    const int dofgid_slave = slavemaps->GID(doflid_slave);
    if (dofgid_slave < 0) dserror("Local ID not found!");

    // apply pseudo Dirichlet conditions to filled matrix, i.e., to local row and column indices
    if (systemmatrix_structure.Filled())
    {
      const int rowlid_slave = systemmatrix_structure.RowMap().LID(dofgid_slave);
      if (rowlid_slave < 0) dserror("Global ID not found!");
      if (systemmatrix_structure.EpetraMatrix()->ReplaceMyValues(
              rowlid_slave, 1, &one, &rowlid_slave))
        dserror("ReplaceMyValues failed!");
    }

    // apply pseudo Dirichlet conditions to unfilled matrix, i.e., to global row and column indices
    else if (systemmatrix_structure.EpetraMatrix()->InsertGlobalValues(
                 dofgid_slave, 1, &one, &dofgid_slave))
      dserror("InsertGlobalValues failed!");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlock::ApplyStructuralDBCSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  // locsys manager of structure
  const auto& locsysmanager_structure = ssi_mono_->StructureField()->LocsysManager();

  // map of structural Dirichlet BCs
  const auto dbcmap_structure = ssi_mono_->StructureField()->GetDBCMapExtractor()->CondMap();

  if (locsysmanager_structure == Teuchos::null)
    systemmatrix->ApplyDirichlet(*dbcmap_structure);
  else
  {
    auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

    // apply structural Dirichlet conditions
    for (int iblock = 0; iblock < systemmatrix_block->Cols(); ++iblock)
    {
      locsysmanager_structure->RotateGlobalToLocal(
          Teuchos::rcp(&systemmatrix_block->Matrix(PositionStructure(), iblock), false));
      systemmatrix_block->Matrix(PositionStructure(), iblock)
          .ApplyDirichletWithTrafo(
              locsysmanager_structure->Trafo(), *dbcmap_structure, (iblock == PositionStructure()));
      locsysmanager_structure->RotateLocalToGlobal(
          Teuchos::rcp(&systemmatrix_block->Matrix(PositionStructure(), iblock), false));
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::ApplyStructuralDBCSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  // locsys manager of structure
  const auto& locsysmanager_structure = ssi_mono_->StructureField()->LocsysManager();

  // map of structural Dirichlet BCs
  const auto& dbcmap_structure = ssi_mono_->StructureField()->GetDBCMapExtractor()->CondMap();

  // structural dof row map
  const auto& dofrowmap_structure = ssi_mono_->StructureField()->DofRowMap();

  if (locsysmanager_structure == Teuchos::null)
    systemmatrix->ApplyDirichlet(*dbcmap_structure);
  else
  {
    auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);

    // extract structural rows of global system matrix
    const Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_structure =
        Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_structure, 27, false, true));
    LINALG::MatrixLogicalSplitAndTransform()(*systemmatrix_sparse, *dofrowmap_structure,
        systemmatrix->DomainMap(), 1.0, nullptr, nullptr, *systemmatrix_structure);
    systemmatrix_structure->Complete(systemmatrix->DomainMap(), *dofrowmap_structure);

    // apply structural Dirichlet conditions
    locsysmanager_structure->RotateGlobalToLocal(systemmatrix_structure);
    systemmatrix_structure->ApplyDirichletWithTrafo(
        locsysmanager_structure->Trafo(), *dbcmap_structure);
    locsysmanager_structure->RotateLocalToGlobal(systemmatrix_structure);

    // assemble structural rows of global system matrix back into global system matrix
    systemmatrix_sparse->Put(*systemmatrix_structure, 1.0, dofrowmap_structure);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::AssembleRHS(Teuchos::RCP<Epetra_Vector>& RHS,
    Teuchos::RCP<Epetra_Vector> RHSscatra, Teuchos::RCP<const Epetra_Vector> RHSstructure,
    Teuchos::RCP<Epetra_Vector> RHSmanifold)
{
  RHS->PutScalar(0.0);

  // assemble scalar transport right-hand side vector into monolithic right-hand side vector
  ssi_mono_->MapsSubProblems()->InsertVector(
      RHSscatra, ssi_mono_->GetProblemPosition(SSI::Subproblem::scalar_transport), RHS);

  if (ssi_mono_->IsScaTraManifold())
  {
    ssi_mono_->MapsSubProblems()->InsertVector(
        RHSmanifold, ssi_mono_->GetProblemPosition(SSI::Subproblem::manifold), RHS);
  }

  if (ssi_mono_->SSIInterfaceMeshtying())
  {
    // perform structural meshtying before assembling structural right-hand side vector into
    // monolithic right-hand side vector

    // make copy of structural right-hand side vector
    Epetra_Vector residual_structure(*RHSstructure);

    // transform slave-side part of structural right-hand side vector to master side
    Teuchos::RCP<Epetra_Vector> slavetomaster = ssi_mono_->MapsCoupStruct()->InsertVector(
        ssi_mono_->InterfaceCouplingAdapterStructure()->SlaveToMaster(
            ssi_mono_->MapsCoupStruct()->ExtractVector(residual_structure, 1)),
        2);

    if (Meshtying3DomainIntersection())
    {
      slavetomaster->Update(1.0,
          *ssi_mono_->MapsCoupStruct3DomainIntersection()->InsertVector(
              ssi_mono_->InterfaceCouplingAdapterStructure3DomainIntersection()->SlaveToMaster(
                  ssi_mono_->MapsCoupStruct3DomainIntersection()->ExtractVector(
                      residual_structure, 1)),
              2),
          1.0);
    }

    // locsys manager of structure
    const auto& locsysmanager_structure = ssi_mono_->StructureField()->LocsysManager();

    // apply pseudo Dirichlet conditions to transformed slave-side part of structural right-hand
    // side vector
    const auto zeros = Teuchos::rcp(new Epetra_Vector(slavetomaster->Map()));
    if (locsysmanager_structure != Teuchos::null)
      locsysmanager_structure->RotateGlobalToLocal(slavetomaster);
    LINALG::ApplyDirichlettoSystem(
        slavetomaster, zeros, *ssi_mono_->StructureField()->GetDBCMapExtractor()->CondMap());
    if (locsysmanager_structure != Teuchos::null)
      locsysmanager_structure->RotateLocalToGlobal(slavetomaster);

    // assemble transformed slave-side part of structural right-hand side vector
    residual_structure.Update(1.0, *slavetomaster, 1.0);

    // zero out slave-side part of structural right-hand side vector
    ssi_mono_->MapsCoupStruct()->PutScalar(residual_structure, 1, 0.0);
    if (Meshtying3DomainIntersection())
      ssi_mono_->MapsCoupStruct3DomainIntersection()->PutScalar(residual_structure, 1, 0.0);

    // assemble final structural right-hand side vector into monolithic right-hand side vector
    ssi_mono_->MapsSubProblems()->AddVector(
        residual_structure, ssi_mono_->GetProblemPosition(SSI::Subproblem::structure), *RHS, -1.0);
  }
  else
  {
    ssi_mono_->MapsSubProblems()->AddVector(
        RHSstructure, ssi_mono_->GetProblemPosition(SSI::Subproblem::structure), RHS, -1.0);
  }

  if (ssi_mono_->SSIInterfaceContact())
  {
    // add the scatra contact contribution
    ssi_mono_->MapsSubProblems()->AddVector(
        ssi_mono_->CoNitscheStrategySsi()->GetRhsBlockPtr(DRT::UTILS::VecBlockType::scatra),
        ssi_mono_->GetProblemPosition(SSI::Subproblem::scalar_transport), RHS);

    // apply the dirichlet boundary conditions
    const auto zeros = Teuchos::rcp(new Epetra_Vector(RHS->Map()));
    LINALG::ApplyDirichlettoSystem(RHS, zeros, *(ssi_mono_->CombinedDBCMap()));
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
Teuchos::RCP<SSI::AssembleStrategyBase> SSI::BuildAssembleStrategy(
    Teuchos::RCP<const SSI::SSIMono> ssi_mono, LINALG::MatrixType matrixtype_ssi,
    LINALG::MatrixType matrixtype_scatra)
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