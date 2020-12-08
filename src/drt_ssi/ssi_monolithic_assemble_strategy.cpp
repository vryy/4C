/*----------------------------------------------------------------------*/
/*! \file
\brief Assemble strategy for monolithic SSI
\level 2

 */
/*----------------------------------------------------------------------*/
#include "ssi_monolithic_assemble_strategy.H"

#include "ssi_monolithic.H"

#include "../drt_io/io_control.H"

#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_matrixtransform.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"

#include "../drt_lib/drt_locsys.H"

#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"


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
    : AssembleStrategyBase(ssi_mono), block_position_scatra_(Teuchos::null), position_structure_(-1)
{
  block_position_scatra_ = ssi_mono_->GetBlockPositions(SSI::Subproblem::scalar_transport);
  position_structure_ = ssi_mono_->GetBlockPositions(SSI::Subproblem::structure)->at(0);

  if (block_position_scatra_ == Teuchos::null) dserror("Cannot get position of scatra blocks");
  if (position_structure_ == -1) dserror("Cannot get position of structure block");
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
void SSI::AssembleStrategyBlockBlock::AssembleScatraDomain(
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
      systemmatrix_block->Assign(BlockPositionScaTra()->at(iblock),
          BlockPositionScaTra()->at(jblock), LINALG::View,
          scatradomain_block->Matrix(iblock, jblock));
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatradomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  auto scatradomain_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(scatradomain);

  // add scalar transport system matrix to global system matrix
  systemmatrix_block->Assign(BlockPositionScaTra()->at(0), BlockPositionScaTra()->at(0),
      LINALG::View, *scatradomain_sparse);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatradomain)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);

  auto scatradomain_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(scatradomain);

  // add scalar transport system matrix to global system matrix
  systemmatrix_sparse->Add(*scatradomain_sparse, false, 1.0, 0.0);
}

/*----------------------------------------------------------------------*
 |                         assemble structure domain into system matrix |
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseMatrix> structuredomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  // add entire block or assemble slave side to master side
  if (!ssi_mono_->SSIInterfaceMeshtying())
  {
    systemmatrix_block->Assign(
        PositionStructure(), PositionStructure(), LINALG::View, *structuredomain);
  }
  else
  {
    AssembleStructureDomainMeshtying(
        systemmatrix_block->Matrix(PositionStructure(), PositionStructure()), structuredomain,
        false);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseMatrix> structuredomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  // add entire block or assemble slave side to master side
  if (!ssi_mono_->SSIInterfaceMeshtying())
  {
    systemmatrix_block->Assign(
        PositionStructure(), PositionStructure(), LINALG::View, *structuredomain);
  }
  else
  {
    AssembleStructureDomainMeshtying(
        systemmatrix_block->Matrix(PositionStructure(), PositionStructure()), structuredomain,
        false);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseMatrix> structuredomain)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);

  // add entire block or assemble slave side to master side
  if (!ssi_mono_->SSIInterfaceMeshtying())
    systemmatrix_sparse->Add(*structuredomain, false, 1.0, 1.0);
  else
    AssembleStructureDomainMeshtying(*systemmatrix_sparse, structuredomain, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::AssembleStructureDomainMeshtying(
    LINALG::SparseMatrix& systemmatrix_structure,
    Teuchos::RCP<LINALG::SparseMatrix> structuredomain, bool zero)
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
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *MapStructureCondensed(),
      *MapStructureCondensed(), 1.0, nullptr, nullptr, systemmatrix_structure, true, zero);

  // assemble  derivs. of surface slave dofs w.r.t. master & interior dofs (block d)
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *MapStructureSlave(),
      *MapStructureCondensed(), 1.0, &StructureSlaveConverter(), nullptr, systemmatrix_structure,
      true, true);

  // assemble derivs. of master & interior w.r.t. surface slave dofs (block b)
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *MapStructureCondensed(),
      *MapStructureSlave(), 1.0, nullptr, &StructureSlaveConverter(), systemmatrix_structure, true,
      true);

  // assemble derivs. of surface slave dofs w.r.t. surface slave dofs (block e)
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *MapStructureSlave(),
      *MapStructureSlave(), 1.0, &StructureSlaveConverter(), &StructureSlaveConverter(),
      systemmatrix_structure, true, true);

  if (Meshtying3DomainIntersection())
  {
    // assemble derivs. of line slave dofs w.r.t. master & interior (block g)
    LINALG::MatrixLogicalSplitAndTransform()(*structuredomain,
        *MapStructureSlave3DomainIntersection(), *MapStructureCondensed(), 1.0,
        &StructureSlaveConverter3DomainIntersection(), nullptr, systemmatrix_structure, true, true);

    // assemble derivs. of master & interior w.r.t. line slave dofs (block c)
    LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *MapStructureCondensed(),
        *MapStructureSlave3DomainIntersection(), 1.0, nullptr,
        &StructureSlaveConverter3DomainIntersection(), systemmatrix_structure, true, true);

    // assemble derivs. of line slave dof w.r.t. line slave dofs (block i)
    LINALG::MatrixLogicalSplitAndTransform()(*structuredomain,
        *MapStructureSlave3DomainIntersection(), *MapStructureSlave3DomainIntersection(), 1.0,
        &StructureSlaveConverter3DomainIntersection(),
        &StructureSlaveConverter3DomainIntersection(), systemmatrix_structure, true, true);

    // assemble derivs. of surface slave dofs w.r.t. line slave dofs (block f)
    LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *MapStructureSlave(),
        *MapStructureSlave3DomainIntersection(), 1.0, &StructureSlaveConverter(),
        &StructureSlaveConverter3DomainIntersection(), systemmatrix_structure, true, true);

    // assemble derivs. of line slave dofs w.r.t. surface slave dofs (block h)
    LINALG::MatrixLogicalSplitAndTransform()(*structuredomain,
        *MapStructureSlave3DomainIntersection(), *MapStructureSlave(), 1.0,
        &StructureSlaveConverter3DomainIntersection(), &StructureSlaveConverter(),
        systemmatrix_structure, true, true);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScatraStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  auto scatrastructuredomain_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(scatrastructuredomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
  {
    // add entire block or assemble slave side to master side
    if (!ssi_mono_->SSIInterfaceMeshtying())
    {
      systemmatrix_block->Assign(BlockPositionScaTra()->at(iblock), PositionStructure(),
          LINALG::View, scatrastructuredomain_block->Matrix(iblock, 0));
    }
    else
    {
      AssembleScatraStructureDomainMeshtying(
          systemmatrix_block->Matrix(BlockPositionScaTra()->at(iblock), PositionStructure()),
          scatrastructuredomain_block->Matrix(iblock, 0), false);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleScatraStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  auto scatrastructuredomain_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(scatrastructuredomain);

  // add entire block or assemble slave side to master side
  if (!ssi_mono_->SSIInterfaceMeshtying())
  {
    systemmatrix_block->Assign(BlockPositionScaTra()->at(0), PositionStructure(), LINALG::View,
        *scatrastructuredomain_sparse);
  }
  else
  {
    AssembleScatraStructureDomainMeshtying(
        systemmatrix_block->Matrix(BlockPositionScaTra()->at(0), PositionStructure()),
        *scatrastructuredomain_sparse, false);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleScatraStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);

  auto scatrastructuredomain_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(scatrastructuredomain);

  // add entire block or assemble slave side to master side
  if (!ssi_mono_->SSIInterfaceMeshtying())
    systemmatrix_sparse->Add(*scatrastructuredomain_sparse, false, 1.0, 1.0);
  else
    AssembleScatraStructureDomainMeshtying(
        *systemmatrix_sparse, *scatrastructuredomain_sparse, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::AssembleScatraStructureDomainMeshtying(
    LINALG::SparseMatrix& systemmatrix_scatra_structure,
    const LINALG::SparseMatrix& scatrastructuredomain, bool zero)
{
  // assemble derivs. of scatra w.r.t. structural master & interior dofs
  LINALG::MatrixLogicalSplitAndTransform()(scatrastructuredomain, scatrastructuredomain.RangeMap(),
      *MapStructureCondensed(), 1.0, nullptr, nullptr, systemmatrix_scatra_structure, true, zero);

  // assemble derivs. of scatra w.r.t. structural surface slave dofs
  LINALG::MatrixLogicalSplitAndTransform()(scatrastructuredomain, scatrastructuredomain.RangeMap(),
      *MapStructureSlave(), 1.0, nullptr, &StructureSlaveConverter(), systemmatrix_scatra_structure,
      true, true);

  if (Meshtying3DomainIntersection())
  {
    // assemble derivs. of scatra w.r.t. structural line slave dofs
    LINALG::MatrixLogicalSplitAndTransform()(scatrastructuredomain,
        scatrastructuredomain.RangeMap(), *MapStructureSlave3DomainIntersection(), 1.0, nullptr,
        &StructureSlaveConverter3DomainIntersection(), systemmatrix_scatra_structure, true, true);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScatraStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  auto scatrastructureinterface_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(scatrastructureinterface);

  // derive linearizations of master-side scatra fluxes w.r.t. master-side structural dofs and
  // assemble into auxiliary system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
  {
    // assemble scatra-structure-interface into system matrix
    AssembleScatraStructureDomainMeshtying(
        systemmatrix_block->Matrix(BlockPositionScaTra()->at(iblock), PositionStructure()),
        scatrastructureinterface_block->Matrix(iblock, 0), true);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleScatraStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  auto scatrastructureinterface_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(scatrastructureinterface);

  // assemble scatra-structure-interface into system matrix
  AssembleScatraStructureDomainMeshtying(
      systemmatrix_block->Matrix(BlockPositionScaTra()->at(0), PositionStructure()),
      *scatrastructureinterface_sparse, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleScatraStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);

  auto scatrastructureinterface_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(scatrastructureinterface);

  // assemble scatra-structure-interface into system matrix
  AssembleScatraStructureDomainMeshtying(
      *systemmatrix_sparse, *scatrastructureinterface_sparse, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleStructureScatraDomain(
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
    if (!ssi_mono_->SSIInterfaceMeshtying())
    {
      systemmatrix_block->Assign(PositionStructure(), BlockPositionScaTra()->at(iblock),
          LINALG::View, structurescatradomain_block->Matrix(0, iblock));
    }
    else
    {
      AssembleStructureScatraDomainMeshtying(
          systemmatrix_block->Matrix(PositionStructure(), BlockPositionScaTra()->at(iblock)),
          structurescatradomain_block->Matrix(0, iblock), false);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleStructureScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> structurescatradomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  auto structurescatradomain_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(structurescatradomain);

  // add entire block or assemble slave side to master side
  if (!ssi_mono_->SSIInterfaceMeshtying())
  {
    systemmatrix_block->Assign(PositionStructure(), BlockPositionScaTra()->at(0), LINALG::View,
        *structurescatradomain_sparse);
  }
  else
  {
    AssembleStructureScatraDomainMeshtying(
        systemmatrix_block->Matrix(PositionStructure(), BlockPositionScaTra()->at(0)),
        *structurescatradomain_sparse, false);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleStructureScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> structurescatradomain)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);

  auto structurescatradomain_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(structurescatradomain);

  // add entire block or assemble slave side to master side
  if (!ssi_mono_->SSIInterfaceMeshtying())
    systemmatrix_sparse->Add(*structurescatradomain_sparse, false, 1.0, 1.0);
  else
    AssembleStructureScatraDomainMeshtying(
        *systemmatrix_sparse, *structurescatradomain_sparse, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::AssembleStructureScatraDomainMeshtying(
    LINALG::SparseMatrix& systemmatrix_structure_scatra,
    const LINALG::SparseMatrix& structurescatradomain, bool zero)
{
  // assemble derivs. of structural master & interior dofs w.r.t. scatra dofs
  LINALG::MatrixLogicalSplitAndTransform()(structurescatradomain, *MapStructureCondensed(),
      structurescatradomain.DomainMap(), 1.0, nullptr, nullptr, systemmatrix_structure_scatra, true,
      zero);

  // assemble derivs. of structural surface slave dofs & interior dofs w.r.t. scatra dofs
  LINALG::MatrixLogicalSplitAndTransform()(structurescatradomain, *MapStructureSlave(),
      structurescatradomain.DomainMap(), 1.0, &StructureSlaveConverter(), nullptr,
      systemmatrix_structure_scatra, true, true);

  if (Meshtying3DomainIntersection())
  {
    // assemble derivs. of structural surface line dofs & interior dofs w.r.t. scatra dofs
    LINALG::MatrixLogicalSplitAndTransform()(structurescatradomain,
        *MapStructureSlave3DomainIntersection(), structurescatradomain.DomainMap(), 1.0,
        &StructureSlaveConverter3DomainIntersection(), nullptr, systemmatrix_structure_scatra, true,
        true);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::ApplyMeshtyingSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  if (ssi_mono_->SSIInterfaceMeshtying())
  {
    auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

    ApplyMeshtyingSysMat(systemmatrix_block->Matrix(PositionStructure(), PositionStructure()));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::ApplyMeshtyingSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  if (ssi_mono_->SSIInterfaceMeshtying())
  {
    auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

    ApplyMeshtyingSysMat(systemmatrix_block->Matrix(PositionStructure(), PositionStructure()));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::ApplyMeshtyingSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  if (ssi_mono_->SSIInterfaceMeshtying())
  {
    auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);

    ApplyMeshtyingSysMat(*systemmatrix_sparse);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::ApplyMeshtyingSysMat(LINALG::SparseMatrix& systemmatrix_structure)
{
  // map for slave side structural degrees of freedom
  Teuchos::RCP<const Epetra_Map> slavemaps;
  if (!Meshtying3DomainIntersection())
    slavemaps = ssi_mono_->MapsCoupStruct()->Map(1);
  else
  {
    slavemaps = LINALG::MultiMapExtractor::MergeMaps({ssi_mono_->MapsCoupStruct()->Map(1),
        ssi_mono_->MapsCoupStruct3DomainIntersection()->Map(1)});
  }

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
  // locsys manager of strucutre
  const auto& locsysmanager_structure = ssi_mono_->StructureField()->LocsysManager();

  // map of strucutral Dirichlet BCs
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
  // locsys manager of strucutre
  const auto& locsysmanager_structure = ssi_mono_->StructureField()->LocsysManager();

  // map of strucutral Dirichlet BCs
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
    Teuchos::RCP<Epetra_Vector> RHSscatra, Teuchos::RCP<const Epetra_Vector> RHSstructure)
{
  RHS->PutScalar(0.0);

  // assemble scalar transport right-hand side vector into monolithic right-hand side vector
  ssi_mono_->MapsSubProblems()->InsertVector(
      RHSscatra, ssi_mono_->GetProblemPosition(SSI::Subproblem::scalar_transport), RHS);

  if (!ssi_mono_->SSIInterfaceMeshtying())
  {
    ssi_mono_->MapsSubProblems()->AddVector(
        RHSstructure, ssi_mono_->GetProblemPosition(SSI::Subproblem::structure), RHS, -1.0);
  }
  else
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
    {
      ssi_mono_->MapsCoupStruct3DomainIntersection()->PutScalar(
          residual_structure, ssi_mono_->GetProblemPosition(SSI::Subproblem::structure), 0.0);
    }

    // assemble final structural right-hand side vector into monolithic right-hand side vector
    ssi_mono_->MapsSubProblems()->AddVector(
        residual_structure, ssi_mono_->GetProblemPosition(SSI::Subproblem::structure), *RHS, -1.0);
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