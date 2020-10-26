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
SSI::AssembleStrategyBase::AssembleStrategyBase(const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono)
    : ssi_mono_(ssi_mono)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategyBlock::AssembleStrategyBlock(const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono)
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
    const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono)
    : AssembleStrategyBlock(ssi_mono)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategyBlockSparse::AssembleStrategyBlockSparse(
    const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono)
    : AssembleStrategyBlock(ssi_mono)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategySparse::AssembleStrategySparse(
    const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono)
    : AssembleStrategyBase(ssi_mono)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatradomain)
{
  // cast system matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastMatrixBlock(systemmatrix, systemmatrix_block);

  // cast scatra block
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> scatradomain_block;
  CastMatrixBlock(scatradomain, scatradomain_block);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < BlockPositionScaTra()->size(); ++iblock)
  {
    for (int jblock = 0; jblock < BlockPositionScaTra()->size(); ++jblock)
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
  // cast system matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastMatrixBlock(systemmatrix, systemmatrix_block);

  // cast scatra block
  Teuchos::RCP<LINALG::SparseMatrix> scatradomain_sparse;
  CastMatrixSparse(scatradomain, scatradomain_sparse);

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
  // cast system matrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastMatrixSparse(systemmatrix, systemmatrix_sparse);

  // cast scatra block
  Teuchos::RCP<LINALG::SparseMatrix> scatradomain_sparse;
  CastMatrixSparse(scatradomain, scatradomain_sparse);

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
  // cast system matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastMatrixBlock(systemmatrix, systemmatrix_block);

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
  // cast system matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastMatrixBlock(systemmatrix, systemmatrix_block);

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
  // cast system matrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastMatrixSparse(systemmatrix, systemmatrix_sparse);

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
  // map for interior and master-side structural degrees of freedom
  const auto& mapstructurecondensed = ssi_mono_->MapStructureCondensed();

  // map for slave side structural degrees of freedom
  const auto& mapstructureslave = ssi_mono_->MapsStructure()->Map(1);

  const auto structure_slave_converter =
      ssi_mono_->InterfaceCouplingAdapterStructureSlaveConverter();

  // assemble interior and master-side rows and columns of structural system matrix into
  // global system matrix
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *mapstructurecondensed,
      *mapstructurecondensed, 1.0, nullptr, nullptr, systemmatrix_structure, true, zero);

  // transform and assemble slave-side rows of structural system matrix into global system matrix
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *mapstructureslave,
      *ssi_mono_->MapStructureCondensed(), 1.0, &structure_slave_converter, nullptr,
      systemmatrix_structure, true, true);

  // transform and assemble slave-side columns of structural system matrix into global system matrix
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *mapstructurecondensed,
      *mapstructureslave, 1.0, nullptr, &structure_slave_converter, systemmatrix_structure, true,
      true);

  // transform and assemble slave-side rows and columns of structural system matrix into global
  // system matrix
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *mapstructureslave, *mapstructureslave,
      1.0, &structure_slave_converter, &structure_slave_converter, systemmatrix_structure, true,
      true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScatraStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain)
{
  // cast system matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastMatrixBlock(systemmatrix, systemmatrix_block);

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> scatrastructuredomain_block;
  CastMatrixBlock(scatrastructuredomain, scatrastructuredomain_block);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < BlockPositionScaTra()->size(); ++iblock)
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
  // cast system matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastMatrixBlock(systemmatrix, systemmatrix_block);

  Teuchos::RCP<LINALG::SparseMatrix> scatrastructuredomain_sparse;
  CastMatrixSparse(scatrastructuredomain, scatrastructuredomain_sparse);

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
  // cast system matrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastMatrixSparse(systemmatrix, systemmatrix_sparse);

  Teuchos::RCP<LINALG::SparseMatrix> scatrastructuredomain_sparse;
  CastMatrixSparse(scatrastructuredomain, scatrastructuredomain_sparse);

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
  // map for interior and master-side structural degrees of freedom
  const auto& mapstructurecondensed = ssi_mono_->MapStructureCondensed();

  // map for slave side structural degrees of freedom
  const auto& mapstructureslave = ssi_mono_->MapsStructure()->Map(1);

  // assemble interior and master-side columns of scatra-structure block into global system matrix
  LINALG::MatrixLogicalSplitAndTransform()(scatrastructuredomain, scatrastructuredomain.RangeMap(),
      *mapstructurecondensed, 1.0, nullptr, nullptr, systemmatrix_scatra_structure, true, zero);

  // transform and assemble slave-side columns of scatra-structure block into global system matrix
  LINALG::MatrixLogicalSplitAndTransform()(scatrastructuredomain, scatrastructuredomain.RangeMap(),
      *mapstructureslave, 1.0, nullptr,
      &ssi_mono_->InterfaceCouplingAdapterStructureSlaveConverter(), systemmatrix_scatra_structure,
      true, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScatraStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface)
{
  // cast system matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastMatrixBlock(systemmatrix, systemmatrix_block);

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> scatrastructureinterface_block;
  CastMatrixBlock(scatrastructureinterface, scatrastructureinterface_block);

  // derive linearizations of master-side scatra fluxes w.r.t. master-side structural dofs and
  // assemble into auxiliary system matrix
  for (int iblock = 0; iblock < BlockPositionScaTra()->size(); ++iblock)
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
  // cast system matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastMatrixBlock(systemmatrix, systemmatrix_block);

  Teuchos::RCP<LINALG::SparseMatrix> scatrastructureinterface_sparse;
  CastMatrixSparse(scatrastructureinterface, scatrastructureinterface_sparse);

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
  // cast system matrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastMatrixSparse(systemmatrix, systemmatrix_sparse);

  Teuchos::RCP<LINALG::SparseMatrix> scatrastructureinterface_sparse;
  CastMatrixSparse(scatrastructureinterface, scatrastructureinterface_sparse);

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
  // cast system matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastMatrixBlock(systemmatrix, systemmatrix_block);

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> structurescatradomain_block;
  CastMatrixBlock(structurescatradomain, structurescatradomain_block);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < BlockPositionScaTra()->size(); ++iblock)
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
  // cast system matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastMatrixBlock(systemmatrix, systemmatrix_block);

  // cast structure-scatra block
  Teuchos::RCP<LINALG::SparseMatrix> structurescatradomain_sparse;
  CastMatrixSparse(structurescatradomain, structurescatradomain_sparse);

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
  // cast system matrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastMatrixSparse(systemmatrix, systemmatrix_sparse);

  // cast structure-scatra block
  Teuchos::RCP<LINALG::SparseMatrix> structurescatradomain_sparse;
  CastMatrixSparse(structurescatradomain, structurescatradomain_sparse);

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
  // map for interior and master-side structural degrees of freedom
  const auto& mapstructurecondensed = ssi_mono_->MapStructureCondensed();

  // map for slave side structural degrees of freedom
  const auto& mapstructureslave = ssi_mono_->MapsStructure()->Map(1);

  // assemble interior and master - side rows of structure - scatra block into global system matrix
  LINALG::MatrixLogicalSplitAndTransform()(structurescatradomain, *mapstructurecondensed,
      structurescatradomain.DomainMap(), 1.0, nullptr, nullptr, systemmatrix_structure_scatra, true,
      zero);

  // transform and assemble slave-side rows of structure-scatra block into global system matrix
  LINALG::MatrixLogicalSplitAndTransform()(structurescatradomain, *mapstructureslave,
      structurescatradomain.DomainMap(), 1.0,
      &ssi_mono_->InterfaceCouplingAdapterStructureSlaveConverter(), nullptr,
      systemmatrix_structure_scatra, true, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::ApplyMeshtyingSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  if (ssi_mono_->SSIInterfaceMeshtying())
  {
    // cast system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
    CastMatrixBlock(systemmatrix, systemmatrix_block);

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
    // cast system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
    CastMatrixBlock(systemmatrix, systemmatrix_block);

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
    // cast system matrix
    Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
    CastMatrixSparse(systemmatrix, systemmatrix_sparse);

    ApplyMeshtyingSysMat(*systemmatrix_sparse);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::ApplyMeshtyingSysMat(LINALG::SparseMatrix& systemmatrix_structure)
{
  // map for slave side structural degrees of freedom
  const auto& mapstructureslave = ssi_mono_->MapsStructure()->Map(1);

  // subject slave-side rows of structural system matrix to pseudo Dirichlet conditions to finalize
  // structural meshtying
  const double one(1.0);
  for (int doflid_slave = 0; doflid_slave < mapstructureslave->NumMyElements(); ++doflid_slave)
  {
    // extract global ID of current slave-side row
    const int dofgid_slave = mapstructureslave->GID(doflid_slave);
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
    // cast systemmatrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
    CastMatrixBlock(systemmatrix, systemmatrix_block);

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
    // cast systemmatrix
    Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
    CastMatrixSparse(systemmatrix, systemmatrix_sparse);

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
  // zero out RHS
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
    Teuchos::RCP<Epetra_Vector> slavetomaster = ssi_mono_->MapsStructure()->InsertVector(
        ssi_mono_->InterfaceCouplingAdapterStructure()->SlaveToMaster(
            ssi_mono_->MapsStructure()->ExtractVector(residual_structure, 1)),
        2);

    // locsys manager of strucutre
    const auto& locsysmanager_structure = ssi_mono_->StructureField()->LocsysManager();

    // apply pseudo Dirichlet conditions to transformed slave-side part of structural right-hand
    // side vector
    const Teuchos::RCP<const Epetra_Vector> zeros =
        Teuchos::rcp(new Epetra_Vector(slavetomaster->Map()));
    if (locsysmanager_structure != Teuchos::null)
      locsysmanager_structure->RotateGlobalToLocal(slavetomaster);
    LINALG::ApplyDirichlettoSystem(
        slavetomaster, zeros, *ssi_mono_->StructureField()->GetDBCMapExtractor()->CondMap());
    if (locsysmanager_structure != Teuchos::null)
      locsysmanager_structure->RotateLocalToGlobal(slavetomaster);

    // assemble transformed slave-side part of structural right-hand side vector
    residual_structure.Update(1.0, *slavetomaster, 1.0);

    // zero out slave-side part of structural right-hand side vector
    ssi_mono_->MapsStructure()->PutScalar(residual_structure, 1, 0.0);

    // assemble final structural right-hand side vector into monolithic right-hand side vector
    ssi_mono_->MapsSubProblems()->AddVector(residual_structure, 1, *RHS, -1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::CastMatrixBlock(Teuchos::RCP<LINALG::SparseOperator> input_matrix,
    Teuchos::RCP<LINALG::BlockSparseMatrixBase>& block_matrix)
{
  block_matrix = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(input_matrix);
  if (block_matrix == Teuchos::null) dserror("Matrix is not a block matrix!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::CastMatrixSparse(Teuchos::RCP<LINALG::SparseOperator> input_matrix,
    Teuchos::RCP<LINALG::SparseMatrix>& sparse_matrix)
{
  sparse_matrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(input_matrix);
  if (sparse_matrix == Teuchos::null) dserror("Matrix is not a sparse matrix!");
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
Teuchos::RCP<SSI::AssembleStrategyBase> SSI::BuildAssembleStrategy(
    Teuchos::RCP<const SSI::SSI_Mono> ssi_mono, LINALG::MatrixType matrixtype_ssi,
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