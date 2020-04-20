/*----------------------------------------------------------------------*/
/*! \file
\brief Assemble strategy for monolithic SSI
\level 2

\maintainer Stephan Sinzig

 */
/*----------------------------------------------------------------------*/
#include "ssi_monolithic_assemble_strategy.H"

#include "ssi_monolithic.H"

#include "../drt_io/io_control.H"

#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_matrixtransform.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_ssiwrapper.H"

#include "../drt_lib/drt_locsys.H"

#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"


/*----------------------------------------------------------------------*
 |                                                         constructors |
 *----------------------------------------------------------------------*/
SSI::AssembleStrategyBase::AssembleStrategyBase(
    const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono, ADAPTER::CouplingSlaveConverter converter)
    : converter_(converter), ssi_mono_(ssi_mono)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategyBlock::AssembleStrategyBlock(
    const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono, ADAPTER::CouplingSlaveConverter converter)
    : AssembleStrategyBase(ssi_mono, converter)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategyBlockBlock::AssembleStrategyBlockBlock(
    const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono, ADAPTER::CouplingSlaveConverter converter)
    : AssembleStrategyBlock(ssi_mono, converter)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategyBlockSparse::AssembleStrategyBlockSparse(
    const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono, ADAPTER::CouplingSlaveConverter converter)
    : AssembleStrategyBlock(ssi_mono, converter)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategySparse::AssembleStrategySparse(
    const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono, ADAPTER::CouplingSlaveConverter converter)
    : AssembleStrategyBase(ssi_mono, converter)
{
  return;
}

/*----------------------------------------------------------------------*
 |                            assemble scatra domain into system matrix |
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatradomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // cast scatrablock
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> scatradomain_block =
      Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(scatradomain);
  if (scatradomain_block == Teuchos::null) dserror("Scatra block is not a block matrix!");

  const int numberscatrablocks = ssi_mono_->MeshtyingStrategyS2I()->BlockMaps().NumMaps();

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < numberscatrablocks; ++iblock)
    for (int jblock = 0; jblock < numberscatrablocks; ++jblock)
      systemmatrix_block->Assign(
          iblock, jblock, LINALG::View, scatradomain_block->Matrix(iblock, jblock));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatradomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // cast scatrablock
  Teuchos::RCP<LINALG::SparseMatrix> scatradomain_sparse =
      Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(scatradomain);
  if (scatradomain_sparse == Teuchos::null) dserror("System matrix is not a sparse matrix!");

  // add scalar transport system matrix to global system matrix
  systemmatrix_block->Assign(0, 0, LINALG::View, *scatradomain_sparse);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatradomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

  // cast scatrablock
  Teuchos::RCP<LINALG::SparseMatrix> scatradomain_sparse =
      Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(scatradomain);
  if (scatradomain_sparse == Teuchos::null) dserror("System matrix is not a sparse matrix!");

  // add scalar transport system matrix to global system matrix
  systemmatrix_sparse->Add(*scatradomain_sparse, false, 1.0, 0.);

  return;
}

/*----------------------------------------------------------------------*
 |                         assemble structure domain into system matrix |
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    Teuchos::RCP<LINALG::SparseMatrix> structuredomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  const int numberscatrablocks = ssi_mono_->MeshtyingStrategyS2I()->BlockMaps().NumMaps();

  // add entire block or assemble slave side to master side
  if (!ssi_mono_->SSIInterfaceMeshtying())
    systemmatrix_block->Assign(
        numberscatrablocks, numberscatrablocks, LINALG::View, *structuredomain);
  else
    AssembleStructureDomainMeshtying(
        systemmatrix_block->Matrix(numberscatrablocks, numberscatrablocks), structuredomain, false);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    Teuchos::RCP<LINALG::SparseMatrix> structuredomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // add entire block or assemble slave side to master side
  if (!ssi_mono_->SSIInterfaceMeshtying())
    systemmatrix_block->Assign(1, 1, LINALG::View, *structuredomain);
  else
    AssembleStructureDomainMeshtying(systemmatrix_block->Matrix(1, 1), structuredomain, false);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    Teuchos::RCP<LINALG::SparseMatrix> structuredomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

  // add entire block or assemble slave side to master side
  if (!ssi_mono_->SSIInterfaceMeshtying())
    systemmatrix_sparse->Add(*structuredomain, false, 1.0, 1.0);
  else
    AssembleStructureDomainMeshtying(*systemmatrix_sparse, structuredomain, true);

  return;
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

  // assemble interior and master-side rows and columns of structural system matrix into
  // global system matrix
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *mapstructurecondensed,
      *mapstructurecondensed, 1.0, nullptr, nullptr, systemmatrix_structure, true, zero);

  // transform and assemble slave-side rows of structural system matrix into global system
  // matrix
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *mapstructureslave,
      *ssi_mono_->MapStructureCondensed(), 1.0, &converter_, nullptr, systemmatrix_structure, true,
      true);

  // transform and assemble slave-side columns of structural system matrix into global
  // system matrix
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *mapstructurecondensed,
      *mapstructureslave, 1.0, nullptr, &converter_, systemmatrix_structure, true, true);

  // transform and assemble slave-side rows and columns of structural system matrix into
  // global system matrix
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *mapstructureslave, *mapstructureslave,
      1.0, &converter_, &converter_, systemmatrix_structure, true, true);

  return;
}

/*----------------------------------------------------------------------*
 |                  assemble scatra-structure domain into system matrix |
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScatraStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  const int numberscatrablocks = ssi_mono_->MeshtyingStrategyS2I()->BlockMaps().NumMaps();

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < numberscatrablocks; ++iblock)
  {
    // add entire block or assemble slave side to master side
    if (!ssi_mono_->SSIInterfaceMeshtying())
      systemmatrix_block->Assign(iblock, numberscatrablocks, LINALG::View,
          Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(scatrastructuredomain)
              ->Matrix(iblock, 0));
    else
      AssembleScatraStructureDomainMeshtying(systemmatrix_block->Matrix(iblock, numberscatrablocks),
          Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(scatrastructuredomain)
              ->Matrix(iblock, 0),
          false);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleScatraStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // cast scatrastructureblock
  const LINALG::SparseMatrix& scatrastructuredomain_sparse =
      *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(scatrastructuredomain);

  // add entire block or assemble slave side to master side
  if (!ssi_mono_->SSIInterfaceMeshtying())
    systemmatrix_block->Assign(0, 1, LINALG::View, scatrastructuredomain_sparse);
  else
    AssembleScatraStructureDomainMeshtying(
        systemmatrix_block->Matrix(0, 1), scatrastructuredomain_sparse, false);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleScatraStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

  // cast scatrastructureblock
  const LINALG::SparseMatrix& scatrastructuredomain_sparse =
      *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(scatrastructuredomain);

  // add entire block or assemble slave side to master side
  if (!ssi_mono_->SSIInterfaceMeshtying())
    systemmatrix_sparse->Add(scatrastructuredomain_sparse, false, 1.0, 1.0);
  else
    AssembleScatraStructureDomainMeshtying(
        *systemmatrix_sparse, scatrastructuredomain_sparse, true);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::AssembleScatraStructureDomainMeshtying(
    LINALG::SparseMatrix& systemmatrix_scatra_structure, LINALG::SparseMatrix scatrastructuredomain,
    bool zero)
{
  // map for interior and master-side structural degrees of freedom
  const auto& mapstructurecondensed = ssi_mono_->MapStructureCondensed();

  // map for slave side structural degrees of freedom
  const auto& mapstructureslave = ssi_mono_->MapsStructure()->Map(1);

  // assemble interior and master-side columns of scatra-structure block into global
  // system matrix
  LINALG::MatrixLogicalSplitAndTransform()(scatrastructuredomain, scatrastructuredomain.RangeMap(),
      *mapstructurecondensed, 1.0, nullptr, nullptr, systemmatrix_scatra_structure, true, zero);

  // transform and assemble slave-side columns of scatra-structure block into global
  // system matrix
  LINALG::MatrixLogicalSplitAndTransform()(scatrastructuredomain, scatrastructuredomain.RangeMap(),
      *mapstructureslave, 1.0, nullptr, &converter_, systemmatrix_scatra_structure, true, true);

  return;
};

/*----------------------------------------------------------------------*
 |               assemble scatra-structure interface into system matrix |
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScatraStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterfaceslaveside)
{
  const int numberscatrablocks = ssi_mono_->MeshtyingStrategyS2I()->BlockMaps().NumMaps();

  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // master and slave side of scatrastructureinterface
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> scatrastructureinterface_block = Teuchos::rcp(
      new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*ssi_mono_->MapStructure(),
          ssi_mono_->MeshtyingStrategyS2I()->BlockMaps(), 81, false, true));

  // copy slave side to master side and scale with -1.0. Add slave side
  CopyScatraStructureInterfaceS2Iblock(
      scatrastructureinterface_block, scatrastructureinterfaceslaveside);

  // derive linearizations of master-side scatra fluxes w.r.t. master-side structural dofs
  // and assemble into auxiliary system matrix
  for (int iblock = 0; iblock < numberscatrablocks; ++iblock)
  {
    // assemble scatrastructureinterface_sparse into system matrix
    AssembleScatraStructureDomainMeshtying(systemmatrix_block->Matrix(iblock, numberscatrablocks),
        scatrastructureinterface_block->Matrix(iblock, 0), true);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleScatraStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterfaceslaveside)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // master and slave side of scatrastructureinterface
  Teuchos::RCP<LINALG::SparseMatrix> scatrastructureinterface_sparse =
      Teuchos::rcp(new LINALG::SparseMatrix(*ssi_mono_->Maps()->Map(0), 27, false, true));

  // copy slave side to master side and scale with -1.0. Add slave side
  CopyScatraStructureInterfaceS2Isparse(
      scatrastructureinterface_sparse, scatrastructureinterfaceslaveside);

  // assemble scatrastructureinterface_sparse into system matrix
  AssembleScatraStructureDomainMeshtying(
      systemmatrix_block->Matrix(0, 1), *scatrastructureinterface_sparse, true);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleScatraStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterfaceslaveside)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

  // master and slave side of scatrastructureinterface
  Teuchos::RCP<LINALG::SparseMatrix> scatrastructureinterface_sparse =
      Teuchos::rcp(new LINALG::SparseMatrix(*ssi_mono_->Maps()->Map(0), 27, false, true));

  // copy slave side to master side and scale with -1.0. Add slave side
  CopyScatraStructureInterfaceS2Isparse(
      scatrastructureinterface_sparse, scatrastructureinterfaceslaveside);

  // assemble scatrastructureinterface_sparse into system matrix
  AssembleScatraStructureDomainMeshtying(
      *systemmatrix_sparse, *scatrastructureinterface_sparse, true);

  return;
}

/*----------------------------------------------------------------------*
 |                  assemble structure-scatra domain into system matrix |
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleStructureScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> structurescatradomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  const int numberscatrablocks = ssi_mono_->MeshtyingStrategyS2I()->BlockMaps().NumMaps();

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < numberscatrablocks; ++iblock)
  {
    // add entire block or assemble slave side to master side
    if (!ssi_mono_->SSIInterfaceMeshtying())
      systemmatrix_block->Assign(numberscatrablocks, iblock, LINALG::View,
          Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(structurescatradomain)
              ->Matrix(0, iblock));
    else
      AssembleStructureScatraDomainMeshtying(systemmatrix_block->Matrix(numberscatrablocks, iblock),
          Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(structurescatradomain)
              ->Matrix(0, iblock),
          false);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleStructureScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> structurescatradomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // cast structurescatrablock
  const LINALG::SparseMatrix& structurescatradomain_sparse =
      *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(structurescatradomain);

  // add entire block or assemble slave side to master side
  if (!ssi_mono_->SSIInterfaceMeshtying())
    systemmatrix_block->Assign(1, 0, LINALG::View, structurescatradomain_sparse);
  else
    AssembleStructureScatraDomainMeshtying(
        systemmatrix_block->Matrix(1, 0), structurescatradomain_sparse, false);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleStructureScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> structurescatradomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

  // cast structurescatrablock
  const LINALG::SparseMatrix& structurescatradomain_sparse =
      *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(structurescatradomain);

  // add entire block or assemble slave side to master side
  if (!ssi_mono_->SSIInterfaceMeshtying())
    systemmatrix_sparse->Add(structurescatradomain_sparse, false, 1.0, 1.0);
  else
    AssembleStructureScatraDomainMeshtying(
        *systemmatrix_sparse, structurescatradomain_sparse, true);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::AssembleStructureScatraDomainMeshtying(
    LINALG::SparseMatrix& systemmatrix_structure_scatra, LINALG::SparseMatrix structurescatradomain,
    bool zero)
{
  // map for interior and master-side structural degrees of freedom
  const auto& mapstructurecondensed = ssi_mono_->MapStructureCondensed();

  // map for slave side structural degrees of freedom
  const auto& mapstructureslave = ssi_mono_->MapsStructure()->Map(1);

  // assemble interior and master - side rows of structure - scatra block into global system
  // matrix
  LINALG::MatrixLogicalSplitAndTransform()(structurescatradomain, *mapstructurecondensed,
      structurescatradomain.DomainMap(), 1.0, nullptr, nullptr, systemmatrix_structure_scatra, true,
      zero);

  // transform and assemble slave-side rows of structure-scatra block into global system
  // matrix
  LINALG::MatrixLogicalSplitAndTransform()(structurescatradomain, *mapstructureslave,
      structurescatradomain.DomainMap(), 1.0, &converter_, nullptr, systemmatrix_structure_scatra,
      true, true);

  return;
};

/*----------------------------------------------------------------------*
 |                       apply meshtying to the assembled system matrix |
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::ApplyMeshtyingSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  if (ssi_mono_->SSIInterfaceMeshtying())
  {
    const int numberscatrablocks = ssi_mono_->MeshtyingStrategyS2I()->BlockMaps().NumMaps();

    // cast systemmatrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
    CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

    ApplyMeshtyingSysMat(systemmatrix_block->Matrix(numberscatrablocks, numberscatrablocks));
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::ApplyMeshtyingSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  if (ssi_mono_->SSIInterfaceMeshtying())
  {
    // cast systemmatrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
    CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

    ApplyMeshtyingSysMat(systemmatrix_block->Matrix(1, 1));
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::ApplyMeshtyingSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  if (ssi_mono_->SSIInterfaceMeshtying())
  {
    // cast systemmatrix
    Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
    CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

    ApplyMeshtyingSysMat(*systemmatrix_sparse);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::ApplyMeshtyingSysMat(LINALG::SparseMatrix& systemmatrix_structure)
{
  // map for slave side structural degrees of freedom
  const auto& mapstructureslave = ssi_mono_->MapsStructure()->Map(1);

  // subject slave-side rows of structural system matrix to pseudo Dirichlet conditions to
  // finalize structural meshtying
  const double one(1.0);
  for (int doflid_slave = 0; doflid_slave < mapstructureslave->NumMyElements(); ++doflid_slave)
  {
    // extract global ID of current slave-side row
    const int dofgid_slave = mapstructureslave->GID(doflid_slave);
    if (dofgid_slave < 0) dserror("Local ID not found!");

    // apply pseudo Dirichlet conditions to filled matrix, i.e., to local row and column
    // indices
    if (systemmatrix_structure.Filled())
    {
      const int rowlid_slave = systemmatrix_structure.RowMap().LID(dofgid_slave);
      if (rowlid_slave < 0) dserror("Global ID not found!");
      if (systemmatrix_structure.EpetraMatrix()->ReplaceMyValues(
              rowlid_slave, 1, &one, &rowlid_slave))
        dserror("ReplaceMyValues failed!");
    }

    // apply pseudo Dirichlet conditions to unfilled matrix, i.e., to global row and
    // column indices
    else if (systemmatrix_structure.EpetraMatrix()->InsertGlobalValues(
                 dofgid_slave, 1, &one, &dofgid_slave))
      dserror("InsertGlobalValues failed!");
  }

  return;
}

/*----------------------------------------------------------------------*
 |                                apply structural DBC to system matrix |
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlock::ApplyStructuralDBCSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix)
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
    CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

    // apply structural Dirichlet conditions
    for (int iblock = 0; iblock < systemmatrix_block->Cols(); ++iblock)
    {
      locsysmanager_structure->RotateGlobalToLocal(
          Teuchos::rcp(&systemmatrix_block->Matrix(systemmatrix_block->Cols() - 1, iblock), false));
      systemmatrix_block->Matrix(systemmatrix_block->Cols() - 1, iblock)
          .ApplyDirichletWithTrafo(locsysmanager_structure->Trafo(), *dbcmap_structure,
              iblock == systemmatrix_block->Cols() - 1 ? true : false);
      locsysmanager_structure->RotateLocalToGlobal(
          Teuchos::rcp(&systemmatrix_block->Matrix(systemmatrix_block->Cols() - 1, iblock), false));
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::ApplyStructuralDBCSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix)
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
    CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

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

  return;
}

/*----------------------------------------------------------------------*
 |                                             assemble right-hand-side |
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::AssembleRHS(Teuchos::RCP<Epetra_Vector>& RHS,
    Teuchos::RCP<Epetra_Vector> RHSscatra, Teuchos::RCP<const Epetra_Vector> RHSstructure)
{
  // zero out RHS
  RHS->PutScalar(0.);

  // assemble scalar transport right-hand side vector into monolithic right-hand side vector
  ssi_mono_->Maps()->InsertVector(RHSscatra, 0, RHS);

  if (!ssi_mono_->SSIInterfaceMeshtying())
    ssi_mono_->Maps()->AddVector(RHSstructure, 1, RHS, -1.0);
  else
  {
    // perform structural meshtying before assembling structural right-hand side vector into
    // monolithic right-hand side vector

    // make copy of structural right-hand side vector
    Epetra_Vector residual_structure(*RHSstructure);

    // transform slave-side part of structural right-hand side vector to master side
    Teuchos::RCP<Epetra_Vector> slavetomaster = ssi_mono_->MapsStructure()->InsertVector(
        ssi_mono_->CouplingAdapterStructure()->SlaveToMaster(
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
    ssi_mono_->MapsStructure()->PutScalar(residual_structure, 1, 0.);

    // assemble final structural right-hand side vector into monolithic right-hand side vector
    ssi_mono_->Maps()->AddVector(residual_structure, 1, *RHS, -1.0);
  }

  return;
}

/*----------------------------------------------------------------------*
 |                                                   cast system matrix |
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlock::CastSystemMatrixBlock(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    Teuchos::RCP<LINALG::BlockSparseMatrixBase>& systemmatrix_block)
{
  systemmatrix_block = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(systemmatrix);
  if (systemmatrix_block == Teuchos::null) dserror("System matrix is not a block matrix!");

  return;
};

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::CastSystemMatrixSparse(
    Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    Teuchos::RCP<LINALG::SparseMatrix>& systemmatrix_sparse)
{
  systemmatrix_sparse = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix);
  if (systemmatrix_sparse == Teuchos::null) dserror("System matrix is not a sparse matrix!");

  return;
}

/*----------------------------------------------------------------------*
 |                                copy slave side values on master side |
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::CopyScatraStructureInterfaceS2Isparse(
    Teuchos::RCP<LINALG::SparseMatrix>& scatrastructureinterface_sparse,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterfaceslaveside)
{
  // cast scatrastructureinterfaceslaveside
  const LINALG::SparseMatrix& scatrastructureinterfaceslaveside_sparse =
      *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(scatrastructureinterfaceslaveside);

  // insert slave side values into
  scatrastructureinterface_sparse->Add(scatrastructureinterfaceslaveside_sparse, false, 1.0, 0.0);

  // copy slave side values to master side and scale with minus 1. Insert into
  // scatrastructureinterface_sparse
  LINALG::MatrixRowColTransform()(scatrastructureinterfaceslaveside_sparse, -1.0,
      ADAPTER::CouplingSlaveConverter(*ssi_mono_->MeshtyingStrategyS2I()->CouplingAdapter()),
      ADAPTER::CouplingSlaveConverter(*ssi_mono_->CouplingAdapterStructure()),
      *scatrastructureinterface_sparse, true, true);

  // finalize
  scatrastructureinterface_sparse->Complete(*ssi_mono_->Maps()->Map(1), *ssi_mono_->Maps()->Map(0));

  return;
};

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::CopyScatraStructureInterfaceS2Iblock(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase>& scatrastructureinterface_block,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterfaceslaveside)
{
  const int numberscatrablocks = ssi_mono_->MeshtyingStrategyS2I()->BlockMaps().NumMaps();

  // cast scatrastructureinterfaceslaveside
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> scatrastructureinterfaceslaveside_block =
      Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>>(
          scatrastructureinterfaceslaveside);

  // initialize auxiliary system matrix for linearizations of master-side scatra fluxes w.r.t.
  // master-side structural dofs
  LINALG::SparseMatrix mastermatrix(
      *ssi_mono_->MeshtyingStrategyS2I()->CouplingAdapter()->MasterDofMap(), 27, false, true);

  // derive linearizations of master-side scatra fluxes w.r.t. master-side structural dofs and
  // assemble into auxiliary system matrix
  for (int iblock = 0; iblock < numberscatrablocks; ++iblock)
    LINALG::MatrixRowColTransform()(scatrastructureinterfaceslaveside_block->Matrix(iblock, 0),
        -1.0,
        ADAPTER::CouplingSlaveConverter(*ssi_mono_->MeshtyingStrategyS2I()->CouplingAdapter()),
        ADAPTER::CouplingSlaveConverter(*ssi_mono_->CouplingAdapterStructure()), mastermatrix, true,
        true);

  // finalize auxiliary system matrix
  mastermatrix.Complete(*ssi_mono_->CouplingAdapterStructure()->MasterDofMap(),
      *ssi_mono_->MeshtyingStrategyS2I()->CouplingAdapter()->MasterDofMap());

  // split auxiliary system matrix and assemble into scatra-structure matrix block
  const Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockmastermatrix =
      mastermatrix.Split<LINALG::DefaultBlockMatrixStrategy>(
          *ssi_mono_->MapStructure(), ssi_mono_->MeshtyingStrategyS2I()->BlockMaps());
  blockmastermatrix->Complete();

  // build block matrix with master and slave side contributions
  for (int iblock = 0; iblock < numberscatrablocks; ++iblock)
  {
    // add slave side contributions
    scatrastructureinterface_block->Matrix(iblock, 0).Add(
        scatrastructureinterfaceslaveside_block->Matrix(iblock, 0), false, 1.0, 0.0);
    // add master side contributions
    scatrastructureinterface_block->Matrix(iblock, 0).Add(
        blockmastermatrix->Matrix(iblock, 0), false, 1.0, 1.0);
  }

  scatrastructureinterface_block->Complete();

  return;
};
