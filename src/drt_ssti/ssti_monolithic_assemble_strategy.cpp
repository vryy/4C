/*----------------------------------------------------------------------*/
/*! \file
\brief Assemble strategy for monolithic SSTI

\level 2

*----------------------------------------------------------------------*/
#include "ssti_monolithic.H"
#include "ssti_monolithic_assemble_strategy.H"
#include "ssti_utils.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_locsys.H"

#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../linalg/linalg_matrixtransform.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"

/*----------------------------------------------------------------------*
 |                                                         constructors |
 *----------------------------------------------------------------------*/
SSTI::AssembleStrategyBase::AssembleStrategyBase(
    Teuchos::RCP<const SSTI::SSTIMono> ssti_mono, ADAPTER::CouplingSlaveConverter converter)
    : converter_(std::move(converter)), ssti_mono_(std::move(ssti_mono))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSTI::AssembleStrategyBlock::AssembleStrategyBlock(
    Teuchos::RCP<const SSTI::SSTIMono> ssti_mono, ADAPTER::CouplingSlaveConverter converter)
    : AssembleStrategyBase(std::move(ssti_mono), std::move(converter))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSTI::AssembleStrategyBlockBlock::AssembleStrategyBlockBlock(
    Teuchos::RCP<const SSTI::SSTIMono> ssti_mono, ADAPTER::CouplingSlaveConverter converter)
    : AssembleStrategyBlock(std::move(ssti_mono), std::move(converter))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSTI::AssembleStrategyBlockSparse::AssembleStrategyBlockSparse(
    Teuchos::RCP<const SSTI::SSTIMono> ssti_mono, ADAPTER::CouplingSlaveConverter converter)
    : AssembleStrategyBlock(std::move(ssti_mono), std::move(converter))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSTI::AssembleStrategySparse::AssembleStrategySparse(
    Teuchos::RCP<const SSTI::SSTIMono> ssti_mono, ADAPTER::CouplingSlaveConverter converter)
    : AssembleStrategyBase(std::move(ssti_mono), std::move(converter))
{
}

/*----------------------------------------------------------------------*
 |                            assemble scatra domain into system matrix |
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatradomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // cast scatrablock
  const auto scatradomain_block =
      Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(scatradomain);
  if (scatradomain_block == Teuchos::null) dserror("Matrix is not a block matrix!");

  const int numberscatrablocks = ssti_mono_->AllMaps()->MapsScatra()->NumMaps();

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < numberscatrablocks; ++iblock)
  {
    for (int jblock = 0; jblock < numberscatrablocks; ++jblock)
      systemmatrix_block->Assign(
          iblock, jblock, LINALG::View, scatradomain_block->Matrix(iblock, jblock));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatradomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // cast scatrablock
  const auto scatradomain_sparse =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(scatradomain);
  if (scatradomain_sparse == Teuchos::null) dserror("Matrix is not a sparse matrix!");

  // add scalar transport system matrix to global system matrix
  systemmatrix_block->Assign(0, 0, LINALG::View, *scatradomain_sparse);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatradomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

  // cast scatrablock
  const auto scatradomain_sparse =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(scatradomain);
  if (scatradomain_sparse == Teuchos::null) dserror("Matrix is not a sparse matrix!");

  // add scalar transport system matrix to global system matrix
  systemmatrix_sparse->Add(*scatradomain_sparse, false, 1.0, 0.);
}

/*----------------------------------------------------------------------*
 |                         assemble structure domain into system matrix |
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseMatrix> structuredomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  const int numberscatrablocks = ssti_mono_->AllMaps()->MapsScatra()->NumMaps();

  // add entire block or assemble slave side to master side
  if (!ssti_mono_->InterfaceMeshtying())
  {
    systemmatrix_block->Assign(
        numberscatrablocks, numberscatrablocks, LINALG::View, *structuredomain);
  }
  else
  {
    AssembleStructureDomainMeshtying(
        systemmatrix_block->Matrix(numberscatrablocks, numberscatrablocks), structuredomain, false);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseMatrix> structuredomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // add entire block or assemble slave side to master side
  if (!ssti_mono_->InterfaceMeshtying())
    systemmatrix_block->Assign(1, 1, LINALG::View, *structuredomain);
  else
    AssembleStructureDomainMeshtying(systemmatrix_block->Matrix(1, 1), structuredomain, false);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseMatrix> structuredomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

  // add entire block or assemble slave side to master side
  if (!ssti_mono_->InterfaceMeshtying())
    systemmatrix_sparse->Add(*structuredomain, false, 1.0, 1.0);
  else
    AssembleStructureDomainMeshtying(*systemmatrix_sparse, structuredomain, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBase::AssembleStructureDomainMeshtying(
    LINALG::SparseMatrix& systemmatrix_structure,
    Teuchos::RCP<LINALG::SparseMatrix> structuredomain, bool zero)
{
  // map for interior and master-side structural degrees of freedom
  const auto& mapstructurecondensed = ssti_mono_->AllMaps()->MapStructureCondensed();

  // map for slave side structural degrees of freedom
  const auto& mapstructureslave = ssti_mono_->AllMaps()->MapsInterfaceStructure()->Map(0);

  // assemble interior and master-side rows and columns of structural system matrix into
  // global system matrix
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *mapstructurecondensed,
      *mapstructurecondensed, 1.0, nullptr, nullptr, systemmatrix_structure, true, zero);

  // transform and assemble slave-side rows of structural system matrix into global system
  // matrix
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *mapstructureslave,
      *ssti_mono_->AllMaps()->MapStructureCondensed(), 1.0, &converter_, nullptr,
      systemmatrix_structure, true, true);

  // transform and assemble slave-side columns of structural system matrix into global
  // system matrix
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *mapstructurecondensed,
      *mapstructureslave, 1.0, nullptr, &converter_, systemmatrix_structure, true, true);

  // transform and assemble slave-side rows and columns of structural system matrix into
  // global system matrix
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *mapstructureslave, *mapstructureslave,
      1.0, &converter_, &converter_, systemmatrix_structure, true, true);
}

/*----------------------------------------------------------------------*
 |                            assemble thermo domain into system matrix |
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleThermoDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> thermodomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // cast scatrablock
  const auto thermodomain_block =
      Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(thermodomain);
  if (thermodomain_block == Teuchos::null) dserror("Matrix is not a block matrix!");

  const int numberscatrablocks = ssti_mono_->AllMaps()->MapsScatra()->NumMaps();
  const int numberthermoblocks = ssti_mono_->AllMaps()->MapsThermo()->NumMaps();

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < numberthermoblocks; ++iblock)
  {
    for (int jblock = 0; jblock < numberthermoblocks; ++jblock)
      systemmatrix_block->Assign(numberscatrablocks + 1 + iblock, numberscatrablocks + 1 + jblock,
          LINALG::View, thermodomain_block->Matrix(iblock, jblock));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleThermoDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> thermodomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // cast scatrablock
  const auto thermodomain_sparse =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(thermodomain);
  if (thermodomain_sparse == Teuchos::null) dserror("Matrix is not a sparse matrix!");

  // add scalar transport system matrix to global system matrix
  systemmatrix_block->Assign(2, 2, LINALG::View, *thermodomain_sparse);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleThermoDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> thermodomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

  // cast thermoblock
  const auto thermodomain_sparse =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(thermodomain);
  if (thermodomain_sparse == Teuchos::null) dserror("Matrix is not a sparse matrix!");

  // add scalar transport system matrix to global system matrix
  systemmatrix_sparse->Add(*thermodomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 |                  assemble scatra-structure domain into system matrix |
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleScatraStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  const int numberscatrablocks = ssti_mono_->AllMaps()->MapsScatra()->NumMaps();

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < numberscatrablocks; ++iblock)
  {
    const auto scatrastructuredomain_subblock =
        Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(scatrastructuredomain)
            ->Matrix(iblock, 0);

    // add entire block or assemble slave side to master side
    if (!ssti_mono_->InterfaceMeshtying())
    {
      systemmatrix_block->Assign(
          iblock, numberscatrablocks, LINALG::View, scatrastructuredomain_subblock);
    }
    else
    {
      AssembleScatraStructureDomainMeshtying(systemmatrix_block->Matrix(iblock, numberscatrablocks),
          scatrastructuredomain_subblock, false);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleScatraStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // cast scatrastructureblock
  const auto scatrastructuredomain_sparse =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(scatrastructuredomain);
  if (scatrastructuredomain_sparse == Teuchos::null) dserror("Matrix is not a sparse matrix!");

  // add entire block or assemble slave side to master side
  if (!ssti_mono_->InterfaceMeshtying())
    systemmatrix_block->Assign(0, 1, LINALG::View, *scatrastructuredomain_sparse);
  else
    AssembleScatraStructureDomainMeshtying(
        systemmatrix_block->Matrix(0, 1), *scatrastructuredomain_sparse, false);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleScatraStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

  // cast scatrastructureblock
  const auto scatrastructuredomain_sparse =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(scatrastructuredomain);
  if (scatrastructuredomain_sparse == Teuchos::null) dserror("Matrix is not a sparse matrix!");

  // add entire block or assemble slave side to master side
  if (!ssti_mono_->InterfaceMeshtying())
    systemmatrix_sparse->Add(*scatrastructuredomain_sparse, false, 1.0, 1.0);
  else
    AssembleScatraStructureDomainMeshtying(
        *systemmatrix_sparse, *scatrastructuredomain_sparse, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBase::AssembleScatraStructureDomainMeshtying(
    LINALG::SparseMatrix& systemmatrix_scatra_structure,
    const LINALG::SparseMatrix& scatrastructuredomain, bool zero)
{
  // map for interior and master-side structural degrees of freedom
  const auto& mapstructurecondensed = ssti_mono_->AllMaps()->MapStructureCondensed();

  // map for slave side structural degrees of freedom
  const auto& mapstructureslave = ssti_mono_->AllMaps()->MapsInterfaceStructure()->Map(0);

  // assemble interior and master-side columns of scatra-structure block into global
  // system matrix
  LINALG::MatrixLogicalSplitAndTransform()(scatrastructuredomain, scatrastructuredomain.RangeMap(),
      *mapstructurecondensed, 1.0, nullptr, nullptr, systemmatrix_scatra_structure, true, zero);

  // transform and assemble slave-side columns of scatra-structure block into global
  // system matrix
  LINALG::MatrixLogicalSplitAndTransform()(scatrastructuredomain, scatrastructuredomain.RangeMap(),
      *mapstructureslave, 1.0, nullptr, &converter_, systemmatrix_scatra_structure, true, true);
}

/*----------------------------------------------------------------------*
 |               assemble scatra-structure interface into system matrix |
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleScatraStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface)
{
  const int numberscatrablocks = ssti_mono_->AllMaps()->MapsScatra()->NumMaps();

  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // derive linearizations of master-side scatra fluxes w.r.t. master-side structural dofs
  // and assemble into auxiliary system matrix
  for (int iblock = 0; iblock < numberscatrablocks; ++iblock)
  {
    const auto scatrastructureinterface_subblock =
        Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(scatrastructureinterface)
            ->Matrix(iblock, 0);

    // assemble scatrastructureinterface_sparse into system matrix
    AssembleScatraStructureDomainMeshtying(systemmatrix_block->Matrix(iblock, numberscatrablocks),
        scatrastructureinterface_subblock, true);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleScatraStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // assemble scatrastructureinterface_sparse into system matrix
  AssembleScatraStructureDomainMeshtying(systemmatrix_block->Matrix(0, 1),
      *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(scatrastructureinterface), true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleScatraStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

  // assemble scatrastructureinterface_sparse into system matrix
  AssembleScatraStructureDomainMeshtying(*systemmatrix_sparse,
      *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(scatrastructureinterface), true);
}

/*----------------------------------------------------------------------*
 |                  assemble scatra-thermo interface into system matrix |
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleScatraThermoInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrathermointerface)
{
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  const auto scatrathermointerface_block =
      Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(scatrathermointerface);
  if (scatrathermointerface_block == Teuchos::null) dserror("Matrix is not a block matrix!");

  const int numberscatrablocks = ssti_mono_->AllMaps()->MapsScatra()->NumMaps();
  const int numberthermoblocks = ssti_mono_->AllMaps()->MapsThermo()->NumMaps();

  const auto scatrainterface = ssti_mono_->AllMaps()->MapInterface(ssti_mono_->MeshtyingScatra());

  LINALG::SparseMatrix masterderiv(*scatrainterface, 27, false, true);

  for (int i = 0; i < numberscatrablocks; ++i)
  {
    for (int j = 0; j < numberthermoblocks; ++j)
    {
      const auto scatrathermointerface_subblock = scatrathermointerface_block->Matrix(i, j);

      // assemble linearizations of slave- and master side scatra fluxes w.r.t. slave temperatures
      // into system matrix
      LINALG::MatrixLogicalSplitAndTransform()(scatrathermointerface_subblock,
          scatrathermointerface_subblock.RangeMap(), scatrathermointerface_subblock.DomainMap(),
          1.0, nullptr, nullptr, systemmatrix_block->Matrix(i, numberscatrablocks + 1 + j), true,
          true);

      // assemble linearizations of slave- and master side scatra fluxes w.r.t. master temperatures
      // into system matrix
      ADAPTER::CouplingSlaveConverter thermo_converter(
          *ssti_mono_->MeshtyingThermo()->CouplingAdapter());

      LINALG::MatrixLogicalSplitAndTransform()(scatrathermointerface_subblock,
          scatrathermointerface_subblock.RangeMap(),
          *ssti_mono_->MeshtyingThermo()->CouplingAdapter()->MasterDofMap(), 1.0, nullptr,
          &thermo_converter, masterderiv, true, true);
    }
  }

  masterderiv.Complete(
      *ssti_mono_->MeshtyingThermo()->CouplingAdapter()->MasterDofMap(), *scatrainterface);

  const auto blockmasterderiv = masterderiv.Split<LINALG::DefaultBlockMatrixStrategy>(
      *ssti_mono_->AllMaps()->MapsThermo(), *ssti_mono_->AllMaps()->MapsScatra());

  blockmasterderiv->Complete();

  for (int i = 0; i < numberscatrablocks; ++i)
  {
    for (int j = 0; j < numberthermoblocks; ++j)
    {
      const auto masterderiv_subblock = blockmasterderiv->Matrix(i, j);

      // assemble linearizations of slave side scatra fluxes w.r.t. slave and master side elch
      // into system matrix
      LINALG::MatrixLogicalSplitAndTransform()(masterderiv_subblock,
          masterderiv_subblock.RangeMap(), masterderiv_subblock.DomainMap(), 1.0, nullptr, nullptr,
          systemmatrix_block->Matrix(i, numberscatrablocks + 1 + j), true, true);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleScatraThermoInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrathermointerface)
{
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  const auto scatrathermointerface_sparse =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(scatrathermointerface);
  if (scatrathermointerface_sparse == Teuchos::null) dserror("Matrix is not a sparse matrix!");

  // assemble linearizations of slave- and master side scatra fluxes w.r.t. slave temperatures into
  // system matrix
  LINALG::MatrixLogicalSplitAndTransform()(*scatrathermointerface_sparse,
      scatrathermointerface_sparse->RangeMap(), scatrathermointerface_sparse->DomainMap(), 1.0,
      nullptr, nullptr, systemmatrix_block->Matrix(0, 2), true, true);

  // assemble linearizations of slave- and master side scatra fluxes w.r.t. master temperatures into
  // system matrix
  ADAPTER::CouplingSlaveConverter thermo_converter(
      *ssti_mono_->MeshtyingThermo()->CouplingAdapter());

  LINALG::MatrixLogicalSplitAndTransform()(*scatrathermointerface_sparse,
      scatrathermointerface_sparse->RangeMap(),
      *ssti_mono_->MeshtyingThermo()->CouplingAdapter()->MasterDofMap(), 1.0, nullptr,
      &thermo_converter, systemmatrix_block->Matrix(0, 2), true, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleScatraThermoInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrathermointerface)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

  // cast thermoblock
  const auto scatrathermointerface_sparse =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(scatrathermointerface);
  if (scatrathermointerface_sparse == Teuchos::null) dserror("Matrix is not a sparse matrix!");

  // assemble linearizations of slave- and master side scatra fluxes w.r.t. slave temperatures into
  // system matrix
  LINALG::MatrixLogicalSplitAndTransform()(*scatrathermointerface_sparse,
      scatrathermointerface_sparse->RangeMap(), scatrathermointerface_sparse->DomainMap(), 1.0,
      nullptr, nullptr, *systemmatrix_sparse, true, true);

  // assemble linearizations of slave- and master side scatra fluxes w.r.t. master temperatures into
  // system matrix
  ADAPTER::CouplingSlaveConverter thermo_converter(
      *ssti_mono_->MeshtyingThermo()->CouplingAdapter());

  LINALG::MatrixLogicalSplitAndTransform()(*scatrathermointerface_sparse,
      scatrathermointerface_sparse->RangeMap(),
      *ssti_mono_->MeshtyingThermo()->CouplingAdapter()->MasterDofMap(), 1.0, nullptr,
      &thermo_converter, *systemmatrix_sparse, true, true);
}

/*----------------------------------------------------------------------*
 |                  assemble structure-scatra domain into system matrix |
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleStructureScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> structurescatradomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  const int numberscatrablocks = ssti_mono_->AllMaps()->MapsScatra()->NumMaps();

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < numberscatrablocks; ++iblock)
  {
    const auto structurescatradomain_subblock =
        Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(structurescatradomain)
            ->Matrix(0, iblock);

    // add entire block or assemble slave side to master side
    if (!ssti_mono_->InterfaceMeshtying())
    {
      systemmatrix_block->Assign(
          numberscatrablocks, iblock, LINALG::View, structurescatradomain_subblock);
    }
    else
    {
      AssembleStructureScatraDomainMeshtying(systemmatrix_block->Matrix(numberscatrablocks, iblock),
          structurescatradomain_subblock, false);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleStructureScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> structurescatradomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // cast structurescatrablock
  const auto structurescatradomain_sparse =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(structurescatradomain);
  if (structurescatradomain_sparse == Teuchos::null) dserror("Matrix is not a sparse matrix!");

  // add entire block or assemble slave side to master side
  if (!ssti_mono_->InterfaceMeshtying())
    systemmatrix_block->Assign(1, 0, LINALG::View, *structurescatradomain_sparse);
  else
    AssembleStructureScatraDomainMeshtying(
        systemmatrix_block->Matrix(1, 0), *structurescatradomain_sparse, false);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleStructureScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> structurescatradomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

  // cast structurescatrablock
  const auto structurescatradomain_sparse =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(structurescatradomain);
  if (structurescatradomain_sparse == Teuchos::null) dserror("Matrix is not a sparse matrix!");

  // add entire block or assemble slave side to master side
  if (!ssti_mono_->InterfaceMeshtying())
    systemmatrix_sparse->Add(*structurescatradomain_sparse, false, 1.0, 1.0);
  else
    AssembleStructureScatraDomainMeshtying(
        *systemmatrix_sparse, *structurescatradomain_sparse, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBase::AssembleStructureScatraDomainMeshtying(
    LINALG::SparseMatrix& systemmatrix_structure_scatra,
    const LINALG::SparseMatrix& structurescatradomain, bool zero)
{
  // map for interior and master-side structural degrees of freedom
  const auto& mapstructurecondensed = ssti_mono_->AllMaps()->MapStructureCondensed();

  // map for slave side structural degrees of freedom
  const auto& mapstructureslave = ssti_mono_->AllMaps()->MapsInterfaceStructure()->Map(0);

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
}
/*----------------------------------------------------------------------*
 |                     assemble thermo-scatra domain into system matrix |
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleThermoScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> thermoscatradomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  const int numberscatrablocks = ssti_mono_->AllMaps()->MapsScatra()->NumMaps();
  const int numberthermoblocks = ssti_mono_->AllMaps()->MapsThermo()->NumMaps();

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < numberthermoblocks; ++iblock)
  {
    for (int jblock = 0; jblock < numberscatrablocks; ++jblock)
    {
      const auto thermoscatradomain_subblock =
          Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(thermoscatradomain)
              ->Matrix(iblock, jblock);
      LINALG::MatrixLogicalSplitAndTransform()(thermoscatradomain_subblock,
          thermoscatradomain_subblock.RangeMap(), thermoscatradomain_subblock.DomainMap(), 1.0,
          nullptr, nullptr, systemmatrix_block->Matrix(numberscatrablocks + 1 + iblock, jblock),
          true, true);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleThermoScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> thermoscatradomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // cast scatrablock
  const auto thermoscatradomain_sparse =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(thermoscatradomain);
  if (thermoscatradomain_sparse == Teuchos::null) dserror("System matrix is not a sparse matrix!");

  // add scalar transport system matrix to global system matrix
  LINALG::MatrixLogicalSplitAndTransform()(*thermoscatradomain_sparse,
      thermoscatradomain_sparse->RangeMap(), thermoscatradomain_sparse->DomainMap(), 1.0, nullptr,
      nullptr, systemmatrix_block->Matrix(2, 0), true, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleThermoScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> thermoscatradomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

  // cast thermoblock
  const auto thermoscatradomain_sparse =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(thermoscatradomain);
  if (thermoscatradomain_sparse == Teuchos::null) dserror("System matrix is not a sparse matrix!");

  // add scalar transport system matrix to global system matrix
  systemmatrix_sparse->Add(*thermoscatradomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 |                     assemble thermo-scatra domain into system matrix |
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleThermoScatraInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> thermoscatrainterface)
{
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  const int numberscatrablocks = ssti_mono_->AllMaps()->MapsScatra()->NumMaps();
  const int numberthermoblocks = ssti_mono_->AllMaps()->MapsThermo()->NumMaps();

  LINALG::SparseMatrix masterflux(
      *ssti_mono_->MeshtyingThermo()->CouplingAdapter()->MasterDofMap(), 27, false, true);

  for (int i = 0; i < numberthermoblocks; ++i)
  {
    for (int j = 0; j < numberscatrablocks; ++j)
    {
      const auto scatrathermointerface_subblock =
          Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(thermoscatrainterface)
              ->Matrix(i, j);

      // assemble linearizations of slave side scatra fluxes w.r.t. slave and master side elch
      // into system matrix
      LINALG::MatrixLogicalSplitAndTransform()(scatrathermointerface_subblock,
          scatrathermointerface_subblock.RangeMap(), scatrathermointerface_subblock.DomainMap(),
          1.0, nullptr, nullptr, systemmatrix_block->Matrix(numberscatrablocks + 1 + i, j), true,
          true);

      // assemble linearizations of master side scatra fluxes w.r.t. slave and master side elch
      // into system matrix
      ADAPTER::CouplingSlaveConverter thermo_converter(
          *ssti_mono_->MeshtyingThermo()->CouplingAdapter());

      LINALG::SparseMatrix slaveflux(
          *ssti_mono_->MeshtyingThermo()->BlockMapsSlave().Map(i), 27, false, true);

      SCATRA::MeshtyingStrategyS2I::ExtractMatrixRows(scatrathermointerface_subblock, slaveflux,
          *ssti_mono_->MeshtyingThermo()->BlockMapsSlave().Map(i));

      slaveflux.Complete(*ssti_mono_->AllMaps()->MapsScatra()->FullMap(),
          *ssti_mono_->MeshtyingThermo()->BlockMapsSlave().Map(i));

      LINALG::MatrixLogicalSplitAndTransform()(scatrathermointerface_subblock,
          *ssti_mono_->MeshtyingThermo()->CouplingAdapter()->MasterDofMap(),
          scatrathermointerface_subblock.DomainMap(), -1.0, &thermo_converter, nullptr, masterflux,
          true, true);
    }
  }

  masterflux.Complete(*ssti_mono_->AllMaps()->MapsScatra()->FullMap(),
      *ssti_mono_->MeshtyingThermo()->CouplingAdapter()->MasterDofMap());

  const auto blockmasterflux = masterflux.Split<LINALG::DefaultBlockMatrixStrategy>(
      *ssti_mono_->AllMaps()->MapsScatra(), *ssti_mono_->AllMaps()->MapsThermo());

  blockmasterflux->Complete();

  for (int i = 0; i < numberthermoblocks; ++i)
  {
    for (int j = 0; j < numberscatrablocks; ++j)
    {
      auto blockmasterflux_subblock = blockmasterflux->Matrix(i, j);

      // assemble linearizations of slave side scatra fluxes w.r.t. slave and master side elch
      // into system matrix
      LINALG::MatrixLogicalSplitAndTransform()(blockmasterflux_subblock,
          blockmasterflux_subblock.RangeMap(), blockmasterflux_subblock.DomainMap(), 1.0, nullptr,
          nullptr, systemmatrix_block->Matrix(numberscatrablocks + 1 + i, j), true, true);
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleThermoScatraInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> thermoscatrainterface)
{
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);
  const auto thermoscatrainterface_sparse =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(thermoscatrainterface);
  if (thermoscatrainterface_sparse == Teuchos::null) dserror("Matrix is not a sparse matrix!");

  // assemble linearizations of slave side scatra fluxes w.r.t. slave and master side elch
  // into system matrix
  LINALG::MatrixLogicalSplitAndTransform()(*thermoscatrainterface_sparse,
      thermoscatrainterface_sparse->RangeMap(), thermoscatrainterface_sparse->DomainMap(), 1.0,
      nullptr, nullptr, systemmatrix_block->Matrix(2, 0), true, true);

  // assemble linearizations of master side thermo fluxes w.r.t. slave and master side elch
  // into system matrix
  ADAPTER::CouplingSlaveConverter thermo_converter(
      *ssti_mono_->MeshtyingThermo()->CouplingAdapter());

  LINALG::MatrixLogicalSplitAndTransform()(*thermoscatrainterface_sparse,
      *ssti_mono_->MeshtyingThermo()->CouplingAdapter()->MasterDofMap(),
      thermoscatrainterface_sparse->DomainMap(), -1.0, &thermo_converter, nullptr,
      systemmatrix_block->Matrix(2, 0), true, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleThermoScatraInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> thermoscatrainterface)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

  // cast thermoblock
  const auto thermoscatrainterface_sparse =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(thermoscatrainterface);
  if (thermoscatrainterface_sparse == Teuchos::null) dserror("Matrix is not a sparse matrix!");

  // assemble linearizations of slave side scatra fluxes w.r.t. slave and master side elch
  // into system matrix
  LINALG::MatrixLogicalSplitAndTransform()(*thermoscatrainterface_sparse,
      thermoscatrainterface_sparse->RangeMap(), thermoscatrainterface_sparse->DomainMap(), 1.0,
      nullptr, nullptr, *systemmatrix_sparse, true, true);

  // assemble linearizations of master side scatra fluxes w.r.t. slave and master side elch
  // into system matrix
  ADAPTER::CouplingSlaveConverter thermo_converter(
      *ssti_mono_->MeshtyingThermo()->CouplingAdapter());

  LINALG::MatrixLogicalSplitAndTransform()(*thermoscatrainterface_sparse,
      *ssti_mono_->MeshtyingThermo()->CouplingAdapter()->MasterDofMap(),
      thermoscatrainterface_sparse->DomainMap(), -1.0, &thermo_converter, nullptr,
      *systemmatrix_sparse, true, true);
}

/*----------------------------------------------------------------------*
 |                  assemble thermo-structure domain into system matrix |
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleThermoStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> thermostructuredomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  const int numberscatrablocks = ssti_mono_->AllMaps()->MapsScatra()->NumMaps();
  const int numberthermoblocks = ssti_mono_->AllMaps()->MapsThermo()->NumMaps();

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < numberthermoblocks; ++iblock)
  {
    const auto thermostructuredomain_subblock =
        Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(thermostructuredomain)
            ->Matrix(iblock, 0);

    // add entire block or assemble slave side to master side
    if (!ssti_mono_->InterfaceMeshtying())
    {
      systemmatrix_block->Assign(numberscatrablocks + 1 + iblock, numberscatrablocks, LINALG::View,
          thermostructuredomain_subblock);
    }
    else
    {
      AssembleThermoStructureDomainMeshtying(
          systemmatrix_block->Matrix(numberscatrablocks + 1 + iblock, numberscatrablocks),
          thermostructuredomain_subblock, false);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleThermoStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> thermostructuredomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // cast scatrastructureblock
  const auto thermostructuredomain_sparse =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(thermostructuredomain);
  if (thermostructuredomain_sparse == Teuchos::null) dserror("Matrix is not a sparse matrix!");

  // add entire block or assemble slave side to master side
  if (!ssti_mono_->InterfaceMeshtying())
    systemmatrix_block->Assign(2, 1, LINALG::View, *thermostructuredomain_sparse);
  else
    AssembleThermoStructureDomainMeshtying(
        systemmatrix_block->Matrix(2, 1), *thermostructuredomain_sparse, false);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleThermoStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> thermostructuredomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

  // cast scatrablock
  const auto thermostructureblock_sparse =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(thermostructuredomain);
  if (thermostructureblock_sparse == Teuchos::null) dserror("Matrix is not a sparse matrix!");

  if (!ssti_mono_->InterfaceMeshtying())
    systemmatrix_sparse->Add(*thermostructureblock_sparse, false, 1.0, 1.0);
  else
    AssembleThermoStructureDomainMeshtying(
        *systemmatrix_sparse, *thermostructureblock_sparse, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBase::AssembleThermoStructureDomainMeshtying(
    LINALG::SparseMatrix& systemmatrix_thermo_structure,
    const LINALG::SparseMatrix& thermostructuredomain, bool zero)
{
  // map for interior and master-side structural degrees of freedom
  const auto& mapstructurecondensed = ssti_mono_->AllMaps()->MapStructureCondensed();

  // map for slave side structural degrees of freedom
  const auto& mapstructureslave = ssti_mono_->AllMaps()->MapsInterfaceStructure()->Map(0);

  // assemble interior and master-side columns of scatra-structure block into global
  // system matrix
  LINALG::MatrixLogicalSplitAndTransform()(thermostructuredomain, thermostructuredomain.RangeMap(),
      *mapstructurecondensed, 1.0, nullptr, nullptr, systemmatrix_thermo_structure, true, zero);

  // transform and assemble slave-side columns of scatra-structure block into global
  // system matrix
  LINALG::MatrixLogicalSplitAndTransform()(thermostructuredomain, thermostructuredomain.RangeMap(),
      *mapstructureslave, 1.0, nullptr, &converter_, systemmatrix_thermo_structure, true, true);
}

/*----------------------------------------------------------------------*
 |               assemble thermo-structure interface into system matrix |
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleThermoStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> thermostructureinterface)
{
  const int numberthermoblocks = ssti_mono_->AllMaps()->MapsScatra()->NumMaps();
  const int numberscatrablocks = ssti_mono_->AllMaps()->MapsScatra()->NumMaps();

  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // derive linearizations of master-side scatra fluxes w.r.t. master-side structural dofs
  // and assemble into auxiliary system matrix
  for (int iblock = 0; iblock < numberthermoblocks; ++iblock)
  {
    const auto thermostructureinterface_subblock =
        Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(thermostructureinterface)
            ->Matrix(iblock, 0);
    // assemble scatrastructureinterface_sparse into system matrix
    AssembleThermoStructureDomainMeshtying(
        systemmatrix_block->Matrix(numberscatrablocks + 1 + iblock, numberscatrablocks),
        thermostructureinterface_subblock, true);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleThermoStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> thermostructureinterface)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // assemble scatrastructureinterface_sparse into system matrix
  AssembleThermoStructureDomainMeshtying(systemmatrix_block->Matrix(2, 1),
      *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(thermostructureinterface), true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleThermoStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> thermostructureinterface)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

  // assemble scatrastructureinterface_sparse into system matrix
  AssembleThermoStructureDomainMeshtying(*systemmatrix_sparse,
      *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(thermostructureinterface), true);
}

/*----------------------------------------------------------------------*
 |                  assemble structure-thermo domain into system matrix |
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleStructureThermoDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> structurethermodomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  const int numberthermoblocks = ssti_mono_->AllMaps()->MapsThermo()->NumMaps();
  const int numberscatrablocks = ssti_mono_->AllMaps()->MapsScatra()->NumMaps();

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < numberthermoblocks; ++iblock)
  {
    const auto structurethermodomain_subblock =
        Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(structurethermodomain)
            ->Matrix(0, iblock);

    // add entire block or assemble slave side to master side
    if (!ssti_mono_->InterfaceMeshtying())
    {
      systemmatrix_block->Assign(numberscatrablocks, numberscatrablocks + 1 + iblock, LINALG::View,
          structurethermodomain_subblock);
    }
    else
    {
      AssembleStructureThermoDomainMeshtying(
          systemmatrix_block->Matrix(numberscatrablocks, numberscatrablocks + 1 + iblock),
          structurethermodomain_subblock, false);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleStructureThermoDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> structurethermodomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
  CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

  // cast structurescatrablock
  const auto structurethermodomain_sparse =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(structurethermodomain);
  if (structurethermodomain_sparse == Teuchos::null) dserror("Matrix is not a sparse matrix!");

  // add entire block or assemble slave side to master side
  if (!ssti_mono_->InterfaceMeshtying())
    systemmatrix_block->Assign(1, 2, LINALG::View, *structurethermodomain_sparse);
  else
    AssembleStructureThermoDomainMeshtying(
        systemmatrix_block->Matrix(1, 2), *structurethermodomain_sparse, false);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleStructureThermoDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> structurethermodomain)
{
  // cast systemmatrix
  Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
  CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

  // cast structurescatrablock
  const auto structurethermodomain_sparse =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(structurethermodomain);
  if (structurethermodomain_sparse == Teuchos::null) dserror("Matrix is not a sparse matrix!");

  // add entire block or assemble slave side to master side
  if (!ssti_mono_->InterfaceMeshtying())
    systemmatrix_sparse->Add(*structurethermodomain_sparse, false, 1.0, 1.0);
  else
    AssembleStructureScatraDomainMeshtying(
        *systemmatrix_sparse, *structurethermodomain_sparse, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBase::AssembleStructureThermoDomainMeshtying(
    LINALG::SparseMatrix& systemmatrix_structure_thermo,
    const LINALG::SparseMatrix& structurethermodomain, bool zero)
{
  // map for interior and master-side structural degrees of freedom
  const auto& mapstructurecondensed = ssti_mono_->AllMaps()->MapStructureCondensed();

  // map for slave side structural degrees of freedom
  const auto& mapstructureslave = ssti_mono_->AllMaps()->MapsInterfaceStructure()->Map(0);

  // assemble interior and master - side rows of structure - scatra block into global system
  // matrix
  LINALG::MatrixLogicalSplitAndTransform()(structurethermodomain, *mapstructurecondensed,
      structurethermodomain.DomainMap(), 1.0, nullptr, nullptr, systemmatrix_structure_thermo, true,
      zero);

  // transform and assemble slave-side rows of structure-scatra block into global system
  // matrix
  LINALG::MatrixLogicalSplitAndTransform()(structurethermodomain, *mapstructureslave,
      structurethermodomain.DomainMap(), 1.0, &converter_, nullptr, systemmatrix_structure_thermo,
      true, true);
};

/*----------------------------------------------------------------------*
 |                       apply meshtying to the assembled system matrix |
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::ApplyMeshtyingSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  if (ssti_mono_->InterfaceMeshtying())
  {
    const int numberscatrablocks = ssti_mono_->AllMaps()->MapsScatra()->NumMaps();

    // cast systemmatrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
    CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

    ApplyMeshtyingSysMat(systemmatrix_block->Matrix(numberscatrablocks, numberscatrablocks));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::ApplyMeshtyingSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  if (ssti_mono_->InterfaceMeshtying())
  {
    // cast systemmatrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
    CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

    ApplyMeshtyingSysMat(systemmatrix_block->Matrix(1, 1));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::ApplyMeshtyingSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  if (ssti_mono_->InterfaceMeshtying())
  {
    // cast systemmatrix
    Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_sparse;
    CastSystemMatrixSparse(systemmatrix, systemmatrix_sparse);

    ApplyMeshtyingSysMat(*systemmatrix_sparse);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBase::ApplyMeshtyingSysMat(LINALG::SparseMatrix& systemmatrix_structure)
{
  // map for slave side structural degrees of freedom
  const auto& mapstructureslave = ssti_mono_->AllMaps()->MapsInterfaceStructure()->Map(0);

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
}

/*----------------------------------------------------------------------*
 |                                apply structural DBC to system matrix |
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlock::ApplyStructuralDBCSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  // locsys manager of strucutre
  const auto& locsysmanager_structure = ssti_mono_->StructureField()->LocsysManager();

  // map of strucutral Dirichlet BCs
  const auto dbcmap_structure = ssti_mono_->StructureField()->GetDBCMapExtractor()->CondMap();

  if (locsysmanager_structure == Teuchos::null)
    systemmatrix->ApplyDirichlet(*dbcmap_structure);
  else
  {
    // cast systemmatrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_block;
    CastSystemMatrixBlock(systemmatrix, systemmatrix_block);

    const int numscatrablocks = ssti_mono_->AllMaps()->MapsScatra()->NumMaps();

    // apply structural Dirichlet conditions
    for (int iblock = 0; iblock < systemmatrix_block->Cols(); ++iblock)
    {
      locsysmanager_structure->RotateGlobalToLocal(
          Teuchos::rcp(&systemmatrix_block->Matrix(numscatrablocks, iblock), false));
      systemmatrix_block->Matrix(numscatrablocks, iblock)
          .ApplyDirichletWithTrafo(
              locsysmanager_structure->Trafo(), *dbcmap_structure, iblock == numscatrablocks);
      locsysmanager_structure->RotateLocalToGlobal(
          Teuchos::rcp(&systemmatrix_block->Matrix(numscatrablocks, iblock), false));
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::ApplyStructuralDBCSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  // locsys manager of strucutre
  const auto& locsysmanager_structure = ssti_mono_->StructureField()->LocsysManager();

  // map of strucutral Dirichlet BCs
  const auto& dbcmap_structure = ssti_mono_->StructureField()->GetDBCMapExtractor()->CondMap();

  // structural dof row map
  const auto& dofrowmap_structure = ssti_mono_->StructureField()->DofRowMap();

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
}

/*----------------------------------------------------------------------*
 |                                             assemble right-hand-side |
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBase::AssembleRHS(Teuchos::RCP<Epetra_Vector> RHS,
    Teuchos::RCP<Epetra_Vector> RHSscatra, Teuchos::RCP<const Epetra_Vector> RHSstructure,
    Teuchos::RCP<const Epetra_Vector> RHSthermo)
{
  // zero out RHS
  RHS->PutScalar(0.);

  // assemble scalar transport right-hand side vector into monolithic right-hand side vector
  ssti_mono_->AllMaps()->MapsSubproblems()->InsertVector(RHSscatra, 0, RHS);
  ssti_mono_->AllMaps()->MapsSubproblems()->InsertVector(RHSthermo, 2, RHS);

  if (!ssti_mono_->InterfaceMeshtying())
    ssti_mono_->AllMaps()->MapsSubproblems()->AddVector(RHSstructure, 1, RHS, -1.0);
  else
  {
    // perform structural meshtying before assembling structural right-hand side vector into
    // monolithic right-hand side vector

    // make copy of structural right-hand side vector
    Epetra_Vector residual_structure(*RHSstructure);

    // transform slave-side part of structural right-hand side vector to master side
    Teuchos::RCP<Epetra_Vector> slavetomaster =
        ssti_mono_->AllMaps()->MapsInterfaceStructure()->InsertVector(
            ssti_mono_->CouplingAdapterStructure()->SlaveToMaster(
                ssti_mono_->AllMaps()->MapsInterfaceStructure()->ExtractVector(
                    residual_structure, 0)),
            1);

    // locsys manager of strucutre
    const auto& locsysmanager_structure = ssti_mono_->StructureField()->LocsysManager();

    // apply pseudo Dirichlet conditions to transformed slave-side part of structural right-hand
    // side vector
    const Teuchos::RCP<const Epetra_Vector> zeros =
        Teuchos::rcp(new Epetra_Vector(slavetomaster->Map()));
    if (locsysmanager_structure != Teuchos::null)
      locsysmanager_structure->RotateGlobalToLocal(slavetomaster);
    LINALG::ApplyDirichlettoSystem(
        slavetomaster, zeros, *ssti_mono_->StructureField()->GetDBCMapExtractor()->CondMap());
    if (locsysmanager_structure != Teuchos::null)
      locsysmanager_structure->RotateLocalToGlobal(slavetomaster);

    // assemble transformed slave-side part of structural right-hand side vector
    residual_structure.Update(1.0, *slavetomaster, 1.0);

    // zero out slave-side part of structural right-hand side vector
    ssti_mono_->AllMaps()->MapsInterfaceStructure()->PutScalar(residual_structure, 0, 0.0);

    // assemble final structural right-hand side vector into monolithic right-hand side vector
    ssti_mono_->AllMaps()->MapsSubproblems()->AddVector(residual_structure, 1, *RHS, -1.0);
  }
}

/*----------------------------------------------------------------------*
 |                                                   cast system matrix |
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlock::CastSystemMatrixBlock(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::BlockSparseMatrixBase>& systemmatrix_block)
{
  systemmatrix_block = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(systemmatrix);
  if (systemmatrix_block == Teuchos::null) dserror("System matrix is not a block matrix!");
};

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::CastSystemMatrixSparse(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseMatrix>& systemmatrix_sparse)
{
  systemmatrix_sparse = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix);
  if (systemmatrix_sparse == Teuchos::null) dserror("System matrix is not a sparse matrix!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<SSTI::AssembleStrategyBase> SSTI::BuildAssembleStrategy(
    Teuchos::RCP<const SSTI::SSTIMono> ssti_mono, const ADAPTER::CouplingSlaveConverter& converter,
    LINALG::MatrixType matrixtype_ssti, LINALG::MatrixType matrixtype_scatra)
{
  Teuchos::RCP<SSTI::AssembleStrategyBase> assemblestrategy = Teuchos::null;

  switch (matrixtype_ssti)
  {
    case LINALG::MatrixType::block_field:
    {
      switch (matrixtype_scatra)
      {
        case LINALG::MatrixType::block_condition:
        {
          assemblestrategy =
              Teuchos::rcp(new SSTI::AssembleStrategyBlockBlock(ssti_mono, converter));
          break;
        }
        case LINALG::MatrixType::sparse:
        {
          assemblestrategy =
              Teuchos::rcp(new SSTI::AssembleStrategyBlockSparse(ssti_mono, converter));
          break;
        }

        default:
        {
          dserror("unknown matrix type");
          break;
        }
      }
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      assemblestrategy = Teuchos::rcp(new SSTI::AssembleStrategySparse(ssti_mono, converter));
      break;
    }
    default:
    {
      dserror("unknown matrix type");
      break;
    }
  }

  return assemblestrategy;
}
