/*----------------------------------------------------------------------*/
/*! \file
\brief Assemble strategy for monolithic SSTI

\level 2

*----------------------------------------------------------------------*/
#include "ssti_monolithic.H"
#include "ssti_monolithic_assemble_strategy.H"
#include "ssti_utils.H"

#include "ad_str_ssiwrapper.H"

#include "io_control.H"

#include "drt_locsys.H"

#include "scatra_timint_meshtying_strategy_s2i.H"

#include "ssi_utils.H"

#include "linalg_matrixtransform.H"
#include "linalg_utils_sparse_algebra_create.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSTI::AssembleStrategyBase::AssembleStrategyBase(Teuchos::RCP<const SSTI::SSTIMono> ssti_mono)
    : ssti_mono_(std::move(ssti_mono))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSTI::AssembleStrategyBlock::AssembleStrategyBlock(Teuchos::RCP<const SSTI::SSTIMono> ssti_mono)
    : AssembleStrategyBase(ssti_mono),
      block_position_scatra_(ssti_mono->GetBlockPositions(SSTI::Subproblem::scalar_transport)),
      block_position_thermo_(ssti_mono->GetBlockPositions(SSTI::Subproblem::thermo)),
      position_structure_(ssti_mono->GetBlockPositions(SSTI::Subproblem::structure)->at(0))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSTI::AssembleStrategyBlockBlock::AssembleStrategyBlockBlock(
    Teuchos::RCP<const SSTI::SSTIMono> ssti_mono)
    : AssembleStrategyBlock(std::move(ssti_mono))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSTI::AssembleStrategyBlockSparse::AssembleStrategyBlockSparse(
    Teuchos::RCP<const SSTI::SSTIMono> ssti_mono)
    : AssembleStrategyBlock(std::move(ssti_mono))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSTI::AssembleStrategySparse::AssembleStrategySparse(Teuchos::RCP<const SSTI::SSTIMono> ssti_mono)
    : AssembleStrategyBase(std::move(ssti_mono))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatradomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatradomain_block = LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatradomain);

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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatradomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatradomain_sparse = LINALG::CastToConstSparseMatrixAndCheckSuccess(scatradomain);

  auto& systemmatrix_block_scatra_scatra =
      systemmatrix_block->Matrix(BlockPositionScaTra()->at(0), BlockPositionScaTra()->at(0));

  systemmatrix_block_scatra_scatra.Add(*scatradomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleScatra(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatradomain)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatradomain_sparse = LINALG::CastToConstSparseMatrixAndCheckSuccess(scatradomain);

  // add scalar transport system matrix to global system matrix
  systemmatrix_sparse->Add(*scatradomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseMatrix> structuredomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  // add entire block or assemble slave side to master side
  if (InterfaceMeshtying())
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
void SSTI::AssembleStrategyBlockSparse::AssembleStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseMatrix> structuredomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  // add entire block or assemble slave side to master side
  if (InterfaceMeshtying())
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
void SSTI::AssembleStrategySparse::AssembleStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseMatrix> structuredomain)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);

  // add entire block or assemble slave side to master side
  if (InterfaceMeshtying())
    AssembleStructureMeshtying(*systemmatrix_sparse, structuredomain);
  else
    systemmatrix_sparse->Add(*structuredomain, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBase::AssembleStructureMeshtying(
    LINALG::SparseMatrix& systemmatrix_structure,
    Teuchos::RCP<const LINALG::SparseMatrix> structuredomain)
{
  /* Transform and assemble the structure matrix into the global system matrix block by block:
   * S_i:  structure interior dofs
   * S_m:  structure master side dofs
   * S_ss: structure slave surface dofs
   * S_sl: structure slave line dofs
   *
   *      | S_i | S_m | S_ss| S_sl|
   *      |-----|-----|-----|-----|
   * S_i  |  a  |  b  |  c  |  d  |
   * S_m  |  e  |  f  |  g  |  -  |
   * S_ss |  h  |  i  |  j  |  k  |
   * S_sl |  l  |  -  |  m  |  n  |
   *      -------------------------
   */

  auto map_structure_interior = SSTIStructureMeshtying()->InteriorMap();
  auto master_dof_map = SSTIStructureMeshtying()->FullMasterSideMap();

  // assemble derivs. of interior dofs w.r.t. interior dofs (block a)
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *map_structure_interior,
      *map_structure_interior, 1.0, nullptr, nullptr, systemmatrix_structure, true, true);

  // assemble derivs. of interior dofs w.r.t. master dofs (block b)
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *map_structure_interior,
      *master_dof_map, 1.0, nullptr, nullptr, systemmatrix_structure, true, true);

  // assemble derivs. of master dofs w.r.t. interior dofs (block e)
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *master_dof_map,
      *map_structure_interior, 1.0, nullptr, nullptr, systemmatrix_structure, true, true);

  // assemble derivs. of master dofs w.r.t. master dofs (block f)
  LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *master_dof_map, *master_dof_map, 1.0,
      nullptr, nullptr, systemmatrix_structure, true, true);

  for (const auto& meshtying : SSTIStructureMeshtying()->MeshtyingHandlers())
  {
    auto cond_slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
    auto converter = meshtying->SlaveSideConverter();

    // assemble derivs. of surface slave dofs w.r.t. interior dofs (block h)
    LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *cond_slave_dof_map,
        *map_structure_interior, 1.0, &(*converter), nullptr, systemmatrix_structure, true, true);

    // assemble derivs. of surface slave dofs w.r.t. master dofs (block i)
    LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *cond_slave_dof_map, *master_dof_map,
        1.0, &(*converter), nullptr, systemmatrix_structure, true, true);

    // assemble derivs. of interior dofs w.r.t. surface slave dofs (block c)
    LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *map_structure_interior,
        *cond_slave_dof_map, 1.0, nullptr, &(*converter), systemmatrix_structure, true, true);

    // assemble derivs. of master dofs w.r.t. surface slave dofs (block g)
    LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *master_dof_map, *cond_slave_dof_map,
        1.0, nullptr, &(*converter), systemmatrix_structure, true, true);

    // assemble derivs. of surface slave dofs w.r.t. surface slave dofs (block j)
    LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *cond_slave_dof_map,
        *cond_slave_dof_map, 1.0, &(*converter), &(*converter), systemmatrix_structure, true, true);

    for (const auto& meshtying2 : SSTIStructureMeshtying()->MeshtyingHandlers())
    {
      if (meshtying2 != meshtying)
      {
        auto cond_slave_dof_map2 = meshtying2->SlaveMasterCoupling()->SlaveDofMap();
        auto converter2 = meshtying2->SlaveSideConverter();

        // assemble derivatives of surface slave dofs w.r.t. line slave dofs (block l)
        LINALG::MatrixLogicalSplitAndTransform()(*structuredomain, *cond_slave_dof_map,
            *cond_slave_dof_map2, 1.0, &(*converter), &(*converter2), systemmatrix_structure, true,
            true);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleThermo(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> thermodomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto thermodomain_block = LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(thermodomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionThermo()->size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(BlockPositionThermo()->size()); ++jblock)
    {
      auto& systemmatrix_block_ithermo_jthermo = systemmatrix_block->Matrix(
          BlockPositionThermo()->at(iblock), BlockPositionThermo()->at(jblock));

      systemmatrix_block_ithermo_jthermo.Add(
          thermodomain_block->Matrix(iblock, jblock), false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleThermo(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> thermodomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto thermodomain_sparse = LINALG::CastToConstSparseMatrixAndCheckSuccess(thermodomain);

  auto& systemmatrix_block_thermo_thermo =
      systemmatrix_block->Matrix(BlockPositionThermo()->at(0), BlockPositionThermo()->at(0));

  systemmatrix_block_thermo_thermo.Add(*thermodomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleThermo(Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> thermodomain)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto thermodomain_sparse = LINALG::CastToConstSparseMatrixAndCheckSuccess(thermodomain);

  // add scalar transport system matrix to global system matrix
  systemmatrix_sparse->Add(*thermodomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleScatraStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatrastructuredomain,
    Teuchos::RCP<const LINALG::SparseOperator> scatrastructureinterface)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatrastructuredomain_block =
      LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatrastructuredomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
  {
    const auto scatrastructuredomain_subblock = scatrastructuredomain_block->Matrix(iblock, 0);

    // add entire block or assemble slave side to master side
    if (InterfaceMeshtying())
    {
      AssembleXXXStructureMeshtying(
          systemmatrix_block->Matrix(BlockPositionScaTra()->at(iblock), PositionStructure()),
          scatrastructuredomain_block->Matrix(iblock, 0));

      auto scatrastructureinterface_block =
          LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatrastructureinterface);

      AssembleXXXStructureMeshtying(
          systemmatrix_block->Matrix(BlockPositionScaTra()->at(iblock), PositionStructure()),
          scatrastructureinterface_block->Matrix(iblock, 0));
    }
    else
    {
      auto& systemmatrix_block_iscatra_struct =
          systemmatrix_block->Matrix(BlockPositionScaTra()->at(iblock), PositionStructure());

      systemmatrix_block_iscatra_struct.Add(scatrastructuredomain_subblock, false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleScatraStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatrastructuredomain,
    Teuchos::RCP<const LINALG::SparseOperator> scatrastructureinterface)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  auto scatrastructuredomain_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(scatrastructuredomain);

  // add entire block or assemble slave side to master side
  if (InterfaceMeshtying())
  {
    AssembleXXXStructureMeshtying(
        systemmatrix_block->Matrix(BlockPositionScaTra()->at(0), PositionStructure()),
        *scatrastructuredomain_sparse);

    auto scatrastructureinterface_sparse =
        LINALG::CastToConstSparseMatrixAndCheckSuccess(scatrastructureinterface);

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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleScatraStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatrastructuredomain,
    Teuchos::RCP<const LINALG::SparseOperator> scatrastructureinterface)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatrastructuredomain_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(scatrastructuredomain);

  // add entire block or assemble slave side to master side
  if (InterfaceMeshtying())
  {
    AssembleXXXStructureMeshtying(*systemmatrix_sparse, *scatrastructuredomain_sparse);

    auto scatrastructureinterface_sparse =
        LINALG::CastToConstSparseMatrixAndCheckSuccess(scatrastructureinterface);

    AssembleXXXStructureMeshtying(*systemmatrix_sparse, *scatrastructureinterface_sparse);
  }

  else
    systemmatrix_sparse->Add(*scatrastructuredomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBase::AssembleXXXStructureMeshtying(
    LINALG::SparseMatrix& systemmatrix_x_structure, const LINALG::SparseMatrix& x_structurematrix)
{
  auto map_structure_interior = SSTIStructureMeshtying()->InteriorMap();
  auto master_dof_map = SSTIStructureMeshtying()->FullMasterSideMap();

  // assemble derivs. of x w.r.t. structural interior dofs
  LINALG::MatrixLogicalSplitAndTransform()(x_structurematrix, x_structurematrix.RangeMap(),
      *map_structure_interior, 1.0, nullptr, nullptr, systemmatrix_x_structure, true, true);

  // assemble derivs. of x w.r.t. structural master dofs
  LINALG::MatrixLogicalSplitAndTransform()(x_structurematrix, x_structurematrix.RangeMap(),
      *master_dof_map, 1.0, nullptr, nullptr, systemmatrix_x_structure, true, true);

  for (const auto& meshtying : SSTIStructureMeshtying()->MeshtyingHandlers())
  {
    auto cond_slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
    auto converter = meshtying->SlaveSideConverter();

    // assemble derivs. of x w.r.t. structural surface slave dofs
    LINALG::MatrixLogicalSplitAndTransform()(x_structurematrix, x_structurematrix.RangeMap(),
        *cond_slave_dof_map(), 1.0, nullptr, &(*converter), systemmatrix_x_structure, true, true);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleScatraThermoDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrathermodomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatrathermodomain_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(scatrathermodomain);

  // assemble blocks of scalar transport-thermo matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(BlockPositionThermo()->size()); ++jblock)
    {
      auto scatrathermodomain_subblock = scatrathermodomain_block->Matrix(iblock, jblock);
      auto systemmatrix_subblock = systemmatrix_block->Matrix(
          BlockPositionScaTra()->at(iblock), BlockPositionThermo()->at(jblock));
      systemmatrix_subblock.UnComplete();
      systemmatrix_subblock.Add(scatrathermodomain_subblock, false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleScatraThermoDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrathermodomain)
{
  auto systemmatrix_subblock =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix)
          ->Matrix(BlockPositionScaTra()->at(0), BlockPositionThermo()->at(0));
  auto scatrathermodomain_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(scatrathermodomain);

  systemmatrix_subblock.UnComplete();
  systemmatrix_subblock.Add(*scatrathermodomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleScatraThermoDomain(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrathermodomain)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatrathermodomain_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(scatrathermodomain);

  // add scalar transport-thermo matrix into global system matrix
  systemmatrix_sparse->Add(*scatrathermodomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleScatraThermoInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatrathermointerface)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatrathermointerface_block =
      LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatrathermointerface);

  LINALG::SparseMatrix masterderiv(*AllMaps()->BlockMapScatra()->FullMap(), 27, false, true);

  for (int i = 0; i < static_cast<int>(BlockPositionScaTra()->size()); ++i)
  {
    for (int j = 0; j < static_cast<int>(BlockPositionThermo()->size()); ++j)
    {
      const auto scatrathermointerface_subblock = scatrathermointerface_block->Matrix(i, j);

      // assemble linearizations of slave- and master side scatra fluxes w.r.t. slave temperatures
      // into system matrix
      auto& systemmatrix_block_iscatra_jthermo =
          systemmatrix_block->Matrix(BlockPositionScaTra()->at(i), BlockPositionThermo()->at(j));
      systemmatrix_block_iscatra_jthermo.Add(scatrathermointerface_subblock, false, 1.0, 1.0);

      // assemble linearizations of slave- and master side scatra fluxes w.r.t. master temperatures
      // into system matrix
      ADAPTER::CouplingSlaveConverter thermo_converter(*MeshtyingThermo()->CouplingAdapter());

      LINALG::MatrixLogicalSplitAndTransform()(scatrathermointerface_subblock,
          scatrathermointerface_subblock.RangeMap(),
          *MeshtyingThermo()->CouplingAdapter()->MasterDofMap(), 1.0, nullptr, &thermo_converter,
          masterderiv, true, true);
    }
  }

  masterderiv.Complete(*MeshtyingThermo()->CouplingAdapter()->MasterDofMap(),
      *AllMaps()->BlockMapScatra()->FullMap());

  const auto blockmasterderiv = masterderiv.Split<LINALG::DefaultBlockMatrixStrategy>(
      *AllMaps()->BlockMapThermo(), *AllMaps()->BlockMapScatra());

  blockmasterderiv->Complete();

  for (int i = 0; i < static_cast<int>(BlockPositionScaTra()->size()); ++i)
  {
    for (int j = 0; j < static_cast<int>(BlockPositionThermo()->size()); ++j)
    {
      // assemble linearizations of slave side scatra fluxes w.r.t. slave and master side elch
      // into system matrix
      auto& systemmatrix_block_iscatra_jthermo =
          systemmatrix_block->Matrix(BlockPositionScaTra()->at(i), BlockPositionThermo()->at(j));
      systemmatrix_block_iscatra_jthermo.Add(blockmasterderiv->Matrix(i, j), false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleScatraThermoInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatrathermointerface)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatrathermointerface_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(scatrathermointerface);

  // assemble linearizations of slave- and master side scatra fluxes w.r.t. slave temperatures into
  // system matrix
  auto& systemmatrix_block_scatra_thermo =
      systemmatrix_block->Matrix(BlockPositionScaTra()->at(0), BlockPositionThermo()->at(0));
  systemmatrix_block_scatra_thermo.Add(*scatrathermointerface_sparse, false, 1.0, 1.0);

  // assemble linearizations of slave- and master side scatra fluxes w.r.t. master temperatures into
  // system matrix
  ADAPTER::CouplingSlaveConverter thermo_converter(*MeshtyingThermo()->CouplingAdapter());

  LINALG::MatrixLogicalSplitAndTransform()(*scatrathermointerface_sparse,
      scatrathermointerface_sparse->RangeMap(),
      *MeshtyingThermo()->CouplingAdapter()->MasterDofMap(), 1.0, nullptr, &thermo_converter,
      systemmatrix_block->Matrix(BlockPositionScaTra()->at(0), BlockPositionThermo()->at(0)), true,
      true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleScatraThermoInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatrathermointerface)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatrathermointerface_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(scatrathermointerface);

  // assemble linearizations of slave- and master side scatra fluxes w.r.t. slave temperatures into
  // system matrix
  systemmatrix_sparse->Add(*scatrathermointerface_sparse, false, 1.0, 1.0);

  // assemble linearizations of slave- and master side scatra fluxes w.r.t. master temperatures into
  // system matrix
  ADAPTER::CouplingSlaveConverter thermo_converter(*MeshtyingThermo()->CouplingAdapter());

  LINALG::MatrixLogicalSplitAndTransform()(*scatrathermointerface_sparse,
      scatrathermointerface_sparse->RangeMap(),
      *MeshtyingThermo()->CouplingAdapter()->MasterDofMap(), 1.0, nullptr, &thermo_converter,
      *systemmatrix_sparse, true, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleStructureScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> structurescatradomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto structurescatradomain_block =
      LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(structurescatradomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
  {
    const auto structurescatradomain_subblock = structurescatradomain_block->Matrix(0, iblock);

    // add entire block or assemble slave side to master side
    if (InterfaceMeshtying())
    {
      AssembleStructureXXXMeshtying(
          systemmatrix_block->Matrix(PositionStructure(), BlockPositionScaTra()->at(iblock)),
          structurescatradomain_subblock);
    }
    else
    {
      auto& systemmatrix_block_struct_iscatra =
          systemmatrix_block->Matrix(PositionStructure(), BlockPositionScaTra()->at(iblock));

      systemmatrix_block_struct_iscatra.Add(structurescatradomain_subblock, false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleStructureScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> structurescatradomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto structurescatradomain_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(structurescatradomain);

  // add entire block or assemble slave side to master side
  if (InterfaceMeshtying())
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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleStructureScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> structurescatradomain)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto structurescatradomain_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(structurescatradomain);

  // add entire block or assemble slave side to master side
  if (InterfaceMeshtying())
    AssembleStructureXXXMeshtying(*systemmatrix_sparse, *structurescatradomain_sparse);
  else
    systemmatrix_sparse->Add(*structurescatradomain_sparse, false, 1.0, 1.0);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBase::AssembleStructureXXXMeshtying(
    LINALG::SparseMatrix& systemmatrix_structure_x, const LINALG::SparseMatrix& structures_x_matrix)
{
  auto map_structure_interior = SSTIStructureMeshtying()->InteriorMap();
  auto master_dof_map = SSTIStructureMeshtying()->FullMasterSideMap();

  // assemble derivs. of structural interior dofs w.r.t. scatra dofs
  LINALG::MatrixLogicalSplitAndTransform()(structures_x_matrix, *map_structure_interior,
      structures_x_matrix.DomainMap(), 1.0, nullptr, nullptr, systemmatrix_structure_x, true, true);

  // assemble derivs. of structural master dofs w.r.t. scatra dofs
  LINALG::MatrixLogicalSplitAndTransform()(structures_x_matrix, *master_dof_map,
      structures_x_matrix.DomainMap(), 1.0, nullptr, nullptr, systemmatrix_structure_x, true, true);

  for (const auto& meshtying : SSTIStructureMeshtying()->MeshtyingHandlers())
  {
    auto cond_slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
    auto converter = meshtying->SlaveSideConverter();

    // assemble derivs. of structural surface slave dofs w.r.t. scatra dofs
    LINALG::MatrixLogicalSplitAndTransform()(structures_x_matrix, *cond_slave_dof_map(),
        structures_x_matrix.DomainMap(), 1.0, &(*converter), nullptr, systemmatrix_structure_x,
        true, true);
  }
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleThermoScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> thermoscatradomain,
    Teuchos::RCP<const LINALG::SparseOperator> thermoscatrainterface)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto thermoscatradomain_block =
      LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(thermoscatradomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionThermo()->size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(BlockPositionScaTra()->size()); ++jblock)
    {
      auto systemmatrix_block_ithermo_jscatra = systemmatrix_block->Matrix(
          BlockPositionThermo()->at(iblock), BlockPositionScaTra()->at(jblock));
      systemmatrix_block_ithermo_jscatra.UnComplete();
      systemmatrix_block_ithermo_jscatra.Add(
          thermoscatradomain_block->Matrix(iblock, jblock), false, 1.0, 1.0);
    }
  }

  if (InterfaceMeshtying()) AssembleThermoScatraInterface(systemmatrix, thermoscatrainterface);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleThermoScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> thermoscatradomain,
    Teuchos::RCP<const LINALG::SparseOperator> thermoscatrainterface)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto thermoscatradomain_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(thermoscatradomain);

  auto& systemmatrix_block_thermo_scatra =
      systemmatrix_block->Matrix(BlockPositionThermo()->at(0), BlockPositionScaTra()->at(0));
  systemmatrix_block_thermo_scatra.UnComplete();
  systemmatrix_block_thermo_scatra.Add(*thermoscatradomain_sparse, false, 1.0, 1.0);

  if (InterfaceMeshtying()) AssembleThermoScatraInterface(systemmatrix, thermoscatrainterface);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleThermoScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> thermoscatradomain,
    Teuchos::RCP<const LINALG::SparseOperator> thermoscatrainterface)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto thermoscatradomain_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(thermoscatradomain);

  // add scalar transport system matrix to global system matrix
  systemmatrix_sparse->Add(*thermoscatradomain_sparse, false, 1.0, 1.0);

  if (InterfaceMeshtying()) AssembleThermoScatraInterface(systemmatrix, thermoscatrainterface);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleThermoScatraInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> thermoscatrainterface)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  auto thermoscatrainterface_block =
      LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(thermoscatrainterface);

  LINALG::SparseMatrix masterflux(*AllMaps()->BlockMapThermo()->FullMap(), 27, false, true);

  for (int i = 0; i < static_cast<int>(BlockPositionThermo()->size()); ++i)
  {
    for (int j = 0; j < static_cast<int>(BlockPositionScaTra()->size()); ++j)
    {
      const auto thermoscatrainterface_subblock = thermoscatrainterface_block->Matrix(i, j);

      // assemble linearizations of slave side scatra fluxes w.r.t. slave and master side elch
      // into system matrix
      auto& systemmatrix_block_ithermo_jscatra =
          systemmatrix_block->Matrix(BlockPositionThermo()->at(i), BlockPositionScaTra()->at(j));
      systemmatrix_block_ithermo_jscatra.Add(thermoscatrainterface_subblock, false, 1.0, 1.0);

      // assemble linearizations of master side thermo fluxes w.r.t. slave and master side elch
      // into system matrix
      ADAPTER::CouplingSlaveConverter thermo_converter(*MeshtyingThermo()->CouplingAdapter());

      LINALG::SparseMatrix slaveflux(*AllMaps()->BlockMapThermo()->FullMap(), 27, false, true);

      SCATRA::MeshtyingStrategyS2I::ExtractMatrixRows(
          thermoscatrainterface_subblock, slaveflux, *MeshtyingThermo()->BlockMapsSlave().Map(i));

      slaveflux.Complete(
          *AllMaps()->BlockMapScatra()->FullMap(), *AllMaps()->BlockMapThermo()->FullMap());

      LINALG::MatrixLogicalSplitAndTransform()(thermoscatrainterface_subblock,
          *MeshtyingThermo()->CouplingAdapter()->MasterDofMap(),
          thermoscatrainterface_subblock.DomainMap(), -1.0, &thermo_converter, nullptr, masterflux,
          true, true);
    }
  }

  masterflux.Complete(
      *AllMaps()->BlockMapScatra()->FullMap(), *AllMaps()->BlockMapThermo()->FullMap());

  const auto blockmasterflux = masterflux.Split<LINALG::DefaultBlockMatrixStrategy>(
      *AllMaps()->BlockMapScatra(), *AllMaps()->BlockMapThermo());

  blockmasterflux->Complete();

  // assemble linearizations of slave side thermo fluxes w.r.t. slave and master side elch
  // into system matrix
  for (int i = 0; i < static_cast<int>(BlockPositionThermo()->size()); ++i)
  {
    for (int j = 0; j < static_cast<int>(BlockPositionScaTra()->size()); ++j)
    {
      auto& systemmatrix_block_ithermo_jscatra =
          systemmatrix_block->Matrix(BlockPositionThermo()->at(i), BlockPositionScaTra()->at(j));

      systemmatrix_block_ithermo_jscatra.Add(blockmasterflux->Matrix(i, j), false, 1.0, 1.0);
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleThermoScatraInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> thermoscatrainterface)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto thermoscatrainterface_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(thermoscatrainterface);

  // assemble linearizations of slave side scatra fluxes w.r.t. slave and master side elch
  // into system matrix
  auto& systemmatrix_block_thermo_scatra =
      systemmatrix_block->Matrix(BlockPositionThermo()->at(0), BlockPositionScaTra()->at(0));
  systemmatrix_block_thermo_scatra.Add(*thermoscatrainterface_sparse, false, 1.0, 1.0);

  // assemble linearizations of master side thermo fluxes w.r.t. slave and master side elch
  // into system matrix
  ADAPTER::CouplingSlaveConverter thermo_converter(*MeshtyingThermo()->CouplingAdapter());

  LINALG::MatrixLogicalSplitAndTransform()(*thermoscatrainterface_sparse,
      *MeshtyingThermo()->CouplingAdapter()->MasterDofMap(),
      thermoscatrainterface_sparse->DomainMap(), -1.0, &thermo_converter, nullptr,
      systemmatrix_block->Matrix(BlockPositionThermo()->at(0), BlockPositionScaTra()->at(0)), true,
      true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleThermoScatraInterface(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> thermoscatrainterface)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto thermoscatrainterface_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(thermoscatrainterface);

  // assemble linearizations of slave side scatra fluxes w.r.t. slave and master side elch
  // into system matrix
  systemmatrix_sparse->Add(*thermoscatrainterface_sparse, false, 1.0, 1.0);

  // assemble linearizations of master side scatra fluxes w.r.t. slave and master side elch
  // into system matrix
  ADAPTER::CouplingSlaveConverter thermo_converter(*MeshtyingThermo()->CouplingAdapter());

  LINALG::MatrixLogicalSplitAndTransform()(*thermoscatrainterface_sparse,
      *MeshtyingThermo()->CouplingAdapter()->MasterDofMap(),
      thermoscatrainterface_sparse->DomainMap(), -1.0, &thermo_converter, nullptr,
      *systemmatrix_sparse, true, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleThermoStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> thermostructuredomain,
    Teuchos::RCP<const LINALG::SparseOperator> thermostructureinterface)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto thermostructuredomain_block =
      LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(thermostructuredomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionThermo()->size()); ++iblock)
  {
    const auto thermostructuredomain_subblock = thermostructuredomain_block->Matrix(iblock, 0);

    // add entire block or assemble slave side to master side
    if (InterfaceMeshtying())
    {
      AssembleXXXStructureMeshtying(
          systemmatrix_block->Matrix(BlockPositionThermo()->at(iblock), PositionStructure()),
          thermostructuredomain_subblock);

      auto thermostructureinterface_block =
          LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(thermostructureinterface);

      AssembleXXXStructureMeshtying(
          systemmatrix_block->Matrix(BlockPositionThermo()->at(iblock), PositionStructure()),
          thermostructureinterface_block->Matrix(iblock, 0));
    }
    else
    {
      auto& systemmatrix_block_ithermo_struct =
          systemmatrix_block->Matrix(BlockPositionThermo()->at(iblock), PositionStructure());

      systemmatrix_block_ithermo_struct.Add(thermostructuredomain_subblock, false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleThermoStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> thermostructuredomain,
    Teuchos::RCP<const LINALG::SparseOperator> thermostructureinterface)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto thermostructuredomain_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(thermostructuredomain);

  // add entire block or assemble slave side to master side
  if (InterfaceMeshtying())
  {
    AssembleXXXStructureMeshtying(
        systemmatrix_block->Matrix(BlockPositionThermo()->at(0), PositionStructure()),
        *thermostructuredomain_sparse);

    auto thermostructureinterface_sparse =
        LINALG::CastToConstSparseMatrixAndCheckSuccess(thermostructureinterface);

    AssembleXXXStructureMeshtying(
        systemmatrix_block->Matrix(BlockPositionThermo()->at(0), PositionStructure()),
        *thermostructureinterface_sparse);
  }
  else
  {
    auto& sytemmatrix_block_thermo_structure =
        systemmatrix_block->Matrix(BlockPositionThermo()->at(0), PositionStructure());

    sytemmatrix_block_thermo_structure.Add(*thermostructuredomain_sparse, false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleThermoStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> thermostructuredomain,
    Teuchos::RCP<const LINALG::SparseOperator> thermostructureinterface)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto thermostructuredomain_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(thermostructuredomain);

  if (InterfaceMeshtying())
  {
    AssembleXXXStructureMeshtying(*systemmatrix_sparse, *thermostructuredomain_sparse);

    auto thermostructureinterface_sparse =
        LINALG::CastToConstSparseMatrixAndCheckSuccess(thermostructureinterface);

    AssembleXXXStructureMeshtying(*systemmatrix_sparse, *thermostructureinterface_sparse);
  }
  else
    systemmatrix_sparse->Add(*thermostructuredomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleStructureThermo(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> structurethermodomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto structurethermodomain_block =
      LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(structurethermodomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionThermo()->size()); ++iblock)
  {
    const auto structurethermodomain_subblock = structurethermodomain_block->Matrix(0, iblock);

    // add entire block or assemble slave side to master side
    if (InterfaceMeshtying())
    {
      AssembleStructureXXXMeshtying(
          systemmatrix_block->Matrix(PositionStructure(), BlockPositionThermo()->at(iblock)),
          structurethermodomain_subblock);
    }
    else
    {
      auto& systemmatrix_block_struct_ithermo =
          systemmatrix_block->Matrix(PositionStructure(), BlockPositionThermo()->at(iblock));

      systemmatrix_block_struct_ithermo.Add(structurethermodomain_subblock, false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleStructureThermo(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> structurethermodomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto structurethermodomain_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(structurethermodomain);

  // add entire block or assemble slave side to master side
  if (InterfaceMeshtying())
  {
    AssembleStructureXXXMeshtying(
        systemmatrix_block->Matrix(PositionStructure(), BlockPositionThermo()->at(0)),
        *structurethermodomain_sparse);
  }
  else
  {
    auto& systemmatrix_block_struct_thermo =
        systemmatrix_block->Matrix(PositionStructure(), BlockPositionThermo()->at(0));

    systemmatrix_block_struct_thermo.Add(*structurethermodomain_sparse, false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleStructureThermo(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const LINALG::SparseOperator> structurethermodomain)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto structurethermodomain_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(structurethermodomain);

  // add entire block or assemble slave side to master side
  if (InterfaceMeshtying())
    AssembleStructureXXXMeshtying(*systemmatrix_sparse, *structurethermodomain_sparse);
  else
    systemmatrix_sparse->Add(*structurethermodomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::ApplyMeshtyingSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  if (InterfaceMeshtying())
  {
    auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

    ApplyMeshtyingSysMat(systemmatrix_block->Matrix(PositionStructure(), PositionStructure()));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::ApplyMeshtyingSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  if (InterfaceMeshtying())
  {
    // cast systemmatrix
    auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

    ApplyMeshtyingSysMat(systemmatrix_block->Matrix(PositionStructure(), PositionStructure()));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::ApplyMeshtyingSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  if (InterfaceMeshtying())
  {
    // cast systemmatrix
    auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);

    ApplyMeshtyingSysMat(*systemmatrix_sparse);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBase::ApplyMeshtyingSysMat(LINALG::SparseMatrix& systemmatrix_structure)
{
  // map for slave side structural degrees of freedom
  auto slavemaps = SSTIStructureMeshtying()->FullSlaveSideMap();

  // subject slave-side rows of structural system matrix to pseudo Dirichlet conditions to
  // finalize structural meshtying
  const double one(1.0);
  for (int doflid_slave = 0; doflid_slave < slavemaps->NumMyElements(); ++doflid_slave)
  {
    // extract global ID of current slave-side row
    const int dofgid_slave = slavemaps->GID(doflid_slave);
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
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlock::ApplyStructuralDBCSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  // locsys manager of strucutre
  const auto& locsysmanager_structure = StructureField()->LocsysManager();

  // map of strucutral Dirichlet BCs
  const auto dbcmap_structure = StructureField()->GetDBCMapExtractor()->CondMap();

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
              locsysmanager_structure->Trafo(), *dbcmap_structure, iblock == PositionStructure());
      locsysmanager_structure->RotateLocalToGlobal(
          Teuchos::rcp(&systemmatrix_block->Matrix(PositionStructure(), iblock), false));
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::ApplyStructuralDBCSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  // locsys manager of strucutre
  const auto& locsysmanager_structure = StructureField()->LocsysManager();

  // map of strucutral Dirichlet BCs
  const auto& dbcmap_structure = StructureField()->GetDBCMapExtractor()->CondMap();

  // structural dof row map
  const auto& dofrowmap_structure = StructureField()->DofRowMap();

  if (locsysmanager_structure == Teuchos::null)
    systemmatrix->ApplyDirichlet(*dbcmap_structure);
  else
  {
    auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);

    // extract structural rows of global system matrix
    auto systemmatrix_structure =
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
void SSTI::AssembleStrategyBase::AssembleRHS(Teuchos::RCP<Epetra_Vector> RHS,
    Teuchos::RCP<const Epetra_Vector> RHSscatra, Teuchos::RCP<const Epetra_Vector> RHSstructure,
    Teuchos::RCP<const Epetra_Vector> RHSthermo)
{
  // zero out RHS
  RHS->PutScalar(0.0);

  // assemble scalar transport right-hand side vector into monolithic right-hand side vector
  AllMaps()->MapsSubProblems()->InsertVector(
      RHSscatra, ssti_mono_->GetProblemPosition(Subproblem::scalar_transport), RHS);
  AllMaps()->MapsSubProblems()->InsertVector(
      RHSthermo, ssti_mono_->GetProblemPosition(Subproblem::thermo), RHS);

  if (InterfaceMeshtying())
  {
    // perform structural meshtying before assembling structural right-hand side vector into
    // monolithic right-hand side vector

    // make copy of structural right-hand side vector
    Epetra_Vector residual_structure(*RHSstructure);

    auto rhs_structure_master = LINALG::CreateVector(*StructureField()->DofRowMap(), true);

    for (const auto& meshtying : SSTIStructureMeshtying()->MeshtyingHandlers())
    {
      auto coupling_adapter = meshtying->SlaveMasterCoupling();
      auto coupling_map_extractor = meshtying->SlaveMasterExtractor();

      // transform slave-side part of structure right-hand side vector to master side
      const auto rhs_structure_only_slave_dofs =
          coupling_map_extractor->ExtractVector(residual_structure, 1);

      const auto rhs_structure_only_master_dofs =
          coupling_adapter->SlaveToMaster(rhs_structure_only_slave_dofs);

      coupling_map_extractor->AddVector(rhs_structure_only_master_dofs, 2, rhs_structure_master);

      // zero out slave-side part of structure right-hand side vector
      coupling_map_extractor->PutScalar(residual_structure, 1, 0.0);
    }

    // locsys manager of strucutre
    const auto& locsysmanager_structure = ssti_mono_->StructureField()->LocsysManager();

    // apply pseudo Dirichlet conditions to transformed slave-side part of structural right-hand
    // side vector
    auto zeros = Teuchos::rcp(new Epetra_Vector(residual_structure.Map()));

    if (locsysmanager_structure != Teuchos::null)
      locsysmanager_structure->RotateGlobalToLocal(rhs_structure_master);
    LINALG::ApplyDirichlettoSystem(rhs_structure_master, zeros,
        *ssti_mono_->StructureField()->GetDBCMapExtractor()->CondMap());
    if (locsysmanager_structure != Teuchos::null)
      locsysmanager_structure->RotateLocalToGlobal(rhs_structure_master);

    // assemble master-side part of structure right-hand side vector
    residual_structure.Update(1.0, *rhs_structure_master, 1.0);

    // assemble final structural right-hand side vector into monolithic right-hand side vector
    AllMaps()->MapsSubProblems()->AddVector(
        residual_structure, ssti_mono_->GetProblemPosition(Subproblem::structure), *RHS, -1.0);
  }
  else
  {
    AllMaps()->MapsSubProblems()->AddVector(
        RHSstructure, ssti_mono_->GetProblemPosition(Subproblem::structure), RHS, -1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<SSTI::AssembleStrategyBase> SSTI::BuildAssembleStrategy(
    Teuchos::RCP<const SSTI::SSTIMono> ssti_mono, LINALG::MatrixType matrixtype_ssti,
    LINALG::MatrixType matrixtype_scatra)
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
          assemblestrategy = Teuchos::rcp(new SSTI::AssembleStrategyBlockBlock(ssti_mono));
          break;
        }
        case LINALG::MatrixType::sparse:
        {
          assemblestrategy = Teuchos::rcp(new SSTI::AssembleStrategyBlockSparse(ssti_mono));
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
      assemblestrategy = Teuchos::rcp(new SSTI::AssembleStrategySparse(ssti_mono));
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
