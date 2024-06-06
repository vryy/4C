/*----------------------------------------------------------------------*/
/*! \file
\brief Assemble strategy for monolithic SSTI

\level 2

*----------------------------------------------------------------------*/
#include "4C_ssti_monolithic_assemble_strategy.hpp"

#include "4C_adapter_str_ssiwrapper.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_discretization_condition_locsys.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_scatra_timint_meshtying_strategy_s2i.hpp"
#include "4C_ssi_utils.hpp"
#include "4C_ssti_monolithic.hpp"
#include "4C_ssti_utils.hpp"

FOUR_C_NAMESPACE_OPEN

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
      position_structure_(ssti_mono->GetBlockPositions(SSTI::Subproblem::structure).at(0))
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
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatradomain)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatradomain_block =
      Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatradomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(block_position_sca_tra().size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(block_position_sca_tra().size()); ++jblock)
    {
      auto& systemmatrix_block_iscatra_jscatra = systemmatrix_block->Matrix(
          block_position_sca_tra().at(iblock), block_position_sca_tra().at(jblock));

      systemmatrix_block_iscatra_jscatra.Add(
          scatradomain_block->Matrix(iblock, jblock), false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleScatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatradomain)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatradomain_sparse = Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatradomain);

  auto& systemmatrix_block_scatra_scatra =
      systemmatrix_block->Matrix(block_position_sca_tra().at(0), block_position_sca_tra().at(0));

  systemmatrix_block_scatra_scatra.Add(*scatradomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleScatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatradomain)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatradomain_sparse = Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatradomain);

  // add scalar transport system matrix to global system matrix
  systemmatrix_sparse->Add(*scatradomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleStructure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> structuredomain)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  // add entire block or assemble slave side to master side
  if (interface_meshtying())
  {
    assemble_structure_meshtying(
        systemmatrix_block->Matrix(position_structure(), position_structure()), structuredomain);
  }
  else
  {
    auto& systemmatrix_block_struct_struct =
        systemmatrix_block->Matrix(position_structure(), position_structure());

    systemmatrix_block_struct_struct.Add(*structuredomain, false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleStructure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> structuredomain)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  // add entire block or assemble slave side to master side
  if (interface_meshtying())
  {
    assemble_structure_meshtying(
        systemmatrix_block->Matrix(position_structure(), position_structure()), structuredomain);
  }
  else
  {
    auto& systemmatrix_block_struct_struct =
        systemmatrix_block->Matrix(position_structure(), position_structure());

    systemmatrix_block_struct_struct.Add(*structuredomain, false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleStructure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> structuredomain)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);

  // add entire block or assemble slave side to master side
  if (interface_meshtying())
    assemble_structure_meshtying(*systemmatrix_sparse, structuredomain);
  else
    systemmatrix_sparse->Add(*structuredomain, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBase::assemble_structure_meshtying(
    Core::LinAlg::SparseMatrix& systemmatrix_structure,
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> structuredomain)
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

  auto map_structure_interior = ssti_structure_meshtying()->InteriorMap();
  auto master_dof_map = ssti_structure_meshtying()->FullMasterSideMap();

  // assemble derivs. of interior dofs w.r.t. interior dofs (block a)
  Core::LinAlg::MatrixLogicalSplitAndTransform()(*structuredomain, *map_structure_interior,
      *map_structure_interior, 1.0, nullptr, nullptr, systemmatrix_structure, true, true);

  // assemble derivs. of interior dofs w.r.t. master dofs (block b)
  Core::LinAlg::MatrixLogicalSplitAndTransform()(*structuredomain, *map_structure_interior,
      *master_dof_map, 1.0, nullptr, nullptr, systemmatrix_structure, true, true);

  // assemble derivs. of master dofs w.r.t. interior dofs (block e)
  Core::LinAlg::MatrixLogicalSplitAndTransform()(*structuredomain, *master_dof_map,
      *map_structure_interior, 1.0, nullptr, nullptr, systemmatrix_structure, true, true);

  // assemble derivs. of master dofs w.r.t. master dofs (block f)
  Core::LinAlg::MatrixLogicalSplitAndTransform()(*structuredomain, *master_dof_map, *master_dof_map,
      1.0, nullptr, nullptr, systemmatrix_structure, true, true);

  for (const auto& meshtying : ssti_structure_meshtying()->MeshTyingHandlers())
  {
    auto cond_slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
    auto converter = meshtying->SlaveSideConverter();

    // assemble derivs. of surface slave dofs w.r.t. interior dofs (block h)
    Core::LinAlg::MatrixLogicalSplitAndTransform()(*structuredomain, *cond_slave_dof_map,
        *map_structure_interior, 1.0, &(*converter), nullptr, systemmatrix_structure, true, true);

    // assemble derivs. of surface slave dofs w.r.t. master dofs (block i)
    Core::LinAlg::MatrixLogicalSplitAndTransform()(*structuredomain, *cond_slave_dof_map,
        *master_dof_map, 1.0, &(*converter), nullptr, systemmatrix_structure, true, true);

    // assemble derivs. of interior dofs w.r.t. surface slave dofs (block c)
    Core::LinAlg::MatrixLogicalSplitAndTransform()(*structuredomain, *map_structure_interior,
        *cond_slave_dof_map, 1.0, nullptr, &(*converter), systemmatrix_structure, true, true);

    // assemble derivs. of master dofs w.r.t. surface slave dofs (block g)
    Core::LinAlg::MatrixLogicalSplitAndTransform()(*structuredomain, *master_dof_map,
        *cond_slave_dof_map, 1.0, nullptr, &(*converter), systemmatrix_structure, true, true);

    // assemble derivs. of surface slave dofs w.r.t. surface slave dofs (block j)
    Core::LinAlg::MatrixLogicalSplitAndTransform()(*structuredomain, *cond_slave_dof_map,
        *cond_slave_dof_map, 1.0, &(*converter), &(*converter), systemmatrix_structure, true, true);

    for (const auto& meshtying2 : ssti_structure_meshtying()->MeshTyingHandlers())
    {
      if (meshtying2 != meshtying)
      {
        auto cond_slave_dof_map2 = meshtying2->SlaveMasterCoupling()->SlaveDofMap();
        auto converter2 = meshtying2->SlaveSideConverter();

        // assemble derivatives of surface slave dofs w.r.t. line slave dofs (block l)
        Core::LinAlg::MatrixLogicalSplitAndTransform()(*structuredomain, *cond_slave_dof_map,
            *cond_slave_dof_map2, 1.0, &(*converter), &(*converter2), systemmatrix_structure, true,
            true);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::AssembleThermo(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> thermodomain)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto thermodomain_block =
      Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(thermodomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(block_position_thermo().size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(block_position_thermo().size()); ++jblock)
    {
      auto& systemmatrix_block_ithermo_jthermo = systemmatrix_block->Matrix(
          block_position_thermo().at(iblock), block_position_thermo().at(jblock));

      systemmatrix_block_ithermo_jthermo.Add(
          thermodomain_block->Matrix(iblock, jblock), false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::AssembleThermo(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> thermodomain)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto thermodomain_sparse = Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(thermodomain);

  auto& systemmatrix_block_thermo_thermo =
      systemmatrix_block->Matrix(block_position_thermo().at(0), block_position_thermo().at(0));

  systemmatrix_block_thermo_thermo.Add(*thermodomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::AssembleThermo(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> thermodomain)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto thermodomain_sparse = Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(thermodomain);

  // add scalar transport system matrix to global system matrix
  systemmatrix_sparse->Add(*thermodomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::assemble_scatra_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrastructuredomain,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrastructureinterface)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatrastructuredomain_block =
      Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatrastructuredomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(block_position_sca_tra().size()); ++iblock)
  {
    const auto scatrastructuredomain_subblock = scatrastructuredomain_block->Matrix(iblock, 0);

    // add entire block or assemble slave side to master side
    if (interface_meshtying())
    {
      assemble_xxx_structure_meshtying(
          systemmatrix_block->Matrix(block_position_sca_tra().at(iblock), position_structure()),
          scatrastructuredomain_block->Matrix(iblock, 0));

      auto scatrastructureinterface_block =
          Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatrastructureinterface);

      assemble_xxx_structure_meshtying(
          systemmatrix_block->Matrix(block_position_sca_tra().at(iblock), position_structure()),
          scatrastructureinterface_block->Matrix(iblock, 0));
    }
    else
    {
      auto& systemmatrix_block_iscatra_struct =
          systemmatrix_block->Matrix(block_position_sca_tra().at(iblock), position_structure());

      systemmatrix_block_iscatra_struct.Add(scatrastructuredomain_subblock, false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::assemble_scatra_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrastructuredomain,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrastructureinterface)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  auto scatrastructuredomain_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatrastructuredomain);

  // add entire block or assemble slave side to master side
  if (interface_meshtying())
  {
    assemble_xxx_structure_meshtying(
        systemmatrix_block->Matrix(block_position_sca_tra().at(0), position_structure()),
        *scatrastructuredomain_sparse);

    auto scatrastructureinterface_sparse =
        Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatrastructureinterface);

    assemble_xxx_structure_meshtying(
        systemmatrix_block->Matrix(block_position_sca_tra().at(0), position_structure()),
        *scatrastructureinterface_sparse);
  }
  else
  {
    auto& systemmatrix_block_scatra_struct =
        systemmatrix_block->Matrix(block_position_sca_tra().at(0), position_structure());

    systemmatrix_block_scatra_struct.Add(*scatrastructuredomain_sparse, false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::assemble_scatra_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrastructuredomain,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrastructureinterface)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatrastructuredomain_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatrastructuredomain);

  // add entire block or assemble slave side to master side
  if (interface_meshtying())
  {
    assemble_xxx_structure_meshtying(*systemmatrix_sparse, *scatrastructuredomain_sparse);

    auto scatrastructureinterface_sparse =
        Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatrastructureinterface);

    assemble_xxx_structure_meshtying(*systemmatrix_sparse, *scatrastructureinterface_sparse);
  }

  else
    systemmatrix_sparse->Add(*scatrastructuredomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBase::assemble_xxx_structure_meshtying(
    Core::LinAlg::SparseMatrix& systemmatrix_x_structure,
    const Core::LinAlg::SparseMatrix& x_structurematrix)
{
  auto map_structure_interior = ssti_structure_meshtying()->InteriorMap();
  auto master_dof_map = ssti_structure_meshtying()->FullMasterSideMap();

  // assemble derivs. of x w.r.t. structural interior dofs
  Core::LinAlg::MatrixLogicalSplitAndTransform()(x_structurematrix, x_structurematrix.RangeMap(),
      *map_structure_interior, 1.0, nullptr, nullptr, systemmatrix_x_structure, true, true);

  // assemble derivs. of x w.r.t. structural master dofs
  Core::LinAlg::MatrixLogicalSplitAndTransform()(x_structurematrix, x_structurematrix.RangeMap(),
      *master_dof_map, 1.0, nullptr, nullptr, systemmatrix_x_structure, true, true);

  for (const auto& meshtying : ssti_structure_meshtying()->MeshTyingHandlers())
  {
    auto cond_slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
    auto converter = meshtying->SlaveSideConverter();

    // assemble derivs. of x w.r.t. structural surface slave dofs
    Core::LinAlg::MatrixLogicalSplitAndTransform()(x_structurematrix, x_structurematrix.RangeMap(),
        *cond_slave_dof_map(), 1.0, nullptr, &(*converter), systemmatrix_x_structure, true, true);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::assemble_scatra_thermo_domain(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatrathermodomain)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatrathermodomain_block =
      Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(scatrathermodomain);

  // assemble blocks of scalar transport-thermo matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(block_position_sca_tra().size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(block_position_thermo().size()); ++jblock)
    {
      auto scatrathermodomain_subblock = scatrathermodomain_block->Matrix(iblock, jblock);
      auto systemmatrix_subblock = systemmatrix_block->Matrix(
          block_position_sca_tra().at(iblock), block_position_thermo().at(jblock));
      systemmatrix_subblock.UnComplete();
      systemmatrix_subblock.Add(scatrathermodomain_subblock, false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::assemble_scatra_thermo_domain(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatrathermodomain)
{
  auto systemmatrix_subblock =
      Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix)
          ->Matrix(block_position_sca_tra().at(0), block_position_thermo().at(0));
  auto scatrathermodomain_sparse =
      Core::LinAlg::CastToSparseMatrixAndCheckSuccess(scatrathermodomain);

  systemmatrix_subblock.UnComplete();
  systemmatrix_subblock.Add(*scatrathermodomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::assemble_scatra_thermo_domain(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatrathermodomain)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatrathermodomain_sparse =
      Core::LinAlg::CastToSparseMatrixAndCheckSuccess(scatrathermodomain);

  // add scalar transport-thermo matrix into global system matrix
  systemmatrix_sparse->Add(*scatrathermodomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::assemble_scatra_thermo_interface(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrathermointerface)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatrathermointerface_block =
      Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatrathermointerface);

  Core::LinAlg::SparseMatrix masterderiv(*all_maps()->BlockMapScatra()->FullMap(), 27, false, true);

  for (int i = 0; i < static_cast<int>(block_position_sca_tra().size()); ++i)
  {
    for (int j = 0; j < static_cast<int>(block_position_thermo().size()); ++j)
    {
      const auto scatrathermointerface_subblock = scatrathermointerface_block->Matrix(i, j);

      // assemble linearizations of slave- and master side scatra fluxes w.r.t. slave temperatures
      // into system matrix
      auto& systemmatrix_block_iscatra_jthermo =
          systemmatrix_block->Matrix(block_position_sca_tra().at(i), block_position_thermo().at(j));
      systemmatrix_block_iscatra_jthermo.Add(scatrathermointerface_subblock, false, 1.0, 1.0);

      // assemble linearizations of slave- and master side scatra fluxes w.r.t. master temperatures
      // into system matrix
      Core::Adapter::CouplingSlaveConverter thermo_converter(
          *meshtying_thermo()->CouplingAdapter());

      Core::LinAlg::MatrixLogicalSplitAndTransform()(scatrathermointerface_subblock,
          scatrathermointerface_subblock.RangeMap(),
          *meshtying_thermo()->CouplingAdapter()->MasterDofMap(), 1.0, nullptr, &thermo_converter,
          masterderiv, true, true);
    }
  }

  masterderiv.Complete(*meshtying_thermo()->CouplingAdapter()->MasterDofMap(),
      *all_maps()->BlockMapScatra()->FullMap());

  const auto blockmasterderiv = masterderiv.Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
      *all_maps()->block_map_thermo(), *all_maps()->BlockMapScatra());

  blockmasterderiv->Complete();

  for (int i = 0; i < static_cast<int>(block_position_sca_tra().size()); ++i)
  {
    for (int j = 0; j < static_cast<int>(block_position_thermo().size()); ++j)
    {
      // assemble linearizations of slave side scatra fluxes w.r.t. slave and master side elch
      // into system matrix
      auto& systemmatrix_block_iscatra_jthermo =
          systemmatrix_block->Matrix(block_position_sca_tra().at(i), block_position_thermo().at(j));
      systemmatrix_block_iscatra_jthermo.Add(blockmasterderiv->Matrix(i, j), false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::assemble_scatra_thermo_interface(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrathermointerface)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatrathermointerface_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatrathermointerface);

  // assemble linearizations of slave- and master side scatra fluxes w.r.t. slave temperatures into
  // system matrix
  auto& systemmatrix_block_scatra_thermo =
      systemmatrix_block->Matrix(block_position_sca_tra().at(0), block_position_thermo().at(0));
  systemmatrix_block_scatra_thermo.Add(*scatrathermointerface_sparse, false, 1.0, 1.0);

  // assemble linearizations of slave- and master side scatra fluxes w.r.t. master temperatures into
  // system matrix
  Core::Adapter::CouplingSlaveConverter thermo_converter(*meshtying_thermo()->CouplingAdapter());

  Core::LinAlg::MatrixLogicalSplitAndTransform()(*scatrathermointerface_sparse,
      scatrathermointerface_sparse->RangeMap(),
      *meshtying_thermo()->CouplingAdapter()->MasterDofMap(), 1.0, nullptr, &thermo_converter,
      systemmatrix_block->Matrix(block_position_sca_tra().at(0), block_position_thermo().at(0)),
      true, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::assemble_scatra_thermo_interface(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatrathermointerface)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatrathermointerface_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatrathermointerface);

  // assemble linearizations of slave- and master side scatra fluxes w.r.t. slave temperatures into
  // system matrix
  systemmatrix_sparse->Add(*scatrathermointerface_sparse, false, 1.0, 1.0);

  // assemble linearizations of slave- and master side scatra fluxes w.r.t. master temperatures into
  // system matrix
  Core::Adapter::CouplingSlaveConverter thermo_converter(*meshtying_thermo()->CouplingAdapter());

  Core::LinAlg::MatrixLogicalSplitAndTransform()(*scatrathermointerface_sparse,
      scatrathermointerface_sparse->RangeMap(),
      *meshtying_thermo()->CouplingAdapter()->MasterDofMap(), 1.0, nullptr, &thermo_converter,
      *systemmatrix_sparse, true, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::assemble_structure_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> structurescatradomain)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto structurescatradomain_block =
      Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(structurescatradomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(block_position_sca_tra().size()); ++iblock)
  {
    const auto structurescatradomain_subblock = structurescatradomain_block->Matrix(0, iblock);

    // add entire block or assemble slave side to master side
    if (interface_meshtying())
    {
      assemble_structure_xxx_meshtying(
          systemmatrix_block->Matrix(position_structure(), block_position_sca_tra().at(iblock)),
          structurescatradomain_subblock);
    }
    else
    {
      auto& systemmatrix_block_struct_iscatra =
          systemmatrix_block->Matrix(position_structure(), block_position_sca_tra().at(iblock));

      systemmatrix_block_struct_iscatra.Add(structurescatradomain_subblock, false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::assemble_structure_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> structurescatradomain)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto structurescatradomain_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(structurescatradomain);

  // add entire block or assemble slave side to master side
  if (interface_meshtying())
  {
    assemble_structure_xxx_meshtying(
        systemmatrix_block->Matrix(position_structure(), block_position_sca_tra().at(0)),
        *structurescatradomain_sparse);
  }
  else
  {
    auto& systemmatrix_block_struct_scatra =
        systemmatrix_block->Matrix(position_structure(), block_position_sca_tra().at(0));

    systemmatrix_block_struct_scatra.Add(*structurescatradomain_sparse, false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::assemble_structure_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> structurescatradomain)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto structurescatradomain_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(structurescatradomain);

  // add entire block or assemble slave side to master side
  if (interface_meshtying())
    assemble_structure_xxx_meshtying(*systemmatrix_sparse, *structurescatradomain_sparse);
  else
    systemmatrix_sparse->Add(*structurescatradomain_sparse, false, 1.0, 1.0);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBase::assemble_structure_xxx_meshtying(
    Core::LinAlg::SparseMatrix& systemmatrix_structure_x,
    const Core::LinAlg::SparseMatrix& structures_x_matrix)
{
  auto map_structure_interior = ssti_structure_meshtying()->InteriorMap();
  auto master_dof_map = ssti_structure_meshtying()->FullMasterSideMap();

  // assemble derivs. of structural interior dofs w.r.t. scatra dofs
  Core::LinAlg::MatrixLogicalSplitAndTransform()(structures_x_matrix, *map_structure_interior,
      structures_x_matrix.DomainMap(), 1.0, nullptr, nullptr, systemmatrix_structure_x, true, true);

  // assemble derivs. of structural master dofs w.r.t. scatra dofs
  Core::LinAlg::MatrixLogicalSplitAndTransform()(structures_x_matrix, *master_dof_map,
      structures_x_matrix.DomainMap(), 1.0, nullptr, nullptr, systemmatrix_structure_x, true, true);

  for (const auto& meshtying : ssti_structure_meshtying()->MeshTyingHandlers())
  {
    auto cond_slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
    auto converter = meshtying->SlaveSideConverter();

    // assemble derivs. of structural surface slave dofs w.r.t. scatra dofs
    Core::LinAlg::MatrixLogicalSplitAndTransform()(structures_x_matrix, *cond_slave_dof_map(),
        structures_x_matrix.DomainMap(), 1.0, &(*converter), nullptr, systemmatrix_structure_x,
        true, true);
  }
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::assemble_thermo_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatradomain,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatrainterface)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto thermoscatradomain_block =
      Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(thermoscatradomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(block_position_thermo().size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(block_position_sca_tra().size()); ++jblock)
    {
      auto systemmatrix_block_ithermo_jscatra = systemmatrix_block->Matrix(
          block_position_thermo().at(iblock), block_position_sca_tra().at(jblock));
      systemmatrix_block_ithermo_jscatra.UnComplete();
      systemmatrix_block_ithermo_jscatra.Add(
          thermoscatradomain_block->Matrix(iblock, jblock), false, 1.0, 1.0);
    }
  }

  if (interface_meshtying()) assemble_thermo_scatra_interface(systemmatrix, thermoscatrainterface);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::assemble_thermo_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatradomain,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatrainterface)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto thermoscatradomain_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(thermoscatradomain);

  auto& systemmatrix_block_thermo_scatra =
      systemmatrix_block->Matrix(block_position_thermo().at(0), block_position_sca_tra().at(0));
  systemmatrix_block_thermo_scatra.UnComplete();
  systemmatrix_block_thermo_scatra.Add(*thermoscatradomain_sparse, false, 1.0, 1.0);

  if (interface_meshtying()) assemble_thermo_scatra_interface(systemmatrix, thermoscatrainterface);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::assemble_thermo_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatradomain,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatrainterface)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto thermoscatradomain_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(thermoscatradomain);

  // add scalar transport system matrix to global system matrix
  systemmatrix_sparse->Add(*thermoscatradomain_sparse, false, 1.0, 1.0);

  if (interface_meshtying()) assemble_thermo_scatra_interface(systemmatrix, thermoscatrainterface);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::assemble_thermo_scatra_interface(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatrainterface)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  auto thermoscatrainterface_block =
      Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(thermoscatrainterface);

  Core::LinAlg::SparseMatrix masterflux(
      *all_maps()->block_map_thermo()->FullMap(), 27, false, true);

  for (int i = 0; i < static_cast<int>(block_position_thermo().size()); ++i)
  {
    for (int j = 0; j < static_cast<int>(block_position_sca_tra().size()); ++j)
    {
      const auto thermoscatrainterface_subblock = thermoscatrainterface_block->Matrix(i, j);

      // assemble linearizations of slave side scatra fluxes w.r.t. slave and master side elch
      // into system matrix
      auto& systemmatrix_block_ithermo_jscatra =
          systemmatrix_block->Matrix(block_position_thermo().at(i), block_position_sca_tra().at(j));
      systemmatrix_block_ithermo_jscatra.Add(thermoscatrainterface_subblock, false, 1.0, 1.0);

      // assemble linearizations of master side thermo fluxes w.r.t. slave and master side elch
      // into system matrix
      Core::Adapter::CouplingSlaveConverter thermo_converter(
          *meshtying_thermo()->CouplingAdapter());

      Core::LinAlg::SparseMatrix slaveflux(
          *all_maps()->block_map_thermo()->FullMap(), 27, false, true);

      ScaTra::MeshtyingStrategyS2I::ExtractMatrixRows(
          thermoscatrainterface_subblock, slaveflux, *meshtying_thermo()->BlockMapsSlave().Map(i));

      slaveflux.Complete(
          *all_maps()->BlockMapScatra()->FullMap(), *all_maps()->block_map_thermo()->FullMap());

      Core::LinAlg::MatrixLogicalSplitAndTransform()(thermoscatrainterface_subblock,
          *meshtying_thermo()->CouplingAdapter()->MasterDofMap(),
          thermoscatrainterface_subblock.DomainMap(), -1.0, &thermo_converter, nullptr, masterflux,
          true, true);
    }
  }

  masterflux.Complete(
      *all_maps()->BlockMapScatra()->FullMap(), *all_maps()->block_map_thermo()->FullMap());

  const auto blockmasterflux = masterflux.Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
      *all_maps()->BlockMapScatra(), *all_maps()->block_map_thermo());

  blockmasterflux->Complete();

  // assemble linearizations of slave side thermo fluxes w.r.t. slave and master side elch
  // into system matrix
  for (int i = 0; i < static_cast<int>(block_position_thermo().size()); ++i)
  {
    for (int j = 0; j < static_cast<int>(block_position_sca_tra().size()); ++j)
    {
      auto& systemmatrix_block_ithermo_jscatra =
          systemmatrix_block->Matrix(block_position_thermo().at(i), block_position_sca_tra().at(j));

      systemmatrix_block_ithermo_jscatra.Add(blockmasterflux->Matrix(i, j), false, 1.0, 1.0);
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::assemble_thermo_scatra_interface(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatrainterface)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto thermoscatrainterface_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(thermoscatrainterface);

  // assemble linearizations of slave side scatra fluxes w.r.t. slave and master side elch
  // into system matrix
  auto& systemmatrix_block_thermo_scatra =
      systemmatrix_block->Matrix(block_position_thermo().at(0), block_position_sca_tra().at(0));
  systemmatrix_block_thermo_scatra.Add(*thermoscatrainterface_sparse, false, 1.0, 1.0);

  // assemble linearizations of master side thermo fluxes w.r.t. slave and master side elch
  // into system matrix
  Core::Adapter::CouplingSlaveConverter thermo_converter(*meshtying_thermo()->CouplingAdapter());

  Core::LinAlg::MatrixLogicalSplitAndTransform()(*thermoscatrainterface_sparse,
      *meshtying_thermo()->CouplingAdapter()->MasterDofMap(),
      thermoscatrainterface_sparse->DomainMap(), -1.0, &thermo_converter, nullptr,
      systemmatrix_block->Matrix(block_position_thermo().at(0), block_position_sca_tra().at(0)),
      true, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::assemble_thermo_scatra_interface(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> thermoscatrainterface)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto thermoscatrainterface_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(thermoscatrainterface);

  // assemble linearizations of slave side scatra fluxes w.r.t. slave and master side elch
  // into system matrix
  systemmatrix_sparse->Add(*thermoscatrainterface_sparse, false, 1.0, 1.0);

  // assemble linearizations of master side scatra fluxes w.r.t. slave and master side elch
  // into system matrix
  Core::Adapter::CouplingSlaveConverter thermo_converter(*meshtying_thermo()->CouplingAdapter());

  Core::LinAlg::MatrixLogicalSplitAndTransform()(*thermoscatrainterface_sparse,
      *meshtying_thermo()->CouplingAdapter()->MasterDofMap(),
      thermoscatrainterface_sparse->DomainMap(), -1.0, &thermo_converter, nullptr,
      *systemmatrix_sparse, true, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::assemble_thermo_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> thermostructuredomain,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> thermostructureinterface)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto thermostructuredomain_block =
      Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(thermostructuredomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(block_position_thermo().size()); ++iblock)
  {
    const auto thermostructuredomain_subblock = thermostructuredomain_block->Matrix(iblock, 0);

    // add entire block or assemble slave side to master side
    if (interface_meshtying())
    {
      assemble_xxx_structure_meshtying(
          systemmatrix_block->Matrix(block_position_thermo().at(iblock), position_structure()),
          thermostructuredomain_subblock);

      auto thermostructureinterface_block =
          Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(thermostructureinterface);

      assemble_xxx_structure_meshtying(
          systemmatrix_block->Matrix(block_position_thermo().at(iblock), position_structure()),
          thermostructureinterface_block->Matrix(iblock, 0));
    }
    else
    {
      auto& systemmatrix_block_ithermo_struct =
          systemmatrix_block->Matrix(block_position_thermo().at(iblock), position_structure());

      systemmatrix_block_ithermo_struct.Add(thermostructuredomain_subblock, false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::assemble_thermo_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> thermostructuredomain,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> thermostructureinterface)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto thermostructuredomain_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(thermostructuredomain);

  // add entire block or assemble slave side to master side
  if (interface_meshtying())
  {
    assemble_xxx_structure_meshtying(
        systemmatrix_block->Matrix(block_position_thermo().at(0), position_structure()),
        *thermostructuredomain_sparse);

    auto thermostructureinterface_sparse =
        Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(thermostructureinterface);

    assemble_xxx_structure_meshtying(
        systemmatrix_block->Matrix(block_position_thermo().at(0), position_structure()),
        *thermostructureinterface_sparse);
  }
  else
  {
    auto& sytemmatrix_block_thermo_structure =
        systemmatrix_block->Matrix(block_position_thermo().at(0), position_structure());

    sytemmatrix_block_thermo_structure.Add(*thermostructuredomain_sparse, false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::assemble_thermo_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> thermostructuredomain,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> thermostructureinterface)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto thermostructuredomain_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(thermostructuredomain);

  if (interface_meshtying())
  {
    assemble_xxx_structure_meshtying(*systemmatrix_sparse, *thermostructuredomain_sparse);

    auto thermostructureinterface_sparse =
        Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(thermostructureinterface);

    assemble_xxx_structure_meshtying(*systemmatrix_sparse, *thermostructureinterface_sparse);
  }
  else
    systemmatrix_sparse->Add(*thermostructuredomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::assemble_structure_thermo(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> structurethermodomain)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto structurethermodomain_block =
      Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(structurethermodomain);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(block_position_thermo().size()); ++iblock)
  {
    const auto structurethermodomain_subblock = structurethermodomain_block->Matrix(0, iblock);

    // add entire block or assemble slave side to master side
    if (interface_meshtying())
    {
      assemble_structure_xxx_meshtying(
          systemmatrix_block->Matrix(position_structure(), block_position_thermo().at(iblock)),
          structurethermodomain_subblock);
    }
    else
    {
      auto& systemmatrix_block_struct_ithermo =
          systemmatrix_block->Matrix(position_structure(), block_position_thermo().at(iblock));

      systemmatrix_block_struct_ithermo.Add(structurethermodomain_subblock, false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::assemble_structure_thermo(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> structurethermodomain)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto structurethermodomain_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(structurethermodomain);

  // add entire block or assemble slave side to master side
  if (interface_meshtying())
  {
    assemble_structure_xxx_meshtying(
        systemmatrix_block->Matrix(position_structure(), block_position_thermo().at(0)),
        *structurethermodomain_sparse);
  }
  else
  {
    auto& systemmatrix_block_struct_thermo =
        systemmatrix_block->Matrix(position_structure(), block_position_thermo().at(0));

    systemmatrix_block_struct_thermo.Add(*structurethermodomain_sparse, false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::assemble_structure_thermo(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> structurethermodomain)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto structurethermodomain_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(structurethermodomain);

  // add entire block or assemble slave side to master side
  if (interface_meshtying())
    assemble_structure_xxx_meshtying(*systemmatrix_sparse, *structurethermodomain_sparse);
  else
    systemmatrix_sparse->Add(*structurethermodomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockBlock::apply_meshtying_system_matrix(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix)
{
  if (interface_meshtying())
  {
    auto systemmatrix_block =
        Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

    apply_meshtying_sys_mat(systemmatrix_block->Matrix(position_structure(), position_structure()));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlockSparse::apply_meshtying_system_matrix(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix)
{
  if (interface_meshtying())
  {
    // cast systemmatrix
    auto systemmatrix_block =
        Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

    apply_meshtying_sys_mat(systemmatrix_block->Matrix(position_structure(), position_structure()));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::apply_meshtying_system_matrix(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix)
{
  if (interface_meshtying())
  {
    // cast systemmatrix
    auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);

    apply_meshtying_sys_mat(*systemmatrix_sparse);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBase::apply_meshtying_sys_mat(
    Core::LinAlg::SparseMatrix& systemmatrix_structure)
{
  // map for slave side structural degrees of freedom
  auto slavemaps = ssti_structure_meshtying()->FullSlaveSideMap();

  // subject slave-side rows of structural system matrix to pseudo Dirichlet conditions to
  // finalize structural meshtying
  const double one(1.0);
  for (int doflid_slave = 0; doflid_slave < slavemaps->NumMyElements(); ++doflid_slave)
  {
    // extract global ID of current slave-side row
    const int dofgid_slave = slavemaps->GID(doflid_slave);
    if (dofgid_slave < 0) FOUR_C_THROW("Local ID not found!");

    // apply pseudo Dirichlet conditions to filled matrix, i.e., to local row and column
    // indices
    if (systemmatrix_structure.Filled())
    {
      const int rowlid_slave = systemmatrix_structure.RowMap().LID(dofgid_slave);
      if (rowlid_slave < 0) FOUR_C_THROW("Global ID not found!");
      if (systemmatrix_structure.EpetraMatrix()->ReplaceMyValues(
              rowlid_slave, 1, &one, &rowlid_slave))
        FOUR_C_THROW("ReplaceMyValues failed!");
    }

    // apply pseudo Dirichlet conditions to unfilled matrix, i.e., to global row and
    // column indices
    else if (systemmatrix_structure.EpetraMatrix()->InsertGlobalValues(
                 dofgid_slave, 1, &one, &dofgid_slave))
      FOUR_C_THROW("InsertGlobalValues failed!");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategyBlock::apply_structural_dbc_system_matrix(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix)
{
  // locsys manager of strucutre
  const auto& locsysmanager_structure = structure_field()->LocsysManager();

  // map of strucutral Dirichlet BCs
  const auto dbcmap_structure = structure_field()->GetDBCMapExtractor()->CondMap();

  if (locsysmanager_structure == Teuchos::null)
    systemmatrix->ApplyDirichlet(*dbcmap_structure);
  else
  {
    auto systemmatrix_block =
        Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

    // apply structural Dirichlet conditions
    for (int iblock = 0; iblock < systemmatrix_block->Cols(); ++iblock)
    {
      locsysmanager_structure->RotateGlobalToLocal(
          Teuchos::rcp(&systemmatrix_block->Matrix(position_structure(), iblock), false));
      systemmatrix_block->Matrix(position_structure(), iblock)
          .apply_dirichlet_with_trafo(
              *locsysmanager_structure->Trafo(), *dbcmap_structure, iblock == position_structure());
      locsysmanager_structure->RotateLocalToGlobal(
          Teuchos::rcp(&systemmatrix_block->Matrix(position_structure(), iblock), false));
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::AssembleStrategySparse::apply_structural_dbc_system_matrix(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix)
{
  // locsys manager of strucutre
  const auto& locsysmanager_structure = structure_field()->LocsysManager();

  // map of strucutral Dirichlet BCs
  const auto& dbcmap_structure = structure_field()->GetDBCMapExtractor()->CondMap();

  // structural dof row map
  const auto& dofrowmap_structure = structure_field()->dof_row_map();

  if (locsysmanager_structure == Teuchos::null)
    systemmatrix->ApplyDirichlet(*dbcmap_structure);
  else
  {
    auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);

    // extract structural rows of global system matrix
    auto systemmatrix_structure =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dofrowmap_structure, 27, false, true));
    Core::LinAlg::MatrixLogicalSplitAndTransform()(*systemmatrix_sparse, *dofrowmap_structure,
        systemmatrix->DomainMap(), 1.0, nullptr, nullptr, *systemmatrix_structure);
    systemmatrix_structure->Complete(systemmatrix->DomainMap(), *dofrowmap_structure);

    // apply structural Dirichlet conditions
    locsysmanager_structure->RotateGlobalToLocal(systemmatrix_structure);
    systemmatrix_structure->apply_dirichlet_with_trafo(
        *locsysmanager_structure->Trafo(), *dbcmap_structure);
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
  all_maps()->MapsSubProblems()->InsertVector(
      RHSscatra, ssti_mono_->GetProblemPosition(Subproblem::scalar_transport), RHS);
  all_maps()->MapsSubProblems()->InsertVector(
      RHSthermo, ssti_mono_->GetProblemPosition(Subproblem::thermo), RHS);

  if (interface_meshtying())
  {
    // perform structural meshtying before assembling structural right-hand side vector into
    // monolithic right-hand side vector

    // make copy of structural right-hand side vector
    Epetra_Vector residual_structure(*RHSstructure);

    auto rhs_structure_master = Core::LinAlg::CreateVector(*structure_field()->dof_row_map(), true);

    for (const auto& meshtying : ssti_structure_meshtying()->MeshTyingHandlers())
    {
      auto coupling_adapter = meshtying->SlaveMasterCoupling();
      auto coupling_map_extractor = meshtying->slave_master_extractor();

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
    const auto& locsysmanager_structure = ssti_mono_->structure_field()->LocsysManager();

    // apply pseudo Dirichlet conditions to transformed slave-side part of structural right-hand
    // side vector
    auto zeros = Teuchos::rcp(new Epetra_Vector(residual_structure.Map()));

    if (locsysmanager_structure != Teuchos::null)
      locsysmanager_structure->RotateGlobalToLocal(rhs_structure_master);
    Core::LinAlg::apply_dirichlet_to_system(*rhs_structure_master, *zeros,
        *ssti_mono_->structure_field()->GetDBCMapExtractor()->CondMap());
    if (locsysmanager_structure != Teuchos::null)
      locsysmanager_structure->RotateLocalToGlobal(rhs_structure_master);

    // assemble master-side part of structure right-hand side vector
    residual_structure.Update(1.0, *rhs_structure_master, 1.0);

    // assemble final structural right-hand side vector into monolithic right-hand side vector
    all_maps()->MapsSubProblems()->AddVector(
        residual_structure, ssti_mono_->GetProblemPosition(Subproblem::structure), *RHS, -1.0);
  }
  else
  {
    all_maps()->MapsSubProblems()->AddVector(
        RHSstructure, ssti_mono_->GetProblemPosition(Subproblem::structure), RHS, -1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<SSTI::AssembleStrategyBase> SSTI::BuildAssembleStrategy(
    Teuchos::RCP<const SSTI::SSTIMono> ssti_mono, Core::LinAlg::MatrixType matrixtype_ssti,
    Core::LinAlg::MatrixType matrixtype_scatra)
{
  Teuchos::RCP<SSTI::AssembleStrategyBase> assemblestrategy = Teuchos::null;

  switch (matrixtype_ssti)
  {
    case Core::LinAlg::MatrixType::block_field:
    {
      switch (matrixtype_scatra)
      {
        case Core::LinAlg::MatrixType::block_condition:
        {
          assemblestrategy = Teuchos::rcp(new SSTI::AssembleStrategyBlockBlock(ssti_mono));
          break;
        }
        case Core::LinAlg::MatrixType::sparse:
        {
          assemblestrategy = Teuchos::rcp(new SSTI::AssembleStrategyBlockSparse(ssti_mono));
          break;
        }
        default:
        {
          FOUR_C_THROW("unknown matrix type");
          break;
        }
      }
      break;
    }
    case Core::LinAlg::MatrixType::sparse:
    {
      assemblestrategy = Teuchos::rcp(new SSTI::AssembleStrategySparse(ssti_mono));
      break;
    }
    default:
    {
      FOUR_C_THROW("unknown matrix type");
      break;
    }
  }

  return assemblestrategy;
}

FOUR_C_NAMESPACE_CLOSE
