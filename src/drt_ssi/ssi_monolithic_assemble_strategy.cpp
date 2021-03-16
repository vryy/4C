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

  if (SSIMono().SSIInterfaceContact())
  {
    // get scatra-scatra-interface block matrix
    const auto& scatra_scatra_interface_blockmatrix =
        SSIMono()
            .CoNitscheStrategySsi()
            ->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::scatra_scatra)
            ->Split<LINALG::DefaultBlockMatrixStrategy>(
                *SSIMono().MapsScatra(), *SSIMono().MapsScatra());
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

  if (SSIMono().SSIInterfaceContact())
  {
    systemmatrix_block_scatra_scatra.Add(*(SSIMono().CoNitscheStrategySsi()->GetMatrixBlockPtr(
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

  if (SSIMono().SSIInterfaceContact())
  {
    systemmatrix_sparse->Add(*(SSIMono().CoNitscheStrategySsi()->GetMatrixBlockPtr(
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
  auto& systemmatrix_block_struct_struct =
      systemmatrix_block->Matrix(PositionStructure(), PositionStructure());

  systemmatrix_block_struct_struct.Add(*structuredomain, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseMatrix> structuredomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto& systemmatrix_block_struct_struct =
      systemmatrix_block->Matrix(PositionStructure(), PositionStructure());

  systemmatrix_block_struct_struct.Add(*structuredomain, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseMatrix> structuredomain)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  systemmatrix_sparse->Add(*structuredomain, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScatraStructure(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain,
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);


  Teuchos::RCP<LINALG::BlockSparseMatrixBase> scatrastructuredomain_block(Teuchos::null);
  if (scatrastructuredomain != Teuchos::null)
  {
    scatrastructuredomain_block =
        LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(scatrastructuredomain);
  }

  // cast interface matrix if necessary
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> scatrastructureinterface_block = Teuchos::null;
  if (SSIMono().SSIInterfaceMeshtying())
    scatrastructureinterface_block =
        LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(scatrastructureinterface);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
  {
    auto& systemmatrix_block_iscatra_struct =
        systemmatrix_block->Matrix(BlockPositionScaTra()->at(iblock), PositionStructure());

    if (scatrastructuredomain != Teuchos::null)
    {
      systemmatrix_block_iscatra_struct.Add(
          scatrastructuredomain_block->Matrix(iblock, 0), false, 1.0, 1.0);
    }

    if (SSIMono().SSIInterfaceMeshtying())
    {
      systemmatrix_block_iscatra_struct.Add(
          scatrastructureinterface_block->Matrix(iblock, 0), false, 1.0, 1.0);
    }
  }

  if (SSIMono().SSIInterfaceContact())
  {
    // get scatra-structure-interface block matrix
    const auto& scatra_struct_interface_blockmatrix =
        SSIMono()
            .CoNitscheStrategySsi()
            ->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::scatra_displ)
            ->Split<LINALG::DefaultBlockMatrixStrategy>(
                *SSIMono().MapStructure(), *SSIMono().MapsScatra());
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

  auto& systemmatrix_block_scatra_struct =
      systemmatrix_block->Matrix(BlockPositionScaTra()->at(0), PositionStructure());

  if (scatrastructuredomain != Teuchos::null)
  {
    auto scatrastructuredomain_sparse =
        LINALG::CastToSparseMatrixAndCheckSuccess(scatrastructuredomain);

    systemmatrix_block_scatra_struct.Add(*scatrastructuredomain_sparse, false, 1.0, 1.0);
  }

  if (SSIMono().SSIInterfaceMeshtying())
  {
    auto scatrastructureinterface_sparse =
        LINALG::CastToSparseMatrixAndCheckSuccess(scatrastructureinterface);
    systemmatrix_block_scatra_struct.Add(*scatrastructureinterface_sparse, false, 1.0, 1.0);
  }

  if (SSIMono().SSIInterfaceContact())
  {
    const auto& scatra_struct_interface_matrix =
        SSIMono().CoNitscheStrategySsi()->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::scatra_displ);

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
  if (scatrastructuredomain != Teuchos::null)
  {
    auto scatrastructuredomain_sparse =
        LINALG::CastToSparseMatrixAndCheckSuccess(scatrastructuredomain);

    systemmatrix_sparse->Add(*scatrastructuredomain_sparse, false, 1.0, 1.0);
  }

  if (SSIMono().SSIInterfaceMeshtying())
  {
    auto scatrastructureinterface_sparse =
        LINALG::CastToSparseMatrixAndCheckSuccess(scatrastructureinterface);
    systemmatrix_sparse->Add(*scatrastructureinterface_sparse, false, 1.0, 1.0);
  }

  if (SSIMono().SSIInterfaceContact())
  {
    systemmatrix_sparse->Add(*(SSIMono().CoNitscheStrategySsi()->GetMatrixBlockPtr(
                                 DRT::UTILS::MatBlockType::scatra_displ)),
        false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScatraScaTraManifold(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatramanifolddomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatramanifolddomain_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(scatramanifolddomain);

  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTra()->size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(BlockPositionScaTraManifold()->size()); ++jblock)
    {
      systemmatrix_block
          ->Matrix(BlockPositionScaTra()->at(iblock), BlockPositionScaTraManifold()->at(jblock))
          .Add(scatramanifolddomain_block->Matrix(iblock, jblock), false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleScatraScaTraManifold(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatramanifolddomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatramanifolddomain_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(scatramanifolddomain);

  systemmatrix_block->Matrix(BlockPositionScaTra()->at(0), BlockPositionScaTraManifold()->at(0))
      .Add(*scatramanifolddomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleScatraScaTraManifold(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> scatramanifolddomain)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatramanifolddomain_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(scatramanifolddomain);

  systemmatrix_sparse->Add(*scatramanifolddomain_sparse, false, 1.0, 1.0);
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
    auto& systemmatrix_block_struct_iscatra =
        systemmatrix_block->Matrix(PositionStructure(), BlockPositionScaTra()->at(iblock));
    systemmatrix_block_struct_iscatra.Add(
        structurescatradomain_block->Matrix(0, iblock), false, 1.0, 1.0);
  }

  if (SSIMono().SSIInterfaceContact())
  {
    // get structure-scatra-interface block matrix
    const auto& struct_scatra_interface_blockmatrix =
        SSIMono()
            .CoNitscheStrategySsi()
            ->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::displ_scatra)
            ->Split<LINALG::DefaultBlockMatrixStrategy>(
                *SSIMono().MapsScatra(), *SSIMono().MapStructure());
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

  auto& systemmatrix_block_struct_scatra =
      systemmatrix_block->Matrix(PositionStructure(), BlockPositionScaTra()->at(0));
  systemmatrix_block_struct_scatra.Add(*structurescatradomain_sparse, false, 1.0, 1.0);

  if (SSIMono().SSIInterfaceContact())
  {
    const auto& struct_scatra_interface_matrix =
        SSIMono().CoNitscheStrategySsi()->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::displ_scatra);

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

  systemmatrix_sparse->Add(*structurescatradomain_sparse, false, 1.0, 1.0);

  // add contact contributions
  if (SSIMono().SSIInterfaceContact())
  {
    systemmatrix_sparse->Add(*(SSIMono().CoNitscheStrategySsi()->GetMatrixBlockPtr(
                                 DRT::UTILS::MatBlockType::displ_scatra)),
        false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::AssembleScaTraManifoldScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> manifoldscatradomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto manifoldscatradomain_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(manifoldscatradomain);

  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTraManifold()->size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(BlockPositionScaTra()->size()); ++jblock)
    {
      systemmatrix_block
          ->Matrix(BlockPositionScaTraManifold()->at(iblock), BlockPositionScaTra()->at(jblock))
          .Add(manifoldscatradomain_block->Matrix(iblock, jblock), false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::AssembleScaTraManifoldScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> manifoldscatradomain)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto manifoldscatradomain_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(manifoldscatradomain);

  systemmatrix_block->Matrix(BlockPositionScaTraManifold()->at(0), BlockPositionScaTra()->at(0))
      .Add(*manifoldscatradomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::AssembleScaTraManifoldScatra(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<LINALG::SparseOperator> manifoldscatradomain)
{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto manifoldscatradomain_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(manifoldscatradomain);

  systemmatrix_sparse->Add(*manifoldscatradomain_sparse, false, 1.0, 1.0);
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

  for (int iblock = 0; iblock < static_cast<int>(BlockPositionScaTraManifold()->size()); ++iblock)
  {
    auto& systemmatrix_block_iscatramanifold_struct =
        systemmatrix_block->Matrix(BlockPositionScaTraManifold()->at(iblock), PositionStructure());
    systemmatrix_block_iscatramanifold_struct.Add(
        manifoldstructuredomain_block->Matrix(iblock, 0), false, 1.0, 1.0);
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

  auto& systemmatrix_block_scatramanifold_struct =
      systemmatrix_block->Matrix(BlockPositionScaTraManifold()->at(0), PositionStructure());
  systemmatrix_block_scatramanifold_struct.Add(*manifoldstructuredomain_sparse, false, 1.0, 1.0);
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

  systemmatrix_sparse->Add(*manifoldstructuredomain_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::AssembleRHS(Teuchos::RCP<Epetra_Vector> rhs,
    Teuchos::RCP<const Epetra_Vector> rhs_scatra, Teuchos::RCP<const Epetra_Vector> rhs_structure,
    Teuchos::RCP<const Epetra_Vector> rhs_manifold,
    Teuchos::RCP<const Epetra_Vector> flux_manifold_scatra_m,
    Teuchos::RCP<const Epetra_Vector> flux_manifold_scatra_s)
{
  SSIMono().MapsSubProblems()->InsertVector(
      rhs_scatra, SSIMono().GetProblemPosition(SSI::Subproblem::scalar_transport), rhs);

  if (SSIMono().IsScaTraManifold())
  {
    SSIMono().MapsSubProblems()->InsertVector(
        rhs_manifold, SSIMono().GetProblemPosition(SSI::Subproblem::manifold), rhs);

    SSIMono().MapsSubProblems()->InsertVector(
        rhs_manifold, SSIMono().GetProblemPosition(SSI::Subproblem::manifold), rhs);

    SSIMono().MapsSubProblems()->AddVector(
        flux_manifold_scatra_m, SSIMono().GetProblemPosition(Subproblem::manifold), rhs);

    SSIMono().MapsSubProblems()->AddVector(
        flux_manifold_scatra_s, SSIMono().GetProblemPosition(Subproblem::scalar_transport), rhs);
  }

  SSIMono().MapsSubProblems()->AddVector(
      rhs_structure, SSIMono().GetProblemPosition(SSI::Subproblem::structure), rhs, -1.0);

  if (SSIMono().SSIInterfaceContact())
  {
    // add the scatra contact contribution
    SSIMono().MapsSubProblems()->AddVector(
        SSIMono().CoNitscheStrategySsi()->GetRhsBlockPtr(DRT::UTILS::VecBlockType::scatra),
        SSIMono().GetProblemPosition(SSI::Subproblem::scalar_transport), rhs);
  }
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