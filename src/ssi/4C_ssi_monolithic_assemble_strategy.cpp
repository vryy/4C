/*----------------------------------------------------------------------*/
/*! \file
\brief Assemble strategy for monolithic SSI
\level 2

 */
/*----------------------------------------------------------------------*/
#include "4C_ssi_monolithic_assemble_strategy.hpp"

#include "4C_adapter_str_ssiwrapper.hpp"
#include "4C_contact_nitsche_strategy_ssi.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_fem_condition_locsys.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_scatra_timint_meshtying_strategy_s2i.hpp"
#include "4C_ssi_monolithic.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategyBase::AssembleStrategyBase(
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, const bool is_scatra_manifold)
    : is_scatra_manifold_(is_scatra_manifold), ssi_maps_(std::move(ssi_maps))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategyBlock::AssembleStrategyBlock(
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, const bool is_scatra_manifold)
    : AssembleStrategyBase(ssi_maps, is_scatra_manifold),
      block_position_scatra_(ssi_maps()->GetBlockPositions(SSI::Subproblem::scalar_transport)),
      position_structure_(ssi_maps()->GetBlockPositions(SSI::Subproblem::structure).at(0))
{
  if (is_sca_tra_manifold())
    block_position_scatra_manifold_ = ssi_maps()->GetBlockPositions(SSI::Subproblem::manifold);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategyBlockBlock::AssembleStrategyBlockBlock(
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, const bool is_scatra_manifold)
    : AssembleStrategyBlock(ssi_maps, is_scatra_manifold)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategyBlockSparse::AssembleStrategyBlockSparse(
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, const bool is_scatra_manifold)
    : AssembleStrategyBlock(ssi_maps, is_scatra_manifold)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::AssembleStrategySparse::AssembleStrategySparse(
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, const bool is_scatra_manifold)
    : AssembleStrategyBase(ssi_maps, is_scatra_manifold)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::assemble_scatra_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatra_scatra_matrix)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatra_scatra_matrix_block =
      Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatra_scatra_matrix);
  systemmatrix_block->UnComplete();

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(block_position_sca_tra().size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(block_position_sca_tra().size()); ++jblock)
    {
      auto& systemmatrix_block_iscatra_jscatra = systemmatrix_block->Matrix(
          block_position_sca_tra().at(iblock), block_position_sca_tra().at(jblock));

      systemmatrix_block_iscatra_jscatra.Add(
          scatra_scatra_matrix_block->Matrix(iblock, jblock), false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::assemble_scatra_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatra_scatra_matrix)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatra_scatra_matrix_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatra_scatra_matrix);

  auto& systemmatrix_block_scatra_scatra =
      systemmatrix_block->Matrix(block_position_sca_tra().at(0), block_position_sca_tra().at(0));

  systemmatrix_block_scatra_scatra.Add(*scatra_scatra_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::assemble_scatra_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatra_scatra_matrix)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatra_scatra_matrix_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatra_scatra_matrix);

  // add scalar transport system matrix to global system matrix
  systemmatrix_sparse->Add(*scatra_scatra_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::assemble_structure_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> structure_structure_matrix)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto& systemmatrix_block_struct_struct =
      systemmatrix_block->Matrix(position_structure(), position_structure());

  systemmatrix_block_struct_struct.Add(*structure_structure_matrix, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::assemble_structure_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> structure_structure_matrix)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto& systemmatrix_block_struct_struct =
      systemmatrix_block->Matrix(position_structure(), position_structure());

  systemmatrix_block_struct_struct.Add(*structure_structure_matrix, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::assemble_structure_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> structure_structure_matrix)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  systemmatrix_sparse->Add(*structure_structure_matrix, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::assemble_scatra_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatra_structure_matrix)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatra_structure_matrix_block =
      Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatra_structure_matrix);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(block_position_sca_tra().size()); ++iblock)
  {
    auto& systemmatrix_block_iscatra_struct =
        systemmatrix_block->Matrix(block_position_sca_tra().at(iblock), position_structure());

    systemmatrix_block_iscatra_struct.Add(
        scatra_structure_matrix_block->Matrix(iblock, 0), false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::assemble_scatra_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatra_structure_matrix)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatra_structure_matrix_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatra_structure_matrix);

  auto& systemmatrix_block_scatra_struct =
      systemmatrix_block->Matrix(block_position_sca_tra().at(0), position_structure());
  systemmatrix_block_scatra_struct.UnComplete();

  systemmatrix_block_scatra_struct.Add(*scatra_structure_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::assemble_scatra_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatra_structure_matrix)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatra_structure_matrix_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatra_structure_matrix);

  systemmatrix_sparse->Add(*scatra_structure_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::assemble_scatra_scatramanifold(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatra_scatramanifold_matrix)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatra_scatramanifold_matrix_block =
      Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatra_scatramanifold_matrix);

  for (int iblock = 0; iblock < static_cast<int>(block_position_sca_tra().size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(block_position_sca_tra_manifold().size());
         ++jblock)
    {
      systemmatrix_block
          ->Matrix(
              block_position_sca_tra().at(iblock), block_position_sca_tra_manifold().at(jblock))
          .Add(scatra_scatramanifold_matrix_block->Matrix(iblock, jblock), false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::assemble_scatra_scatramanifold(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatra_scatramanifold_matrix)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatra_scatramanifold_matrix_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatra_scatramanifold_matrix);

  systemmatrix_block->Matrix(
                        block_position_sca_tra().at(0), block_position_sca_tra_manifold().at(0))
      .Add(*scatra_scatramanifold_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::assemble_scatra_scatramanifold(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatra_scatramanifold_matrix)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatra_scatramanifold_matrix_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatra_scatramanifold_matrix);

  systemmatrix_sparse->Add(*scatra_scatramanifold_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::assemble_structure_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> structure_scatra_matrix)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto structure_scatra_matrix_block =
      Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(structure_scatra_matrix);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(block_position_sca_tra().size()); ++iblock)
  {
    auto& systemmatrix_block_struct_iscatra =
        systemmatrix_block->Matrix(position_structure(), block_position_sca_tra().at(iblock));
    systemmatrix_block_struct_iscatra.Add(
        structure_scatra_matrix_block->Matrix(0, iblock), false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::assemble_structure_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> structure_scatra_matrix)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto structure_scatra_matrix_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(structure_scatra_matrix);

  auto& systemmatrix_block_struct_scatra =
      systemmatrix_block->Matrix(position_structure(), block_position_sca_tra().at(0));
  systemmatrix_block_struct_scatra.UnComplete();
  systemmatrix_block_struct_scatra.Add(*structure_scatra_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::assemble_structure_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> structure_scatra_matrix)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto structure_scatra_matrix_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(structure_scatra_matrix);

  systemmatrix_sparse->Add(*structure_scatra_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::assemble_scatramanifold_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatramanifold_scatra_matrix)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatramanifold_scatra_matrix_block =
      Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatramanifold_scatra_matrix);

  for (int iblock = 0; iblock < static_cast<int>(block_position_sca_tra_manifold().size());
       ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(block_position_sca_tra().size()); ++jblock)
    {
      systemmatrix_block
          ->Matrix(
              block_position_sca_tra_manifold().at(iblock), block_position_sca_tra().at(jblock))
          .Add(scatramanifold_scatra_matrix_block->Matrix(iblock, jblock), false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::assemble_scatramanifold_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatramanifold_scatra_matrix)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatramanifold_scatra_matrix_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatramanifold_scatra_matrix);

  systemmatrix_block->Matrix(
                        block_position_sca_tra_manifold().at(0), block_position_sca_tra().at(0))
      .Add(*scatramanifold_scatra_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::assemble_scatramanifold_scatra(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatramanifold_scatra_matrix)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatramanifold_scatra_matrix_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatramanifold_scatra_matrix);

  systemmatrix_sparse->Add(*scatramanifold_scatra_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::assemble_scatramanifold_scatramanifold(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatramanifold_scatramanifold_matrix)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatramanifold_scatramanifold_matrix_block =
      Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(
          scatramanifold_scatramanifold_matrix);

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(block_position_sca_tra_manifold().size());
       ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(block_position_sca_tra_manifold().size());
         ++jblock)
    {
      auto& systemmatrix_block_iscatramanifold_jscatramanifold =
          systemmatrix_block->Matrix(block_position_sca_tra_manifold().at(iblock),
              block_position_sca_tra_manifold().at(jblock));

      systemmatrix_block_iscatramanifold_jscatramanifold.Add(
          scatramanifold_scatramanifold_matrix_block->Matrix(iblock, jblock), false, 1.0, 1.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::assemble_scatramanifold_scatramanifold(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatramanifold_scatramanifold_matrix)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatramanifold_scatramanifold_matrix_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatramanifold_scatramanifold_matrix);

  auto& systemmatrix_block_scatramanifold_scatramanifold = systemmatrix_block->Matrix(
      block_position_sca_tra_manifold().at(0), block_position_sca_tra_manifold().at(0));

  systemmatrix_block_scatramanifold_scatramanifold.Add(
      *scatramanifold_scatramanifold_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::assemble_scatramanifold_scatramanifold(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatramanifold_scatramanifold_matrix)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatramanifold_scatramanifold_matrix_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatramanifold_scatramanifold_matrix);

  systemmatrix_sparse->Add(*scatramanifold_scatramanifold_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::assemble_scatramanifold_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatramanifold_structure_matrix)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatramanifold_structure_matrix_block =
      Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(
          scatramanifold_structure_matrix);

  for (int iblock = 0; iblock < static_cast<int>(block_position_sca_tra_manifold().size());
       ++iblock)
  {
    auto& systemmatrix_block_iscatramanifold_struct = systemmatrix_block->Matrix(
        block_position_sca_tra_manifold().at(iblock), position_structure());
    systemmatrix_block_iscatramanifold_struct.Add(
        scatramanifold_structure_matrix_block->Matrix(iblock, 0), false, 1.0, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::assemble_scatramanifold_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatramanifold_structure_matrix)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto scatramanifold_structure_matrix_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatramanifold_structure_matrix);

  auto& systemmatrix_block_scatramanifold_struct =
      systemmatrix_block->Matrix(block_position_sca_tra_manifold().at(0), position_structure());
  systemmatrix_block_scatramanifold_struct.Add(
      *scatramanifold_structure_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::assemble_scatramanifold_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseOperator> scatramanifold_structure_matrix)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  auto scatramanifold_structure_matrix_sparse =
      Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(scatramanifold_structure_matrix);

  systemmatrix_sparse->Add(*scatramanifold_structure_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::AssembleRHS(Teuchos::RCP<Epetra_Vector> rhs,
    Teuchos::RCP<const Epetra_Vector> rhs_scatra, Teuchos::RCP<const Epetra_Vector> rhs_structure,
    Teuchos::RCP<const Epetra_Vector> rhs_manifold)
{
  ssi_maps()->MapsSubProblems()->InsertVector(
      rhs_scatra, UTILS::SSIMaps::GetProblemPosition(SSI::Subproblem::scalar_transport), rhs);

  if (is_sca_tra_manifold())
  {
    ssi_maps()->MapsSubProblems()->InsertVector(
        rhs_manifold, UTILS::SSIMaps::GetProblemPosition(SSI::Subproblem::manifold), rhs);
  }

  ssi_maps()->MapsSubProblems()->AddVector(
      rhs_structure, UTILS::SSIMaps::GetProblemPosition(SSI::Subproblem::structure), rhs, -1.0);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
Teuchos::RCP<SSI::AssembleStrategyBase> SSI::BuildAssembleStrategy(
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, const bool is_scatra_manifold,
    Core::LinAlg::MatrixType matrixtype_ssi, Core::LinAlg::MatrixType matrixtype_scatra)
{
  Teuchos::RCP<SSI::AssembleStrategyBase> assemblestrategy = Teuchos::null;

  switch (matrixtype_ssi)
  {
    case Core::LinAlg::MatrixType::block_field:
    {
      switch (matrixtype_scatra)
      {
        case Core::LinAlg::MatrixType::block_condition:
        case Core::LinAlg::MatrixType::block_condition_dof:
        {
          assemblestrategy =
              Teuchos::rcp(new SSI::AssembleStrategyBlockBlock(ssi_maps, is_scatra_manifold));
          break;
        }
        case Core::LinAlg::MatrixType::sparse:
        {
          assemblestrategy =
              Teuchos::rcp(new SSI::AssembleStrategyBlockSparse(ssi_maps, is_scatra_manifold));
          break;
        }

        default:
        {
          FOUR_C_THROW("unknown matrix type of ScaTra field");
          break;
        }
      }
      break;
    }
    case Core::LinAlg::MatrixType::sparse:
    {
      assemblestrategy =
          Teuchos::rcp(new SSI::AssembleStrategySparse(ssi_maps, is_scatra_manifold));
      break;
    }
    default:
    {
      FOUR_C_THROW("unknown matrix type of SSI problem");
      break;
    }
  }

  return assemblestrategy;
}
FOUR_C_NAMESPACE_CLOSE
