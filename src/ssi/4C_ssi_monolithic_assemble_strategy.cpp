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
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, const bool is_scatra_manifold_value)
    : AssembleStrategyBase(ssi_maps, is_scatra_manifold_value),
      block_position_scatra_(ssi_maps()->get_block_positions(SSI::Subproblem::scalar_transport)),
      position_structure_(ssi_maps()->get_block_positions(SSI::Subproblem::structure).at(0))
{
  if (is_scatra_manifold())
    block_position_scatra_manifold_ = ssi_maps()->get_block_positions(SSI::Subproblem::manifold);
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
  systemmatrix_block->un_complete();

  // assemble blocks of scalar transport system matrix into global system matrix
  for (int iblock = 0; iblock < static_cast<int>(block_position_scatra().size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(block_position_scatra().size()); ++jblock)
    {
      auto& systemmatrix_block_iscatra_jscatra = systemmatrix_block->matrix(
          block_position_scatra().at(iblock), block_position_scatra().at(jblock));

      systemmatrix_block_iscatra_jscatra.add(
          scatra_scatra_matrix_block->matrix(iblock, jblock), false, 1.0, 1.0);
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
      systemmatrix_block->matrix(block_position_scatra().at(0), block_position_scatra().at(0));

  systemmatrix_block_scatra_scatra.add(*scatra_scatra_matrix_sparse, false, 1.0, 1.0);
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
  systemmatrix_sparse->add(*scatra_scatra_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockBlock::assemble_structure_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> structure_structure_matrix)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto& systemmatrix_block_struct_struct =
      systemmatrix_block->matrix(position_structure(), position_structure());

  systemmatrix_block_struct_struct.add(*structure_structure_matrix, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBlockSparse::assemble_structure_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> structure_structure_matrix)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);
  auto& systemmatrix_block_struct_struct =
      systemmatrix_block->matrix(position_structure(), position_structure());

  systemmatrix_block_struct_struct.add(*structure_structure_matrix, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategySparse::assemble_structure_structure(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> structure_structure_matrix)
{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(systemmatrix);
  systemmatrix_sparse->add(*structure_structure_matrix, false, 1.0, 1.0);
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
  for (int iblock = 0; iblock < static_cast<int>(block_position_scatra().size()); ++iblock)
  {
    auto& systemmatrix_block_iscatra_struct =
        systemmatrix_block->matrix(block_position_scatra().at(iblock), position_structure());

    systemmatrix_block_iscatra_struct.add(
        scatra_structure_matrix_block->matrix(iblock, 0), false, 1.0, 1.0);
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
      systemmatrix_block->matrix(block_position_scatra().at(0), position_structure());
  systemmatrix_block_scatra_struct.un_complete();

  systemmatrix_block_scatra_struct.add(*scatra_structure_matrix_sparse, false, 1.0, 1.0);
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

  systemmatrix_sparse->add(*scatra_structure_matrix_sparse, false, 1.0, 1.0);
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

  for (int iblock = 0; iblock < static_cast<int>(block_position_scatra().size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(block_position_scatra_manifold().size());
         ++jblock)
    {
      systemmatrix_block
          ->matrix(block_position_scatra().at(iblock), block_position_scatra_manifold().at(jblock))
          .add(scatra_scatramanifold_matrix_block->matrix(iblock, jblock), false, 1.0, 1.0);
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

  systemmatrix_block->matrix(block_position_scatra().at(0), block_position_scatra_manifold().at(0))
      .add(*scatra_scatramanifold_matrix_sparse, false, 1.0, 1.0);
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

  systemmatrix_sparse->add(*scatra_scatramanifold_matrix_sparse, false, 1.0, 1.0);
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
  for (int iblock = 0; iblock < static_cast<int>(block_position_scatra().size()); ++iblock)
  {
    auto& systemmatrix_block_struct_iscatra =
        systemmatrix_block->matrix(position_structure(), block_position_scatra().at(iblock));
    systemmatrix_block_struct_iscatra.add(
        structure_scatra_matrix_block->matrix(0, iblock), false, 1.0, 1.0);
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
      systemmatrix_block->matrix(position_structure(), block_position_scatra().at(0));
  systemmatrix_block_struct_scatra.un_complete();
  systemmatrix_block_struct_scatra.add(*structure_scatra_matrix_sparse, false, 1.0, 1.0);
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

  systemmatrix_sparse->add(*structure_scatra_matrix_sparse, false, 1.0, 1.0);
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

  for (int iblock = 0; iblock < static_cast<int>(block_position_scatra_manifold().size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(block_position_scatra().size()); ++jblock)
    {
      systemmatrix_block
          ->matrix(block_position_scatra_manifold().at(iblock), block_position_scatra().at(jblock))
          .add(scatramanifold_scatra_matrix_block->matrix(iblock, jblock), false, 1.0, 1.0);
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

  systemmatrix_block->matrix(block_position_scatra_manifold().at(0), block_position_scatra().at(0))
      .add(*scatramanifold_scatra_matrix_sparse, false, 1.0, 1.0);
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

  systemmatrix_sparse->add(*scatramanifold_scatra_matrix_sparse, false, 1.0, 1.0);
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
  for (int iblock = 0; iblock < static_cast<int>(block_position_scatra_manifold().size()); ++iblock)
  {
    for (int jblock = 0; jblock < static_cast<int>(block_position_scatra_manifold().size());
         ++jblock)
    {
      auto& systemmatrix_block_iscatramanifold_jscatramanifold = systemmatrix_block->matrix(
          block_position_scatra_manifold().at(iblock), block_position_scatra_manifold().at(jblock));

      systemmatrix_block_iscatramanifold_jscatramanifold.add(
          scatramanifold_scatramanifold_matrix_block->matrix(iblock, jblock), false, 1.0, 1.0);
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

  auto& systemmatrix_block_scatramanifold_scatramanifold = systemmatrix_block->matrix(
      block_position_scatra_manifold().at(0), block_position_scatra_manifold().at(0));

  systemmatrix_block_scatramanifold_scatramanifold.add(
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

  systemmatrix_sparse->add(*scatramanifold_scatramanifold_matrix_sparse, false, 1.0, 1.0);
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

  for (int iblock = 0; iblock < static_cast<int>(block_position_scatra_manifold().size()); ++iblock)
  {
    auto& systemmatrix_block_iscatramanifold_struct = systemmatrix_block->matrix(
        block_position_scatra_manifold().at(iblock), position_structure());
    systemmatrix_block_iscatramanifold_struct.add(
        scatramanifold_structure_matrix_block->matrix(iblock, 0), false, 1.0, 1.0);
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
      systemmatrix_block->matrix(block_position_scatra_manifold().at(0), position_structure());
  systemmatrix_block_scatramanifold_struct.add(
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

  systemmatrix_sparse->add(*scatramanifold_structure_matrix_sparse, false, 1.0, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::AssembleStrategyBase::assemble_rhs(Teuchos::RCP<Epetra_Vector> rhs,
    Teuchos::RCP<const Epetra_Vector> rhs_scatra, Teuchos::RCP<const Epetra_Vector> rhs_structure,
    Teuchos::RCP<const Epetra_Vector> rhs_manifold)
{
  ssi_maps()->maps_sub_problems()->insert_vector(
      rhs_scatra, UTILS::SSIMaps::get_problem_position(SSI::Subproblem::scalar_transport), rhs);

  if (is_scatra_manifold())
  {
    ssi_maps()->maps_sub_problems()->insert_vector(
        rhs_manifold, UTILS::SSIMaps::get_problem_position(SSI::Subproblem::manifold), rhs);
  }

  ssi_maps()->maps_sub_problems()->add_vector(
      rhs_structure, UTILS::SSIMaps::get_problem_position(SSI::Subproblem::structure), rhs, -1.0);
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
