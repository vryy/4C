/*----------------------------------------------------------------------*/
/*! \file
\brief Dirichlet boundary condition handler for monolithic SSI

\level 2

 */
/*----------------------------------------------------------------------*/

#include "4C_ssi_monolithic_dbc_handler.hpp"

#include "4C_adapter_str_ssiwrapper.hpp"
#include "4C_fem_condition_locsys.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_ssi_monolithic.hpp"
#include "4C_ssi_utils.hpp"

FOUR_C_NAMESPACE_OPEN

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::DBCHandlerBase::DBCHandlerBase(const bool is_scatra_manifold,
    Teuchos::RCP<ScaTra::ScaTraTimIntImpl> scatra,
    Teuchos::RCP<ScaTra::ScaTraTimIntImpl> scatra_manifold,
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
    Teuchos::RCP<Adapter::SSIStructureWrapper> structure)
    : is_scatra_manifold_(is_scatra_manifold),
      scatra_(std::move(scatra)),
      scatra_manifold_(std::move(scatra_manifold)),
      ssi_maps_(std::move(ssi_maps)),
      structure_(std::move(structure))
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::DBCHandlerSparse::DBCHandlerSparse(const bool is_scatra_manifold,
    Teuchos::RCP<ScaTra::ScaTraTimIntImpl> scatra,
    Teuchos::RCP<ScaTra::ScaTraTimIntImpl> scatra_manifold,
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
    Teuchos::RCP<Adapter::SSIStructureWrapper> structure)
    : DBCHandlerBase(is_scatra_manifold, scatra, scatra_manifold, ssi_maps, structure)
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::DBCHandlerBlock::DBCHandlerBlock(const bool is_scatra_manifold,
    Teuchos::RCP<ScaTra::ScaTraTimIntImpl> scatra,
    Teuchos::RCP<ScaTra::ScaTraTimIntImpl> scatra_manifold,
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
    Teuchos::RCP<Adapter::SSIStructureWrapper> structure)
    : DBCHandlerBase(is_scatra_manifold, scatra, scatra_manifold, ssi_maps, structure),
      position_structure_(ssi_maps()->get_block_positions(SSI::Subproblem::structure).at(0))
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::DBCHandlerBase::apply_dbc_to_rhs(Teuchos::RCP<Epetra_Vector> rhs)
{
  // apply Dirichlet boundary conditions to the structure part of the right hand side
  const auto& locsysmanager_structure = structure_field()->locsys_manager();
  auto rhs_struct = ssi_maps()->maps_sub_problems()->extract_vector(
      rhs, UTILS::SSIMaps::get_problem_position(SSI::Subproblem::structure));
  const auto zeros_struct =
      Teuchos::rcp(new Epetra_Vector(*structure_field()->get_dbc_map_extractor()->cond_map()));

  if (locsysmanager_structure != Teuchos::null)
    locsysmanager_structure->rotate_global_to_local(rhs_struct);
  Core::LinAlg::apply_dirichlet_to_system(
      *rhs_struct, *zeros_struct, *structure_field()->get_dbc_map_extractor()->cond_map());
  if (locsysmanager_structure != Teuchos::null)
    locsysmanager_structure->rotate_local_to_global(rhs_struct);

  ssi_maps()->maps_sub_problems()->insert_vector(
      rhs_struct, UTILS::SSIMaps::get_problem_position(SSI::Subproblem::structure), rhs);

  // apply Dirichlet boundary conditions to the scatra part of the right hand side
  const auto zeros_scatra =
      Teuchos::rcp(new Epetra_Vector(*sca_tra_field()->dirich_maps()->cond_map()));
  Core::LinAlg::apply_dirichlet_to_system(
      *rhs, *zeros_scatra, *sca_tra_field()->dirich_maps()->cond_map());

  // apply Dirichlet boundary conditions to the scatra manifold part of the right hand side
  if (is_sca_tra_manifold())
  {
    const auto zeros_scatramanifold =
        Teuchos::rcp(new Epetra_Vector(*sca_tra_manifold_field()->dirich_maps()->cond_map()));
    Core::LinAlg::apply_dirichlet_to_system(
        *rhs, *zeros_scatramanifold, *sca_tra_manifold_field()->dirich_maps()->cond_map());
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::DBCHandlerBase::apply_dbc_to_system_matrix(
    Teuchos::RCP<Core::LinAlg::SparseOperator> system_matrix)
{
  // apply the scalar transport Dirichlet boundary conditions to the global system matrix
  system_matrix->apply_dirichlet(*sca_tra_field()->dirich_maps()->cond_map(), true);

  // apply the scalar transport on manifolds Dirichlet boundary conditions to the global system
  // matrix
  if (is_sca_tra_manifold())
    system_matrix->apply_dirichlet(*sca_tra_manifold_field()->dirich_maps()->cond_map(), true);

  // apply the structure Dirichlet boundary conditions to the global system matrix
  apply_structure_dbc_to_system_matrix(system_matrix);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::DBCHandlerBase::apply_structure_dbc_to_system_matrix(
    Teuchos::RCP<Core::LinAlg::SparseOperator> system_matrix)
{
  // locsys manager of structure
  const auto& locsysmanager_structure = structure_field()->locsys_manager();

  // map of structure Dirichlet BCs
  const auto& dbcmap_structure = structure_field()->get_dbc_map_extractor()->cond_map();

  if (locsysmanager_structure == Teuchos::null)
    system_matrix->apply_dirichlet(*dbcmap_structure);
  else
    apply_structure_dbc_with_loc_sys_rotation_to_system_matrix(
        system_matrix, dbcmap_structure, locsysmanager_structure);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::DBCHandlerSparse::apply_structure_dbc_with_loc_sys_rotation_to_system_matrix(
    Teuchos::RCP<Core::LinAlg::SparseOperator> system_matrix,
    const Teuchos::RCP<const Epetra_Map>& dbcmap_structure,
    Teuchos::RCP<const Core::Conditions::LocsysManager> locsysmanager_structure)

{
  auto systemmatrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(system_matrix);

  // structure dof row map
  const auto& dofrowmap_structure = structure_field()->dof_row_map();

  // extract structure rows of global system matrix
  const auto systemmatrix_structure =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dofrowmap_structure, 27, false, true));
  Core::LinAlg::MatrixLogicalSplitAndTransform()(*systemmatrix_sparse, *dofrowmap_structure,
      system_matrix->domain_map(), 1.0, nullptr, nullptr, *systemmatrix_structure);
  systemmatrix_structure->complete(system_matrix->domain_map(), *dofrowmap_structure);

  // apply structure Dirichlet conditions
  locsysmanager_structure->rotate_global_to_local(systemmatrix_structure);
  systemmatrix_structure->apply_dirichlet_with_trafo(
      *locsysmanager_structure->trafo(), *dbcmap_structure);
  locsysmanager_structure->rotate_local_to_global(systemmatrix_structure);

  // assemble structure rows of global system matrix back into global system matrix
  systemmatrix_sparse->put(*systemmatrix_structure, 1.0, dofrowmap_structure);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::DBCHandlerBlock::apply_structure_dbc_with_loc_sys_rotation_to_system_matrix(
    Teuchos::RCP<Core::LinAlg::SparseOperator> system_matrix,
    const Teuchos::RCP<const Epetra_Map>& dbcmap_structure,
    Teuchos::RCP<const Core::Conditions::LocsysManager> locsysmanager_structure)
{
  auto systemmatrix_block = Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(system_matrix);

  // apply structure Dirichlet conditions
  for (int iblock = 0; iblock < systemmatrix_block->cols(); ++iblock)
  {
    locsysmanager_structure->rotate_global_to_local(
        Teuchos::rcp(&systemmatrix_block->matrix(position_structure(), iblock), false));
    systemmatrix_block->matrix(position_structure(), iblock)
        .apply_dirichlet_with_trafo(
            *locsysmanager_structure->trafo(), *dbcmap_structure, (iblock == position_structure()));
    locsysmanager_structure->rotate_local_to_global(
        Teuchos::rcp(&systemmatrix_block->matrix(position_structure(), iblock), false));
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
Teuchos::RCP<SSI::DBCHandlerBase> SSI::BuildDBCHandler(const bool is_scatra_manifold,
    Core::LinAlg::MatrixType matrixtype_ssi, Teuchos::RCP<ScaTra::ScaTraTimIntImpl> scatra,
    Teuchos::RCP<ScaTra::ScaTraTimIntImpl> scatra_manifold,
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
    Teuchos::RCP<Adapter::SSIStructureWrapper> structure)
{
  Teuchos::RCP<SSI::DBCHandlerBase> dbc_handler = Teuchos::null;

  switch (matrixtype_ssi)
  {
    case Core::LinAlg::MatrixType::block_field:
    {
      dbc_handler = Teuchos::rcp(new SSI::DBCHandlerBlock(
          is_scatra_manifold, scatra, scatra_manifold, ssi_maps, structure));
      break;
    }
    case Core::LinAlg::MatrixType::sparse:
    {
      dbc_handler = Teuchos::rcp(new SSI::DBCHandlerSparse(
          is_scatra_manifold, scatra, scatra_manifold, ssi_maps, structure));
      break;
    }
    default:
    {
      FOUR_C_THROW("unknown matrix type of SSI problem");
      break;
    }
  }

  return dbc_handler;
}
FOUR_C_NAMESPACE_CLOSE
