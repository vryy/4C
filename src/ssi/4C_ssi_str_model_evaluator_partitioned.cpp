/*-----------------------------------------------------------*/
/*! \file
\brief Model evaluator for structure part of partitioned ssi.

\level 3


*/
/*-----------------------------------------------------------*/


#include "4C_ssi_str_model_evaluator_partitioned.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_adapter_str_ssiwrapper.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_solver_nonlin_nox_group.hpp"
#include "4C_ssi_partitioned.hpp"
#include "4C_ssi_utils.hpp"
#include "4C_structure_new_dbc.hpp"
#include "4C_structure_new_impl_generic.hpp"
#include "4C_structure_new_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::PartitionedSSI::PartitionedSSI(const Teuchos::RCP<const SSI::SSIPart> ssi_part)
    : ssi_part_(ssi_part)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::PartitionedSSI::assemble_jacobian(
    CORE::LINALG::SparseOperator& jac, const double& timefac_np) const
{
  // perform structural meshtying
  if (ssi_part_->ssi_interface_meshtying())
  {
    // cast old Jacobian
    auto& jac_sparse = dynamic_cast<CORE::LINALG::SparseMatrix&>(jac);

    auto map_structure_interior = ssi_part_->ssi_structure_mesh_tying()->InteriorMap();
    auto cond_master_dof_map = ssi_part_->ssi_structure_mesh_tying()->FullMasterSideMap();

    // initialize new Jacobian
    CORE::LINALG::SparseMatrix jac_new(*GState().dof_row_map(), 81, true, true);

    // assemble interior rows and columns of original Jacobian into new Jacobian
    CORE::LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *map_structure_interior,
        *map_structure_interior, 1.0, nullptr, nullptr, jac_new, true, true);

    // assemble interior rows and master-side columns of original Jacobian into new Jacobian
    CORE::LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *map_structure_interior,
        *cond_master_dof_map, 1.0, nullptr, nullptr, jac_new, true, true);

    // assemble master-side rows and interior columns of original Jacobian into new Jacobian
    CORE::LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *cond_master_dof_map,
        *map_structure_interior, 1.0, nullptr, nullptr, jac_new, true, true);

    // assemble master-side rows and columns of original Jacobian into new Jacobian
    CORE::LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *cond_master_dof_map,
        *cond_master_dof_map, 1.0, nullptr, nullptr, jac_new, true, true);

    for (const auto& meshtying : ssi_part_->ssi_structure_mesh_tying()->MeshTyingHandlers())
    {
      auto cond_slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
      auto converter = meshtying->SlaveSideConverter();

      // transform and assemble slave-side rows of original Jacobian into new Jacobian (interior
      // columns)
      CORE::LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *cond_slave_dof_map,
          *map_structure_interior, 1.0, &(*converter), nullptr, jac_new, true, true);

      // transform and assemble slave-side rows of original Jacobian into new Jacobian (master-side
      // columns)
      CORE::LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *cond_slave_dof_map,
          *cond_master_dof_map, 1.0, &(*converter), nullptr, jac_new, true, true);

      // transform and assemble slave-side columns of original Jacobian into new Jacobian (interior
      // rows)
      CORE::LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *map_structure_interior,
          *cond_slave_dof_map, 1.0, nullptr, &(*converter), jac_new, true, true);

      // transform and assemble slave-side columns of original Jacobian into new Jacobian
      // (master-side rows)
      CORE::LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *cond_master_dof_map,
          *cond_slave_dof_map, 1.0, nullptr, &(*converter), jac_new, true, true);

      for (const auto& meshtying2 : ssi_part_->ssi_structure_mesh_tying()->MeshTyingHandlers())
      {
        auto cond_slave_dof_map2 = meshtying2->SlaveMasterCoupling()->SlaveDofMap();
        auto converter2 = meshtying2->SlaveSideConverter();

        // assemble derivatives of surface slave dofs w.r.t. line slave dofs (block l)
        CORE::LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *cond_slave_dof_map,
            *cond_slave_dof_map2, 1.0, &(*converter), &(*converter2), jac_new, true, true);
      }
    }

    auto slavemaps = ssi_part_->ssi_structure_mesh_tying()->FullSlaveSideMap();

    // subject slave-side rows of new Jacobian to pseudo Dirichlet conditions to finalize
    // structural meshtying
    jac_new.Complete();
    jac_new.ApplyDirichlet(*slavemaps);

    // replace old Jacobian by new one
    jac_sparse.Assign(CORE::LINALG::View, jac_new);
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::PartitionedSSI::run_pre_compute_x(
    const Epetra_Vector& xold, Epetra_Vector& dir_mutable, const NOX::NLN::Group& curr_grp)
{
  // perform structural meshtying
  if (ssi_part_->ssi_interface_meshtying())
  {
    for (const auto& meshtying : ssi_part_->ssi_structure_mesh_tying()->MeshTyingHandlers())
    {
      auto coupling_map_extractor = meshtying->slave_master_extractor();

      // transform and assemble master-side part of structural increment vector to slave side
      coupling_map_extractor->InsertVector(
          *meshtying->SlaveMasterCoupling()->MasterToSlave(
              coupling_map_extractor->ExtractVector(dir_mutable, 2)),
          1, dir_mutable);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::PartitionedSSI::Setup()
{
  check_init();

  STR::MODELEVALUATOR::BaseSSI::Setup();

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::PartitionedSSI::assemble_force(
    Epetra_Vector& f, const double& timefac_np) const
{
  // perform structural meshtying
  if (ssi_part_->ssi_interface_meshtying() and ssi_part_->is_setup())
  {
    for (const auto& meshtying : ssi_part_->ssi_structure_mesh_tying()->MeshTyingHandlers())
    {
      auto coupling_map_extractor = meshtying->slave_master_extractor();
      // transform and assemble slave-side part of structural right-hand side vector to master side
      coupling_map_extractor->AddVector(*meshtying->SlaveMasterCoupling()->SlaveToMaster(
                                            coupling_map_extractor->ExtractVector(f, 1)),
          2, f);

      // zero out slave-side part of structural right-hand side vector
      coupling_map_extractor->PutScalar(f, 1, 0.0);
    }
  }

  return true;
}

FOUR_C_NAMESPACE_CLOSE
