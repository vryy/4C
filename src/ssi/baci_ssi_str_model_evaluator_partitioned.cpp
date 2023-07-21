/*-----------------------------------------------------------*/
/*! \file
\brief Model evaluator for structure part of partitioned ssi.

\level 3


*/
/*-----------------------------------------------------------*/


#include "baci_ssi_str_model_evaluator_partitioned.H"

#include "baci_ssi_utils.H"

#include "baci_adapter_str_ssiwrapper.H"
#include "baci_coupling_adapter.H"
#include "baci_coupling_adapter_converter.H"
#include "baci_adapter_scatra_base_algorithm.H"

#include "baci_scatra_timint_implicit.H"

#include "baci_ssi_partitioned.H"

#include "baci_structure_new_dbc.H"
#include "baci_structure_new_impl_generic.H"
#include "baci_structure_new_timint_implicit.H"

#include "baci_solver_nonlin_nox_group.H"

#include "baci_linalg_utils_sparse_algebra_math.H"

#include "baci_linalg_matrixtransform.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::PartitionedSSI::PartitionedSSI(const Teuchos::RCP<const SSI::SSIPart> ssi_part)
    : ssi_part_(ssi_part)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::PartitionedSSI::AssembleJacobian(
    CORE::LINALG::SparseOperator& jac, const double& timefac_np) const
{
  // perform structural meshtying
  if (ssi_part_->SSIInterfaceMeshtying())
  {
    // cast old Jacobian
    auto& jac_sparse = dynamic_cast<CORE::LINALG::SparseMatrix&>(jac);

    auto map_structure_interior = ssi_part_->SSIStructureMeshTying()->InteriorMap();
    auto cond_master_dof_map = ssi_part_->SSIStructureMeshTying()->FullMasterSideMap();

    // initialize new Jacobian
    CORE::LINALG::SparseMatrix jac_new(*GState().DofRowMap(), 81, true, true);

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

    for (const auto& meshtying : ssi_part_->SSIStructureMeshTying()->MeshTyingHandlers())
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

      for (const auto& meshtying2 : ssi_part_->SSIStructureMeshTying()->MeshTyingHandlers())
      {
        auto cond_slave_dof_map2 = meshtying2->SlaveMasterCoupling()->SlaveDofMap();
        auto converter2 = meshtying2->SlaveSideConverter();

        // assemble derivatives of surface slave dofs w.r.t. line slave dofs (block l)
        CORE::LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *cond_slave_dof_map,
            *cond_slave_dof_map2, 1.0, &(*converter), &(*converter2), jac_new, true, true);
      }
    }

    auto slavemaps = ssi_part_->SSIStructureMeshTying()->FullSlaveSideMap();

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
void STR::MODELEVALUATOR::PartitionedSSI::RunPreComputeX(
    const Epetra_Vector& xold, Epetra_Vector& dir_mutable, const NOX::NLN::Group& curr_grp)
{
  // perform structural meshtying
  if (ssi_part_->SSIInterfaceMeshtying())
  {
    for (const auto& meshtying : ssi_part_->SSIStructureMeshTying()->MeshTyingHandlers())
    {
      auto coupling_map_extractor = meshtying->SlaveMasterExtractor();

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
  CheckInit();

  STR::MODELEVALUATOR::BaseSSI::Setup();

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::PartitionedSSI::AssembleForce(
    Epetra_Vector& f, const double& timefac_np) const
{
  // perform structural meshtying
  if (ssi_part_->SSIInterfaceMeshtying() and ssi_part_->IsSetup())
  {
    for (const auto& meshtying : ssi_part_->SSIStructureMeshTying()->MeshTyingHandlers())
    {
      auto coupling_map_extractor = meshtying->SlaveMasterExtractor();
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
