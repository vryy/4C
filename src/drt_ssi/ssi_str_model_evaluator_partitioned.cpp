/*-----------------------------------------------------------*/
/*! \file
\brief Model evaluator for structure part of partitioned ssi.

\level 3


*/
/*-----------------------------------------------------------*/


#include "ssi_str_model_evaluator_partitioned.H"

#include "ssi_utils.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_ssi/ssi_partitioned.H"

#include "../drt_structure_new/str_dbc.H"
#include "../drt_structure_new/str_impl_generic.H"
#include "../drt_structure_new/str_timint_implicit.H"

#include "../solver_nonlin_nox/nox_nln_group.H"

#include "../linalg/linalg_utils_sparse_algebra_math.H"

#include "../linalg/linalg_matrixtransform.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 01/18 |
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::PartitionedSSI::PartitionedSSI(const Teuchos::RCP<const SSI::SSIPart>
        ssi_part  //!< partitioned algorithm for scalar-structure interaction
    )
    : ssi_part_(ssi_part)
{
}

/*----------------------------------------------------------------------*
 | assemble Jacobian                                         fang 01/18 |
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::PartitionedSSI::AssembleJacobian(
    LINALG::SparseOperator& jac, const double& timefac_np) const
{
  // perform structural meshtying
  if (ssi_part_->SSIInterfaceMeshtying())
  {
    // cast old Jacobian
    auto& jac_sparse = dynamic_cast<LINALG::SparseMatrix&>(jac);

    auto map_structure_interior = ssi_part_->SSIStructureMeshTying()->InteriorMap();
    auto cond_master_dof_map = ssi_part_->SSIStructureMeshTying()->FullMasterSideMap();

    // initialize new Jacobian
    LINALG::SparseMatrix jac_new(*GState().DofRowMap(), 81, true, true);

    // assemble interior rows and columns of original Jacobian into new Jacobian
    LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *map_structure_interior,
        *map_structure_interior, 1.0, nullptr, nullptr, jac_new, true, true);

    // assemble interior rows and master-side columns of original Jacobian into new Jacobian
    LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *map_structure_interior,
        *cond_master_dof_map, 1.0, nullptr, nullptr, jac_new, true, true);

    // assemble master-side rows and interior columns of original Jacobian into new Jacobian
    LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *cond_master_dof_map,
        *map_structure_interior, 1.0, nullptr, nullptr, jac_new, true, true);

    // assemble master-side rows and columns of original Jacobian into new Jacobian
    LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *cond_master_dof_map, *cond_master_dof_map,
        1.0, nullptr, nullptr, jac_new, true, true);

    for (const auto& meshtying : ssi_part_->SSIStructureMeshTying()->MeshtyingHandlers())
    {
      auto cond_slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
      auto converter = meshtying->SlaveSideConverter();

      // transform and assemble slave-side rows of original Jacobian into new Jacobian (interior
      // columns)
      LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *cond_slave_dof_map,
          *map_structure_interior, 1.0, &(*converter), nullptr, jac_new, true, true);

      // transform and assemble slave-side rows of original Jacobian into new Jacobian (master-side
      // columns)
      LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *cond_slave_dof_map,
          *cond_master_dof_map, 1.0, &(*converter), nullptr, jac_new, true, true);

      // transform and assemble slave-side columns of original Jacobian into new Jacobian (interior
      // rows)
      LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *map_structure_interior,
          *cond_slave_dof_map, 1.0, nullptr, &(*converter), jac_new, true, true);

      // transform and assemble slave-side columns of original Jacobian into new Jacobian
      // (master-side rows)
      LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *cond_master_dof_map,
          *cond_slave_dof_map, 1.0, nullptr, &(*converter), jac_new, true, true);

      for (const auto& meshtying2 : ssi_part_->SSIStructureMeshTying()->MeshtyingHandlers())
      {
        auto cond_slave_dof_map2 = meshtying2->SlaveMasterCoupling()->SlaveDofMap();
        auto converter2 = meshtying2->SlaveSideConverter();

        // assemble derivatives of surface slave dofs w.r.t. line slave dofs (block l)
        LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *cond_slave_dof_map,
            *cond_slave_dof_map2, 1.0, &(*converter), &(*converter2), jac_new, true, true);
      }
    }

    auto slavemaps = ssi_part_->SSIStructureMeshTying()->FullSlaveSideMap();

    // subject slave-side rows of new Jacobian to pseudo Dirichlet conditions to finalize
    // structural meshtying
    jac_new.Complete();
    jac_new.ApplyDirichlet(*slavemaps);

    // replace old Jacobian by new one
    jac_sparse.Assign(LINALG::View, jac_new);
  }

  return true;
}

/*----------------------------------------------------------------------*
 | pre-compute solution vector                               fang 01/18 |
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::PartitionedSSI::RunPreComputeX(
    const Epetra_Vector& xold, Epetra_Vector& dir_mutable, const NOX::NLN::Group& curr_grp)
{
  // perform structural meshtying
  if (ssi_part_->SSIInterfaceMeshtying())
  {
    for (const auto& meshtying : ssi_part_->SSIStructureMeshTying()->MeshtyingHandlers())
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
  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::PartitionedSSI::GetBlockDofRowMapPtr() const
{
  CheckInitSetup();
  return GState().DofRowMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::PartitionedSSI::GetCurrentSolutionPtr() const
{
  CheckInit();
  return GState().GetDisNp();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::PartitionedSSI::GetLastTimeStepSolutionPtr()
    const
{
  CheckInit();
  return GState().GetDisN();
}

/*----------------------------------------------------------------------*
 | assemble right-hand side vector                           fang 01/18 |
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::PartitionedSSI::AssembleForce(
    Epetra_Vector& f, const double& timefac_np) const
{
  // perform structural meshtying
  if (ssi_part_->SSIInterfaceMeshtying() and ssi_part_->IsSetup())
  {
    for (const auto& meshtying : ssi_part_->SSIStructureMeshTying()->MeshtyingHandlers())
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

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::PartitionedSSI::UpdateStepState(const double& timefac_n) {}
