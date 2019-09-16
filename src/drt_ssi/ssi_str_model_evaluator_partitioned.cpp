/*-----------------------------------------------------------*/
/*! \file
\brief Model evaluator for structure part of partitioned ssi.

\level 3

\maintainer Christoph Schmidt

*/
/*-----------------------------------------------------------*/


#include "ssi_str_model_evaluator_partitioned.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_fsi/fsi_matrixtransform.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_ssi/ssi_partitioned.H"

#include "../drt_structure_new/str_dbc.H"
#include "../drt_structure_new/str_impl_generic.H"
#include "../drt_structure_new/str_timint_implicit.H"
#include "../drt_structure_new/str_nln_solver_generic.H"
#include "../drt_structure_new/str_timint_basedataglobalstate.H"

#include "../solver_nonlin_nox/nox_nln_group.H"

#include "../linalg/linalg_utils.H"

#include "Epetra_Comm.h"

/*----------------------------------------------------------------------*
 | constructor                                               fang 01/18 |
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::PartitionedSSI::PartitionedSSI(const Teuchos::RCP<const SSI::SSI_Part>
        ssi_part  //!< partitioned algorithm for scalar-structure interaction
    )
    : ssi_part_(ssi_part)
{
  return;
}

/*----------------------------------------------------------------------*
 | assemble Jacobian                                         fang 01/18 |
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::PartitionedSSI::AssembleJacobian(
    LINALG::SparseOperator& jac, const double& timefac_np) const
{
  // perform structural meshtying for scatra-scatra interface coupling
  if (ssi_part_->ScaTraField()->ScaTraField()->S2ICoupling())
  {
    // cast old Jacobian
    LINALG::SparseMatrix& jac_sparse = dynamic_cast<LINALG::SparseMatrix&>(jac);

    // initialize new Jacobian
    LINALG::SparseMatrix jac_new(*GState().DofRowMap(), 81, true, true);

    // assemble interior and master-side rows and columns of original Jacobian into new Jacobian
    FSI::UTILS::MatrixLogicalSplitAndTransform()(jac_sparse, *ssi_part_->MapStructureCondensed(),
        *ssi_part_->MapStructureCondensed(), 1., NULL, NULL, jac_new);

    // transform and assemble slave-side rows of original Jacobian into new Jacobian
    ADAPTER::CouplingSlaveConverter converter(*ssi_part_->CouplingAdapterStructure());
    FSI::UTILS::MatrixLogicalSplitAndTransform()(jac_sparse,
        *ssi_part_->CouplingAdapterStructure()->SlaveDofMap(), *ssi_part_->MapStructureCondensed(),
        1., &converter, NULL, jac_new, true, true);

    // transform and assemble slave-side columns of original Jacobian into new Jacobian
    FSI::UTILS::MatrixLogicalSplitAndTransform()(jac_sparse, *ssi_part_->MapStructureCondensed(),
        *ssi_part_->CouplingAdapterStructure()->SlaveDofMap(), 1., NULL, &converter, jac_new, true,
        true);

    // transform and assemble slave-side rows and columns of original Jacobian into new Jacobian
    FSI::UTILS::MatrixLogicalSplitAndTransform()(jac_sparse,
        *ssi_part_->CouplingAdapterStructure()->SlaveDofMap(),
        *ssi_part_->CouplingAdapterStructure()->SlaveDofMap(), 1., &converter, &converter, jac_new,
        true, true);

    // subject slave-side rows of new Jacobian to pseudo Dirichlet conditions to finalize structural
    // meshtying
    jac_new.Complete();
    jac_new.ApplyDirichlet(*ssi_part_->CouplingAdapterStructure()->SlaveDofMap());

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
  // perform structural meshtying for scatra-scatra interface coupling
  if (ssi_part_->ScaTraField()->ScaTraField()->S2ICoupling())
    // transform and assemble master-side part of structural increment vector to slave side
    ssi_part_->MapsStructure()->InsertVector(
        *ssi_part_->CouplingAdapterStructure()->MasterToSlave(
            ssi_part_->MapsStructure()->ExtractVector(dir_mutable, 2)),
        1, dir_mutable);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::PartitionedSSI::Setup()
{
  CheckInit();
  // set flag
  issetup_ = true;

  return;
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
  // perform structural meshtying for scatra-scatra interface coupling
  if (ssi_part_->ScaTraField()->ScaTraField()->S2ICoupling() and ssi_part_->IsSetup())
  {
    // transform and assemble slave-side part of structural right-hand side vector to master side
    ssi_part_->MapsStructure()->AddVector(*ssi_part_->CouplingAdapterStructure()->SlaveToMaster(
                                              ssi_part_->MapsStructure()->ExtractVector(f, 1)),
        2, f);

    // zero out slave-side part of structural right-hand side vector
    ssi_part_->MapsStructure()->PutScalar(f, 1, 0.);
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::PartitionedSSI::UpdateStepState(const double& timefac_n) { return; }
