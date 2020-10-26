/*-----------------------------------------------------------*/
/*! \file
\brief Model evaluator for structure part of partitioned ssi.

\level 3


*/
/*-----------------------------------------------------------*/


#include "ssi_str_model_evaluator_partitioned.H"

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

    // initialize new Jacobian
    LINALG::SparseMatrix jac_new(*GState().DofRowMap(), 81, true, true);

    // assemble interior and master-side rows and columns of original Jacobian into new Jacobian
    LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *ssi_part_->MapStructureCondensed(),
        *ssi_part_->MapStructureCondensed(), 1., nullptr, nullptr, jac_new);

    // transform and assemble slave-side rows of original Jacobian into new Jacobian
    LINALG::MatrixLogicalSplitAndTransform()(jac_sparse,
        *ssi_part_->InterfaceCouplingAdapterStructure()->SlaveDofMap(),
        *ssi_part_->MapStructureCondensed(), 1.0,
        &ssi_part_->InterfaceCouplingAdapterStructureSlaveConverter(), nullptr, jac_new, true,
        true);

    // transform and assemble slave-side columns of original Jacobian into new Jacobian
    LINALG::MatrixLogicalSplitAndTransform()(jac_sparse, *ssi_part_->MapStructureCondensed(),
        *ssi_part_->InterfaceCouplingAdapterStructure()->SlaveDofMap(), 1.0, nullptr,
        &ssi_part_->InterfaceCouplingAdapterStructureSlaveConverter(), jac_new, true, true);

    // transform and assemble slave-side rows and columns of original Jacobian into new Jacobian
    LINALG::MatrixLogicalSplitAndTransform()(jac_sparse,
        *ssi_part_->InterfaceCouplingAdapterStructure()->SlaveDofMap(),
        *ssi_part_->InterfaceCouplingAdapterStructure()->SlaveDofMap(), 1.0,
        &ssi_part_->InterfaceCouplingAdapterStructureSlaveConverter(),
        &ssi_part_->InterfaceCouplingAdapterStructureSlaveConverter(), jac_new, true, true);

    // subject slave-side rows of new Jacobian to pseudo Dirichlet conditions to finalize structural
    // meshtying
    jac_new.Complete();
    jac_new.ApplyDirichlet(*ssi_part_->InterfaceCouplingAdapterStructure()->SlaveDofMap());

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
    // transform and assemble master-side part of structural increment vector to slave side
    ssi_part_->MapsStructure()->InsertVector(
        *ssi_part_->InterfaceCouplingAdapterStructure()->MasterToSlave(
            ssi_part_->MapsStructure()->ExtractVector(dir_mutable, 2)),
        1, dir_mutable);
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
    // transform and assemble slave-side part of structural right-hand side vector to master side
    ssi_part_->MapsStructure()->AddVector(
        *ssi_part_->InterfaceCouplingAdapterStructure()->SlaveToMaster(
            ssi_part_->MapsStructure()->ExtractVector(f, 1)),
        2, f);

    // zero out slave-side part of structural right-hand side vector
    ssi_part_->MapsStructure()->PutScalar(f, 1, 0.);
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::PartitionedSSI::UpdateStepState(const double& timefac_n) {}
