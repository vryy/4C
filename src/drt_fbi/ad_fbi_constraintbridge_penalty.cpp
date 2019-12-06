/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation connecting the penalty constraint enforcement technique with a discretization
approach for Fluid-beam interaction.

\level 2

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/
#include "ad_fbi_constraintbridge_penalty.H"
#include "partitioned_penaltycoupling_assembly_manager_direct.H"
#include "beam_to_fluid_meshtying_params.H"
#include <Epetra_FEVector.h>
#include "../linalg/linalg_sparsematrix.H"

void ADAPTER::FBIConstraintBridgePenalty::Setup(
    const Epetra_Map* beam_map, const Epetra_Map* fluid_map)
{
  // Initialize all necessary vectors and matrices
  FBIConstraintBridge::Setup(beam_map, fluid_map);
  fs_ = Teuchos::rcp(new Epetra_FEVector(*beam_map));
  ff_ = Teuchos::rcp(new Epetra_FEVector(*fluid_map));
  Cff_ = Teuchos::rcp(new LINALG::SparseMatrix(*fluid_map, 30, true, true,
      LINALG::SparseMatrix::FE_MATRIX));  // todo Is there a better estimator?
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIConstraintBridgePenalty::Evaluate(
    const DRT::Discretization& discretization1, const DRT::Discretization& discretization2)
{
  // Create assembly manager.. todo this will have to change as soon as we add mortar pairs.. hand
  // in by dependency injection?
  BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManagerDirect
      assembly_manager =
          BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManagerDirect(
              *(GetPairs()));
  // compute and assembly the coupling matrices and vectors
  assembly_manager.EvaluateForceStiff(
      discretization1, discretization2, ff_, fs_, Cff_, Css_, Csf_, Cfs_);
  Cff_->Complete();

  // Unset the dirichlet flag in case we were doing a fluid solve
  UnsetWeakDirichletFlag();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIConstraintBridgePenalty::ResetBridge()
{
  ADAPTER::FBIConstraintBridge::ResetBridge();
  fs_->PutScalar(0.0);
  Cff_->Reset();
  ff_->PutScalar(0.0);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIConstraintBridgePenalty::SetWeakDirichletFlag()
{
  beam_interaction_params_->SetWeakDirichletFlag();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIConstraintBridgePenalty::UnsetWeakDirichletFlag()
{
  beam_interaction_params_->UnsetWeakDirichletFlag();
}
