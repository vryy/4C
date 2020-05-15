/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation connecting the penalty constraint enforcement technique with a discretization
approach for Fluid-beam interaction.

\level 2

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/
#include "ad_fbi_constraintbridge_penalty.H"
#include "partitioned_penaltycoupling_assembly_manager_direct.H"
#include "partitioned_penaltycoupling_assembly_manager_indirect.H"
#include "beam_to_fluid_assembly_manager_factory.H"
#include "beam_to_fluid_meshtying_params.H"
#include "../drt_beaminteraction/beam_contact_pair.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_sparsematrix.H"
#include <Epetra_FEVector.h>

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
    Teuchos::RCP<const DRT::Discretization> discretization1,
    Teuchos::RCP<const DRT::Discretization> discretization2,
    Teuchos::RCP<const Epetra_Vector> fluid_vel, Teuchos::RCP<const Epetra_Vector> beam_vel)
{
  // Create assembly manager.. todo this will have to change as soon as we add mortar pairs.. hand
  // in by dependency injection?
  Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManager>
      assembly_manager = BEAMINTERACTION::BeamToFluidAssemblyManagerFactory::CreateAssemblyManager(
          discretization1, discretization2, *(GetPairs()), GetParams());
  // compute and assembly the coupling matrices and vectors
  assembly_manager->EvaluateForceStiff(
      *discretization1, *discretization2, ff_, fs_, Cff_, Css_, Csf_, Cfs_, fluid_vel, beam_vel);
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

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIConstraintBridgePenalty::ScalePenaltyStructureContributions()
{
  if (fs_->Scale(GetParams()->GetPenaltyParameter()))
    dserror("Scaling of the penalty force was unsuccessful!\n");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIConstraintBridgePenalty::ScalePenaltyFluidContributions()
{
  if (Cff_->Scale(GetParams()->GetPenaltyParameter()) ||
      ff_->Scale(GetParams()->GetPenaltyParameter()))
    dserror("Scaling of the penalty force was unsuccessful!\n");
}
