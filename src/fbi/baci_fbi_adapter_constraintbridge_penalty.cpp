/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation connecting the penalty constraint enforcement technique with a discretization
approach for Fluid-beam interaction.

\level 2

*----------------------------------------------------------------------*/
#include "baci_fbi_adapter_constraintbridge_penalty.hpp"

#include "baci_beaminteraction_contact_pair.hpp"
#include "baci_fbi_beam_to_fluid_assembly_manager_factory.hpp"
#include "baci_fbi_beam_to_fluid_meshtying_params.hpp"
#include "baci_fbi_partitioned_penaltycoupling_assembly_manager_direct.hpp"
#include "baci_fbi_partitioned_penaltycoupling_assembly_manager_indirect.hpp"
#include "baci_lib_discret.hpp"
#include "baci_linalg_sparsematrix.hpp"
#include "baci_linalg_sparseoperator.hpp"

#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN

void ADAPTER::FBIConstraintBridgePenalty::Setup(const Epetra_Map* beam_map,
    const Epetra_Map* fluid_map, Teuchos::RCP<CORE::LINALG::SparseOperator> fluidmatrix,
    bool fluidmeshtying)
{
  // Initialize all necessary vectors and matrices
  FBIConstraintBridge::Setup(beam_map, fluid_map, fluidmatrix, fluidmeshtying);
  fs_ = Teuchos::rcp(new Epetra_FEVector(*beam_map));
  ff_ = Teuchos::rcp(new Epetra_FEVector(*fluid_map));
  Cff_ = fluidmatrix;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIConstraintBridgePenalty::Evaluate(
    Teuchos::RCP<const DRT::Discretization> discretization1,
    Teuchos::RCP<const DRT::Discretization> discretization2,
    Teuchos::RCP<const Epetra_Vector> fluid_vel, Teuchos::RCP<const Epetra_Vector> beam_vel)
{
  // Create assembly manager..
  Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManager>
      assembly_manager = BEAMINTERACTION::BeamToFluidAssemblyManagerFactory::CreateAssemblyManager(
          discretization1, discretization2, *(GetPairs()), GetParams(), assemblystrategy_);
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
  fs_->PutScalar(0.0);
  Cff_->Reset();
  ff_->PutScalar(0.0);
  fluid_scaled_ = false;
  structure_scaled_ = false;
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
  if (!structure_scaled_)
  {
    if (fs_->Scale(GetParams()->GetPenaltyParameter()))
      dserror("Scaling of the penalty force was unsuccessful!\n");
    structure_scaled_ = true;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIConstraintBridgePenalty::ScalePenaltyFluidContributions()
{
  if (!fluid_scaled_)
  {
    if (Cff_->Scale(GetParams()->GetPenaltyParameter()) ||
        ff_->Scale(GetParams()->GetPenaltyParameter()))
      dserror("Scaling of the penalty force was unsuccessful!\n");
    fluid_scaled_ = true;
  }
}

FOUR_C_NAMESPACE_CLOSE
