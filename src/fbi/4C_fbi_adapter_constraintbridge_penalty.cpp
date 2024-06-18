/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation connecting the penalty constraint enforcement technique with a discretization
approach for Fluid-beam interaction.

\level 2

*----------------------------------------------------------------------*/
#include "4C_fbi_adapter_constraintbridge_penalty.hpp"

#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_fbi_beam_to_fluid_assembly_manager_factory.hpp"
#include "4C_fbi_beam_to_fluid_meshtying_params.hpp"
#include "4C_fbi_partitioned_penaltycoupling_assembly_manager_direct.hpp"
#include "4C_fbi_partitioned_penaltycoupling_assembly_manager_indirect.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_sparseoperator.hpp"

#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN

void Adapter::FBIConstraintBridgePenalty::setup(const Epetra_Map* beam_map,
    const Epetra_Map* fluid_map, Teuchos::RCP<Core::LinAlg::SparseOperator> fluidmatrix,
    bool fluidmeshtying)
{
  // Initialize all necessary vectors and matrices
  FBIConstraintBridge::setup(beam_map, fluid_map, fluidmatrix, fluidmeshtying);
  fs_ = Teuchos::rcp(new Epetra_FEVector(*beam_map));
  ff_ = Teuchos::rcp(new Epetra_FEVector(*fluid_map));
  cff_ = fluidmatrix;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIConstraintBridgePenalty::evaluate(
    Teuchos::RCP<const Core::FE::Discretization> discretization1,
    Teuchos::RCP<const Core::FE::Discretization> discretization2,
    Teuchos::RCP<const Epetra_Vector> fluid_vel, Teuchos::RCP<const Epetra_Vector> beam_vel)
{
  // Create assembly manager..
  Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManager>
      assembly_manager =
          BEAMINTERACTION::BeamToFluidAssemblyManagerFactory::create_assembly_manager(
              discretization1, discretization2, *(GetPairs()), GetParams(), assemblystrategy_);
  // compute and assembly the coupling matrices and vectors
  assembly_manager->evaluate_force_stiff(
      *discretization1, *discretization2, ff_, fs_, cff_, css_, csf_, cfs_, fluid_vel, beam_vel);
  cff_->Complete();

  // Unset the dirichlet flag in case we were doing a fluid solve
  unset_weak_dirichlet_flag();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIConstraintBridgePenalty::ResetBridge()
{
  fs_->PutScalar(0.0);
  cff_->Reset();
  ff_->PutScalar(0.0);
  fluid_scaled_ = false;
  structure_scaled_ = false;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIConstraintBridgePenalty::set_weak_dirichlet_flag()
{
  beam_interaction_params_->set_weak_dirichlet_flag();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIConstraintBridgePenalty::unset_weak_dirichlet_flag()
{
  beam_interaction_params_->unset_weak_dirichlet_flag();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIConstraintBridgePenalty::scale_penalty_structure_contributions()
{
  if (!structure_scaled_)
  {
    if (fs_->Scale(GetParams()->GetPenaltyParameter()))
      FOUR_C_THROW("Scaling of the penalty force was unsuccessful!\n");
    structure_scaled_ = true;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIConstraintBridgePenalty::scale_penalty_fluid_contributions()
{
  if (!fluid_scaled_)
  {
    if (cff_->Scale(GetParams()->GetPenaltyParameter()) ||
        ff_->Scale(GetParams()->GetPenaltyParameter()))
      FOUR_C_THROW("Scaling of the penalty force was unsuccessful!\n");
    fluid_scaled_ = true;
  }
}

FOUR_C_NAMESPACE_CLOSE
