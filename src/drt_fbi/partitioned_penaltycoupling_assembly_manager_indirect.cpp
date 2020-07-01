/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble pair based contributions into global matrices. The pairs in this class can
not be directly assembled into the global matrices. They have to be assembled into the global
coupling matrices M and D first.


\level 3

*/


#include "partitioned_penaltycoupling_assembly_manager_indirect.H"
#include "beam_to_fluid_mortar_manager.H"
#include "beam_to_fluid_meshtying_params.H"
#include "../drt_fbi/fbi_calc_utils.H"

#include "../drt_beaminteraction/beam_contact_pair.H"
#include "../drt_beaminteraction/beaminteraction_calc_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"


/**
 *
 */
BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManagerIndirect::
    PartitionedBeamInteractionAssemblyManagerIndirect(
        std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> assembly_contact_elepairs,
        Teuchos::RCP<const DRT::Discretization>& discretization1,
        Teuchos::RCP<const DRT::Discretization>& discretization2,
        Teuchos::RCP<FBI::BeamToFluidMeshtyingParams> beam_contact_params_ptr)
    : PartitionedBeamInteractionAssemblyManager(assembly_contact_elepairs)
{
  // Create the mortar manager.
  mortar_manager_ = Teuchos::rcp<BEAMINTERACTION::BeamToFluidMortarManager>(
      new BEAMINTERACTION::BeamToFluidMortarManager(discretization1, discretization2,
          beam_contact_params_ptr, discretization1->DofRowMap()->MaxAllGID()));

  // Setup the mortar manager.
  mortar_manager_->Setup();
  mortar_manager_->SetLocalMaps(assembly_contact_elepairs_);
}


/**
 *
 */
void BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManagerIndirect::
    EvaluateForceStiff(const DRT::Discretization& discretization1,
        const DRT::Discretization& discretization2, Teuchos::RCP<Epetra_FEVector>& ff,
        Teuchos::RCP<Epetra_FEVector>& fb, Teuchos::RCP<LINALG::SparseMatrix>& cff,
        Teuchos::RCP<LINALG::SparseMatrix>& cbb, Teuchos::RCP<LINALG::SparseMatrix>& cfb,
        Teuchos::RCP<LINALG::SparseMatrix>& cbf, Teuchos::RCP<const Epetra_Vector> fluid_vel,
        Teuchos::RCP<const Epetra_Vector> beam_vel)
{
  for (auto& elepairptr : assembly_contact_elepairs_)
  {
    // PreEvaluate the pair
    elepairptr->PreEvaluate();
  }
  // Evaluate the global mortar matrices.
  mortar_manager_->EvaluateGlobalDM(assembly_contact_elepairs_);

  // Add the global mortar matrices to the force vector and stiffness matrix.
  mortar_manager_->AddGlobalForceStiffnessContributions(
      ff, fb, cbb, cbf, cff, cfb, beam_vel, fluid_vel);
}
