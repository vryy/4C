/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble pair based contributions into global matrices. The pairs in this class can
not be directly assembled into the global matrices. They have to be assembled into the global
coupling matrices M and D first.


\level 3

*/


#include "fbi_partitioned_penaltycoupling_assembly_manager_indirect.H"
#include "fbi_beam_to_fluid_mortar_manager.H"
#include "fbi_beam_to_fluid_meshtying_params.H"
#include "fbi_calc_utils.H"

#include "beaminteraction_contact_pair.H"
#include "beaminteraction_calc_utils.H"

#include "lib_discret.H"
#include "lib_element.H"
#include "linalg_serialdensematrix.H"
#include "linalg_serialdensevector.H"
#include "linalg_sparsematrix.H"
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

/**
 *
 */
BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManagerIndirect::
    PartitionedBeamInteractionAssemblyManagerIndirect(
        std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>& assembly_contact_elepairs,
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
        Teuchos::RCP<Epetra_FEVector>& fb, Teuchos::RCP<LINALG::SparseOperator> cff,
        Teuchos::RCP<LINALG::SparseMatrix>& cbb, Teuchos::RCP<LINALG::SparseMatrix>& cfb,
        Teuchos::RCP<LINALG::SparseMatrix>& cbf, Teuchos::RCP<const Epetra_Vector> fluid_vel,
        Teuchos::RCP<const Epetra_Vector> beam_vel)
{
  Teuchos::RCP<Teuchos::Time> t =
      Teuchos::TimeMonitor::getNewTimer("FBI::PartitionedAssemblyManagerIndirect");
  Teuchos::TimeMonitor monitor(*t);

  for (auto& elepairptr : assembly_contact_elepairs_)
  {
    // PreEvaluate the pair
    elepairptr->PreEvaluate();
  }
  // Evaluate the global mortar matrices.
  mortar_manager_->EvaluateGlobalDM(assembly_contact_elepairs_);

  // Add the global mortar matrices to the force vector and stiffness matrix.
  mortar_manager_->AddGlobalForceStiffnessContributions(ff, fb, cbb, cbf,
      Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(cff, true), cfb, beam_vel, fluid_vel);
}
