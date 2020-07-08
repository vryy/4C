/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble pair based contributions into global matrices. The pairs in this class can
not be directly assembled into the global matrices. They have to be assembled into the global
coupling matrices M and D first.


\level 3

*/


#include "beaminteraction_submodel_evaluator_beamcontact_assembly_manager_indirect.H"

#include "beam_contact_pair.H"
#include "beaminteraction_calc_utils.H"
#include "beam_to_solid_mortar_manager.H"
#include "str_model_evaluator_beaminteraction_datastate.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"


/**
 *
 */
BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerInDirect::
    BeamContactAssemblyManagerInDirect(
        const std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>&
            assembly_contact_elepairs,
        const Teuchos::RCP<const DRT::Discretization>& discret,
        const Teuchos::RCP<const BEAMINTERACTION::BeamToSolidParamsBase>& beam_to_solid_params)
    : BeamContactAssemblyManager()
{
  // Create the mortar manager.
  mortar_manager_ = Teuchos::rcp<BEAMINTERACTION::BeamToSolidMortarManager>(
      new BEAMINTERACTION::BeamToSolidMortarManager(
          discret, beam_to_solid_params, discret->DofRowMap()->MaxAllGID()));

  // Setup the mortar manager.
  mortar_manager_->Setup();
  mortar_manager_->SetLocalMaps(assembly_contact_elepairs);
}


/**
 *
 */
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerInDirect::EvaluateForceStiff(
    Teuchos::RCP<DRT::Discretization> discret,
    const Teuchos::RCP<const STR::MODELEVALUATOR::BeamInteractionDataState>& data_state,
    Teuchos::RCP<Epetra_FEVector> fe_sysvec, Teuchos::RCP<LINALG::SparseMatrix> fe_sysmat)
{
  // Evaluate the global mortar matrices.
  mortar_manager_->EvaluateGlobalDM();

  // Add the global mortar matrices to the force vector and stiffness matrix.
  mortar_manager_->AddGlobalForceStiffnessPenaltyContributions(data_state, fe_sysmat, fe_sysvec);
}


double BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerInDirect::GetEnergy(
    const Teuchos::RCP<const Epetra_Vector>& disp) const
{
  return mortar_manager_->GetEnergy(disp);
}
