/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble pair based contributions into global matrices. The pairs in this class can
not be directly assembled into the global matrices. They have to be assembled into the global
coupling matrices M and D first.


\level 3

*/


#include "beaminteraction_submodel_evaluator_beamcontact_assembly_manager_indirect.H"

#include "beaminteraction_contact_pair.H"
#include "beaminteraction_calc_utils.H"
#include "beaminteraction_beam_to_solid_mortar_manager.H"
#include "beaminteraction_str_model_evaluator_datastate.H"

#include "lib_discret.H"
#include "lib_element.H"
#include "linalg_serialdensematrix.H"
#include "linalg_serialdensevector.H"


/**
 *
 */
BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerInDirect::
    BeamContactAssemblyManagerInDirect(
        const std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>&
            assembly_contact_elepairs,
        const Teuchos::RCP<const ::DRT::Discretization>& discret,
        const Teuchos::RCP<const BEAMINTERACTION::BeamToSolidParamsBase>& beam_to_solid_params)
    : BeamContactAssemblyManager()
{
  // Create the mortar manager. We add 1 to the MaxAllGID since this gives the maximum GID and NOT
  // the length of the GIDs.
  mortar_manager_ = Teuchos::rcp<BEAMINTERACTION::BeamToSolidMortarManager>(
      new BEAMINTERACTION::BeamToSolidMortarManager(
          discret, beam_to_solid_params, discret->DofRowMap()->MaxAllGID() + 1));

  // Setup the mortar manager.
  mortar_manager_->Setup();
  mortar_manager_->SetLocalMaps(assembly_contact_elepairs);
}


/**
 *
 */
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerInDirect::EvaluateForceStiff(
    Teuchos::RCP<::DRT::Discretization> discret,
    const Teuchos::RCP<const STR::MODELEVALUATOR::BeamInteractionDataState>& data_state,
    Teuchos::RCP<Epetra_FEVector> fe_sysvec, Teuchos::RCP<LINALG::SparseMatrix> fe_sysmat)
{
  // Evaluate the global mortar matrices.
  mortar_manager_->EvaluateGlobalCouplingContributions(data_state->GetDisColNp());

  // Add the global mortar matrices to the force vector and stiffness matrix.
  mortar_manager_->AddGlobalForceStiffnessPenaltyContributions(data_state, fe_sysmat, fe_sysvec);
}


double BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerInDirect::GetEnergy(
    const Teuchos::RCP<const Epetra_Vector>& disp) const
{
  return mortar_manager_->GetEnergy();
}
