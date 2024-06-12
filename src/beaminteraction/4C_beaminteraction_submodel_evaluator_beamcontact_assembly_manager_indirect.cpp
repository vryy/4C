/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble pair based contributions into global matrices. The pairs in this class can
not be directly assembled into the global matrices. They have to be assembled into the global
coupling matrices M and D first.


\level 3

*/


#include "4C_beaminteraction_submodel_evaluator_beamcontact_assembly_manager_indirect.hpp"

#include "4C_beaminteraction_beam_to_solid_mortar_manager.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_beaminteraction_str_model_evaluator_datastate.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
void BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerInDirect::evaluate_force_stiff(
    Teuchos::RCP<Core::FE::Discretization> discret,
    const Teuchos::RCP<const STR::MODELEVALUATOR::BeamInteractionDataState>& data_state,
    Teuchos::RCP<Epetra_FEVector> fe_sysvec, Teuchos::RCP<Core::LinAlg::SparseMatrix> fe_sysmat)
{
  // Evaluate the global coupling terms
  mortar_manager_->evaluate_global_coupling_contributions(data_state->GetDisColNp());

  // Add the global mortar matrices to the force vector and stiffness matrix
  mortar_manager_->add_global_force_stiffness_penalty_contributions(
      data_state, fe_sysmat, fe_sysvec);
}


double BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerInDirect::get_energy(
    const Teuchos::RCP<const Epetra_Vector>& disp) const
{
  return mortar_manager_->get_energy();
}

FOUR_C_NAMESPACE_CLOSE
