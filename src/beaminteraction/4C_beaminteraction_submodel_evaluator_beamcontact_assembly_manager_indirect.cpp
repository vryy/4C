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
  mortar_manager_->evaluate_force_stiff_penalty_regularization(data_state, fe_sysmat, fe_sysvec);
}


double BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerInDirect::get_energy(
    const Teuchos::RCP<const Epetra_Vector>& disp) const
{
  const double global_mortar_energy = mortar_manager_->get_energy();

  // The value we returned here is summed up over all processors. Since we already have the global
  // energy here, we only return it on rank 0.
  if (disp->Comm().MyPID() == 0)
    return global_mortar_energy;
  else
    return 0.0;
}

FOUR_C_NAMESPACE_CLOSE
