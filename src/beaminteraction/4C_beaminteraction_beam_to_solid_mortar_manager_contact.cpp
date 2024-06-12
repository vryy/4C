/*----------------------------------------------------------------------*/
/*! \file

\brief Manage the creation of additional DOFs for mortar couplings between beams and solids in
contact simulations

\level 3
*/


#include "4C_beaminteraction_beam_to_solid_mortar_manager_contact.hpp"

#include "4C_beaminteraction_beam_to_solid_mortar_manager.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_contact_params.hpp"
#include "4C_beaminteraction_beam_to_solid_utils.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_beaminteraction_str_model_evaluator_datastate.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_utils_fad.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCPDecl.hpp>

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
BEAMINTERACTION::BeamToSolidMortarManagerContact::BeamToSolidMortarManagerContact(
    const Teuchos::RCP<const Core::FE::Discretization>& discret,
    const Teuchos::RCP<const BEAMINTERACTION::BeamToSolidParamsBase>& params,
    int start_value_lambda_gid)
    : BeamToSolidMortarManager(discret, params, start_value_lambda_gid)
{
}

/**
 *
 */
std::tuple<Teuchos::RCP<Epetra_Vector>, Teuchos::RCP<Epetra_Vector>, Teuchos::RCP<Epetra_Vector>>
BEAMINTERACTION::BeamToSolidMortarManagerContact::get_penalty_regularization(
    const bool compute_linearization) const
{
  using fad_type = fad_type_1st_order_2_variables;
  const auto beam_to_solid_conact_params =
      Teuchos::rcp_dynamic_cast<const BEAMINTERACTION::BeamToSolidSurfaceContactParams>(
          beam_to_solid_params_, true);

  // Get the penalty regularized Lagrange multipliers and the derivative w.r.t. the constraint
  // vector (averaged gap) and the scaling vector (kappa)
  auto create_lambda_row_vector_with_zeros = [this]()
  {
    auto row_vector = Teuchos::rcp(new Epetra_Vector(*lambda_dof_rowmap_));
    row_vector->PutScalar(0.0);
    return row_vector;
  };
  auto lambda = create_lambda_row_vector_with_zeros();
  auto lambda_lin_constraint = create_lambda_row_vector_with_zeros();
  auto lambda_lin_kappa = create_lambda_row_vector_with_zeros();
  for (int lid = 0; lid < lambda_dof_rowmap_->NumMyElements(); lid++)
  {
    if (lambda_active_->Values()[lid] > 0.1)
    {
      const fad_type weighted_gap =
          Core::FADUtils::HigherOrderFadValue<fad_type>::apply(2, 0, constraint_->Values()[lid]);
      const fad_type kappa =
          Core::FADUtils::HigherOrderFadValue<fad_type>::apply(2, 1, kappa_->Values()[lid]);
      const fad_type scaled_gap = weighted_gap / kappa;

      // The -1 here is due to the way the lagrange multipliers are defined in the coupling
      // constraints.
      const fad_type local_lambda = -1.0 * PenaltyForce(scaled_gap, beam_to_solid_conact_params);
      lambda->ReplaceMyValue(lid, 0, Core::FADUtils::CastToDouble(local_lambda));
      lambda_lin_constraint->ReplaceMyValue(lid, 0, local_lambda.dx(0));
      lambda_lin_kappa->ReplaceMyValue(lid, 0, local_lambda.dx(1));
    }
  }
  return {lambda, lambda_lin_constraint, lambda_lin_kappa};
}

FOUR_C_NAMESPACE_CLOSE