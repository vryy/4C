/*----------------------------------------------------------------------*/
/*!
\file mortar_paramsinterface.H

\brief Mortar parameter interface. Necessary for the communication
       between the structural time integration framework and the
       mortar strategies.

\level 2

\maintainer Michael Hiermeier
\date Apr 18, 2016

*/
/*-----------------------------------------------------------------------*/

#ifndef SRC_DRT_MORTAR_MORTAR_PARAMSINTERFACE_H_
#define SRC_DRT_MORTAR_MORTAR_PARAMSINTERFACE_H_

#include "../headers/any_data_container.H"

namespace MORTAR {
  enum ActionType : int {
    eval_none,                 /*!< No evaluation type has been chosen */
    eval_force_stiff,          /*!< Evaluation of the contact/meshtying right-hand-side and the
                                    contact/meshtying jacobian. We call this method also, when we
                                    are only interested in the jacobian, since the created overhead
                                    is negligible. */
    eval_force,                 /*!< Evaluation of the contact/meshtying right-hand-side only. Necessary
                                     and meaningful for line search strategies for example. */
    eval_run_post_evaluate,     /*!< Run in the end of a call to STR::ModelEvaluator::EvaluteForce/Stiff/ForceStiff */
    eval_weighted_gap,          /*!< Evaluate the weighted gap only. */
    eval_run_post_compute_x,    /*!< recover internal quantities, e.g. Lagrange multipliers */
    eval_reset,                 /*!< reset internal quantities, e.g. displacement state and/or Lagrange multipliers */
    eval_run_pre_compute_x,     /*!< augment the solution direction at the very beginning of a ComputeX call */
    eval_run_post_iterate,      /*!< run in the end of a NOX::Solver::Step() call */
    eval_contact_potential,     /*!< Evaluate the contact potential */
    eval_wgap_gradient_error,   /*!< Evaluate the error of the weighted gap gradient */
    eval_static_constraint_rhs, /*!< Evaluate only the contributions to the constraint rhs.
                                     The active set is not updated during the evaluation. */
    eval_run_post_apply_jacobian_inverse, /*!< run in the end of a
                                     NOX::NLN::LinearSystem::applyJacobianInverse call */
    eval_correct_parameters
  };

  static inline const std::string ActionType2String(const enum ActionType& act)
  {
    switch (act)
    {
      case eval_none:
        return "eval_none";
      case eval_force_stiff:
        return "eval_force_stiff";
      case eval_force:
        return "eval_force";
      case eval_run_post_evaluate:
        return "eval_run_post_evaluate";
      case eval_weighted_gap:
        return "eval_weighted_gap";
      case MORTAR::eval_reset:
        return "eval_reset";
      case MORTAR::eval_run_post_compute_x:
        return "eval_run_post_compute_x";
      case MORTAR::eval_run_pre_compute_x:
        return "eval_run_pre_compute_x";
      case MORTAR::eval_run_post_iterate:
        return "eval_run_post_iterate";
      case MORTAR::eval_contact_potential:
        return "eval_contact_potential";
      case MORTAR::eval_wgap_gradient_error:
        return "eval_wgap_gradient_error";
      case MORTAR::eval_static_constraint_rhs:
        return "eval_static_constraint_rhs";
      case MORTAR::eval_run_post_apply_jacobian_inverse:
        return "eval_run_post_apply_jacobian_inverse";
      case MORTAR::eval_correct_parameters:
        return "eval_correct_parameters";
      default:
        return "unknown";
    }
    return "";
  };

  class ParamsInterface : public GEN::AnyDataContainer
  {
  public:
    //! destructor
    virtual ~ParamsInterface() {};

    virtual enum ActionType GetActionType() const = 0;

    virtual int GetNlnIter() const = 0;

    virtual int GetStepNp() const = 0;

    virtual int GetRestartStep() const = 0;

  };  // class ParamsInterface
} // namespace MORTAR


#endif /* SRC_DRT_MORTAR_MORTAR_PARAMSINTERFACE_H_ */
