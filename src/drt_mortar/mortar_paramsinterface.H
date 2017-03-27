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


namespace MORTAR {
  enum ActionType {
    eval_none,               /*!< No evaluation type has been chosen */
    eval_force_stiff,        /*!< Evaluation of the contact/meshtying right-hand-side and the
                                  contact/meshtying jacobian. We call this method also, when we
                                  are only interested in the jacobian, since the created overhead
                                  is negligible. */
    eval_force,              /*!< Evaluation of the contact/meshtying right-hand-side only. Necessary
                                  and meaningful for line search strategies for example. */
    eval_run_post_evaluate,  /*!< Run in the end of a call to STR::ModelEvaluator::EvaluteForce/Stiff/ForceStiff */
    eval_weighted_gap,       /*!< Evaluate the weighted gap only. */
    eval_recover,            /*!< recover internal quantities, e.g. Lagrange multipliers */
    eval_reset,              /*!< reset internal quantities, e.g. displacement state and/or Lagrange multipliers */
    eval_run_pre_compute_x,  /*!< augment the solution direction at the very beginning of a ComputeX call */
    eval_run_post_iterate,   /*!< run in the end of a NOX::Solver::Step() call */
    eval_contact_potential   /*!< Evaluate the contact potential */
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
      case MORTAR::eval_recover:
        return "eval_recover";
      case MORTAR::eval_run_pre_compute_x:
        return "eval_run_pre_compute_x";
      case MORTAR::eval_run_post_iterate:
        return "eval_run_post_iterate";
      case MORTAR::eval_contact_potential:
        return "eval_contact_potential";
      default:
        return "unknown";
    }
    return "";
  };

  class ParamsInterface
  {
  public:
    //! destructor
    virtual ~ParamsInterface() {};

    virtual enum ActionType GetActionType() const = 0;

    virtual int GetNlnIter() const = 0;

    virtual int GetStepNp() const = 0;

  };  // class ParamsInterface
} // namespace MORTAR


#endif /* SRC_DRT_MORTAR_MORTAR_PARAMSINTERFACE_H_ */
