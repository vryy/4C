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
    eval_recover,            /*!< recover internal quantities, e.g. Lagrange multipliers */
    eval_reset               /*!< reset internal quantities, e.g. displacement state and/or Lagrange multipliers */
  };

  static inline const std::string ActionType2String(const enum ActionType& act)
  {
    switch (act)
    {
      case eval_none:
        return "eval_none";
        break;
      case eval_force_stiff:
        return "eval_force_stiff";
        break;
      case eval_force:
        return "eval_force";
        break;
      default:
        return "unknown";
        break;
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
