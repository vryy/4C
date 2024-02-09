/*-----------------------------------------------------------*/
/*! \file

\brief Generic class for all constraint submodel evaluators.


\level 3
*/
#ifndef BACI_CONSTRAINT_FRAMEWORK_SUBMODELEVALUATOR_BASE_HPP
#define BACI_CONSTRAINT_FRAMEWORK_SUBMODELEVALUATOR_BASE_HPP

#include "baci_config.hpp"

#include "baci_structure_new_model_evaluator_generic.hpp"

BACI_NAMESPACE_OPEN

namespace CONSTRAINTS::SUBMODELEVALUATOR
{
  /*! \brief Interface class of all submodel evaluators managing
   *  constraint terms
   */
  class Constraint
  {
   public:
    /*! \brief Constructor
     *
     */
    Constraint() = default;

    /**
     * \brief Destructor
     */
    virtual ~Constraint() = default;

    //! \brief Evaluate the current tangential stiffness matrices at \f$t_{n+1}\f$
    virtual void EvaluateStiff() = 0;

    //! \brief Evaluate the current right-hand-sides at \f$t_{n+1}\f$
    virtual void EvaluateForce() = 0;
  };
}  // namespace CONSTRAINTS::SUBMODELEVALUATOR


BACI_NAMESPACE_CLOSE

#endif
