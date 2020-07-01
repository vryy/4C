/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_LINESEARCH_PREPOSTOPERATOR_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_LINESEARCH_PREPOSTOPERATOR_H_

#include <Teuchos_RCP.hpp>
#include "nox_nln_abstract_prepostoperator.H"

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos
namespace NOX
{
  namespace NLN
  {
    namespace LineSearch
    {
      enum PrePostOpType
      {
        prepost_output_every_iter
      };


      class PrePostOperator
      {
       public:
        typedef std::map<enum PrePostOpType, Teuchos::RCP<NOX::NLN::Abstract::PrePostOperator>> map;

        /// disallow the following
        PrePostOperator() = delete;
        PrePostOperator(const PrePostOperator& ppo) = delete;
        PrePostOperator& operator=(const PrePostOperator& ppo) = delete;

        /// allowed constructor
        PrePostOperator(Teuchos::ParameterList& linesearchSublist);

        virtual ~PrePostOperator(){};

        //! Resets the pre/post operator.
        virtual void reset(Teuchos::ParameterList& linesearchSublist);

        /*! Returns the inherent pre/post operator std::map of the "Line Search"
         *  sublist. If the corresponding parameter called "User Defined Pre/Post Operator"
         *  is not yet defined, a empty std::map is generated and set into the
         *  parameter list first. */
        static map& GetMutableMap(Teuchos::ParameterList& p_ls_list);

        /** User defined method that will be executed before the step is modified in
         *  the line search routine.
         *
         * \param solver     : Access to the underlying solver object.
         * \param linesearch : Access to the line search object. */
        virtual void runPreModifyStepLength(
            const NOX::Solver::Generic& solver, const NOX::LineSearch::Generic& linesearch);

       protected:
        //! Flag that determines if a pre/post operator has been supplied by user.
        bool havePrePostOperator_ = false;

        //! Points to user defined pre/post operator for the linear system.
        Teuchos::RCP<map> prePostOperatorMapPtr_ = Teuchos::null;
      };

      inline void PrePostOperator::runPreModifyStepLength(
          const NOX::Solver::Generic& solver, const NOX::LineSearch::Generic& linesearch)
      {
        if (havePrePostOperator_)
        {
          map::iterator it;
          for (auto& it : *prePostOperatorMapPtr_)
            it.second->runPreModifyStepLength(solver, linesearch);
        }
      }
    }  // namespace LineSearch
  }    // namespace NLN
}  // namespace NOX


#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_LINESEARCH_PREPOSTOPERATOR_H_ */
