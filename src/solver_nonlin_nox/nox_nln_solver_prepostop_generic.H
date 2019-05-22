/*-----------------------------------------------------------*/
/*!
\file nox_nln_solver_prepostop_generic.H

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_PREPOSTOPERATOR_H_
#define NOX_NLN_PREPOSTOPERATOR_H_

#include <NOX_Abstract_PrePostOperator.H>  // base class

namespace NOX
{
  namespace NLN
  {
    namespace Solver
    {
      namespace PrePostOp
      {
        /*! \brief Non-linear solver helper class
         *
         * The NOX::NLN::Solver s are supposed to use this helper class to save
         * solver dependent things, such as the number of nonlinear solver steps.
         *
         * ToDo Extend the functionality to non-linear solvers which are not
         * derived from the LineSearchBased class.
         *
         * \author Michael Hiermeier */
        class Generic : public NOX::Abstract::PrePostOperator
        {
         public:
          //! constructor
          Generic();

          //! destructor
          virtual ~Generic(){};

          //! derived
          virtual void runPreIterate(const NOX::Solver::Generic& solver);

          //! derived
          virtual void runPreSolve(const NOX::Solver::Generic& nlnSolver);
        };
      }  // namespace PrePostOp
    }    // namespace Solver
  }      // namespace NLN
}  // namespace NOX


#endif /* NOX_NLN_PREPOSTOPERATOR_H_ */
