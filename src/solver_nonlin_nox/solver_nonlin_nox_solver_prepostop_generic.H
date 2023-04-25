/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SOLVER_NONLIN_NOX_SOLVER_PREPOSTOP_GENERIC_H
#define SOLVER_NONLIN_NOX_SOLVER_PREPOSTOP_GENERIC_H

#include <NOX_Observer.hpp>

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

        class Generic : public NOX::Observer
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


#endif  // SOLVER_NONLIN_NOX_SOLVER_PREPOSTOP_GENERIC_H
