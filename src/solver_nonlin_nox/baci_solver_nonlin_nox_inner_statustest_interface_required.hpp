/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_INNER_STATUSTEST_INTERFACE_REQUIRED_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_INNER_STATUSTEST_INTERFACE_REQUIRED_HPP

#include "baci_config.hpp"

#include "baci_solver_nonlin_nox_forward_decl.hpp"

#include <NOX_Abstract_Group.H>
#include <NOX_StatusTest_Generic.H>

BACI_NAMESPACE_OPEN

namespace NOX
{
  namespace NLN
  {
    namespace INNER
    {
      namespace StatusTest
      {
        enum StatusType : int;
        namespace Interface
        {
          class Required
          {
           public:
            //! constructor
            Required(){};

            //! destructor
            virtual ~Required() = default;

            //! Get the number of inner-loop nonlinear iterations (e.g. line search iterations)
            virtual int GetNumIterations() const = 0;

            //! Get the objective or meritfunction
            virtual const ::NOX::MeritFunction::Generic& GetMeritFunction() const = 0;

            //! Execute the inner status test
            virtual NOX::NLN::INNER::StatusTest::StatusType CheckInnerStatus(
                const ::NOX::Solver::Generic& solver, const ::NOX::Abstract::Group& grp,
                ::NOX::StatusTest::CheckType checkType) const = 0;
          };  // class Required
        }     // namespace Interface
      }       // namespace StatusTest
    }         // namespace INNER
  }           // namespace NLN
}  // namespace NOX

BACI_NAMESPACE_CLOSE

#endif
