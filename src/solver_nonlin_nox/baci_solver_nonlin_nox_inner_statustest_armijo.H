/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_SOLVER_NONLIN_NOX_INNER_STATUSTEST_ARMIJO_H
#define BACI_SOLVER_NONLIN_NOX_INNER_STATUSTEST_ARMIJO_H

#include "baci_solver_nonlin_nox_inner_statustest_generic.H"  // base class

#include <NOX_StatusTest_Generic.H>

#include <deque>

namespace NOX
{
  namespace Abstract
  {
    class Group;
  }  // namespace Abstract
}  // namespace NOX

BACI_NAMESPACE_OPEN

namespace NOX
{
  namespace NLN
  {
    namespace LineSearch
    {
      class Generic;
    }  // namespace LineSearch
    namespace INNER
    {
      namespace StatusTest
      {
        class Armijo : public Generic
        {
         public:
          //! constructor
          Armijo(
              const double& c_1, const bool& isMonotone = true, const std::size_t& maxHistSize = 1);

          //! %Test the line search stopping criterion
          /*!
            The test can (and should, if possible) be skipped if
            checkType is NOX::StatusType::None. If the test is skipped, then
            the status should be set to ::NOX::StatusTest::Unevaluated.
          */
          NOX::NLN::INNER::StatusTest::StatusType CheckStatus(
              const NOX::NLN::INNER::StatusTest::Interface::Required& interface,
              const ::NOX::Solver::Generic& solver, const ::NOX::Abstract::Group& grp,
              ::NOX::StatusTest::CheckType checkType) override;

          //! Return the result of the most recent checkStatus call
          NOX::NLN::INNER::StatusTest::StatusType GetStatus() const override;

          ///! Output formatted description of stopping test to output stream.
          std::ostream& Print(std::ostream& stream, int indent = 0) const override;

         protected:
          bool Setup(
              const NOX::NLN::LineSearch::Generic& linesearch, const ::NOX::Abstract::Group& grp);

         private:
          void throwError(const std::string& functionName, const std::string& errorMsg) const;

         protected:
          //! Status
          NOX::NLN::INNER::StatusTest::StatusType status_;

          //! slope scaling parameter
          double c_1_;

          //! reference function value of the last accepted state, or if we use a non-monotone
          //! behavior the maximum of the last maxHistSize accepted steps.
          double fref_;

          //! current function value
          double fcurr_;

          //! slope in the current search direction
          double slope_;

          //! current step length
          double step_;

          //! boolean which indicates if we use a monotone or a non-monotone behavior
          bool isMonotone_;

          //! maximal size of the history vector
          std::size_t maxHistSize_;

          //! history vector
          std::deque<double> histVector_;
        };
      }  // namespace StatusTest
    }    // namespace INNER
  }      // namespace NLN
}  // namespace NOX


BACI_NAMESPACE_CLOSE

#endif  // SOLVER_NONLIN_NOX_INNER_STATUSTEST_ARMIJO_H
