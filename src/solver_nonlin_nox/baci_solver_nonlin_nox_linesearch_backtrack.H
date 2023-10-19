/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN backtracking line search implementation.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_SOLVER_NONLIN_NOX_LINESEARCH_BACKTRACK_H
#define BACI_SOLVER_NONLIN_NOX_LINESEARCH_BACKTRACK_H

#include "baci_solver_nonlin_nox_floating_point_exception.H"
#include "baci_solver_nonlin_nox_inner_statustest_generic.H"
#include "baci_solver_nonlin_nox_linesearch_generic.H"  // base class

#include <NOX_StatusTest_Generic.H>

// forward declaration
namespace NOX
{
  class Utils;
  namespace MeritFunction
  {
    class Generic;
  }
  namespace StatusTest
  {
    class Generic;
  }

  namespace NLN
  {
    namespace INNER
    {
      namespace StatusTest
      {
        class Generic;
      }  // namespace StatusTest
    }    // namespace INNER
    namespace LineSearch
    {
      class Backtrack : public Generic
      {
       public:
        //! Constructor
        Backtrack(const Teuchos::RCP<NOX::GlobalData>& gd,
            const Teuchos::RCP<NOX::StatusTest::Generic> outerTests,
            const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> innerTests,
            Teuchos::ParameterList& params);

        //! Destructor
        ~Backtrack() override{};

        /// hard reset
        bool reset(const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params);

        /// weak reset
        void reset();

        bool compute(NOX::Abstract::Group& newgrp, double& step, const NOX::Abstract::Vector& dir,
            const NOX::Solver::Generic& s) override;

        NOX::NLN::INNER::StatusTest::StatusType CheckInnerStatus(const NOX::Solver::Generic& solver,
            const NOX::Abstract::Group& grp, NOX::StatusTest::CheckType checkType) const override;

        //! @name Access functionality
        //@{
        //! get the number of line search iterations
        int GetNumIterations() const override;

        //! get the merit function
        const NOX::MeritFunction::Generic& GetMeritFunction() const override;

        //! get the current search direction
        const NOX::Abstract::Vector& GetSearchDirection() const override;

        //! get current step length
        double GetStepLength() const override;

        //!@}

        //! @name Mutator functionality
        //! @{
        //! set current step length
        void SetStepLength(double step) override;
        //! @}

       protected:
        //! print inner status test results
        void PrintUpdate(std::ostream& os) const;

       private:
        //! throw NOX error
        void throwError(const std::string& functionName, const std::string& errorMsg) const;

       private:
        //! handle floating point exceptions
        FloatingPointException fp_except_;

        //! inner iteration counter
        int lsIters_;

        //! Printing Utilities
        Teuchos::RCP<NOX::Utils> utils_;

        //! Merit function
        Teuchos::RCP<NOX::MeritFunction::Generic> meritFunctionPtr_;

        //! Current step pointer, points to the step variable of the
        //! NOX::NLN::Solver::LineSearchBased class
        double* stepPtr_;

        //! Default step
        double defaultStep_;

        //! line search reduction factor
        double reductionFactor_;

        //! inner status test checktype
        NOX::StatusTest::CheckType checkType_;

        //! inner status type
        NOX::NLN::INNER::StatusTest::StatusType status_;

        //! search direction
        Teuchos::RCP<const NOX::Abstract::Vector> searchDirectionPtr_;

        //! outer convergence tests
        Teuchos::RCP<NOX::StatusTest::Generic> outerTestsPtr_;

        //! line search stopping tests
        Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> innerTestsPtr_;
      };
    }  // namespace LineSearch
  }    // namespace NLN
}  // namespace NOX


#endif  // SOLVER_NONLIN_NOX_LINESEARCH_BACKTRACK_H
