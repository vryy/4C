/*-----------------------------------------------------------*/
/*!
\file nox_nln_linesearch_backtrack.H

\brief %NOX::NLN backtracking line search implementation.

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_LINESEARCH_BACKTRACK_H_
#define NOX_NLN_LINESEARCH_BACKTRACK_H_

#include "nox_nln_linesearch_generic.H"  // base class
#include "nox_nln_inner_statustest_generic.H"

#include "nox_nln_floating_point_exception.H"
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
        virtual ~Backtrack(){};

        /// hard reset
        bool reset(const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params);

        /// weak reset
        void reset();

        // derived
        bool compute(NOX::Abstract::Group& newgrp, double& step, const NOX::Abstract::Vector& dir,
            const NOX::Solver::Generic& s);

        /// derived
        NOX::NLN::INNER::StatusTest::StatusType CheckInnerStatus(const NOX::Solver::Generic& solver,
            const NOX::Abstract::Group& grp, NOX::StatusTest::CheckType checkType) const override;

        //! @name Access functionality
        //@{
        //! get the number of line search iterations
        virtual const int& GetNumIterations() const;

        //! get the merit function
        virtual const NOX::MeritFunction::Generic& GetMeritFunction() const;

        //! get the current search direction
        virtual const NOX::Abstract::Vector& GetSearchDirection() const;

        //! get current step length
        const double& GetStepLength() const;

        //!@}

        //! @name Mutator functionality
        //! @{
        //! set current step length
        virtual void SetStepLength(double step);
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


#endif /* NOX_NLN_LINESEARCH_BACKTRACK_H_ */
