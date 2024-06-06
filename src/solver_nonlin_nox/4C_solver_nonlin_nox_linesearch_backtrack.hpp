/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN backtracking line search implementation.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_LINESEARCH_BACKTRACK_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_LINESEARCH_BACKTRACK_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_floating_point_exception.hpp"
#include "4C_solver_nonlin_nox_forward_decl.hpp"
#include "4C_solver_nonlin_nox_inner_statustest_generic.hpp"
#include "4C_solver_nonlin_nox_linesearch_generic.hpp"  // base class

#include <NOX_StatusTest_Generic.H>

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace Nln
  {
    namespace Inner
    {
      namespace StatusTest
      {
        class Generic;
      }  // namespace StatusTest
    }    // namespace Inner
    namespace LineSearch
    {
      class Backtrack : public Generic
      {
       public:
        //! Constructor
        Backtrack(const Teuchos::RCP<::NOX::GlobalData>& gd,
            const Teuchos::RCP<::NOX::StatusTest::Generic> outerTests,
            const Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic> innerTests,
            Teuchos::ParameterList& params);

        /// hard reset
        bool reset(const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params);

        /// weak reset
        void reset();

        bool compute(::NOX::Abstract::Group& newgrp, double& step,
            const ::NOX::Abstract::Vector& dir, const ::NOX::Solver::Generic& s) override;

        NOX::Nln::Inner::StatusTest::StatusType CheckInnerStatus(
            const ::NOX::Solver::Generic& solver, const ::NOX::Abstract::Group& grp,
            ::NOX::StatusTest::CheckType checkType) const override;

        //! @name Access functionality
        //@{
        //! get the number of line search iterations
        int GetNumIterations() const override;

        //! get the merit function
        const ::NOX::MeritFunction::Generic& GetMeritFunction() const override;

        //! get the current search direction
        const ::NOX::Abstract::Vector& GetSearchDirection() const override;

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
        void print_update(std::ostream& os) const;

       private:
        //! throw NOX error
        void throw_error(const std::string& functionName, const std::string& errorMsg) const;

       private:
        //! handle floating point exceptions
        FloatingPointException fp_except_;

        //! inner iteration counter
        int ls_iters_;

        //! Printing Utilities
        Teuchos::RCP<::NOX::Utils> utils_;

        //! Merit function
        Teuchos::RCP<::NOX::MeritFunction::Generic> merit_function_ptr_;

        //! Current step pointer, points to the step variable of the
        //! NOX::Nln::Solver::LineSearchBased class
        double* step_ptr_;

        //! Default step
        double default_step_;

        //! line search reduction factor
        double reduction_factor_;

        //! inner status test checktype
        ::NOX::StatusTest::CheckType check_type_;

        //! inner status type
        NOX::Nln::Inner::StatusTest::StatusType status_;

        //! search direction
        Teuchos::RCP<const ::NOX::Abstract::Vector> search_direction_ptr_;

        //! outer convergence tests
        Teuchos::RCP<::NOX::StatusTest::Generic> outer_tests_ptr_;

        //! line search stopping tests
        Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic> inner_tests_ptr_;
      };
    }  // namespace LineSearch
  }    // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
