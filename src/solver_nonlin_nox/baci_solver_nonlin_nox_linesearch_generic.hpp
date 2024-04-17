/*-----------------------------------------------------------*/
/*! \file

\brief generic class for %NOX::NLN backtracking line search



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_LINESEARCH_GENERIC_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_LINESEARCH_GENERIC_HPP

#include "baci_config.hpp"

#include "baci_solver_nonlin_nox_inner_statustest_interface_required.hpp"

#include <NOX_LineSearch_Generic.H>

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace NLN
  {
    namespace LineSearch
    {
      class PrePostOperator;

      class Generic : public ::NOX::LineSearch::Generic,
                      public NOX::NLN::INNER::StatusTest::Interface::Required
      {
       public:
        //! constructor
        Generic() = default;


        //! @name NOX::NLN::LineSearch::Generic
        //! @{
        //! returns the slope in the current search direction
        virtual const ::NOX::Abstract::Vector& GetSearchDirection() const = 0;

        //! returns the stepSize
        virtual double GetStepLength() const = 0;

        //! sets the stepSize
        virtual void SetStepLength(double step) = 0;
        //! @}

        //! @name ::NOX::LineSearch::Generic
        //! @{
        bool compute(::NOX::Abstract::Group& grp, double& step, const ::NOX::Abstract::Vector& dir,
            const ::NOX::Solver::Generic& s) override = 0;
        //! @}

        //! @name NOX::NLN::INNER::StatusTest::Interface::Required
        //! @{
        //! get the number of line search iterations
        int GetNumIterations() const override = 0;

        //! get the merit function
        const ::NOX::MeritFunction::Generic& GetMeritFunction() const override = 0;
        //! @}

       protected:
        Teuchos::RCP<PrePostOperator> prePostOperatorPtr_ = Teuchos::null;
      };
    }  // namespace LineSearch
  }    // namespace NLN
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
