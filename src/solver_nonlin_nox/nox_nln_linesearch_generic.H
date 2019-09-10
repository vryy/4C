/*-----------------------------------------------------------*/
/*! \file

\brief generic class for %NOX::NLN backtracking line search

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_LINESEARCH_GENERIC_H_
#define NOX_NLN_LINESEARCH_GENERIC_H_

#include <NOX_LineSearch_Generic.H>
#include "nox_nln_inner_statustest_interface_required.H"

namespace NOX
{
  namespace NLN
  {
    namespace LineSearch
    {
      class PrePostOperator;

      class Generic : public NOX::LineSearch::Generic,
                      public NOX::NLN::INNER::StatusTest::Interface::Required
      {
       public:
        //! constructor
        Generic() = default;

        //! destructor
        virtual ~Generic() = default;

        //! @name NOX::NLN::LineSearch::Generic
        //! @{
        //! returns the slope in the current search direction
        virtual const NOX::Abstract::Vector& GetSearchDirection() const = 0;

        //! returns the stepSize
        virtual const double& GetStepLength() const = 0;

        //! sets the stepSize
        virtual void SetStepLength(double step) = 0;
        //! @}

        //! @name NOX::LineSearch::Generic
        //! @{
        //! [derived]
        virtual bool compute(NOX::Abstract::Group& grp, double& step,
            const NOX::Abstract::Vector& dir, const NOX::Solver::Generic& s) = 0;
        //! @}

        //! @name NOX::NLN::INNER::StatusTest::Interface::Required
        //! @{
        //! get the number of line search iterations
        virtual const int& GetNumIterations() const = 0;

        //! get the merit function
        virtual const NOX::MeritFunction::Generic& GetMeritFunction() const = 0;
        //! @}

       protected:
        Teuchos::RCP<PrePostOperator> prePostOperatorPtr_ = Teuchos::null;
      };
    }  // namespace LineSearch
  }    // namespace NLN
}  // namespace NOX


#endif /* NOX_NLN_LINESEARCH_GENERIC_H_ */
