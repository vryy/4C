/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN Weighted root mean square test of the solution
       increment. A detailed description can be found in the NOX
       documentation.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_SOLVER_NONLIN_NOX_STATUSTEST_NORMWRMS_H
#define BACI_SOLVER_NONLIN_NOX_STATUSTEST_NORMWRMS_H

#include "baci_solver_nonlin_nox_enum_lists.H"

#include <NOX_StatusTest_Generic.H>  // base class
#include <Teuchos_RCP.hpp>

#include <vector>

BACI_NAMESPACE_OPEN

namespace NOX
{
  namespace NLN
  {
    namespace StatusTest
    {
      class NormWRMS : public ::NOX::StatusTest::Generic
      {
       public:
        /*! \brief Constructor
         *  At the moment we support only scalar \c ATOL input values for each quantity.
         *  Extensions to a vector-based ATOL are possible. See \c ::NOX::StatusTest::NormWRMS
         *  for more information. */
        NormWRMS(const std::vector<NOX::NLN::StatusTest::QuantityType>& checkList,
            const std::vector<double>& rtol, const std::vector<double>& atol,
            const std::vector<double>& BDFMultiplier, const std::vector<double>& tolerance,
            const double& alpha, const double& beta,
            const std::vector<bool>& disable_implicit_weighting);

        ::NOX::StatusTest::StatusType checkStatus(
            const ::NOX::Solver::Generic& problem, ::NOX::StatusTest::CheckType checkType) override;

        //! Check for the given quantity
        bool IsQuantity(const NOX::NLN::StatusTest::QuantityType& qType) const;

        ::NOX::StatusTest::StatusType getStatus() const override;

        //! returns the absolute tolerance of the given quantity
        double GetAbsoluteTolerance(const NOX::NLN::StatusTest::QuantityType& qType) const;

        //! returns the relative tolerance of the given quantity
        double GetRelativeTolerance(const NOX::NLN::StatusTest::QuantityType& qType) const;

        std::ostream& print(std::ostream& stream, int indent) const override;

       private:
        //! calculated norm for the different quantities
        Teuchos::RCP<std::vector<double>> normWRMS_;

        //! number of quantities to check
        std::size_t nChecks_;

        //! nox_nln_statustest quantities which are checked
        std::vector<NOX::NLN::StatusTest::QuantityType> checkList_;

        //! relative tolerance
        std::vector<double> rtol_;

        //! absolute tolerance
        std::vector<double> atol_;

        //! Time integration method multiplier (BDF Multiplier)
        std::vector<double> factor_;

        //! Required tolerance for the NormWRMS to be declared converged.
        std::vector<double> tol_;

        //! Minimum step size allowed during a line search for WRMS norm to be flagged as converged.
        double alpha_;

        //! Actual step size used during line search.
        double computedStepSize_;

        //! Maximum linear solve tolerance allowed for WRMS norm to be flagged as converged.
        double beta_;

        //! Actual tolerance achieved by the linear solver during the last linear solve.
        double achievedTol_;

        //! Global status
        ::NOX::StatusTest::StatusType gStatus_;

        //! Status of each quantity
        std::vector<::NOX::StatusTest::StatusType> status_;

        //! Flag that tells the print method whether to print the criteria 2 information.
        bool printCriteria2Info_;

        //! Flag that tells the print method whether to print the criteria 3 information.
        bool printCriteria3Info_;

        //! If true, the implicit weighting will be disabled during norm calculation
        std::vector<bool> disable_implicit_weighting_;
      };
    }  // namespace StatusTest
  }    // namespace NLN
}  // namespace NOX


BACI_NAMESPACE_CLOSE

#endif  // SOLVER_NONLIN_NOX_STATUSTEST_NORMWRMS_H
