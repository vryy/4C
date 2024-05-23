/*-----------------------------------------------------------*/
/*! \file

\brief inner status test that restricts value of update vector



\level 3
*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_INNER_STATUSTEST_UPPERBOUND_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_INNER_STATUSTEST_UPPERBOUND_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_forward_decl.hpp"
#include "4C_solver_nonlin_nox_inner_statustest_generic.hpp"  // base class

#include <NOX_Abstract_Vector.H>
#include <NOX_StatusTest_Generic.H>

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace NLN
  {
    namespace StatusTest
    {
      enum QuantityType : int;
    }  // namespace StatusTest
    namespace LineSearch
    {
      class Generic;
    }  // namespace LineSearch
    namespace INNER
    {
      namespace StatusTest
      {
        class UpperBound : public Generic
        {
         public:
          //! constructor
          UpperBound(const double& upperboundval, const ::NOX::Abstract::Vector::NormType normtype,
              const NOX::NLN::StatusTest::QuantityType qtype);

          //! Test the line search stopping criterion
          NOX::NLN::INNER::StatusTest::StatusType CheckStatus(
              const NOX::NLN::INNER::StatusTest::Interface::Required& interface,
              const ::NOX::Solver::Generic& solver, const ::NOX::Abstract::Group& grp,
              ::NOX::StatusTest::CheckType checkType) override;

          //! Return the result of the most recent checkStatus call
          NOX::NLN::INNER::StatusTest::StatusType GetStatus() const override;

          ///! Output formatted description of stopping test to output stream.
          std::ostream& Print(std::ostream& stream, int indent = 0) const override;

         protected:
          double get_search_direction_length(const NOX::NLN::LineSearch::Generic& linesearch,
              const ::NOX::Solver::Generic& solver, const ::NOX::Abstract::Group& grp) const;

         private:
          void throwError(const std::string& functionName, const std::string& errorMsg) const;

         protected:
          //! Status
          NOX::NLN::INNER::StatusTest::StatusType status_;

          //! norm type
          ::NOX::Abstract::Vector::NormType normtype_;

          //! degree of freedom type
          NOX::NLN::StatusTest::QuantityType qtype_;

          //! value for upper bound of Newton step size
          double upperboundval_;

          double reduction_fac_;

          //! current maximal value of step size
          double stepmaxval_;
        };
      }  // namespace StatusTest
    }    // namespace INNER
  }      // namespace NLN
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
