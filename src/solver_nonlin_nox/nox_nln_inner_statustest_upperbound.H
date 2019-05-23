/*-----------------------------------------------------------*/
/*!

\brief inner status test that restricts value of update vector

\maintainer Maximilian Grill


\level 3
*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_INNER_STATUSTEST_UPPERBOUND_H_
#define NOX_NLN_INNER_STATUSTEST_UPPERBOUND_H_

#include "nox_nln_inner_statustest_generic.H"  // base class
#include <NOX_StatusTest_Generic.H>
#include <NOX_Abstract_Vector.H>

namespace NOX
{
  namespace Abstract
  {
    class Group;
    class Vector;
  }  // namespace Abstract
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
          UpperBound(const double& upperboundval, const NOX::Abstract::Vector::NormType normtype,
              const NOX::NLN::StatusTest::QuantityType qtype);

          //! destructor
          virtual ~UpperBound(){};

          //! Test the line search stopping criterion
          virtual NOX::NLN::INNER::StatusTest::StatusType CheckStatus(
              const NOX::NLN::INNER::StatusTest::Interface::Required& interface,
              const NOX::Solver::Generic& solver, const NOX::Abstract::Group& grp,
              NOX::StatusTest::CheckType checkType);

          //! Return the result of the most recent checkStatus call
          virtual NOX::NLN::INNER::StatusTest::StatusType GetStatus() const;

          ///! Output formatted description of stopping test to output stream.
          virtual std::ostream& Print(std::ostream& stream, int indent = 0) const;

         protected:
          double GetSearchDirectionLength(const NOX::NLN::LineSearch::Generic& linesearch,
              const NOX::Solver::Generic& solver, const NOX::Abstract::Group& grp) const;

         private:
          void throwError(const std::string& functionName, const std::string& errorMsg) const;

         protected:
          //! Status
          NOX::NLN::INNER::StatusTest::StatusType status_;

          //! norm type
          NOX::Abstract::Vector::NormType normtype_;

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

#endif /* NOX_NLN_INNER_STATUSTEST_UPPERBOUND_H_ */
