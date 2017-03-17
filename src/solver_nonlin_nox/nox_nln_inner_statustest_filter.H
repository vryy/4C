/*----------------------------------------------------------------------------*/
/*!
\file nox_nln_inner_statustest_filter.H

\brief Inner status test class for constraint problems. Filter
       techniques are based on ideas from multi-objective optimization:

       - Control of the two distinct goals of minimization of the objective
         function and satisfaction of the constraints.

       - Unlike merit functions, filter methods keep these two goals separate

\maintainer Michael Hiermeier

\date Mar 6, 2017

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_INNER_STATUSTEST_FILTER_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_INNER_STATUSTEST_FILTER_H_

#include "nox_nln_inner_statustest_generic.H"

namespace NOX {
namespace NLN {
namespace INNER {
namespace StatusTest {

  class Filter : public Generic
  {
  public:
    /// constructor
    Filter();

    /// destructor
    virtual ~Filter() {};

    //! %Test the inner stopping criterion
    /*!
      The test can (and should, if possible) be skipped if
      checkType is NOX::StatusType::None.  If the test is skipped, then
      the status should be set to NOX::StatusTest::Unevaluated.
    */
    virtual StatusType CheckStatus(
        const Interface::Required &  interface,
        const NOX::Abstract::Group & grp,
        NOX::StatusTest::CheckType   checkType);

    //! Return the result of the most recent inner checkStatus call
    virtual StatusType GetStatus() const;

    ///! Output formatted description of inner stopping test to output stream.
    virtual std::ostream& Print(std::ostream& stream, int indent = 0) const;

  protected:
    //! Status
    NOX::NLN::INNER::StatusTest::StatusType status_;

  };

} // namespace StatusTest
} // namespace INNER
} // namespace NLN
} // NOX


#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_INNER_STATUSTEST_FILTER_H_ */
