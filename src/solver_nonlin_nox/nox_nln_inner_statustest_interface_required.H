/*-----------------------------------------------------------*/
/*!

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_INNER_STATUSTEST_INTERFACE_REQUIRED_H_
#define NOX_NLN_INNER_STATUSTEST_INTERFACE_REQUIRED_H_

#include <NOX_StatusTest_Generic.H>

// forward declaration
namespace NOX
{
  namespace MeritFunction
  {
    class Generic;
  }  // namespace MeritFunction
  namespace NLN
  {
    namespace INNER
    {
      namespace StatusTest
      {
        enum StatusType : int;
        namespace Interface
        {
          class Required
          {
           public:
            //! constructor
            Required(){};

            //! destructor
            virtual ~Required(){};

            //! Get the number of inner-loop nonlinear iterations (e.g. line search iterations)
            virtual const int& GetNumIterations() const = 0;

            //! Get the objective or meritfunction
            virtual const NOX::MeritFunction::Generic& GetMeritFunction() const = 0;

            //! Execute the inner status test
            virtual NOX::NLN::INNER::StatusTest::StatusType CheckInnerStatus(
                const NOX::Solver::Generic& solver, const NOX::Abstract::Group& grp,
                NOX::StatusTest::CheckType checkType) const = 0;
          };  // class Required
        }     // namespace Interface
      }       // namespace StatusTest
    }         // namespace INNER
  }           // namespace NLN
}  // namespace NOX


#endif /* NOX_NLN_INNER_STATUSTEST_INTERFACE_REQUIRED_H_ */
