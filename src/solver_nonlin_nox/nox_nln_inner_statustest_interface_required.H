/*-----------------------------------------------------------*/
/*!
\file nox_nln_inner_statustest_interface_required.H

\maintainer Michael Hiermeier

\date Aug 17, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_INNER_STATUSTEST_INTERFACE_REQUIRED_H_
#define NOX_NLN_INNER_STATUSTEST_INTERFACE_REQUIRED_H_

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
          };  // class Required
        }     // namespace Interface
      }       // namespace StatusTest
    }         // namespace INNER
  }           // namespace NLN
}  // namespace NOX


#endif /* NOX_NLN_INNER_STATUSTEST_INTERFACE_REQUIRED_H_ */
