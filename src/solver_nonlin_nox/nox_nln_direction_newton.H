/*-----------------------------------------------------------*/
/*!
\file nox_nln_direction_newton.H

\maintainer Michael Hiermeier

\date Aug 6, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_DIRECTION_NEWTON_H_
#define NOX_NLN_DIRECTION_NEWTON_H_

#include <NOX_Direction_Newton.H>  // base class

// forward declarations...
namespace NOX
{
  class Utils;
  namespace Abstract
  {
    class Vector;
    class Group;
  }  // namespace Abstract
  namespace Solver
  {
    class Generic;
  }  // namespace Solver

  namespace NLN
  {
    namespace Direction
    {
      class Newton : public NOX::Direction::Newton
      {
       public:
        //! Constructor
        Newton(const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params);

        //! Destructor
        virtual ~Newton(){};

        // derived
        virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
            const NOX::Solver::Generic& solver);

       private:
        // throw NOX error
        void throwError(const std::string& functionName, const std::string& errorMsg);

       private:
        //! NOX_Utils pointer
        Teuchos::RCP<NOX::Utils> utils_;
      };
    }  // namespace Direction
  }    // namespace NLN
}  // namespace NOX


#endif /* NOX_NLN_DIRECTION_NEWTON_H_ */
