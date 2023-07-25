/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_SOLVER_NONLIN_NOX_DIRECTION_NEWTON_H
#define BACI_SOLVER_NONLIN_NOX_DIRECTION_NEWTON_H

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


#endif  // SOLVER_NONLIN_NOX_DIRECTION_NEWTON_H
