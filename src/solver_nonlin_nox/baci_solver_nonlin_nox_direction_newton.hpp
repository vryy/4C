/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_DIRECTION_NEWTON_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_DIRECTION_NEWTON_HPP

#include "baci_config.hpp"

#include "baci_solver_nonlin_nox_forward_decl.hpp"

#include <NOX_Direction_Newton.H>  // base class

BACI_NAMESPACE_OPEN

namespace NOX
{
  namespace NLN
  {
    namespace Direction
    {
      class Newton : public ::NOX::Direction::Newton
      {
       public:
        //! Constructor
        Newton(const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params);


        bool compute(::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& group,
            const ::NOX::Solver::Generic& solver) override;

       private:
        // throw NOX error
        void throwError(const std::string& functionName, const std::string& errorMsg);

       private:
        //! NOX_Utils pointer
        Teuchos::RCP<::NOX::Utils> utils_;
      };
    }  // namespace Direction
  }    // namespace NLN
}  // namespace NOX

BACI_NAMESPACE_CLOSE

#endif
