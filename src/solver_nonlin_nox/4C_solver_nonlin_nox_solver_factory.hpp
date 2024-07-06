/*-----------------------------------------------------------*/
/*! \file

\brief Factory to create the desired non-linear solver object.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_SOLVER_FACTORY_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_SOLVER_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_forward_decl.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace Nln
  {
    class GlobalData;
    namespace Inner
    {
      namespace StatusTest
      {
        class Generic;
      }  // namespace StatusTest
    }    // namespace Inner

    namespace Solver
    {
      /*! \brief %Modified Factory class to control the creation of solvers derived from the
      ::NOX::Solver::Generic object.

      \author Michael Hiermeier
      */

      class Factory
      {
       public:
        //! Constructor.
        Factory();

        //! Destructor.
        virtual ~Factory() = default;

        Teuchos::RCP<::NOX::Solver::Generic> build_solver(
            const Teuchos::RCP<::NOX::Abstract::Group>& grp,
            const Teuchos::RCP<::NOX::StatusTest::Generic>& outerTests,
            const Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>& innerTests,
            const Teuchos::RCP<NOX::Nln::GlobalData>& nlnGlobalData);
      };

      /*! \brief Nonmember helper function for the NOX::Constraint::Solver::Factory.

      \relates NOX::NlnSol::Constraint::Solver::Factory

      */

      Teuchos::RCP<::NOX::Solver::Generic> BuildSolver(
          const Teuchos::RCP<::NOX::Abstract::Group>& grp,
          const Teuchos::RCP<::NOX::StatusTest::Generic>& outerTests,
          const Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>& innerTests,
          const Teuchos::RCP<NOX::Nln::GlobalData>& nlnGlobalData);

    }  // namespace Solver
  }    // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
