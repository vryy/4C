/*-----------------------------------------------------------*/
/*! \file

\brief Factory to create the desired non-linear solver object.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_SOLVER_FACTORY_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_SOLVER_FACTORY_HPP

#include "baci_config.hpp"

#include "baci_solver_nonlin_nox_forward_decl.hpp"

#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

namespace NOX
{
  namespace NLN
  {
    class GlobalData;
    namespace INNER
    {
      namespace StatusTest
      {
        class Generic;
      }  // namespace StatusTest
    }    // namespace INNER

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

        Teuchos::RCP<::NOX::Solver::Generic> BuildSolver(
            const Teuchos::RCP<::NOX::Abstract::Group>& grp,
            const Teuchos::RCP<::NOX::StatusTest::Generic>& outerTests,
            const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests,
            const Teuchos::RCP<NOX::NLN::GlobalData>& nlnGlobalData);
      };

      /*! \brief Nonmember helper function for the NOX::Constraint::Solver::Factory.

      \relates NOX::NLNSOL::Constraint::Solver::Factory

      */

      Teuchos::RCP<::NOX::Solver::Generic> BuildSolver(
          const Teuchos::RCP<::NOX::Abstract::Group>& grp,
          const Teuchos::RCP<::NOX::StatusTest::Generic>& outerTests,
          const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests,
          const Teuchos::RCP<NOX::NLN::GlobalData>& nlnGlobalData);

    }  // namespace Solver
  }    // namespace NLN
}  // namespace NOX

BACI_NAMESPACE_CLOSE

#endif
