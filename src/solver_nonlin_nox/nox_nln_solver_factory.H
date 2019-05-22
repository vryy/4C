/*-----------------------------------------------------------*/
/*!
\file nox_nln_solver_factory.H

\brief Factory to create the desired non-linear solver object.

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_SOLVER_FACTORY_H_
#define NOX_NLN_SOLVER_FACTORY_H_

#include <Teuchos_RCP.hpp>

// forward declaration
namespace NOX
{
  namespace Abstract
  {
    class Group;
  }  // namespace Abstract

  namespace StatusTest
  {
    class Generic;
  }  // namespace StatusTest

  namespace Solver
  {
    class Generic;
  }  // namespace Solver

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
      NOX::Solver::Generic object.

      \author Michael Hiermeier
      */

      class Factory
      {
       public:
        //! Constructor.
        Factory();

        //! Destructor.
        virtual ~Factory(){};

        Teuchos::RCP<NOX::Solver::Generic> BuildSolver(
            const Teuchos::RCP<NOX::Abstract::Group>& grp,
            const Teuchos::RCP<NOX::StatusTest::Generic>& outerTests,
            const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests,
            const Teuchos::RCP<NOX::NLN::GlobalData>& params);
      };

      /*! \brief Nonmember helper function for the NOX::Constraint::Solver::Factory.

      \relates NOX::NLNSOL::Constraint::Solver::Factory

      */

      Teuchos::RCP<NOX::Solver::Generic> BuildSolver(const Teuchos::RCP<NOX::Abstract::Group>& grp,
          const Teuchos::RCP<NOX::StatusTest::Generic>& outerTests,
          const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests,
          const Teuchos::RCP<NOX::NLN::GlobalData>& params);

    }  // namespace Solver
  }    // namespace NLN
}  // namespace NOX


#endif /* NOX_NLN_SOLVER_FACTORY_H_ */
