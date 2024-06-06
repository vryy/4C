/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN factory to create a %::NOX::Epetra::LinearSystem.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_LINEARSYSTEM_FACTORY_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_LINEARSYSTEM_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_enum_lists.hpp"
#include "4C_solver_nonlin_nox_forward_decl.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Core::LinAlg
{
  class Solver;
  class SparseOperator;
}  // namespace Core::LinAlg

namespace NOX
{
  namespace Nln
  {
    class GlobalData;
    namespace LinSystem
    {
      class Factory
      {
       public:
        //! Constructor.
        Factory();


        Teuchos::RCP<::NOX::Epetra::LinearSystem> BuildLinearSystem(
            const NOX::Nln::LinSystem::LinearSystemType& linsystype,
            NOX::Nln::GlobalData& noxNlnGlobalData,
            const Teuchos::RCP<Core::LinAlg::SparseOperator>& jac,
            const Teuchos::RCP<::NOX::Epetra::Vector>& cloneVector,
            const Teuchos::RCP<Core::LinAlg::SparseOperator>& precMat,
            const Teuchos::RCP<::NOX::Epetra::Scaling>& scalingObject) const;
      };

      /*! \brief Nonmember helper function for the NOX::Nln::LinearSystem::Factory.

      \relates NOX::Nln::LinearSystem::Factory

      */
      Teuchos::RCP<::NOX::Epetra::LinearSystem> BuildLinearSystem(
          const NOX::Nln::LinSystem::LinearSystemType& linsystype,
          NOX::Nln::GlobalData& noxNlnGlobalData,
          const Teuchos::RCP<Core::LinAlg::SparseOperator>& jac,
          const Teuchos::RCP<::NOX::Epetra::Vector>& cloneVector,
          const Teuchos::RCP<Core::LinAlg::SparseOperator>& precMat,
          const Teuchos::RCP<::NOX::Epetra::Scaling>& scalingObject);
    }  // namespace LinSystem
  }    // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
