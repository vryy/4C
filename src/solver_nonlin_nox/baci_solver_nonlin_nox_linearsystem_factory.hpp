/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN factory to create a %::NOX::Epetra::LinearSystem.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_LINEARSYSTEM_FACTORY_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_LINEARSYSTEM_FACTORY_HPP

#include "baci_config.hpp"

#include "baci_solver_nonlin_nox_enum_lists.hpp"
#include "baci_solver_nonlin_nox_forward_decl.hpp"

#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

// forward declaration
namespace CORE::LINALG
{
  class Solver;
  class SparseOperator;
}  // namespace CORE::LINALG

namespace NOX
{
  namespace NLN
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
            const NOX::NLN::LinSystem::LinearSystemType& linsystype,
            NOX::NLN::GlobalData& noxNlnGlobalData,
            const Teuchos::RCP<CORE::LINALG::SparseOperator>& jac,
            const Teuchos::RCP<::NOX::Epetra::Vector>& cloneVector,
            const Teuchos::RCP<CORE::LINALG::SparseOperator>& precMat,
            const Teuchos::RCP<::NOX::Epetra::Scaling>& scalingObject) const;
      };

      /*! \brief Nonmember helper function for the NOX::NLN::LinearSystem::Factory.

      \relates NOX::NLN::LinearSystem::Factory

      */
      Teuchos::RCP<::NOX::Epetra::LinearSystem> BuildLinearSystem(
          const NOX::NLN::LinSystem::LinearSystemType& linsystype,
          NOX::NLN::GlobalData& noxNlnGlobalData,
          const Teuchos::RCP<CORE::LINALG::SparseOperator>& jac,
          const Teuchos::RCP<::NOX::Epetra::Vector>& cloneVector,
          const Teuchos::RCP<CORE::LINALG::SparseOperator>& precMat,
          const Teuchos::RCP<::NOX::Epetra::Scaling>& scalingObject);
    }  // namespace LinSystem
  }    // namespace NLN
}  // namespace NOX

BACI_NAMESPACE_CLOSE

#endif
