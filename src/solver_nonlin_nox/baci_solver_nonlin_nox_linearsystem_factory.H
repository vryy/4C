/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN factory to create a %NOX::Epetra::LinearSystem.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_SOLVER_NONLIN_NOX_LINEARSYSTEM_FACTORY_H
#define BACI_SOLVER_NONLIN_NOX_LINEARSYSTEM_FACTORY_H

#include "baci_solver_nonlin_nox_enum_lists.H"

#include <Teuchos_RCP.hpp>

// forward declaration
namespace CORE::LINALG
{
  class Solver;
  class SparseOperator;
}  // namespace CORE::LINALG

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

namespace NOX
{
  namespace Epetra
  {
    class Vector;
    class Scaling;
    class LinearSystem;
    namespace Interface
    {
      class Required;
      class Jacobian;
      class Preconditioner;
    }  // namespace Interface
  }    // namespace Epetra

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


        Teuchos::RCP<NOX::Epetra::LinearSystem> BuildLinearSystem(
            const NOX::NLN::LinSystem::LinearSystemType& linsystype,
            NOX::NLN::GlobalData& noxNlnGlobalData,
            const Teuchos::RCP<CORE::LINALG::SparseOperator>& jac,
            const Teuchos::RCP<NOX::Epetra::Vector>& cloneVector,
            const Teuchos::RCP<CORE::LINALG::SparseOperator>& precMat,
            const Teuchos::RCP<NOX::Epetra::Scaling>& scalingObject) const;
      };

      /*! \brief Nonmember helper function for the NOX::NLN::LinearSystem::Factory.

      \relates NOX::NLN::LinearSystem::Factory

      */
      Teuchos::RCP<NOX::Epetra::LinearSystem> BuildLinearSystem(
          const NOX::NLN::LinSystem::LinearSystemType& linsystype,
          NOX::NLN::GlobalData& noxNlnGlobalData,
          const Teuchos::RCP<CORE::LINALG::SparseOperator>& jac,
          const Teuchos::RCP<NOX::Epetra::Vector>& cloneVector,
          const Teuchos::RCP<CORE::LINALG::SparseOperator>& precMat,
          const Teuchos::RCP<NOX::Epetra::Scaling>& scalingObject);
    }  // namespace LinSystem
  }    // namespace NLN
}  // namespace NOX


#endif  // SOLVER_NONLIN_NOX_LINEARSYSTEM_FACTORY_H
