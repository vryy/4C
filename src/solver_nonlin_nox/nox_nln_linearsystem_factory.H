/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN factory to create a %NOX::Epetra::LinearSystem.

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_LINEARSYSTEM_FACTORY_H_
#define NOX_NLN_LINEARSYSTEM_FACTORY_H_

#include <Teuchos_RCP.hpp>
#include "nox_nln_enum_lists.H"

// forward declaration
namespace LINALG
{
  class Solver;
  class SparseOperator;
}  // namespace LINALG

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

        //! Destructor.
        ~Factory(){};

        Teuchos::RCP<NOX::Epetra::LinearSystem> BuildLinearSystem(
            const NOX::NLN::LinSystem::LinearSystemType& linsystype,
            NOX::NLN::GlobalData& noxNlnGlobalData, const Teuchos::RCP<LINALG::SparseOperator>& jac,
            const Teuchos::RCP<NOX::Epetra::Vector>& cloneVector,
            const Teuchos::RCP<LINALG::SparseOperator>& precMat,
            const Teuchos::RCP<NOX::Epetra::Scaling>& scalingObject) const;
      };

      /*! \brief Nonmember helper function for the NOX::NLN::LinearSystem::Factory.

      \relates NOX::NLN::LinearSystem::Factory

      */
      Teuchos::RCP<NOX::Epetra::LinearSystem> BuildLinearSystem(
          const NOX::NLN::LinSystem::LinearSystemType& linsystype,
          NOX::NLN::GlobalData& noxNlnGlobalData, const Teuchos::RCP<LINALG::SparseOperator>& jac,
          const Teuchos::RCP<NOX::Epetra::Vector>& cloneVector,
          const Teuchos::RCP<LINALG::SparseOperator>& precMat,
          const Teuchos::RCP<NOX::Epetra::Scaling>& scalingObject);
    }  // namespace LinSystem
  }    // namespace NLN
}  // namespace NOX


#endif /* NOX_NLN_LINEARSYSTEM_FACTORY_H_ */
