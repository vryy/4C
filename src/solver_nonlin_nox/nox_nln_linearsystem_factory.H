/*-----------------------------------------------------------*/
/*!
\file nox_nln_linearsystem_factory.H

\maintainer Michael Hiermeier

\date Jul 14, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_LINEARSYSTEM_FACTORY_H_
#define NOX_NLN_LINEARSYSTEM_FACTORY_H_

#include <Teuchos_RCP.hpp>
#include "nox_nln_globaldata.H"

// forward declaration
namespace LINALG {
  class Solver;
  class SparseOperator;
} // namespace LINALG

namespace Teuchos {
  class ParameterList;
} // namespace Teuchos

namespace NOX {
namespace Epetra {
  class Vector;
  class Scaling;
  class LinearSystem;
namespace Interface {
  class Required;
  class Jacobian;
  class Preconditioner;
} // namespace Interface
} // namespace Epetra

namespace NLN {
namespace LinSystem {
  class Factory
  {
  public:
    //! Constructor.
    Factory();

    //! Destructor.
    ~Factory() {};

    Teuchos::RCP<NOX::Epetra::LinearSystem>
    BuildLinearSystem(
        Teuchos::ParameterList& printParams,
        Teuchos::ParameterList& linearSolverParams,
        const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& linSolvers,
        const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
        const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
        const Teuchos::RCP<LINALG::SparseOperator>& J,
        const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
        const Teuchos::RCP<LINALG::SparseOperator>& M,
        const NOX::Epetra::Vector& cloneVector,
        const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject);
  };

  /*! \brief Nonmember helper function for the NOX::NLN::LinearSystem::Factory.

  \relates NOX::NLN::LinearSystem::Factory

  */
  Teuchos::RCP<NOX::Epetra::LinearSystem>
  BuildLinearSystem(
      Teuchos::ParameterList& printParams,
      Teuchos::ParameterList& linearSolverParams,
      const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& linSolvers,
      const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
      const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
      const Teuchos::RCP<LINALG::SparseOperator>& J,
      const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
      const Teuchos::RCP<LINALG::SparseOperator>& M,
      const NOX::Epetra::Vector& cloneVector,
      const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject);
} // namespace LinearSystem
} // namespace NLN
} // namespace NOX


#endif /* NOX_NLN_LINEARSYSTEM_FACTORY_H_ */
