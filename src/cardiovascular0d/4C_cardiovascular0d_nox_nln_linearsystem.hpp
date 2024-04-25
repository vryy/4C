/*-----------------------------------------------------------*/
/*! \file

\brief Derived class which manages the special requirements to the linear
       solver for 0D cardiovascular-structural problems.


\date Jul 15, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_CARDIOVASCULAR0D_NOX_NLN_LINEARSYSTEM_HPP
#define FOUR_C_CARDIOVASCULAR0D_NOX_NLN_LINEARSYSTEM_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_linearsystem.hpp"

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace NLN
  {
    namespace CARDIOVASCULAR0D
    {
      class LinearSystem : public NOX::NLN::LinearSystem
      {
       public:
        //! Standard constructor with full functionality.
        LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams,
            const std::map<NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>& solvers,
            const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
            const Teuchos::RCP<CORE::LINALG::SparseOperator>& J,
            const Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>& iPrec,
            const Teuchos::RCP<CORE::LINALG::SparseOperator>& M,
            const ::NOX::Epetra::Vector& cloneVector,
            const Teuchos::RCP<::NOX::Epetra::Scaling> scalingObject);

        //! Constructor without scaling object
        LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams,
            const std::map<NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>& solvers,
            const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
            const Teuchos::RCP<CORE::LINALG::SparseOperator>& J,
            const Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>& iPrec,
            const Teuchos::RCP<CORE::LINALG::SparseOperator>& M,
            const ::NOX::Epetra::Vector& cloneVector);

        //! Constructor without preconditioner
        LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams,
            const std::map<NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>& solvers,
            const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
            const Teuchos::RCP<CORE::LINALG::SparseOperator>& J,
            const ::NOX::Epetra::Vector& cloneVector,
            const Teuchos::RCP<::NOX::Epetra::Scaling> scalingObject);

        //! Constructor without preconditioner and scaling object
        LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams,
            const std::map<NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>& solvers,
            const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
            const Teuchos::RCP<CORE::LINALG::SparseOperator>& J,
            const ::NOX::Epetra::Vector& cloneVector);

        //! sets the options of the underlying solver
        CORE::LINALG::SolverParams SetSolverOptions(Teuchos::ParameterList& p,
            Teuchos::RCP<CORE::LINALG::Solver>& solverPtr,
            const NOX::NLN::SolutionType& solverType) override;

        //! Returns a pointer to the linear solver, which has to be used
        NOX::NLN::SolutionType GetActiveLinSolver(
            const std::map<NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>& solvers,
            Teuchos::RCP<CORE::LINALG::Solver>& currSolver) override;

       private:
        //! throws an error message
        void throwError(const std::string& functionName, const std::string& errorMsg) const;

      };  // class LinearSystem
    }     // namespace CARDIOVASCULAR0D
  }       // namespace NLN
}  // namespace NOX


FOUR_C_NAMESPACE_CLOSE

#endif
