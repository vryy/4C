/*-----------------------------------------------------------*/
/*! \file

\brief Derived class which manages the special requirements to the linear
       solver for 0D cardiovascular-structural problems.


\date Jul 15, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef CARDIOVASCULAR0D_NOX_NLN_LINEARSYSTEM_H
#define CARDIOVASCULAR0D_NOX_NLN_LINEARSYSTEM_H

#include "solver_nonlin_nox_linearsystem.H"

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
            const std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& solvers,
            const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
            const Teuchos::RCP<LINALG::SparseOperator>& J,
            const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
            const Teuchos::RCP<LINALG::SparseOperator>& M, const NOX::Epetra::Vector& cloneVector,
            const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject);

        //! Constructor without scaling object
        LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams,
            const std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& solvers,
            const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
            const Teuchos::RCP<LINALG::SparseOperator>& J,
            const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
            const Teuchos::RCP<LINALG::SparseOperator>& M, const NOX::Epetra::Vector& cloneVector);

        //! Constructor without preconditioner
        LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams,
            const std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& solvers,
            const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
            const Teuchos::RCP<LINALG::SparseOperator>& J, const NOX::Epetra::Vector& cloneVector,
            const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject);

        //! Constructor without preconditioner and scaling object
        LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams,
            const std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& solvers,
            const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
            const Teuchos::RCP<LINALG::SparseOperator>& J, const NOX::Epetra::Vector& cloneVector);

        //! Destructor.
        virtual ~LinearSystem(){};

        //! sets the options of the underlying solver
        virtual void SetSolverOptions(Teuchos::ParameterList& p,
            Teuchos::RCP<LINALG::Solver>& solverPtr, const NOX::NLN::SolutionType& solverType);

        //! Returns a pointer to the linear solver, which has to be used
        virtual NOX::NLN::SolutionType GetActiveLinSolver(
            const std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& solvers,
            Teuchos::RCP<LINALG::Solver>& currSolver);

       private:
        //! throws an error message
        void throwError(const std::string& functionName, const std::string& errorMsg) const;

      };  // class LinearSystem
    }     // namespace CARDIOVASCULAR0D
  }       // namespace NLN
}  // namespace NOX


#endif  // CARDIOVASCULAR0D_NOX_NLN_LINEARSYSTEM_H
