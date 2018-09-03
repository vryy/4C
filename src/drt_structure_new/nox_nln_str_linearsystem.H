/*-----------------------------------------------------------*/
/*!
\file nox_nln_str_linearsystem.H

\brief Derived class which manages the special requirements to the linear
       solver for structural problems.

\maintainer Michael Hiermeier

\date Aug 7, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_STR_LINEARSYSTEM_H_
#define NOX_STR_LINEARSYSTEM_H_

#include "../solver_nonlin_nox/nox_nln_linearsystem.H"

namespace NOX
{
  namespace NLN
  {
    namespace STR
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

      };  // class LinearSystem
    }     // namespace STR
  }       // namespace NLN
}  // namespace NOX


#endif /* NOX_STR_LINEARSYSTEM_H_ */
