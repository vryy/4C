/*-----------------------------------------------------------*/
/*!
\file nox_nln_str_linearsystem.H

\maintainer Michael Hiermeier

\date Aug 7, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_STR_LINEARSYSTEM_H_
#define NOX_STR_LINEARSYSTEM_H_

#include "../solver_nonlin_nox/nox_nln_linearsystem.H"

namespace NOX {
namespace NLN {
namespace STR {
  class LinearSystem : public NOX::NLN::LinearSystem
  {
  public:
    //! Constructor.
    LinearSystem(Teuchos::ParameterList& printParams,
        Teuchos::ParameterList& linearSolverParams,
        const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& solvers,
        const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
        const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
        const Teuchos::RCP<LINALG::SparseOperator>& J,
        const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
        const Teuchos::RCP<LINALG::SparseOperator>& M,
        const NOX::Epetra::Vector& cloneVector,
        const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject = Teuchos::null);

    //! Destructor.
    virtual ~LinearSystem(){};

    //! sets the options of the underlying solver
    virtual void SetSolverOptions(Teuchos::RCP<LINALG::Solver>& solverPtr,
        const NOX::NLN::SolutionType& solverType);

    //! Returns a pointer to linear solver, which has to be used
    virtual NOX::NLN::SolutionType GetCurrentLinSolver(
        const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& solvers,
        Teuchos::RCP<LINALG::Solver>& currSolver);

  private:
    //! throws an error message
    void throwError(
        const std::string& functionName,
        const std::string& errorMsg) const;

  };
} // namespace STR
} // namespace NLN
} // namespace NOX


#endif /* NOX_STR_LINEARSYSTEM_H_ */
