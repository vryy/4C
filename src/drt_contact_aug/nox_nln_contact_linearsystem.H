/*-----------------------------------------------------------*/
/*!
\file nox_nln_contact_linearsystem.H

\maintainer Michael Hiermeier

\date Jul 14, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_CONTACT_LINEARSYSTEM_H_
#define NOX_NLN_CONTACT_LINEARSYSTEM_H_

#include "../solver_nonlin_nox/nox_nln_linearsystem.H"

// forward declaration
namespace MORTAR {
  class StrategyBase;
}

namespace NOX {
namespace NLN {
namespace CONTACT {
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

    MORTAR::StrategyBase& GetStrategy();

  private:
    //! throws an error message
    void throwError(
        const std::string& functionName,
        const std::string& errorMsg) const;

  private:
    Teuchos::RCP<MORTAR::StrategyBase> strat_;

    bool isMeshtying_;

    bool isContact_;
  };
} // namespace CONTACT
} // namespace NLN
} // namespace NOX


#endif /* NOX_CONTACT_LINEARSYSTEM_H_ */
