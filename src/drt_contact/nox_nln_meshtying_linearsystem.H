/*----------------------------------------------------------------------*/
/*! \file
\brief Derived class which manages the special requirements to the linear
       solver for mesh tying problems.

\level 3

\maintainer Matthias Mayr

*/
/*----------------------------------------------------------------------*/
#ifndef NOX_NLN_MESHTYING_LINEARSYSTEM_H_
#define NOX_NLN_MESHTYING_LINEARSYSTEM_H_

#include "../solver_nonlin_nox/nox_nln_linearsystem.H"
#include "../solver_nonlin_nox/nox_nln_constraint_interface_required.H"
#include "../solver_nonlin_nox/nox_nln_constraint_interface_preconditioner.H"

// forward declaration
namespace MORTAR
{
  class StrategyBase;
}

namespace NOX
{
  namespace NLN
  {
    namespace MESHTYING
    {
      class LinearSystem : public NOX::NLN::LinearSystem
      {
       public:
        //! Standard constructor with full functionality.
        LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams, const SolverMap& solvers,
            const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
            const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr,
            const Teuchos::RCP<LINALG::SparseOperator>& J,
            const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
            const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec,
            const Teuchos::RCP<LINALG::SparseOperator>& M, const NOX::Epetra::Vector& cloneVector,
            const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject);

        //! Constructor without scaling object
        LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams, const SolverMap& solvers,
            const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
            const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr,
            const Teuchos::RCP<LINALG::SparseOperator>& J,
            const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
            const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec,
            const Teuchos::RCP<LINALG::SparseOperator>& M, const NOX::Epetra::Vector& cloneVector);

        //! Destructor.
        virtual ~LinearSystem(){};

        //! Sets the options of the underlying solver
        virtual void SetSolverOptions(Teuchos::ParameterList& p,
            Teuchos::RCP<LINALG::Solver>& solverPtr, const NOX::NLN::SolutionType& solverType);

        //! Returns a pointer to linear solver, which has to be used
        virtual NOX::NLN::SolutionType GetActiveLinSolver(
            const std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& solvers,
            Teuchos::RCP<LINALG::Solver>& currSolver);

       private:
        //! throws an error message
        void throwError(const std::string& functionName, const std::string& errorMsg) const;

       private:
        //! map of NOX::NLN::CONSTRAINT::Interface::Required objects
        NOX::NLN::CONSTRAINT::ReqInterfaceMap iConstr_;

        //! map of NOX::NLN::CONSTRAINT::Interface::Preconditioner objects
        NOX::NLN::CONSTRAINT::PrecInterfaceMap iConstrPrec_;

      };  // class LinearSystem
    }     // namespace MESHTYING
  }       // namespace NLN
}  // namespace NOX


#endif /* NOX_MESHTYING_LINEARSYSTEM_H_ */
