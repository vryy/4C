/*-----------------------------------------------------------*/
/*! \file

\brief Derived class which manages the special requirements to the linear
       solver for structural-constraint problems.


\date Jul 15, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef CONSTRAINT_NOX_NLN_LAGPENCONSTRAINT_LINEARSYSTEM_H
#define CONSTRAINT_NOX_NLN_LAGPENCONSTRAINT_LINEARSYSTEM_H

#include "solver_nonlin_nox_linearsystem.H"
#include "solver_nonlin_nox_constraint_interface_required.H"
#include "solver_nonlin_nox_constraint_interface_preconditioner.H"

namespace NOX
{
  namespace NLN
  {
    namespace LAGPENCONSTRAINT
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
            const Teuchos::RCP<CORE::LINALG::SparseOperator>& J,
            const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
            const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec,
            const Teuchos::RCP<CORE::LINALG::SparseOperator>& M,
            const NOX::Epetra::Vector& cloneVector,
            const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject);

        //! Constructor without scaling object
        LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams, const SolverMap& solvers,
            const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
            const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr,
            const Teuchos::RCP<CORE::LINALG::SparseOperator>& J,
            const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
            const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec,
            const Teuchos::RCP<CORE::LINALG::SparseOperator>& M,
            const NOX::Epetra::Vector& cloneVector);

        //! Destructor.
        virtual ~LinearSystem(){};

        //! Sets the options of the underlying solver
        virtual void SetSolverOptions(Teuchos::ParameterList& p,
            Teuchos::RCP<CORE::LINALG::Solver>& solverPtr,
            const NOX::NLN::SolutionType& solverType);

        //! Returns a pointer to linear solver, which has to be used
        virtual NOX::NLN::SolutionType GetActiveLinSolver(
            const std::map<NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>& solvers,
            Teuchos::RCP<CORE::LINALG::Solver>& currSolver);

       private:
        //! throws an error message
        void throwError(const std::string& functionName, const std::string& errorMsg) const;

       private:
        //! map of NOX::NLN::CONSTRAINT::Interface::Required objects
        NOX::NLN::CONSTRAINT::ReqInterfaceMap iConstr_;

        //! map of NOX::NLN::CONSTRAINT::Interface::Preconditioner objects
        NOX::NLN::CONSTRAINT::PrecInterfaceMap iConstrPrec_;
      };  // class LinearSystem
    }     // namespace LAGPENCONSTRAINT
  }       // namespace NLN
}  // namespace NOX


#endif  // CONSTRAINT_NOX_NLN_LAGPENCONSTRAINT_LINEARSYSTEM_H
