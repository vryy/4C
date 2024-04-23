/*----------------------------------------------------------------------*/
/*! \file
\brief Derived class which manages the special requirements to the linear
       solver for mesh tying problems.

\level 3


*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_NOX_NLN_MESHTYING_LINEARSYSTEM_HPP
#define FOUR_C_CONTACT_NOX_NLN_MESHTYING_LINEARSYSTEM_HPP

#include "baci_config.hpp"

#include "baci_solver_nonlin_nox_constraint_interface_preconditioner.hpp"
#include "baci_solver_nonlin_nox_constraint_interface_required.hpp"
#include "baci_solver_nonlin_nox_linearsystem.hpp"

FOUR_C_NAMESPACE_OPEN

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
            const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
            const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr,
            const Teuchos::RCP<CORE::LINALG::SparseOperator>& J,
            const Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>& iPrec,
            const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec,
            const Teuchos::RCP<CORE::LINALG::SparseOperator>& M,
            const ::NOX::Epetra::Vector& cloneVector,
            const Teuchos::RCP<::NOX::Epetra::Scaling> scalingObject);

        //! Constructor without scaling object
        LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams, const SolverMap& solvers,
            const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
            const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr,
            const Teuchos::RCP<CORE::LINALG::SparseOperator>& J,
            const Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>& iPrec,
            const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec,
            const Teuchos::RCP<CORE::LINALG::SparseOperator>& M,
            const ::NOX::Epetra::Vector& cloneVector);

        //! Sets the options of the underlying solver
        CORE::LINALG::SolverParams SetSolverOptions(Teuchos::ParameterList& p,
            Teuchos::RCP<CORE::LINALG::Solver>& solverPtr,
            const NOX::NLN::SolutionType& solverType) override;

        //! Returns a pointer to linear solver, which has to be used
        NOX::NLN::SolutionType GetActiveLinSolver(
            const std::map<NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>& solvers,
            Teuchos::RCP<CORE::LINALG::Solver>& currSolver) override;

       private:
        //! throws an error message
        void throwError(const std::string& functionName, const std::string& errorMsg) const;

       private:
        //! map of NOX::NLN::CONSTRAINT::Interface::Required objects
        NOX::NLN::CONSTRAINT::ReqInterfaceMap i_constr_;

        //! map of NOX::NLN::CONSTRAINT::Interface::Preconditioner objects
        NOX::NLN::CONSTRAINT::PrecInterfaceMap i_constr_prec_;

      };  // class LinearSystem
    }     // namespace MESHTYING
  }       // namespace NLN
}  // namespace NOX


FOUR_C_NAMESPACE_CLOSE

#endif
