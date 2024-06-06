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
  namespace Nln
  {
    namespace CARDIOVASCULAR0D
    {
      class LinearSystem : public NOX::Nln::LinearSystem
      {
       public:
        //! Standard constructor with full functionality.
        LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams,
            const std::map<NOX::Nln::SolutionType, Teuchos::RCP<Core::LinAlg::Solver>>& solvers,
            const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
            const Teuchos::RCP<Core::LinAlg::SparseOperator>& J,
            const Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>& iPrec,
            const Teuchos::RCP<Core::LinAlg::SparseOperator>& M,
            const ::NOX::Epetra::Vector& cloneVector,
            const Teuchos::RCP<::NOX::Epetra::Scaling> scalingObject);

        //! Constructor without scaling object
        LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams,
            const std::map<NOX::Nln::SolutionType, Teuchos::RCP<Core::LinAlg::Solver>>& solvers,
            const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
            const Teuchos::RCP<Core::LinAlg::SparseOperator>& J,
            const Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>& iPrec,
            const Teuchos::RCP<Core::LinAlg::SparseOperator>& M,
            const ::NOX::Epetra::Vector& cloneVector);

        //! Constructor without preconditioner
        LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams,
            const std::map<NOX::Nln::SolutionType, Teuchos::RCP<Core::LinAlg::Solver>>& solvers,
            const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
            const Teuchos::RCP<Core::LinAlg::SparseOperator>& J,
            const ::NOX::Epetra::Vector& cloneVector,
            const Teuchos::RCP<::NOX::Epetra::Scaling> scalingObject);

        //! Constructor without preconditioner and scaling object
        LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams,
            const std::map<NOX::Nln::SolutionType, Teuchos::RCP<Core::LinAlg::Solver>>& solvers,
            const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
            const Teuchos::RCP<Core::LinAlg::SparseOperator>& J,
            const ::NOX::Epetra::Vector& cloneVector);

        //! sets the options of the underlying solver
        Core::LinAlg::SolverParams set_solver_options(Teuchos::ParameterList& p,
            Teuchos::RCP<Core::LinAlg::Solver>& solverPtr,
            const NOX::Nln::SolutionType& solverType) override;

        //! Returns a pointer to the linear solver, which has to be used
        NOX::Nln::SolutionType get_active_lin_solver(
            const std::map<NOX::Nln::SolutionType, Teuchos::RCP<Core::LinAlg::Solver>>& solvers,
            Teuchos::RCP<Core::LinAlg::Solver>& currSolver) override;

       private:
        //! throws an error message
        void throw_error(const std::string& functionName, const std::string& errorMsg) const;

      };  // class LinearSystem
    }     // namespace CARDIOVASCULAR0D
  }       // namespace Nln
}  // namespace NOX


FOUR_C_NAMESPACE_CLOSE

#endif
