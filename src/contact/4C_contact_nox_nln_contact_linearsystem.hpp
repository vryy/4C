/*----------------------------------------------------------------------*/
/*! \file
\brief Derived class which manages the special requirements to the linear
       solver for contact problems.

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_NOX_NLN_CONTACT_LINEARSYSTEM_HPP
#define FOUR_C_CONTACT_NOX_NLN_CONTACT_LINEARSYSTEM_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_constraint_interface_preconditioner.hpp"
#include "4C_solver_nonlin_nox_constraint_interface_required.hpp"
#include "4C_solver_nonlin_nox_linearsystem.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Mortar
{
  class StrategyBase;
}

namespace NOX
{
  namespace Nln
  {
    namespace CONTACT
    {
      class LinearSystem : public NOX::Nln::LinearSystem
      {
       public:
        //! Standard constructor with full functionality.
        LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams, const SolverMap& solvers,
            const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
            const NOX::Nln::CONSTRAINT::ReqInterfaceMap& iConstr,
            const Teuchos::RCP<Core::LinAlg::SparseOperator>& J,
            const Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>& iPrec,
            const NOX::Nln::CONSTRAINT::PrecInterfaceMap& iConstrPrec,
            const Teuchos::RCP<Core::LinAlg::SparseOperator>& M,
            const ::NOX::Epetra::Vector& cloneVector,
            const Teuchos::RCP<::NOX::Epetra::Scaling> scalingObject);

        //! Constructor without scaling object
        LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams, const SolverMap& solvers,
            const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
            const NOX::Nln::CONSTRAINT::ReqInterfaceMap& iConstr,
            const Teuchos::RCP<Core::LinAlg::SparseOperator>& J,
            const Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>& iPrec,
            const NOX::Nln::CONSTRAINT::PrecInterfaceMap& iConstrPrec,
            const Teuchos::RCP<Core::LinAlg::SparseOperator>& M,
            const ::NOX::Epetra::Vector& cloneVector);

        //! Sets the options of the underlying solver
        Core::LinAlg::SolverParams set_solver_options(Teuchos::ParameterList& p,
            Teuchos::RCP<Core::LinAlg::Solver>& solverPtr,
            const NOX::Nln::SolutionType& solverType) override;

        //! Returns a pointer to linear solver, which has to be used
        NOX::Nln::SolutionType get_active_lin_solver(
            const std::map<NOX::Nln::SolutionType, Teuchos::RCP<Core::LinAlg::Solver>>& solvers,
            Teuchos::RCP<Core::LinAlg::Solver>& currSolver) override;

        //! derived
        void set_linear_problem_for_solve(Epetra_LinearProblem& linear_problem,
            Core::LinAlg::SparseOperator& jac, Epetra_Vector& lhs,
            Epetra_Vector& rhs) const override;

        //! Combine the linear solution parts [derived]
        void complete_solution_after_solve(
            const Epetra_LinearProblem& linProblem, Epetra_Vector& lhs) const override;

       private:
        //! throws an error message
        void throw_error(const std::string& functionName, const std::string& errorMsg) const;

        //! Solve a linear system containing a diagonal matrix
        void apply_diagonal_inverse(
            Core::LinAlg::SparseMatrix& mat, Epetra_Vector& lhs, const Epetra_Vector& rhs) const;

        /// return a pointer to the currently active linear solver
        Teuchos::RCP<Core::LinAlg::Solver> get_linear_contact_solver(
            const std::map<NOX::Nln::SolutionType, Teuchos::RCP<Core::LinAlg::Solver>>& solvers)
            const;

       private:
        struct LinearSubProblem
        {
          LinearSubProblem(const LinearSystem& linsys)
              : linsys_(linsys),
                p_jac_(Teuchos::null),
                p_lhs_(Teuchos::null),
                p_rhs_(Teuchos::null){/* empty */};

          inline void reset()
          {
            p_jac_ = Teuchos::null;
            p_lhs_ = Teuchos::null;
            p_rhs_ = Teuchos::null;
          }

          /** \brief Extract an active linear sub-problem
           *
           *  This routine checks if there is an inactive set of blocks in this
           *  system of equations. Inactive means, that there are only empty off-diagonal
           *  blocks and a diagonal matrix on the diagonal block for a set of
           *  row and the corresponding column blocks. If this is the case the
           *  "active" problem is extracted as a sub-problem and the very simple
           *  "inactive" problem is solved directly by inverting the diagonal matrix.
           *
           *  \author hiermeier \date 04/17 */
          void extract_active_blocks(
              Core::LinAlg::SparseOperator& mat, Epetra_Vector& lhs, Epetra_Vector& rhs);

          /** \brief Set the original linear problem as sub-problem
           *
           *  This is the default case if no simple pseudo problem can be detected.
           *
           *  \author hiermeier \date 04/17 */
          void set_original_system(
              Core::LinAlg::SparseOperator& mat, Epetra_Vector& lhs, Epetra_Vector& rhs);

          /** \brief insert left hand side of the linear sub-problem into the global
           *  left hand side
           *
           *  \note Nothing is happening if the original system is no block
           *  sparse matrix or, in a more general sense, if no linear sub-problem
           *  could be extracted.
           *
           *  \param[out] glhs: global left hand side vector
           *
           *  \author hiermeier \date 03/18 */
          void insert_into_global_lhs(Epetra_Vector& glhs) const;

          const LinearSystem& linsys_;

          Teuchos::RCP<Core::LinAlg::SparseOperator> p_jac_;
          Teuchos::RCP<Epetra_Vector> p_lhs_;
          Teuchos::RCP<Epetra_Vector> p_rhs_;
        };

        //! map of NOX::Nln::CONSTRAINT::Interface::Required objects
        NOX::Nln::CONSTRAINT::ReqInterfaceMap i_constr_;

        //! map of NOX::Nln::CONSTRAINT::Interface::Preconditioner objects
        NOX::Nln::CONSTRAINT::PrecInterfaceMap i_constr_prec_;

        mutable LinearSubProblem p_lin_prob_;
      };  // class LinearSystem
    }     // namespace CONTACT
  }       // namespace Nln
}  // namespace NOX


FOUR_C_NAMESPACE_CLOSE

#endif
