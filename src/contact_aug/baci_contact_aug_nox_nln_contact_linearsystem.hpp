/*----------------------------------------------------------------------*/
/*! \file
\brief Derived class which manages the special requirements to the linear
       solver for contact problems.

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_AUG_NOX_NLN_CONTACT_LINEARSYSTEM_HPP
#define FOUR_C_CONTACT_AUG_NOX_NLN_CONTACT_LINEARSYSTEM_HPP

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
    namespace CONTACT
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

        //! derived
        void SetLinearProblemForSolve(Epetra_LinearProblem& linear_problem,
            CORE::LINALG::SparseOperator& jac, Epetra_Vector& lhs,
            Epetra_Vector& rhs) const override;

        //! Combine the linear solution parts [derived]
        void CompleteSolutionAfterSolve(
            const Epetra_LinearProblem& linProblem, Epetra_Vector& lhs) const override;

       private:
        //! throws an error message
        void throwError(const std::string& functionName, const std::string& errorMsg) const;

        //! Solve a linear system containing a diagonal matrix
        void applyDiagonalInverse(
            CORE::LINALG::SparseMatrix& mat, Epetra_Vector& lhs, const Epetra_Vector& rhs) const;

        /// return a pointer to the currently active linear solver
        Teuchos::RCP<CORE::LINALG::Solver> GetLinearContactSolver(
            const std::map<NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>& solvers)
            const;

       private:
        struct LinearSubProblem
        {
          LinearSubProblem(const LinearSystem& linsys)
              : linsys_(linsys),
                p_jac_(Teuchos::null),
                p_lhs_(Teuchos::null),
                p_rhs_(Teuchos::null){/* empty */};

          inline void Reset()
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
          void ExtractActiveBlocks(
              CORE::LINALG::SparseOperator& mat, Epetra_Vector& lhs, Epetra_Vector& rhs);

          /** \brief Set the original linear problem as sub-problem
           *
           *  This is the default case if no simple pseudo problem can be detected.
           *
           *  \author hiermeier \date 04/17 */
          void SetOriginalSystem(
              CORE::LINALG::SparseOperator& mat, Epetra_Vector& lhs, Epetra_Vector& rhs);

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
          void InsertIntoGlobalLhs(Epetra_Vector& glhs) const;

          const LinearSystem& linsys_;

          Teuchos::RCP<CORE::LINALG::SparseOperator> p_jac_;
          Teuchos::RCP<Epetra_Vector> p_lhs_;
          Teuchos::RCP<Epetra_Vector> p_rhs_;
        };

        //! map of NOX::NLN::CONSTRAINT::Interface::Required objects
        NOX::NLN::CONSTRAINT::ReqInterfaceMap iConstr_;

        //! map of NOX::NLN::CONSTRAINT::Interface::Preconditioner objects
        NOX::NLN::CONSTRAINT::PrecInterfaceMap iConstrPrec_;

        mutable LinearSubProblem p_lin_prob_;
      };  // class LinearSystem
    }     // namespace CONTACT
  }       // namespace NLN
}  // namespace NOX


FOUR_C_NAMESPACE_CLOSE

#endif
