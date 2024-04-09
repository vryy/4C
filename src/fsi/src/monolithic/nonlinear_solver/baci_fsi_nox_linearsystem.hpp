/*----------------------------------------------------------------------*/
/*! \file

\brief FSI linear system interface to the nonlinear solver NOX

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FSI_NOX_LINEARSYSTEM_HPP
#define FOUR_C_FSI_NOX_LINEARSYSTEM_HPP

#include "baci_config.hpp"

#include <NOX.H>
#include <NOX_Common.H>
#include <NOX_Epetra_Group.H>
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_LinearSystem.H>
#include <NOX_Epetra_Scaling.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Utils.H>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>

#include <vector>

BACI_NAMESPACE_OPEN

namespace CORE::LINALG
{
  class Solver;
}

namespace NOX::FSI
{
  class LinearSystem : public ::NOX::Epetra::LinearSystem
  {
   private:
    enum OperatorType
    {
      EpetraOperator,
      EpetraRowMatrix,
      EpetraVbrMatrix,
      EpetraCrsMatrix,
      SparseMatrix,
      BlockSparseMatrix
    };

   public:
    LinearSystem(Teuchos::ParameterList& printParams,  ///< printing parameters
        Teuchos::ParameterList& linearSolverParams,    ///< parameters for linear solution
        const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>&
            iJac,                                  ///< NOX interface to Jacobian
        const Teuchos::RCP<Epetra_Operator>& J,    ///< the Jacobian or stiffness matrix
        const ::NOX::Epetra::Vector& cloneVector,  ///< initial guess of the solution process
        Teuchos::RCP<CORE::LINALG::Solver>
            structure_solver,  ///< (used-defined) linear algebraic solver
        const Teuchos::RCP<::NOX::Epetra::Scaling> scalingObject =
            Teuchos::null);  ///< scaling of the linear system

    /// provide storage pattern of tangent matrix, i.e. the operator
    OperatorType getOperatorType(const Epetra_Operator& Op);

    ///
    void reset(Teuchos::ParameterList& linearSolverParams);

    /// Applies Jacobian to the given input vector and puts the answer in the result.
    bool applyJacobian(
        const ::NOX::Epetra::Vector& input, ::NOX::Epetra::Vector& result) const override;

    /// Applies Jacobian-Transpose to the given input vector and puts the answer in the result.
    bool applyJacobianTranspose(
        const ::NOX::Epetra::Vector& input, ::NOX::Epetra::Vector& result) const override;

    /// Applies the inverse of the Jacobian matrix to the given input vector and puts the answer
    /// in result.
    bool applyJacobianInverse(Teuchos::ParameterList& params, const ::NOX::Epetra::Vector& input,
        ::NOX::Epetra::Vector& result) override;

    /// Apply right preconditiong to the given input vector.
    bool applyRightPreconditioning(bool useTranspose, Teuchos::ParameterList& params,
        const ::NOX::Epetra::Vector& input, ::NOX::Epetra::Vector& result) const override;

    /// Get the scaling object.
    Teuchos::RCP<::NOX::Epetra::Scaling> getScaling() override;

    /// Sets the diagonal scaling vector(s) used in scaling the linear system.
    void resetScaling(const Teuchos::RCP<::NOX::Epetra::Scaling>& s) override;

    /// Evaluates the Jacobian based on the solution vector x.
    bool computeJacobian(const ::NOX::Epetra::Vector& x) override;

    /// Explicitly constructs a preconditioner based on the solution vector x and the parameter
    /// list p.
    bool createPreconditioner(const ::NOX::Epetra::Vector& x, Teuchos::ParameterList& p,
        bool recomputeGraph) const override;

    /// Deletes the preconditioner.
    bool destroyPreconditioner() const override;

    /// Recalculates the preconditioner using an already allocated graph.
    bool recomputePreconditioner(
        const ::NOX::Epetra::Vector& x, Teuchos::ParameterList& linearSolverParams) const override;

    /// Evaluates the preconditioner policy at the current state.
    PreconditionerReusePolicyType getPreconditionerPolicy(bool advanceReuseCounter = true) override;

    /// Indicates whether a preconditioner has been constructed.
    bool isPreconditionerConstructed() const override;

    /// Indicates whether the linear system has a preconditioner.
    bool hasPreconditioner() const override;

    /// Return Jacobian operator.
    Teuchos::RCP<const Epetra_Operator> getJacobianOperator() const override;

    /// Return Jacobian operator.
    Teuchos::RCP<Epetra_Operator> getJacobianOperator() override;

    /// Return preconditioner operator.
    Teuchos::RCP<const Epetra_Operator> getGeneratedPrecOperator() const override;

    /// Return preconditioner operator.
    Teuchos::RCP<Epetra_Operator> getGeneratedPrecOperator() override;

    /// Set Jacobian operator for solve.
    void setJacobianOperatorForSolve(
        const Teuchos::RCP<const Epetra_Operator>& solveJacOp) override;

    /// Set preconditioner operator for solve.
    void setPrecOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solvePrecOp) override;

   private:
    /// throw an error
    void throwError(const std::string& functionName, const std::string& errorMsg) const;

    ::NOX::Utils utils_;

    Teuchos::RCP<::NOX::Epetra::Interface::Jacobian> jacInterfacePtr_;
    Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner> precInterfacePtr_;
    OperatorType jacType_;
    mutable Teuchos::RCP<Epetra_Operator> jacPtr_;
    mutable Teuchos::RCP<Epetra_Operator> precPtr_;
    Teuchos::RCP<::NOX::Epetra::Scaling> scaling_;
    mutable Teuchos::RCP<::NOX::Epetra::Vector> tmpVectorPtr_;

    bool outputSolveDetails_;
    bool zeroInitialGuess_;
    bool manualScaling_;

    /// index of Newton iteration
    int callcount_;

    /// linear algebraic solver
    Teuchos::RCP<CORE::LINALG::Solver> solver_;

    Teuchos::Time timer_;
  };
}  // namespace NOX::FSI

BACI_NAMESPACE_CLOSE

#endif
