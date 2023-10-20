/*----------------------------------------------------------------------*/
/*! \file
\brief Use #NOX as non-linear solution technique for implicit
       structureal time integration


\level 3
*/

/*----------------------------------------------------------------------*/
/* macros */

#ifndef BACI_STRUCTURE_TIMINT_NOXLINSYS_H
#define BACI_STRUCTURE_TIMINT_NOXLINSYS_H

/*----------------------------------------------------------------------*/
/* headers */
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

/*----------------------------------------------------------------------*/
/* forward declarations */
namespace CORE::LINALG
{
  class Solver;
}

namespace NOX
{
  namespace STR
  {
    /// This class enables the structural time integration #STR::TimIntImpl to use
    /// #NOX as non-linear solution technique while preserving the user-defined
    /// linear algebraic solver #CORE::LINALG::Solver.
    ///
    /// \author bborn
    /// \date 04/09
    class LinearSystem : public NOX::Epetra::LinearSystem
    {
     protected:
      /// kind of storage and access pattern of tangent matrix
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
      /// constructor
      LinearSystem(Teuchos::ParameterList& printParams,  ///< printing parameters
          Teuchos::ParameterList& linearSolverParams,    ///< parameters for linear solution
          const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>&
              iJac,  ///< NOX interface to Jacobian, i.e. #STR::TimIntImpl
          const Teuchos::RCP<Epetra_Operator>& J,  ///< the Jacobian or stiffness matrix
          const NOX::Epetra::Vector& cloneVector,
          Teuchos::RCP<CORE::LINALG::Solver>
              structure_solver,  ///< (used-defined) linear algebraic solver
          const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject = Teuchos::null);


      /// provide storage pattern of tangent matrix, i.e. the operator
      virtual OperatorType getOperatorType(const Epetra_Operator& Op);

      ///
      virtual void reset(Teuchos::ParameterList& linearSolverParams);

      /// Applies Jacobian to the given input vector and puts the answer in the result.
      bool applyJacobian(
          const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const override;

      /// Applies Jacobian-Transpose to the given input vector and puts the answer in the result.
      bool applyJacobianTranspose(
          const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const override;

      /// Applies the inverse of the Jacobian matrix to the given input vector and puts the answer
      /// in result.
      bool applyJacobianInverse(Teuchos::ParameterList& params, const NOX::Epetra::Vector& input,
          NOX::Epetra::Vector& result) override;

      /// Apply right preconditiong to the given input vector.
      bool applyRightPreconditioning(bool useTranspose, Teuchos::ParameterList& params,
          const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const override;

      /// Get the scaling object.
      Teuchos::RCP<NOX::Epetra::Scaling> getScaling() override;

      /// Sets the diagonal scaling vector(s) used in scaling the linear system.
      void resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& s) override;

      /// Evaluates the Jacobian based on the solution vector x.
      bool computeJacobian(const NOX::Epetra::Vector& x) override;

      /// Explicitly constructs a preconditioner based on the solution vector x and the parameter
      /// list p.
      bool createPreconditioner(const NOX::Epetra::Vector& x, Teuchos::ParameterList& p,
          bool recomputeGraph) const override;

      /// Deletes the preconditioner.
      bool destroyPreconditioner() const override;

      /// Recalculates the preconditioner using an already allocated graph.
      bool recomputePreconditioner(
          const NOX::Epetra::Vector& x, Teuchos::ParameterList& linearSolverParams) const override;

      /// Evaluates the preconditioner policy at the current state.
      PreconditionerReusePolicyType getPreconditionerPolicy(
          bool advanceReuseCounter = true) override;

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

     protected:
      /// throw an error
      virtual void throwError(const std::string& functionName, const std::string& errorMsg) const;

     protected:
      NOX::Utils utils_;

      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> jacInterfacePtr_;
      Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> precInterfacePtr_;
      OperatorType jacType_;
      OperatorType precType_;
      mutable Teuchos::RCP<Epetra_Operator> jacPtr_;
      mutable Teuchos::RCP<Epetra_Operator> precPtr_;
      Teuchos::RCP<NOX::Epetra::Scaling> scaling_;
      mutable Teuchos::RCP<NOX::Epetra::Vector> tmpVectorPtr_;
      mutable double conditionNumberEstimate_;

      bool outputSolveDetails_;
      bool zeroInitialGuess_;
      bool manualScaling_;

      /// index of Newton iteration
      int callcount_;

      /// linear algebraic solver
      Teuchos::RCP<CORE::LINALG::Solver> structureSolver_;

      Teuchos::Time timer_;
      mutable double timeApplyJacbianInverse_;
    };

  }  // namespace STR
}  // namespace NOX


#endif  // STRUCTURE_TIMINT_NOXLINSYS_H
