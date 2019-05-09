/*----------------------------------------------------------------------*/
/*!
\brief Use #NOX as non-linear solution technique for implicit
       structureal time integration

\maintainer Alexander Popp

\level 3
*/

/*----------------------------------------------------------------------*/
/* macros */

#ifndef STRTIMINT_NOXLINSYS_H
#define STRTIMINT_NOXLINSYS_H

/*----------------------------------------------------------------------*/
/* headers */
#include <vector>

#include <Epetra_Time.h>

#include <NOX_Common.H>
#include <NOX_Epetra_Group.H>
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Epetra_LinearSystem.H>
#include <NOX_Epetra_Scaling.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Utils.H>
#include <NOX.H>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------*/
/* forward declarations */
namespace LINALG
{
  class Solver;
}

namespace NOX
{
  namespace STR
  {
    /// This class enables the structural time integration #STR::TimIntImpl to use
    /// #NOX as non-linear solution technique while preserving the user-defined
    /// linear algebraic solver #LINALG::Solver.
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
          Teuchos::RCP<LINALG::Solver>
              structure_solver,  ///< (used-defined) linear algebraic solver
          const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject = Teuchos::null);

      /// destructor
      ~LinearSystem();

      /// provide storage pattern of tangent matrix, i.e. the operator
      virtual OperatorType getOperatorType(const Epetra_Operator& Op);

      ///
      virtual void reset(Teuchos::ParameterList& linearSolverParams);

      /// Applies Jacobian to the given input vector and puts the answer in the result.
      virtual bool applyJacobian(
          const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const;

      /// Applies Jacobian-Transpose to the given input vector and puts the answer in the result.
      virtual bool applyJacobianTranspose(
          const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const;

      /// Applies the inverse of the Jacobian matrix to the given input vector and puts the answer
      /// in result.
      virtual bool applyJacobianInverse(Teuchos::ParameterList& params,
          const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result);

      /// Apply right preconditiong to the given input vector.
      virtual bool applyRightPreconditioning(bool useTranspose, Teuchos::ParameterList& params,
          const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const;

      /// Get the scaling object.
      virtual Teuchos::RCP<NOX::Epetra::Scaling> getScaling();

      /// Sets the diagonal scaling vector(s) used in scaling the linear system.
      virtual void resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& s);

      /// Evaluates the Jacobian based on the solution vector x.
      virtual bool computeJacobian(const NOX::Epetra::Vector& x);

      /// Explicitly constructs a preconditioner based on the solution vector x and the parameter
      /// list p.
      virtual bool createPreconditioner(
          const NOX::Epetra::Vector& x, Teuchos::ParameterList& p, bool recomputeGraph) const;

      /// Deletes the preconditioner.
      virtual bool destroyPreconditioner() const;

      /// Recalculates the preconditioner using an already allocated graph.
      virtual bool recomputePreconditioner(
          const NOX::Epetra::Vector& x, Teuchos::ParameterList& linearSolverParams) const;

      /// Evaluates the preconditioner policy at the current state.
      virtual PreconditionerReusePolicyType getPreconditionerPolicy(
          bool advanceReuseCounter = true);

      /// Indicates whether a preconditioner has been constructed.
      virtual bool isPreconditionerConstructed() const;

      /// Indicates whether the linear system has a preconditioner.
      virtual bool hasPreconditioner() const;

      /// Return Jacobian operator.
      virtual Teuchos::RCP<const Epetra_Operator> getJacobianOperator() const;

      /// Return Jacobian operator.
      virtual Teuchos::RCP<Epetra_Operator> getJacobianOperator();

      /// Return preconditioner operator.
      virtual Teuchos::RCP<const Epetra_Operator> getGeneratedPrecOperator() const;

      /// Return preconditioner operator.
      virtual Teuchos::RCP<Epetra_Operator> getGeneratedPrecOperator();

      /// Set Jacobian operator for solve.
      virtual void setJacobianOperatorForSolve(
          const Teuchos::RCP<const Epetra_Operator>& solveJacOp);

      /// Set preconditioner operator for solve.
      virtual void setPrecOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solvePrecOp);

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
      Teuchos::RCP<LINALG::Solver> structureSolver_;

      Epetra_Time timer_;
      mutable double timeApplyJacbianInverse_;
    };

  }  // namespace STR
}  // namespace NOX


#endif  //#ifndef STRTIMINT_NOXLINSYS_H
