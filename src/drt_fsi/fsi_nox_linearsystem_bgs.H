#ifdef CCADISCRET

#ifndef FSI_NOX_LINEARSYSTEM_BGS_H
#define FSI_NOX_LINEARSYSTEM_BGS_H

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
#include <Teuchos_RefCountPtr.hpp>

#include "../drt_lib/linalg_systemmatrix.H"

class SparseMatrix;
class BlockSparseMatrix;
class LINALG::Solver;

namespace NOX {
namespace FSI {

  /// NOX interface to linear block Gauss-Seidel FSI system
  class LinearBGS : public NOX::Epetra::Interface::Required
  {
  public:

    LinearBGS(const LINALG::BlockSparseMatrixBase& A,
              const Epetra_Vector &y,
              Teuchos::RCP<LINALG::Solver> structure_solver,
              Teuchos::RCP<LINALG::Solver> fluid_solver,
              Teuchos::RCP<LINALG::Solver> ale_solver);

    /// compute (full) FSI residual
    bool computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

  private:

    const LINALG::BlockSparseMatrixBase& A_;
    const Epetra_Vector &y_;

    int callcount_;

    Teuchos::RCP<LINALG::Solver> structureSolver_;
    Teuchos::RCP<LINALG::Solver> fluidSolver_;
    Teuchos::RCP<LINALG::Solver> aleSolver_;

    Teuchos::RCP<Epetra_Vector> sy_;
    Teuchos::RCP<Epetra_Vector> fy_;
    Teuchos::RCP<Epetra_Vector> ay_;

    Teuchos::RCP<Epetra_Vector> tmpsx_;
    Teuchos::RCP<Epetra_Vector> tmpfx_;
    Teuchos::RCP<Epetra_Vector> tmpax_;

    Teuchos::RCP<Epetra_Vector> tmpsy_;
    Teuchos::RCP<Epetra_Vector> tmpfy_;
    Teuchos::RCP<Epetra_Vector> tmpay_;
  };


  /// Linear Block-Gauss-Seidel solver build on linear field solvers
  class LinearBGSSolver : public NOX::Epetra::LinearSystem
  {
  protected:

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

    LinearBGSSolver(
      Teuchos::ParameterList& printParams,
      Teuchos::ParameterList& linearSolverParams,
      const Teuchos::RefCountPtr< NOX::Epetra::Interface::Jacobian>& iJac,
      const Teuchos::RefCountPtr<Epetra_Operator>& J,
      const NOX::Epetra::Vector& cloneVector,
      Teuchos::RCP < LINALG::Solver > structure_solver,
      Teuchos::RCP < LINALG::Solver > fluid_solver,
      Teuchos::RCP < LINALG::Solver > ale_solver,
      const Teuchos::RefCountPtr< NOX::Epetra::Scaling> scalingObject = Teuchos::null);

    ~LinearBGSSolver();

    virtual void reset(Teuchos::ParameterList& linearSolverParams);
    virtual OperatorType getOperatorType(const Epetra_Operator& o);

    virtual bool applyJacobian(const NOX::Epetra::Vector& input,
                               NOX::Epetra::Vector& result) const;

    virtual bool applyJacobianTranspose(const NOX::Epetra::Vector& input,
                                        NOX::Epetra::Vector& result) const;

    virtual bool applyJacobianInverse(Teuchos::ParameterList &params,
                                      const NOX::Epetra::Vector &input,
                                      NOX::Epetra::Vector &result);

    virtual bool applyRightPreconditioning(bool useTranspose,
                                           Teuchos::ParameterList& params,
                                           const NOX::Epetra::Vector& input,
                                           NOX::Epetra::Vector& result) const;

    virtual Teuchos::RefCountPtr< NOX::Epetra::Scaling> getScaling();

    virtual void resetScaling(const Teuchos::RefCountPtr< NOX::Epetra::Scaling>& s);

    virtual bool computeJacobian(const NOX::Epetra::Vector& x);

    virtual bool createPreconditioner(const NOX::Epetra::Vector& x,
                                      Teuchos::ParameterList& p,
                                      bool recomputeGraph) const;

    virtual bool destroyPreconditioner() const;

    virtual bool recomputePreconditioner(const NOX::Epetra::Vector& x,
                                         Teuchos::ParameterList& linearSolverParams) const;

    virtual PreconditionerReusePolicyType getPreconditionerPolicy(bool advanceReuseCounter=true);

    virtual bool isPreconditionerConstructed() const;

    virtual bool hasPreconditioner() const;

    virtual Teuchos::RefCountPtr<const Epetra_Operator> getJacobianOperator() const;

    virtual Teuchos::RefCountPtr<Epetra_Operator> getJacobianOperator();


    virtual Teuchos::RefCountPtr<const Epetra_Operator> getGeneratedPrecOperator() const;

    virtual Teuchos::RefCountPtr<Epetra_Operator> getGeneratedPrecOperator();

    virtual void setJacobianOperatorForSolve(const Teuchos::RefCountPtr<const Epetra_Operator>& solveJacOp);


    virtual void setPrecOperatorForSolve(const Teuchos::RefCountPtr<const Epetra_Operator>& solvePrecOp);

  private:

    void bgs(const LINALG::BlockSparseMatrixBase& A,
             NOX::Epetra::Vector& result,
             const NOX::Epetra::Vector& input,
             int& maxit,
             double& tol);

    Teuchos::RCP<NOX::StatusTest::Combo> CreateStatusTest(ParameterList& nlParams,
                                                          Teuchos::RCP<NOX::Epetra::Group> grp);

    void CreateStatusTest(ParameterList& nlParams,
                          Teuchos::RCP<NOX::Epetra::Group> grp,
                          Teuchos::RCP<NOX::StatusTest::Combo> converged);

  protected:

    virtual void throwError(const string& functionName,
                            const string& errorMsg) const;

  protected:

    NOX::Utils utils_;

    Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> jacInterfacePtr_;
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Preconditioner> precInterfacePtr_;
    OperatorType jacType_;
    OperatorType precType_;
    mutable Teuchos::RefCountPtr<Epetra_Operator> jacPtr_;
    mutable Teuchos::RefCountPtr<Epetra_Operator> precPtr_;
    Teuchos::RefCountPtr<NOX::Epetra::Scaling> scaling_;
    mutable Teuchos::RefCountPtr<NOX::Epetra::Vector> tmpVectorPtr_;
    mutable double conditionNumberEstimate_;

    bool outputSolveDetails_;
    bool zeroInitialGuess_;
    bool manualScaling_;

    Teuchos::RCP<LINALG::Solver> structureSolver_;
    Teuchos::RCP<LINALG::Solver> fluidSolver_;
    Teuchos::RCP<LINALG::Solver> aleSolver_;

    Epetra_Time timer_;
    mutable double timeApplyJacbianInverse_;
  };

}
}


#endif
#endif
