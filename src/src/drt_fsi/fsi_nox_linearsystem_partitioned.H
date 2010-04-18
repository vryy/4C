#ifdef CCADISCRET

#ifndef FSI_NOX_LINEARSYSTEM_PARTITIONED_H
#define FSI_NOX_LINEARSYSTEM_PARTITIONED_H

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

#include "../drt_adapter/adapter_structure.H"
#include "../drt_adapter/adapter_fluid.H"
#include "../drt_ale/ale.H"

#include "../drt_inpar/inpar_fsi.H"

namespace FSI
{
  class PartitionedMonolithic;
}

namespace NOX {
namespace FSI {

  /// NOX interface to linear partitioned FSI system
  class LinearPartitioned : public NOX::Epetra::Interface::Required
  {
  public:

    LinearPartitioned(::FSI::PartitionedMonolithic& algorithm,
                      ADAPTER::Structure& structurefield,
                      ADAPTER::Fluid& fluidfield,
                      ALE::Ale& alefield);

    /// compute (full) FSI residual
    bool computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

    Teuchos::RCP<Epetra_Vector> StructOp(Teuchos::RCP<const Epetra_Vector> iforce, const FillType fillFlag);

    Teuchos::RCP<Epetra_Vector> FluidOp(Teuchos::RCP<const Epetra_Vector> idisp, const FillType fillFlag);

    void AleOp(Teuchos::RCP<const Epetra_Vector> idisp);

    /// interface velocity calculation given the interface displacements
    Teuchos::RCP<Epetra_Vector> InterfaceVelocity(Teuchos::RCP<const Epetra_Vector> idispnp) const;

    void ExtractResult(Teuchos::RCP<const Epetra_Vector> idisp, Epetra_Vector& result);

  private:

    int callcount_;

    ::FSI::PartitionedMonolithic& algorithm_;

    ADAPTER::Structure& structurefield_;
    ADAPTER::Fluid& fluidfield_;
    ALE::Ale& alefield_;

    Teuchos::RCP<const Epetra_Vector> srhs_;
    Teuchos::RCP<const Epetra_Vector> frhs_;
    Teuchos::RCP<const Epetra_Vector> arhs_;

    Teuchos::RCP<Epetra_Vector> sx_;
    Teuchos::RCP<Epetra_Vector> fx_;
    Teuchos::RCP<Epetra_Vector> ax_;

    Teuchos::RCP<Epetra_Vector> stmp_;
    Teuchos::RCP<Epetra_Vector> ftmp_;
    Teuchos::RCP<Epetra_Vector> atmp_;

    Teuchos::RCP<LINALG::SparseMatrix> s_;
    Teuchos::RCP<LINALG::SparseMatrix> f_;
    Teuchos::RCP<LINALG::SparseMatrix> a_;

    Teuchos::RCP<LINALG::Solver> slin_;
    Teuchos::RCP<LINALG::Solver> flin_;
    Teuchos::RCP<LINALG::Solver> alin_;

    Teuchos::RCP<Epetra_Vector> idispn_;
    Teuchos::RCP<Epetra_Vector> iveln_;

    Teuchos::RCP<LINALG::SparseMatrix> fluiddirichlet_;
  };

  /// Linear partitioned solver build on linear field solvers
  class LinearPartitionedSolver : public NOX::Epetra::LinearSystem
  {
  public:

    LinearPartitionedSolver(
      Teuchos::ParameterList& printParams,
      Teuchos::ParameterList& linearSolverParams,
      const LINALG::MultiMapExtractor& extractor,
      ::FSI::PartitionedMonolithic& algorithm,
      ADAPTER::Structure& structurefield,
      ADAPTER::Fluid& fluidfield,
      ALE::Ale& alefield,
      INPAR::FSI::LinearBlockSolver linearsolverstrategy);

    ~LinearPartitionedSolver();

    virtual void reset(Teuchos::ParameterList& linearSolverParams);

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

    virtual Teuchos::RCP<NOX::Epetra::Scaling> getScaling();

    virtual void resetScaling(const Teuchos::RCP< NOX::Epetra::Scaling>& s);

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

    virtual Teuchos::RCP<const Epetra_Operator> getJacobianOperator() const;

    virtual Teuchos::RCP<Epetra_Operator> getJacobianOperator();


    virtual Teuchos::RCP<const Epetra_Operator> getGeneratedPrecOperator() const;

    virtual Teuchos::RCP<Epetra_Operator> getGeneratedPrecOperator();

    virtual void setJacobianOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solveJacOp);


    virtual void setPrecOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solvePrecOp);

  private:

    void LinearPartitionedSolve(NOX::Epetra::Vector& result,
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

    bool outputSolveDetails_;
    bool zeroInitialGuess_;

    const LINALG::MultiMapExtractor& extractor_;
    ::FSI::PartitionedMonolithic& algorithm_;

    ADAPTER::Structure& structurefield_;
    ADAPTER::Fluid& fluidfield_;
    ALE::Ale& alefield_;

    INPAR::FSI::LinearBlockSolver linearsolverstrategy_;
  };

}
}


#endif
#endif
