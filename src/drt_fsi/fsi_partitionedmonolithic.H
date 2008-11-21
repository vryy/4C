#ifdef CCADISCRET

#ifndef FSI_PARTITIONEDMONOLITHIC_H
#define FSI_PARTITIONEDMONOLITHIC_H

#include "../drt_inpar/inpar_fsi.H"
#include "fsi_monolithic.H"

namespace NOX
{
  namespace FSI
  {
    class LinearPartitioned;
  }
}

namespace FSI
{

  /// Monolithic FSI based on partitioned (interface based) linear solver
  /*!

    A monolithic FSI solver (a Newton iteration on the fully coupled
    equations) with an partitioned linear solver internally. That is we do not
    need a block matrix here. The linear systems are solved in a
    Schur-complement fashion just like the nonlinear systems in the real
    partitioned schemes are. The same solution schemes apply.

    The rhs is build just like we do it for fluid-split overlapping MFSI.

    \author u.kue
    \date 10/08
   */
  class PartitionedMonolithic : public Monolithic
  {
    friend class NOX::FSI::LinearPartitioned;

  public:
    explicit PartitionedMonolithic(Epetra_Comm& comm);

    virtual void SetupRHS(Epetra_Vector&, bool);

    virtual void SetupSystemMatrix() {}

    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector>);

    virtual Teuchos::RCP<NOX::Epetra::LinearSystem> CreateLinearSystem(Teuchos::ParameterList&, NOX::Epetra::Vector&, Teuchos::RCP<NOX::Utils>);

    virtual Teuchos::RCP<NOX::StatusTest::Combo> CreateStatusTest(Teuchos::ParameterList&, Teuchos::RCP<NOX::Epetra::Group>);

    virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector>, Teuchos::RCP<const Epetra_Vector>&, Teuchos::RCP<const Epetra_Vector>&, Teuchos::RCP<const Epetra_Vector>&);

  private:

    /// build block vector from field vectors
    void SetupVector(Epetra_Vector &f,
                     Teuchos::RCP<const Epetra_Vector> sv,
                     Teuchos::RCP<const Epetra_Vector> fv,
                     Teuchos::RCP<const Epetra_Vector> av,
                     double fluidscale);

    /// coupling of fluid and ale (interface only)
    ADAPTER::Coupling icoupfa_;

    /// preconditioned block Krylov or block Gauss-Seidel linear solver
    INPAR::FSI::LinearBlockSolver linearsolverstrategy_;
  };
}

#endif
#endif
