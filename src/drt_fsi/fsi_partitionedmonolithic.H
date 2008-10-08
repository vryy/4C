#ifdef CCADISCRET

#ifndef FSI_PARTITIONEDMONOLITHIC_H
#define FSI_PARTITIONEDMONOLITHIC_H

#include "fsi_monolithic.H"

namespace FSI
{

  /// Monolithic FSI based on partitioned (interface based) linear solver
  class PartitionedMonolithic : public Monolithic
  {
  public:
    explicit PartitionedMonolithic(Epetra_Comm& comm);

    virtual void SetupRHS(Epetra_Vector&, bool);

    virtual void SetupSystemMatrix() {}

    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector>);

    virtual Teuchos::RCP<NOX::Epetra::LinearSystem> CreateLinearSystem(Teuchos::ParameterList&, NOX::Epetra::Vector&, Teuchos::RCP<NOX::Utils>);

    virtual Teuchos::RCP<NOX::StatusTest::Combo> CreateStatusTest(Teuchos::ParameterList&, Teuchos::RCP<NOX::Epetra::Group>);

    virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector>, Teuchos::RCP<const Epetra_Vector>&, Teuchos::RCP<const Epetra_Vector>&, Teuchos::RCP<const Epetra_Vector>&);

  };
}

#endif
#endif
