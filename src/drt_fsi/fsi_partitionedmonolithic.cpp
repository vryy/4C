#ifdef CCADISCRET

#include "fsi_partitionedmonolithic.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::PartitionedMonolithic::PartitionedMonolithic(Epetra_Comm& comm)
  : Monolithic(comm)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::PartitionedMonolithic::SetupRHS(Epetra_Vector&, bool)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::PartitionedMonolithic::InitialGuess(Teuchos::RCP<Epetra_Vector>)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem> FSI::PartitionedMonolithic::CreateLinearSystem(Teuchos::ParameterList&, NOX::Epetra::Vector&, Teuchos::RCP<NOX::Utils>)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::StatusTest::Combo> FSI::PartitionedMonolithic::CreateStatusTest(Teuchos::ParameterList&, Teuchos::RCP<NOX::Epetra::Group>)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::PartitionedMonolithic::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector>, Teuchos::RCP<const Epetra_Vector>&, Teuchos::RCP<const Epetra_Vector>&, Teuchos::RCP<const Epetra_Vector>&)
{
}

#endif
