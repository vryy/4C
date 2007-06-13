
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fsi_dyn.H"
#include "fsi_dirichletneumann.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


void fsi_ale_drt()
{
#ifdef PARALLEL
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Teuchos::RefCountPtr<FSI::DirichletNeumannCoupling> fsi = rcp(new FSI::DirichletNeumannCoupling(Comm));

  fsi->Setup();
  fsi->Timeloop(fsi);
}

#endif
#endif
