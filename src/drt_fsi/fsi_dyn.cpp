
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fsi_dyn.H"
#include "fsi_dirichletneumann.H"

#include "../drt_lib/drt_resulttest.H"
#include "../drt_fluid/fluidresulttest.H"

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
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  Teuchos::RefCountPtr<FSI::DirichletNeumannCoupling> fsi = rcp(new FSI::DirichletNeumannCoupling(comm));

  fsi->Timeloop(fsi);

#ifdef RESULTTEST
  DRT::ResultTestManager testmanager(comm);
  testmanager.AddFieldTest(rcp(new FluidResultTest(fsi->FluidField())));
  testmanager.TestAll();
#endif
}

#endif
#endif
