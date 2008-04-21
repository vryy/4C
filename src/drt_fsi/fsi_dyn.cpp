
#ifdef CCADISCRET

#include <string>
#include <vector>
#include <set>
#include <functional>

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "fsi_dyn.H"
#include "fsi_dirichletneumann.H"
#include "fsi_monolithicoverlap.H"
#include "fsi_structureale.H"
#include "fsi_fluid_ale.H"
#include "fsi_utils.H"

#include "../drt_lib/drt_resulttest.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "../drt_lib/drt_globalproblem.H"

#include "fsi_create_ale.H"
#include "fsi_create_boundary.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*/
// entry point for Fluid on Ale in DRT
/*----------------------------------------------------------------------*/
void fluid_ale_drt()
{
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  RefCountPtr<DRT::Discretization> aledis = DRT::Problem::Instance()->Dis(genprob.numaf,0);
  if (!aledis->Filled()) aledis->FillComplete();

  // create ale elements if the ale discretization is empty
  if (aledis->NumGlobalNodes()==0)
    CreateAleDiscretization();

  Teuchos::RCP<FSI::FluidAleAlgorithm> fluid = Teuchos::rcp(new FSI::FluidAleAlgorithm(comm));

  fluid->Timeloop();

  DRT::ResultTestManager testmanager(comm);
  testmanager.AddFieldTest(fluid->MBFluidField().CreateFieldTest());
  testmanager.TestAll();
}


/*----------------------------------------------------------------------*/
// entry point for (pure) free surface in DRT
/*----------------------------------------------------------------------*/
void fluid_freesurf_drt()
{
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  RefCountPtr<DRT::Discretization> aledis = DRT::Problem::Instance()->Dis(genprob.numaf,0);
  if (!aledis->Filled()) aledis->FillComplete();

  // create ale elements if the ale discretization is empty
  if (aledis->NumGlobalNodes()==0)
    CreateAleDiscretization();

  Teuchos::RCP<FSI::FluidAleAlgorithm> fluid = Teuchos::rcp(new FSI::FluidAleAlgorithm(comm));

  fluid->Timeloop();

  DRT::ResultTestManager testmanager(comm);
  testmanager.AddFieldTest(fluid->MBFluidField().CreateFieldTest());
  testmanager.TestAll();
}


/*----------------------------------------------------------------------*/
// entry point for FSI in DRT
/*----------------------------------------------------------------------*/
void fsi_ale_drt()
{
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  RefCountPtr<DRT::Discretization> aledis = DRT::Problem::Instance()->Dis(genprob.numaf,0);
  if (!aledis->Filled()) aledis->FillComplete();

  // create ale elements if the ale discretization is empty
  if (aledis->NumGlobalNodes()==0)
    CreateAleDiscretization();

  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  if ((Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO") != fsi_iter_monolithic)&&(Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO") != fsi_pseudo_structureale))
  {
    Teuchos::RCP<FSI::DirichletNeumannCoupling> fsi = Teuchos::rcp(new FSI::DirichletNeumannCoupling(comm));

    if (genprob.restart)
    {
      // read the restart information, set vectors and variables
      fsi->ReadRestart(genprob.restart);
    }

    fsi->Timeloop(fsi);
    DRT::ResultTestManager testmanager(comm);
    testmanager.AddFieldTest(fsi->MBFluidField().CreateFieldTest());
    testmanager.AddFieldTest(fsi->StructureField().CreateFieldTest());
    testmanager.TestAll();
  }
  else if (Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO") != fsi_pseudo_structureale)
  {
    Teuchos::RCP<FSI::MonolithicOverlap> fsi = Teuchos::rcp(new FSI::MonolithicOverlap(comm));

    if (genprob.restart)
    {
      // read the restart information, set vectors and variables
      fsi->ReadRestart(genprob.restart);
    }
    
    fsi->Timeloop(fsi);

    DRT::ResultTestManager testmanager(comm);
    testmanager.AddFieldTest(fsi->FluidField().CreateFieldTest());
    testmanager.AddFieldTest(fsi->StructureField().CreateFieldTest());
    testmanager.TestAll();
  }
  else 
  {
    Teuchos::RCP<FSI::StructureALE> fsi = Teuchos::rcp(new FSI::StructureALE(comm));

    if (genprob.restart)
    {
      // read the restart information, set vectors and variables
      fsi->ReadRestart(genprob.restart);
    }

    fsi->Timeloop();

    DRT::ResultTestManager testmanager(comm);
    testmanager.AddFieldTest(fsi->StructureField().CreateFieldTest());
    testmanager.TestAll();    
  }
  Teuchos::TimeMonitor::summarize();
}

/*----------------------------------------------------------------------*/
// entry point for FSI using XFEM in DRT
/*----------------------------------------------------------------------*/
void xfsi_drt()
{
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  
  Teuchos::RefCountPtr<FSI::DirichletNeumannCoupling> fsi = rcp(new FSI::DirichletNeumannCoupling(comm));

  if (genprob.restart)
  {
      // read the restart information, set vectors and variables
      fsi->ReadRestart(genprob.restart);
  }

  fsi->Timeloop(fsi);

  DRT::ResultTestManager testmanager(comm);
  testmanager.AddFieldTest(fsi->MBFluidField().CreateFieldTest());
  testmanager.TestAll();

  Teuchos::TimeMonitor::summarize();
}

#endif
