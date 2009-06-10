
#ifdef CCADISCRET

#include <string>
#include <vector>
#include <set>
#include <functional>

#include <Teuchos_TimeMonitor.hpp>

#include "fsi_dyn.H"
#include "fsi_dirichletneumann.H"
#include "fsi_robinneumann.H"
#include "fsi_robin.H"
#include "fsi_monolithicoverlap.H"
#include "fsi_monolithiclagrange.H"
#include "fsi_monolithicstructuresplit.H"
#include "fsi_partitionedmonolithic.H"
#include "fsi_structureale.H"
#include "fsi_fluid_ale.H"
#include "fsi_fluid_xfem.H"
#include "fsi_utils.H"

#include "../drt_inpar/inpar_fsi.H"
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

#include "../drt_lib/drt_condition_utils.H"

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
  {
    FSI::UTILS::CreateAleDiscretization();

    if(comm.MyPID()==0)
    {
      cout << "\n\nCreating ALE discretisation ....\n\n";
    }
  }

  Teuchos::RCP<FSI::FluidAleAlgorithm> fluid = Teuchos::rcp(new FSI::FluidAleAlgorithm(comm));
  if (genprob.restart)
  {
    // read the restart information, set vectors and variables
    fluid->ReadRestart(genprob.restart);
  }
  fluid->Timeloop();

  DRT::ResultTestManager testmanager(comm);
  testmanager.AddFieldTest(fluid->MBFluidField().CreateFieldTest());
  testmanager.TestAll();
}


/*----------------------------------------------------------------------*/
// entry point for Fluid on Ale in DRT
/*----------------------------------------------------------------------*/
void fluid_xfem_drt()
{
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif


  // make sure the solid dis is filled
  RCP<DRT::Problem> problem = DRT::Problem::Instance();
  problem->Dis(genprob.numsf,0)->FillComplete();
  
  Teuchos::RCP<FSI::FluidXFEMAlgorithm> xfluid = Teuchos::rcp(new FSI::FluidXFEMAlgorithm(comm));
  if (genprob.restart)
  {
    // read the restart information, set vectors and variables
    xfluid->ReadRestart(genprob.restart);
  }
  xfluid->Timeloop();

  DRT::ResultTestManager testmanager(comm);
  testmanager.AddFieldTest(xfluid->MBFluidField().CreateFieldTest());
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
    FSI::UTILS::CreateAleDiscretization();

  Teuchos::RCP<FSI::FluidAleAlgorithm> fluid = Teuchos::rcp(new FSI::FluidAleAlgorithm(comm));

  fluid->Timeloop();

  DRT::ResultTestManager testmanager(comm);
  testmanager.AddFieldTest(fluid->MBFluidField().CreateFieldTest());
  testmanager.TestAll();
}


/*----------------------------------------------------------------------*/
// entry point for FSI using ALE in DRT
/*----------------------------------------------------------------------*/
void fsi_ale_drt()
{
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  RCP<DRT::Problem> problem = DRT::Problem::Instance();

  // make sure the three discretizations are filled in the right order
  // this creates dof numbers with
  //
  //       structure dof < fluid dof < ale dof
  //
  // We rely on this ordering in certain non-intuitive places!

  problem->Dis(genprob.numsf,0)->FillComplete();
  problem->Dis(genprob.numff,0)->FillComplete();
  problem->Dis(genprob.numaf,0)->FillComplete();

  // create ale elements if the ale discretization is empty
  RCP<DRT::Discretization> aledis = problem->Dis(genprob.numaf,0);
  if (aledis->NumGlobalNodes()==0)
    FSI::UTILS::CreateAleDiscretization();

  const Teuchos::ParameterList& fsidyn   = problem->FSIDynamicParams();

  int coupling = Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO");
  switch (coupling)
  {
  case fsi_pseudo_structureale:
  {
    // pseudo FSI problem used to find starting configuration

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
    break;
  }
  case fsi_iter_monolithic:
  case fsi_iter_monolithicstructuresplit:
  case fsi_iter_monolithiclagrange:
  {
    Teuchos::RCP<FSI::Monolithic> fsi;

    INPAR::FSI::LinearBlockSolver linearsolverstrategy = Teuchos::getIntegralValue<INPAR::FSI::LinearBlockSolver>(fsidyn,"LINEARBLOCKSOLVER");

    // call constructor to initialise the base class
    if (linearsolverstrategy==INPAR::FSI::PartitionedAitken or
        linearsolverstrategy==INPAR::FSI::PartitionedVectorExtrapolation or
        linearsolverstrategy==INPAR::FSI::PartitionedJacobianFreeNewtonKrylov)
    {
      fsi = Teuchos::rcp(new FSI::PartitionedMonolithic(comm));
    }
    else if (coupling==fsi_iter_monolithic)
    {
      fsi = Teuchos::rcp(new FSI::MonolithicOverlap(comm));
    }
    else if (coupling==fsi_iter_monolithicstructuresplit)
    {
      fsi = Teuchos::rcp(new FSI::MonolithicStructureSplit(comm));
    }
    else if (coupling==fsi_iter_monolithiclagrange)
    {
      fsi = Teuchos::rcp(new FSI::MonolithicLagrange(comm));
    }
    else
    {
      dserror("Cannot find appropriate monolithic solver for coupling %d and linear strategy %d",coupling,linearsolverstrategy);
    }

    // read the restart information, set vectors and variables --- 
    // be careful, dofmaps might be changed here in a Redistribute call
    if (genprob.restart)
    {
      fsi->ReadRestart(genprob.restart);
    }

    // now do the coupling setup an create the combined dofmap
    fsi->SetupSystem();

    // here we go...
    fsi->Timeloop(fsi);

    DRT::ResultTestManager testmanager(comm);
    testmanager.AddFieldTest(fsi->FluidField().CreateFieldTest());
    testmanager.AddFieldTest(fsi->StructureField().CreateFieldTest());
    testmanager.TestAll();
    break;
  }
  default:
  {
    // Any partitioned algorithm. Stable of working horses.

    Teuchos::RCP<FSI::Partitioned> fsi;

    INPAR::FSI::PartitionedCouplingMethod method =
      Teuchos::getIntegralValue<INPAR::FSI::PartitionedCouplingMethod>(fsidyn,"PARTITIONED");

    if (method==INPAR::FSI::DirichletNeumann)
    {
      fsi = Teuchos::rcp(new FSI::DirichletNeumann(comm));
    }
    else if (method==INPAR::FSI::RobinNeumann)
    {
      fsi = Teuchos::rcp(new FSI::RobinNeumann(comm));
    }
    else
    {
      fsi = Teuchos::rcp(new FSI::Robin(comm));
    }

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
  }

  Teuchos::TimeMonitor::summarize(std::cout, false, true, false);
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

  if (comm.MyPID() == 0)
  {
    cout << endl;
    cout << YELLOW_LIGHT << "       @..@    " << END_COLOR << endl;
    cout << YELLOW_LIGHT << "      (----)      " << END_COLOR << endl;
    cout << YELLOW_LIGHT << "     ( >__< )   " << END_COLOR << endl;
    cout << YELLOW_LIGHT << "     ^^ ~~ ^^  " << END_COLOR << endl;
    cout << YELLOW_LIGHT << "     _     _ _______ _______ _____" << END_COLOR << endl;
    cout << YELLOW_LIGHT << "      \\\\__/  |______ |______   |  " << END_COLOR << endl;
    cout << YELLOW_LIGHT << "     _/  \\\\_ |       ______| __|__" << END_COLOR << endl;
    cout <<  endl << endl;
  }

  Teuchos::RCP<FSI::DirichletNeumann> fsi = rcp(new FSI::DirichletNeumann(comm));

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

  Teuchos::TimeMonitor::summarize();
}

#endif
