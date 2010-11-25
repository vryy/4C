/*!----------------------------------------------------------------------*/
/*!
\file tsi_dyn.cpp
\brief Control routine for thermo-structure-interaction problems.


<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 12/09 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifdef PARALLEL
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

/*----------------------------------------------------------------------*
 |  headers                                                  dano 12/09 |
 *----------------------------------------------------------------------*/
#include "tsi_dyn.H"
#include "tsi_algorithm.H"
#include "tsi_monolithic.H"
#include "tsi_utils.H"
#include "../drt_inpar/inpar_tsi.H"

#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"

#include <Teuchos_TimeMonitor.hpp>

//#if 0
//#include "../drt_io/io_gmsh.H"
//#endif

/*----------------------------------------------------------------------*
 | general problem data                                     m.gee 06/01 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 | entry point for TSI in DRT                                dano 12/09 |
 *----------------------------------------------------------------------*/
//void tsi_dyn(int disnumsf,int disnumtf,int restart)
void tsi_dyn_drt()
{
  // create a communicator
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // print TSI-Logo to screen
  if (comm.MyPID()==0) TSI::printlogo();

  // access the structure discretization, make sure it is filled
  Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
  structdis = DRT::Problem::Instance()->Dis(genprob.numsf,0);
  // set degrees of freedom in the discretization
  if (!structdis->Filled() or !structdis->HaveDofs()) structdis->FillComplete();

  // access the thermo discretization
  Teuchos::RCP<DRT::Discretization> thermdis = Teuchos::null;
  thermdis = DRT::Problem::Instance()->Dis(genprob.numtf,0);
  if (!thermdis->Filled()) thermdis->FillComplete();

   // we use the structure discretization as layout for the temperature discretization
  if (structdis->NumGlobalNodes()==0) dserror("Structure discretization is empty!");

  // create thermo elements if the temperature discretization is empty
  if (thermdis->NumGlobalNodes()==0)
  {
    Epetra_Time time(comm);

    // fetch the desired material id for the thermo elements
    const int matid = -1;

    // create the thermo discretization
    {
      Teuchos::RCP<DRT::UTILS::DiscretizationCreator<TSI::UTILS::ThermoStructureCloneStrategy> > clonewizard
        = Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<TSI::UTILS::ThermoStructureCloneStrategy>() );

      clonewizard->CreateMatchingDiscretization(structdis,thermdis,matid);
    }

    if (comm.MyPID()==0)
      cout<<"Created thermo discretization from structure field in...."
      <<time.ElapsedTime() << " secs\n\n";
  }
  else
      dserror("Structure AND Thermo discretization present. This is not supported.");

  // access the problem-specific parameter list
  const Teuchos::ParameterList& tsidyn = DRT::Problem::Instance()->TSIDynamicParams();
  const INPAR::TSI::SolutionSchemeOverFields coupling  = Teuchos::getIntegralValue<INPAR::TSI::SolutionSchemeOverFields>(tsidyn,"COUPALGO");

  // choose algorithm depending on solution type
  switch (coupling)
  {
  case INPAR::TSI::Monolithic:
  {
    // create an TSI::Monolithic instance
    Teuchos::RCP<TSI::Monolithic> tsi = Teuchos::rcp(new TSI::Monolithic(comm));

    if (genprob.restart)
    {
      // read the restart information, set vectors and variables
      tsi->ReadRestart(genprob.restart);
    }

    // now do the coupling setup an create the combined dofmap
    tsi->SetupSystem();

    // solve the whole tsi problem
    tsi->TimeLoop();

    // summarize the performance measurements
    Teuchos::TimeMonitor::summarize();

    // perform the result test
    DRT::Problem::Instance()->AddFieldTest(tsi->StructureField().CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(tsi->ThermoField().CreateFieldTest());
    DRT::Problem::Instance()->TestAll(comm);
    break;

  }  // monolithic case
  case INPAR::TSI::OneWay:
  case INPAR::TSI::SequStagg:
  case INPAR::TSI::IterStagg:
  {
    // Any partitioned algorithm. Stable of working horses.
    // create an TSI::Algorithm instance
    Teuchos::RCP<TSI::Algorithm> tsi = Teuchos::rcp(new TSI::Algorithm(comm));

    if (genprob.restart)
    {
      // read the restart information, set vectors and variables
      tsi->ReadRestart(genprob.restart);
    }

    // solve the whole tsi problem
    tsi->TimeLoop();

    // summarize the performance measurements
    Teuchos::TimeMonitor::summarize();

    // perform the result test
    DRT::Problem::Instance()->AddFieldTest(tsi->StructureField().CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(tsi->ThermoField().CreateFieldTest());
    DRT::Problem::Instance()->TestAll(comm);
    break;
  }  // partitioned case
  default:
     dserror("Unknown solutiontype for thermo-structure interaction: %d",coupling);
  }  // end switch

  return;
} // tsi_dyn_drt()

/*----------------------------------------------------------------------*/
#endif  // CCADISCRET
