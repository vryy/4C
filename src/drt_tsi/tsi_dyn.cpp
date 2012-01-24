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
#include "../drt_comm/comm_utils.H"
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
  const Epetra_Comm& comm = DRT::Problem::Instance()->Dis(genprob.numsf,0)->Comm();
#else
  Epetra_SerialComm comm;
#endif

  // print TSI-Logo to screen
  if (comm.MyPID()==0) TSI::printlogo();

  // setup of the discretizations, including clone strategy
  TSI::UTILS::SetupTSI(comm);

  // access the problem-specific parameter list
  const Teuchos::ParameterList& tsidyn
    = DRT::Problem::Instance()->TSIDynamicParams();
  // access the problem-specific parameter list
  const Teuchos::ParameterList& sdynparams
    = DRT::Problem::Instance()->StructuralDynamicParams();
  const INPAR::TSI::SolutionSchemeOverFields coupling
    = DRT::INPUT::IntegralValue<INPAR::TSI::SolutionSchemeOverFields>(tsidyn,"COUPALGO");

  // choose algorithm depending on solution type
  switch (coupling)
  {
  case INPAR::TSI::Monolithic:
  {
    // create an TSI::Monolithic instance
    Teuchos::RCP<TSI::Monolithic> tsi = Teuchos::rcp(new TSI::Monolithic(comm,sdynparams));

    if (genprob.restart)
    {
      // read the restart information, set vectors and variables
      tsi->ReadRestart(genprob.restart);
    }

    // now do the coupling setup and create the combined dofmap
    tsi->SetupSystem();

    // solve the whole tsi problem
    tsi->TimeLoop(sdynparams);

    // summarize the performance measurements
    Teuchos::RCP<const Teuchos::Comm<int> > TeuchosComm = COMM_UTILS::toTeuchosComm(comm);
    Teuchos::TimeMonitor::summarize(TeuchosComm.ptr());

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
    Teuchos::RCP<const Teuchos::Comm<int> > TeuchosComm = COMM_UTILS::toTeuchosComm(comm);
    Teuchos::TimeMonitor::summarize(TeuchosComm.ptr());

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
