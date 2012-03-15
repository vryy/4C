/*!------------------------------------------------------------------------------------------------*
 \file poroelast.cpp

 \brief control routine of poroelasticity problems

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifdef PARALLEL
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

/*----------------------------------------------------------------------*
 |  headers                                                             |
 *----------------------------------------------------------------------*/
#include "poroelast.H"
#include "poroelast_monolithic.H"
#include "poroelast_utils.H"
#include "../drt_inpar/inpar_poroelast.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 | general problem data                                     m.gee 06/01 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB genprob;

/*------------------------------------------------------------------------------------------------*
 | main control routine for poroelasticity problems                                   vuong 01/12 |
 *------------------------------------------------------------------------------------------------*/
void poroelast_drt()
{
  // create a communicator
#ifdef PARALLEL
  const Epetra_Comm& comm =
      DRT::Problem::Instance()->Dis(genprob.numsf, 0)->Comm();
#else
  Epetra_SerialComm comm;
#endif

  // setup of the discretizations, including clone strategy
  POROELAST::UTILS::SetupPoro(comm);

  // access the problem-specific parameter list
  const Teuchos::ParameterList& poroelastdyn =
      DRT::Problem::Instance()->PoroelastDynamicParams();
  // access the problem-specific parameter list
  const Teuchos::ParameterList& sdynparams =
      DRT::Problem::Instance()->StructuralDynamicParams();
  const INPAR::POROELAST::SolutionSchemeOverFields coupling =
      DRT::INPUT::IntegralValue<INPAR::POROELAST::SolutionSchemeOverFields>(
          poroelastdyn, "COUPALGO");

  // choose algorithm depending on solution type (only monolithic type implemented)
  switch (coupling)
  {
    case INPAR::POROELAST::Monolithic:
    {
      // create an POROELAST::Monolithic instance
      Teuchos::RCP<POROELAST::Monolithic> poroelast = Teuchos::rcp(
          new POROELAST::Monolithic(comm, sdynparams));

      if (genprob.restart)
      {
        // read the restart information, set vectors and variables
        poroelast->ReadRestart(genprob.restart);
      }

      // now do the coupling setup and create the combined dofmap
      poroelast->SetupSystem();

      // solve the whole problem
      poroelast->TimeLoop(sdynparams);

      // summarize the performance measurements
      Teuchos::TimeMonitor::summarize();

      // perform the result test
      DRT::Problem::Instance()->AddFieldTest(
          poroelast->StructureField().CreateFieldTest());
      DRT::Problem::Instance()->AddFieldTest(
          poroelast->FluidField().CreateFieldTest());
      DRT::Problem::Instance()->TestAll(comm);
      break;
    } // monolithic case

    default:
      dserror("Unknown solutiontype for poroelasticity: %d",coupling);
    } // end switch

    return;
  }//poroelast_drt()

#endif  // CCADISCRET
