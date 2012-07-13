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


/*----------------------------------------------------------------------*
 |  headers                                                             |
 *----------------------------------------------------------------------*/
#include "poro_dyn.H"
#include "poro_scatra.H"
#include "poroelast_monolithic.H"
#include "poro_monolithicstructuresplit.H"
#include "poro_monolithicfluidsplit.H"
#include "poroelast_utils.H"
#include "../drt_inpar/inpar_poroelast.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include <Teuchos_TimeMonitor.hpp>

/*------------------------------------------------------------------------------------------------*
 | main control routine for poroelasticity problems                                   vuong 01/12 |
 *------------------------------------------------------------------------------------------------*/
void poroelast_drt()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // create a communicator
#ifdef PARALLEL
  const Epetra_Comm& comm =
      problem->GetDis("structure")->Comm();
#else
  Epetra_SerialComm comm;
#endif

  // setup of the discretizations, including clone strategy
  POROELAST::UTILS::SetupPoro(comm);

  // access the problem-specific parameter list
  const Teuchos::ParameterList& poroelastdyn =
      problem->PoroelastDynamicParams();
  // access the problem-specific parameter list
  const Teuchos::ParameterList& sdynparams =
      problem->StructuralDynamicParams();
  const INPAR::POROELAST::SolutionSchemeOverFields coupling =
      DRT::INPUT::IntegralValue<INPAR::POROELAST::SolutionSchemeOverFields>(
          poroelastdyn, "COUPALGO");

  std::string damping = sdynparams.get<std::string>("DAMPING");
  if(damping != "Material")
    dserror("Material damping has to be used for poroelasticity! Set DAMPING to Material in the STRUCTURAL DYNAMIC section.");

  // create an empty Poroelast::Algorithm instance
  Teuchos::RCP<POROELAST::PoroBase> poroelast = Teuchos::null;

  // choose algorithm depending on solution type (only monolithic type implemented)
  switch (coupling)
  {
    case INPAR::POROELAST::Monolithic:
    {
      // create an POROELAST::Monolithic instance
      poroelast = Teuchos::rcp(new POROELAST::Monolithic(comm, poroelastdyn));
      break;
    } // monolithic case

    case INPAR::POROELAST::Monolithic_structuresplit:
    {
      // create an POROELAST::MonolithicStructureSplit instance
      poroelast = Teuchos::rcp(new POROELAST::MonolithicStructureSplit(comm, poroelastdyn));
      break;
    }
    case INPAR::POROELAST::Monolithic_fluidsplit:
    {
      // create an POROELAST::MonolithicFluidSplit instance
      poroelast = Teuchos::rcp(new POROELAST::MonolithicFluidSplit(comm, poroelastdyn));
      break;
    }
    default:
      dserror("Unknown solutiontype for poroelasticity: %d",coupling);
    } // end switch

  // read the restart information, set vectors and variables
  const int restart = problem->Restart();
  poroelast->ReadRestart(restart);

  // now do the coupling setup and create the combined dofmap
  poroelast->SetupSystem();

  // solve the whole problem
  poroelast->TimeLoop();

  // summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // perform the result test
  poroelast->TestResults(comm);

  return;
}//poroelast_drt()

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void poro_scatra_drt()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  //1.- Initialization
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  //2.- Parameter reading
  const Teuchos::ParameterList& poroscatradynparams = problem->PoroScatraControlParams();

  //3.- Creation of Poroelastic + Scalar_Transport problem. (Discretization called inside)
  Teuchos::RCP<PORO_SCATRA::PartPORO_SCATRA> poro_scatra = Teuchos::rcp(
      new PORO_SCATRA::PartPORO_SCATRA(comm, poroscatradynparams));

  //3.1- Read restart if needed. (Discretization called inside)
  poro_scatra->ReadRestart();

  //4.- Run of the actual problem.

  // 4.1.- Some setup needed for the poroelastic subproblem.
  poro_scatra->SetupSystem();

  // 4.2.- Solve the whole problem
  poro_scatra->Timeloop();

  // 4.3.- Summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // 5. - perform the result test
  poro_scatra->TestResults(comm);

}
