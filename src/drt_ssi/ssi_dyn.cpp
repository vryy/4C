/*!------------------------------------------------------------------------------------------------*
 \file ssi_dyn.cpp

 \brief control routine for scalar structure interaction

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *------------------------------------------------------------------------------------------------*/

#include "ssi_dyn.H"
#include "ssi_partitioned_1wc.H"
#include "ssi_partitioned_2wc.H"
#include "../drt_inpar/inpar_ssi.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ssi_drt()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  //1.- Initialization
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  //2.- Parameter reading
  const Teuchos::ParameterList& ssiparams = problem->SSIControlParams();
  const INPAR::SSI::SolutionSchemeOverFields coupling
    = DRT::INPUT::IntegralValue<INPAR::SSI::SolutionSchemeOverFields>(ssiparams,"COUPALGO");

  //3.- Creation of Poroelastic + Scalar_Transport problem. (Discretization called inside)
  Teuchos::RCP<SSI::SSI_Base> ssi = Teuchos::null;

  // choose algorithm depending on solution type
  switch(coupling)
  {
  case INPAR::SSI::Part_SolidToScatra:
    ssi = Teuchos::rcp(new SSI::SSI_Part1WC_SolidToScatra(comm, ssiparams));
    break;
  case INPAR::SSI::Part_ScatraToSolid:
    ssi = Teuchos::rcp(new SSI::SSI_Part1WC_ScatraToSolid(comm, ssiparams));
    break;
  case INPAR::SSI::Part_TwoWay:
    ssi = Teuchos::rcp(new SSI::SSI_Part2WC(comm, ssiparams));
    break;
  default:
    dserror("unknown coupling algorithm for SSI!");
  }

  //3.1- Read restart if needed. (Discretization called inside)
  const int restart = problem->Restart();
  ssi->ReadRestart(restart);

  //4.- Run of the actual problem.

  // 4.1.- Some setup needed for the poroelastic subproblem.
  ssi->SetupSystem();

  // 4.2.- Solve the whole problem
  ssi->Timeloop();

  // 4.3.- Summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // 5. - perform the result test
  ssi->TestResults(comm);

}
