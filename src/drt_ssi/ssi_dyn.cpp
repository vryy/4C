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

#include "../drt_io/io_control.H"

#include <Teuchos_TimeMonitor.hpp>

// forward declaration
namespace SSI{
    namespace Utils{
      void ChangeTimeParameter(const Epetra_Comm&,
          Teuchos::ParameterList&,
          Teuchos::ParameterList&,
          Teuchos::ParameterList&);
    };
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ssi_drt()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  //1.- Initialization
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  //2.- Parameter reading
  Teuchos::ParameterList& ssiparams = const_cast<Teuchos::ParameterList&>( problem->SSIControlParams() );
  // access scatra params list
  Teuchos::ParameterList& scatradyn = const_cast<Teuchos::ParameterList&>( problem->ScalarTransportDynamicParams() );
  // access structural dynamic params list which will be possibly modified while creating the time integrator
  Teuchos::ParameterList& sdyn      = const_cast<Teuchos::ParameterList&>( DRT::Problem::Instance()->StructuralDynamicParams() );

//  //Modification of time parameter list
  SSI::Utils::ChangeTimeParameter(comm, ssiparams, scatradyn, sdyn);

  const INPAR::SSI::SolutionSchemeOverFields coupling
    = DRT::INPUT::IntegralValue<INPAR::SSI::SolutionSchemeOverFields>(ssiparams,"COUPALGO");

  //

  //3.- Creation of Poroelastic + Scalar_Transport problem. (Discretization called inside)
  Teuchos::RCP<SSI::SSI_Base> ssi = Teuchos::null;

  //3.1 choose algorithm depending on solution type
  switch(coupling)
  {
  case INPAR::SSI::Part_SolidToScatra:
    ssi = Teuchos::rcp(new SSI::SSI_Part1WC_SolidToScatra(comm, ssiparams, scatradyn, sdyn));
    break;
  case INPAR::SSI::Part_ScatraToSolid:
    ssi = Teuchos::rcp(new SSI::SSI_Part1WC_ScatraToSolid(comm, ssiparams, scatradyn, sdyn));
    break;
  case INPAR::SSI::Part_TwoWay:
    ssi = Teuchos::rcp(new SSI::SSI_Part2WC(comm, ssiparams, scatradyn, sdyn));
    break;
  default:
    dserror("unknown coupling algorithm for SSI!");
    break;
  }

  //3.2- Read restart if needed. (Discretization called inside)
  const int restart = problem->Restart();

  const double restarttime = problem->RestartTime();
  if (restarttime > 0.0)
    ssi->ReadRestartfromTime(restarttime);

  else
    if (restart)
      ssi->ReadRestart(restart);

  // 3.3 AFTER restart: reset inputfilename of the problem so that results from other runs can be read
  bool flag_readscatra = DRT::INPUT::IntegralValue<bool>(ssiparams,"SCATRA_FROM_RESTART_FILE");
  if(coupling == INPAR::SSI::Part_ScatraToSolid and flag_readscatra ){
       std::string filename = Teuchos::getNumericStringParameter(ssiparams,"SCATRA_FILENAME");
       Teuchos::RCP<IO::InputControl> inputscatra = Teuchos::rcp(new IO::InputControl(filename, comm));
       problem->SetInputControlFile(inputscatra);
  }

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
