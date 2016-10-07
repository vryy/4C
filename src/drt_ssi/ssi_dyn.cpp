/*!------------------------------------------------------------------------------------------------*
 \file ssi_dyn.cpp

 \brief control routine for scalar structure interaction

 \level 1

 \maintainer Anh-Tu Vuong
             vuong@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15264

 *------------------------------------------------------------------------------------------------*/

#include "ssi_dyn.H"
#include "ssi_partitioned_1wc.H"
#include "ssi_partitioned_2wc.H"
#include "ssi_utils.H"
#include "../drt_inpar/inpar_ssi.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_io/io_control.H"

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_utils_parallel.H"


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


  //3.- Creation of Poroelastic + Scalar_Transport problem. (Discretization called inside)
  Teuchos::RCP<SSI::SSI_Base> ssi = Teuchos::null;

  //3.1 choose algorithm depending on solution type
  switch(coupling)
  {
  case INPAR::SSI::ssi_OneWay_ScatraToSolid:
    ssi = Teuchos::rcp(new SSI::SSI_Part1WC_ScatraToSolid(comm, ssiparams));
    break;
  case INPAR::SSI::ssi_OneWay_SolidToScatra:
    ssi = Teuchos::rcp(new SSI::SSI_Part1WC_SolidToScatra(comm, ssiparams));
    break;
  case INPAR::SSI::ssi_IterStagg:
    ssi = Teuchos::rcp(new SSI::SSI_Part2WC(comm, ssiparams));
    break;
  case INPAR::SSI::ssi_IterStaggFixedRel_ScatraToSolid:
    ssi = Teuchos::rcp(new SSI::SSI_Part2WC_ScatraToSolid_Relax(comm, ssiparams));
    break;
  case INPAR::SSI::ssi_IterStaggFixedRel_SolidToScatra:
    ssi = Teuchos::rcp(new SSI::SSI_Part2WC_SolidToScatra_Relax(comm, ssiparams));
    break;
  case INPAR::SSI::ssi_IterStaggAitken_ScatraToSolid:
    ssi = Teuchos::rcp(new SSI::SSI_Part2WC_ScatraToSolid_Relax_Aitken(comm, ssiparams));
        break;
  case INPAR::SSI::ssi_IterStaggAitken_SolidToScatra:
    ssi = Teuchos::rcp(new SSI::SSI_Part2WC_SolidToScatra_Relax_Aitken(comm, ssiparams));
    break;
  default:
    dserror("unknown coupling algorithm for SSI!");
    break;
  }

  problem->GetDis("structure")->FillComplete(true,true,true);
  problem->GetDis("scatra")->FillComplete(true,true,true);

  //3.1.1 init the chosen ssi algorithm
  int redistribute = (int)SSI::none;
  redistribute = ssi -> Init(comm, ssiparams, scatradyn, sdyn, "structure", "scatra");

  if (redistribute == (int)SSI::match)
  {
    // first we bin the scatra discretization
    std::vector<Teuchos::RCP<DRT::Discretization> > dis;
    dis.push_back(problem->GetDis("scatra"));
    DRT::UTILS::RedistributeDiscretizationsByBinning(dis,false);

    // now we redistribute the structure dis to match the scatra dis
    DRT::UTILS::MatchElementDistributionOfMatchingDiscretizations(
        *problem->GetDis("scatra"),
        *problem->GetDis("structure")
    );
  }
  else if(redistribute == (int)SSI::binning)
  {
    // create vector of discr.
    std::vector<Teuchos::RCP<DRT::Discretization> > dis;
    dis.push_back(problem->GetDis("structure"));
    dis.push_back(problem->GetDis("scatra"));

    DRT::UTILS::RedistributeDiscretizationsByBinning(dis,false);
  }

  // now we can finally fill our discretizations
  problem->GetDis("structure")->FillComplete(true,false,true);
  problem->GetDis("scatra")->FillComplete(true,false,true);

  // now as we redistributed our discretizations we can construct all
  // objects relying on the parallel distribution
  ssi -> Setup();

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
  if(coupling == INPAR::SSI::ssi_OneWay_ScatraToSolid and flag_readscatra ){
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
