/*!------------------------------------------------------------------------------------------------*
 \file ssi_base.cpp

 \brief base class for all scalar structure algorithms

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *------------------------------------------------------------------------------------------------*/

#include "../drt_lib/drt_globalproblem.H"
#include "ssi_base.H"
#include "ssi_partitioned.H"
#include "ssi_utils.H"
#include "../drt_inpar/inpar_ssi.H"

#include "../drt_adapter/ad_str_wrapper.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

//for cloning
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_scatra/scatra_utils_clonestrategy.H"

#include "../linalg/linalg_utils.H"

#include "../drt_particle/binning_strategy.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_Time.h>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSI_Base::SSI_Base(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams):
    AlgorithmBase(comm, globaltimeparams)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // get the solver number used for ScalarTransport solver
  const int linsolvernumber = scatraparams.get<int>("LINEAR_SOLVER");

  //2.- Setup discretizations.
  SetupDiscretizations(comm);

  //3.- Create the two uncoupled subproblems.
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");

  // Set isale to false what should be the case in scatratosolid algorithm
  const INPAR::SSI::SolutionSchemeOverFields coupling
      = DRT::INPUT::IntegralValue<INPAR::SSI::SolutionSchemeOverFields>(problem->SSIControlParams(),"COUPALGO");

  bool isale = true;
  if(coupling == INPAR::SSI::ssi_OneWay_ScatraToSolid) isale = false;

  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure =
      Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(structparams, const_cast<Teuchos::ParameterList&>(structparams), structdis));
  structure_ = Teuchos::rcp_dynamic_cast<ADAPTER::Structure>(structure->StructureField());
  scatra_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(scatraparams,isale,"scatra", problem->SolverParams(linsolvernumber)));
  zeros_ = LINALG::CreateVector(*structure_->DofRowMap(), true);

}


/*----------------------------------------------------------------------*
 | read restart information for given time step (public)   vuong 01/12  |
 *----------------------------------------------------------------------*/
void SSI::SSI_Base::ReadRestart( int restart )
{

  if (restart)
  {
    scatra_->ScaTraField()->ReadRestart(restart);
    structure_->ReadRestart(restart);
    SetTimeStep(structure_->TimeOld(), restart);
  }

  return;
}

/*----------------------------------------------------------------------*
 | read restart information for given time (public)        AN, JH 10/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Base::ReadRestartfromTime( double restarttime )
{
  if ( restarttime > 0.0 )
  {

    int restartstructure = SSI::Utils::CheckTimeStepping(structure_->Dt(), restarttime);
    int restartscatra    = SSI::Utils::CheckTimeStepping(scatra_->ScaTraField()->Dt(), restarttime);

    scatra_->ScaTraField()->ReadRestart(restartscatra);
    structure_->ReadRestart(restartstructure);
    SetTimeStep(structure_->TimeOld(), restartstructure);

  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Base::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  problem->AddFieldTest(structure_->CreateFieldTest());
  problem->AddFieldTest(scatra_->CreateScaTraFieldTest());
  problem->TestAll(comm);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Base::SetupDiscretizations(const Epetra_Comm& comm)
{
  // Scheme   : the structure discretization is received from the input. Then, an ale-scatra disc. is cloned.

  DRT::Problem* problem = DRT::Problem::Instance();

  //1.-Initialization.
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");
  if(!scatradis->Filled())
    scatradis->FillComplete();

  if (scatradis->NumGlobalNodes()==0)
  {
    DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(structdis,scatradis);
  }
  else
  {

   std::map<std::string,std::string> conditions_to_copy;
   SCATRA::ScatraFluidCloneStrategy clonestrategy;
   conditions_to_copy = clonestrategy.ConditionsToCopy();
   DRT::UTILS::DiscretizationCreatorBase creator;
   creator.CopyConditions(scatradis,scatradis,conditions_to_copy);

    // redistribute discr. with help of binning strategy
    if(scatradis->Comm().NumProc()>1)
    {
      scatradis->FillComplete();
      structdis->FillComplete();
      // create vector of discr.
      std::vector<Teuchos::RCP<DRT::Discretization> > dis;
      dis.push_back(structdis);
      dis.push_back(scatradis);

      std::vector<Teuchos::RCP<Epetra_Map> > stdelecolmap;
      std::vector<Teuchos::RCP<Epetra_Map> > stdnodecolmap;

      /// binning strategy is created and parallel redistribution is performed
      Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy =
        Teuchos::rcp(new BINSTRATEGY::BinningStrategy(dis,stdelecolmap,stdnodecolmap));
    }
  }

}
