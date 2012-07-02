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

//for cloning
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_scatra/scatra_utils_clonestrategy.H"

#include "../linalg/linalg_utils.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_Time.h>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSI_Base::SSI_Base(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams):
    AlgorithmBase(comm, timeparams)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  const Teuchos::ParameterList& scatradyn  = problem->ScalarTransportDynamicParams();
  // get the solver number used for ScalarTransport solver
  const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");

  //2.- Setup discretizations.
  SetupDiscretizations(comm);

  //3.- Create the two uncoupled subproblems.
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure =
      Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(timeparams));
  structure_ = rcp_dynamic_cast<ADAPTER::Structure>(structure->StructureFieldrcp());
  scatra_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(timeparams,
      true, 0, problem->SolverParams(linsolvernumber)));

  zeros_ = LINALG::CreateVector(*structure_->DofRowMap(), true);
}


/*----------------------------------------------------------------------*
 | read restart information for given time step (public)   vuong 01/12  |
 *----------------------------------------------------------------------*/
void SSI::SSI_Base::ReadRestart( int restart)
{
  if (restart)
  {
    scatra_->ScaTraField().ReadRestart(restart);
    structure_->ReadRestart(restart);

    SetTimeStep(structure_->GetTime(), restart);
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
  Teuchos::RCP<DRT::Discretization> structdis = problem->Dis(genprob.numsf, 0); // Dis(0,0)
  Teuchos::RCP<DRT::Discretization> scatradis = problem->Dis(genprob.numscatra,0); // Dis(1,0)

  if(!scatradis->Filled())
    scatradis->FillComplete();

  if (scatradis->NumGlobalNodes()==0)
  {
    Epetra_Time time(comm);
    // create the fluid scatra discretization
    {
      Teuchos::RCP<DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy> > clonewizard =
      Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy>() );
      std::pair<string,string> key("structure","scatra");
      std::map<int,int> structmatmap = problem->ClonedMaterialMap()[key];
      clonewizard->CreateMatchingDiscretization(structdis,scatradis,structmatmap);
    }
    if (comm.MyPID()==0)
    cout <<"Created scalar transport discretization from structure discretization in... "<< time.ElapsedTime() << " secs\n\n";
  }
  else
  dserror("Structure AND ScaTra discretization present. This is not supported.");
}
