/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_base.cpp

 \brief base algorithm for scalar transport within multiphase porous medium

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#include "poromultiphase_scatra_base.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_poromultiphase/poromultiphase_base.H"
#include "../drt_poromultiphase/poromultiphase_utils.H"

#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_poromulti.H"

/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16  |
 *----------------------------------------------------------------------*/
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::PoroMultiPhaseScaTraBase(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams):
    AlgorithmBase(comm, globaltimeparams),
    poromulti_(Teuchos::null),
    scatra_(Teuchos::null),
    ndsporofluid_scatra_(-1)
{

}


/*----------------------------------------------------------------------*
 | initialize algorithm                                    vuong 08/16  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::Init(
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& algoparams,
    const Teuchos::ParameterList& poroparams,
    const Teuchos::ParameterList& structparams,
    const Teuchos::ParameterList& fluidparams,
    const Teuchos::ParameterList& scatraparams,
    const std::string& struct_disname,
    const std::string& fluid_disname,
    const std::string& scatra_disname,
    bool isale,
    int nds_disp,
    int nds_vel,
    int nds_solidpressure,
    int ndsporofluid_scatra)
{
  //save the dofset number of the scatra on the fluid dis
  ndsporofluid_scatra_ =ndsporofluid_scatra;

  // access the global problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // Create the two uncoupled subproblems.

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // coupling scheme
  // -------------------------------------------------------------------
  INPAR::POROMULTIPHASE::SolutionSchemeOverFields solscheme =
    DRT::INPUT::IntegralValue<INPAR::POROMULTIPHASE::SolutionSchemeOverFields>(poroparams,"COUPALGO");

  poromulti_ = POROMULTIPHASE::UTILS::CreatePoroMultiPhaseAlgorithm(solscheme,globaltimeparams,Comm());

  // initialize
  poromulti_->Init(
                globaltimeparams,
                poroparams,
                structparams,
                fluidparams,
                struct_disname,
                fluid_disname,
                isale,
                nds_disp,
                nds_vel,
                nds_solidpressure,
                ndsporofluid_scatra);

  // get the solver number used for ScalarTransport solver
  const int linsolvernumber = scatraparams.get<int>("LINEAR_SOLVER");

  // scatra problem
  scatra_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());

  // initialize the base algo.
  // scatra time integrator is constructed and initialized inside.
  scatra_->Init(
      globaltimeparams,
      scatraparams,
      problem->SolverParams(linsolvernumber),
      scatra_disname,
      true);

  // only now we must call Setup() on the scatra time integrator.
  // all objects relying on the parallel distribution are
  // created and pointers are set.
  // calls Setup() on the scatra time integrator inside.
  scatra_->ScaTraField()->Setup();

  //done.
  return;
}

/*----------------------------------------------------------------------*
 | read restart information for given time step (public)   vuong 08/16  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::ReadRestart( int restart )
{
  if (restart)
  {
    // read restart data for structure field (will set time and step internally)
    poromulti_->ReadRestart(restart);

    // read restart data for fluid field (will set time and step internally)
    scatra_->ScaTraField()->ReadRestart(restart);

    // reset time and step for the global algorithm
    SetTimeStep(poromulti_->Time(), restart);
  }

  return;
}

/*----------------------------------------------------------------------*
 | Test the results of all subproblems                       vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::CreateFieldTest()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  poromulti_->CreateFieldTest();
  problem->AddFieldTest(scatra_->CreateScaTraFieldTest());
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::SetPoroSolution()
{
  SetMeshDisp();
  SetSolutionFields();
}

/*---------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::SetSolutionFields()
{
  // cast
  Teuchos::RCP<SCATRA::ScaTraTimIntPoroMulti> poroscatra =
      Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntPoroMulti>(scatra_->ScaTraField());

  if(poroscatra==Teuchos::null)
    dserror("cast to ScaTraTimIntPoroMulti failed!");

  // set the solution
  poroscatra->SetSolutionFields(
      poromulti_->FluidFlux(),
      1,
      poromulti_->FluidPressure(),
      2,
      poromulti_->FluidSaturation(),
      2,
      poromulti_->SolidPressure(),
      3
      );
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::SetMeshDisp()
{
  scatra_->ScaTraField()->ApplyMeshMovement(
      poromulti_->StructDispnp(),
      1
      );
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::SetScatraSolution()
{
  poromulti_->SetScatraSolution(ndsporofluid_scatra_,scatra_->ScaTraField()->Phinp());
  return;
}
