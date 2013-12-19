/*!----------------------------------------------------------------------
\file levelset_timint_ost.cpp

\brief one-step theta time integration scheme for level-set problems

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#include "levelset_timint_ost.H"

#include "../drt_scatra_ele/scatra_ele_action.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
//#include "../linalg/linalg_solver.H"

#include "../drt_particle/scatra_particle_coupling.H"


#include "../drt_inpar/drt_validparameters.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                rasthofer 09/13 |
 *----------------------------------------------------------------------*/
SCATRA::LevelSetTimIntOneStepTheta::LevelSetTimIntOneStepTheta(
  Teuchos::RCP<DRT::Discretization>      actdis,
  Teuchos::RCP<LINALG::Solver>           solver,
  Teuchos::RCP<Teuchos::ParameterList>   params,
  Teuchos::RCP<Teuchos::ParameterList>   sctratimintparams,
  Teuchos::RCP<Teuchos::ParameterList>   extraparams,
  Teuchos::RCP<IO::DiscretizationWriter> output)
: ScaTraTimIntImpl(actdis,solver,sctratimintparams,extraparams,output),
  LevelSetAlgorithm(actdis,solver,params,sctratimintparams,extraparams,output),
  TimIntOneStepTheta(actdis,solver,sctratimintparams,extraparams,output)
{
  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                              rasthofer 09/13 |
*-----------------------------------------------------------------------*/
SCATRA::LevelSetTimIntOneStepTheta::~LevelSetTimIntOneStepTheta()
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                         rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetTimIntOneStepTheta::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntOneStepTheta::Init();
  LevelSetAlgorithm::Init();

  return;
}


/*----------------------------------------------------------------------*
| Print information about current time step to screen   rasthofer 09/13 |
*-----------------------------------------------------------------------*/
void SCATRA::LevelSetTimIntOneStepTheta::PrintTimeStepInfo()
{
  if (myrank_==0)
  {
    if (not switchreinit_)
      TimIntOneStepTheta::PrintTimeStepInfo();
    else
      printf("PSEUDOTIMESTEP: %11.4E      %s          THETA = %11.4E   PSEUDOSTEP = %4d/%4d \n",
                  dtau_,MethodTitle().c_str(),thetareinit_,pseudostep_,pseudostepmax_);
  }
  return;
}


/*----------------------------------------------------------------------*
 | Initialization procedure before the first time step  rasthofer 12/13 |
 -----------------------------------------------------------------------*/
void SCATRA::LevelSetTimIntOneStepTheta::PrepareFirstTimeStep()
{
  if (not switchreinit_)
    TimIntOneStepTheta::PrepareFirstTimeStep();
  else
  {
    // set required general parameters
    Teuchos::ParameterList eleparams;

    eleparams.set<int>("action",SCATRA::set_lsreinit_scatra_parameter);

    // set type of scalar transport problem
    eleparams.set<int>("scatratype",scatratype_);

    // reinitialization equation id given in convective form
    // ale is not intended here
    eleparams.set<int>("form of convective term",INPAR::SCATRA::convform_convective);
    eleparams.set("isale",false);

    // parameters for stabilization
    eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

    // set level-set reitialization specific parameters
    eleparams.sublist("REINITIALIZATION") = levelsetparams_->sublist("REINITIALIZATION");
    // turn off stabilization
    Teuchos::setStringToIntegralParameter<int>("STABTYPEREINIT",
          "no_stabilization",
          "type of stabilization (if any)",
          Teuchos::tuple<std::string>("no_stabilization"),
          Teuchos::tuple<std::string>("Do not use any stabilization"),
          Teuchos::tuple<int>(
              INPAR::SCATRA::stabtype_no_stabilization),
              &eleparams.sublist("REINITIALIZATION"));
    // turn off artificial diffusion
    Teuchos::setStringToIntegralParameter<int>("ARTDIFFREINIT",
            "no",
            "potential incorporation of all-scale subgrid diffusivity (a.k.a. discontinuity-capturing) term",
            Teuchos::tuple<std::string>("no"),
            Teuchos::tuple<std::string>("no artificial diffusion"),
            Teuchos::tuple<int>(INPAR::SCATRA::artdiff_none),
            &eleparams.sublist("REINITIALIZATION"));

    // call standard loop over elements
    discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

    // note: time-integration parameter list has not to be overwritten here, since we rely on incremental solve
    //       as already set in PrepareTimeLoopReinit()

    // compute time derivative of phi at pseudo-time tau=0
    CalcInitialPhidt();

    // eventually, undo changes in general parameter list
    SetReinitializationElementParameters();
  }

  return;
}


/*----------------------------------------------------------------------*
 | set part of the residual vector belonging to the old timestep        |
 |                                                      rasthofer 12/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetTimIntOneStepTheta::SetOldPartOfRighthandside()
{
  if (not switchreinit_)
    TimIntOneStepTheta::SetOldPartOfRighthandside();
  else
  // hist_ = phin_ + dt*(1-Theta)*phidtn_
   hist_->Update(1.0, *phin_, dtau_*(1.0-thetareinit_), *phidtn_, 0.0);

  return;
}


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                      rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetTimIntOneStepTheta::Update(const int num)
{
  if (not switchreinit_)
  {
    // compute time derivative at time n+1
    ComputeTimeDerivative();
    
    // compute flux vector field for later output BEFORE time shift of results
    // is performed below !!
    if (writeflux_!=INPAR::SCATRA::flux_no)
    {
      if (DoOutput() or DoBoundaryFluxStatistics())
        flux_ = CalcFlux(true, num);
    }

    // after the next command (time shift of solutions) do NOT call
    // ComputeTimeDerivative() anymore within the current time step!!!

    // solution of this step becomes most recent solution of the last step
    phin_ ->Update(1.0,*phinp_,0.0);

    // time deriv. of this step becomes most recent time derivative of
    // last step
    phidtn_->Update(1.0,*phidtnp_,0.0);
  }
  else
  {
    // solution of this step becomes most recent solution of the last step
    phin_ ->Update(1.0,*phinp_,0.0);

    // reinitialization is done, reset flag
    switchreinit_ = false;

    // compute time derivative at time n (and n+1)

    // we also have reset the time-integration parameter list for two reasons
    // 1: use of reinitialization equation overwrites time-integration parameter list (this is corrected afterwards)
    // 2: incremental solver has to be overwritten if used
    Teuchos::ParameterList eletimeparams;

    eletimeparams.set<int>("action",SCATRA::set_time_parameter);
    // set type of scalar transport problem (after preevaluate evaluate, which need scatratype is called)
    eletimeparams.set<int>("scatratype",scatratype_);

    eletimeparams.set<bool>("using generalized-alpha time integration",false);
    eletimeparams.set<bool>("using stationary formulation",false);
    eletimeparams.set<bool>("incremental solver",true); // this is important to have here

    eletimeparams.set<double>("time-step length",dta_);
    eletimeparams.set<double>("total time",time_);
    eletimeparams.set<double>("time factor",theta_*dta_);
    eletimeparams.set<double>("alpha_F",1.0);

    // call standard loop over elements
    discret_->Evaluate(eletimeparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

    CalcInitialPhidt();
    // reset element time-integration parameters
    SetElementTimeParameter();
  }

  // update also particle field
  if (particle_ != Teuchos::null)
    particle_->TransferAndUpdate();

  return;
}


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 | used within reinitialization loop                    rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetTimIntOneStepTheta::UpdateReinit()
{
  //TODO: Fkt hier raus nehmen
  // compute time derivative at time n+1
  // time derivative of phi:
  // phidt(n+1) = (phi(n+1)-phi(n)) / (theta*dt) + (1-(1/theta))*phidt(n)
  const double fact1 = 1.0/(thetareinit_*dtau_);
  const double fact2 = 1.0 - (1.0/thetareinit_);
  phidtnp_->Update(fact2,*phidtn_,0.0);
  phidtnp_->Update(fact1,*phinp_,-fact1,*phin_,1.0);

  // solution of this step becomes most recent solution of the last step
  phin_ ->Update(1.0,*phinp_,0.0);

  // time deriv. of this step becomes most recent time derivative of
  // last step
  phidtn_->Update(1.0,*phidtnp_,0.0);

  return;
}


/*--------------------------------------------------------------------------------------------*
 | Redistribute the scatra discretization and vectors according to nodegraph  rasthofer 07/11 |
 |                                                                            DA wichmann     |
 *--------------------------------------------------------------------------------------------*/
void SCATRA::LevelSetTimIntOneStepTheta::Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph)
{
  // let the base class do the basic redistribution and transfer of the base class members
  LevelSetAlgorithm::Redistribute(nodegraph);

  // now do all the ost specfic steps
  const Epetra_Map* newdofrowmap = discret_->DofRowMap();
  Teuchos::RCP<Epetra_Vector> old;

  if (fsphinp_ != Teuchos::null)
  {
    old = fsphinp_;
    fsphinp_ = LINALG::CreateVector(*newdofrowmap,true);
    LINALG::Export(*old, *fsphinp_);
  }

  return;
}
