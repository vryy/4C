/*----------------------------------------------------------------------*/
/*!
\file adapter_scatra_base_algorithm.cpp

\brief scalar transport field base algorithm

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "../drt_io/io_control.H"
#include "adapter_scatra_base_algorithm.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_scatra/scatra_timint_stat.H"
#include "../drt_scatra/scatra_timint_ost.H"
#include "../drt_scatra/scatra_timint_bdf2.H"
#include "../drt_scatra/scatra_timint_genalpha.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraBaseAlgorithm::ScaTraBaseAlgorithm(const Teuchos::ParameterList& prbdyn)
{
  /// setup scalar transport algorithm (overriding some dynamic parameters with
  /// values specified in given problem-dependent ParameterList prbdyn)

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numscatra,0);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  RCP<IO::DiscretizationWriter> output =
    rcp(new IO::DiscretizationWriter(actdis));
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------

  const Teuchos::ParameterList& scatradyn =
    DRT::Problem::Instance()->ScalarTransportDynamicParams();

  // print out default parameters of scalar tranport parameter list
  if (actdis->Comm().MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, scatradyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RCP<LINALG::Solver> solver =
    rcp(new LINALG::Solver(DRT::Problem::Instance()->ScalarTransportSolverParams(),
                           actdis->Comm(),
                           DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  RCP<ParameterList> scatratimeparams= rcp(new ParameterList());

  // ----problem type (type of scalar transport problem we want to solve)
  scatratimeparams->set<string>("problem type",DRT::Problem::Instance()->ProblemType());

  // --------------------type of time-integration (or stationary) scheme
  INPUTPARAMS::ScaTraTimeIntegrationScheme timintscheme =
    Teuchos::getIntegralValue<INPUTPARAMS::ScaTraTimeIntegrationScheme>(scatradyn,"TIMEINTEGR");
  scatratimeparams->set<INPUTPARAMS::ScaTraTimeIntegrationScheme>("time int algo",timintscheme);

  // --------------------------------------- time integration parameters
  // the default time step size
  scatratimeparams->set<double>   ("time step size"           ,prbdyn.get<double>("TIMESTEP"));
  // maximum simulation time
  scatratimeparams->set<double>   ("total time"               ,prbdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  scatratimeparams->set<int>      ("max number timesteps"     ,prbdyn.get<int>("NUMSTEP"));

  // ----------------------------------------------- restart and output
  // restart
  scatratimeparams->set           ("write restart every"       ,prbdyn.get<int>("RESTARTEVRY"));
  // solution output
  scatratimeparams->set           ("write solution every"      ,prbdyn.get<int>("UPRES"));
  //scatratimeparams->set           ("write solution every"      ,prbdyn.get<int>("WRITESOLEVRY"));
  // write also flux vectors when solution is written out?
  scatratimeparams->set<string>   ("write flux"   ,scatradyn.get<string>("WRITEFLUX"));

  // ---------------------------------------------------- initial field
  scatratimeparams->set<int>("scalar initial field" ,Teuchos::getIntegralValue<int>(scatradyn,"INITIALFIELD"));
  scatratimeparams->set<int>("scalar initial field func number",scatradyn.get<int>("INITFUNCNO"));

  // ----------------------------------------------------velocity field
  scatratimeparams->set<int>("velocity field" ,Teuchos::getIntegralValue<int>(scatradyn,"VELOCITYFIELD"));
  scatratimeparams->set<int>("velocity function number",scatradyn.get<int>("VELFUNCNO"));

  // -------------------- compute error compared to analytical solution
  scatratimeparams->set<int>("CALCERROR",Teuchos::getIntegralValue<int>(scatradyn,"CALCERROR"));

  // -------------------------------- (fine-scale) subgrid diffusivity?
  scatratimeparams->set<string>("fs subgrid diffusivity",scatradyn.get<string>("FSSUGRVISC"));

  // -------------------- block preconditioning (only supported by ELCH)
  scatratimeparams->set<int>("BLOCKPRECOND",Teuchos::getIntegralValue<int>(scatradyn,"BLOCKPRECOND"));

  // -----------------------sublist containing stabilization parameters
  scatratimeparams->sublist("STABILIZATION")=scatradyn.sublist("STABILIZATION");

  // ----------------sublist containing parameters for newton iteration
  scatratimeparams->sublist("NONLINEAR") = scatradyn.sublist("NONLINEAR");

  // --------------sublist for combustion-specific gfunction parameters
  /* This sublist COMBUSTION DYNAMIC/GFUNCTION contains parameters for the gfunction field
   * which are only relevant for a combustion problem.                         07/08 henke */
  if (genprob.probtyp == prb_combust)
  {
    scatratimeparams->sublist("COMBUSTION GFUNCTION")=prbdyn.sublist("COMBUSTION GFUNCTION");
  }

  // -------------------sublist for electrochemistry-specific parameters
  if (genprob.probtyp == prb_elch)
  {
    scatratimeparams->set<double>("TEMPERATURE",prbdyn.get<double>("TEMPERATURE"));

    // -------------------------------------------------------------------
    // create a 2nd solver for block-preconditioning if chosen from input
    // -------------------------------------------------------------------
    if (scatratimeparams->get<int>("BLOCKPRECOND"))
    {
      // switch to the SIMPLE(R) algorithms
      solver->PutSolverParamsToSubParams("SIMPLER",
         DRT::Problem::Instance()->ScalarTransportElectricPotentialSolverParams());
    }
  }

  // -------------------------------------------------------------------
  // additional parameters and algorithm construction depending on
  // respective time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  if(timintscheme == INPUTPARAMS::timeint_stationary)
  {
    //------------------------------------------------------------------
    // create instance of time integration class (call the constructor)
    //------------------------------------------------------------------
    scatra_ = rcp(new SCATRA::TimIntStationary::TimIntStationary(actdis, solver, scatratimeparams, output));
  }
  else if (timintscheme == INPUTPARAMS::timeint_one_step_theta)
  {
    // -----------------------------------------------------------------
    // set additional parameters in list for one-step-theta scheme
    // -----------------------------------------------------------------

    // parameter theta for time-integration schemes
    scatratimeparams->set<double>("theta",scatradyn.get<double>("THETA"));

    //------------------------------------------------------------------
    // create instance of time integration class (call the constructor)
    //------------------------------------------------------------------
    scatra_ = rcp(new SCATRA::TimIntOneStepTheta::TimIntOneStepTheta(actdis, solver, scatratimeparams, output));
  }
  else if (timintscheme == INPUTPARAMS::timeint_bdf2)
  {
    //------------------------------------------------------------------
    // create instance of time integration class (call the constructor)
    //------------------------------------------------------------------
    scatra_ = rcp(new SCATRA::TimIntBDF2::TimIntBDF2(actdis, solver, scatratimeparams, output));
  }
  else if (timintscheme == INPUTPARAMS::timeint_gen_alpha)
  {
    // -------------------------------------------------------------------
    // set additional parameters in list for generalized-alpha scheme
    // -------------------------------------------------------------------
    // parameter alpha_M for for generalized-alpha scheme
    scatratimeparams->set<double>("alpha_M",scatradyn.get<double>("ALPHA_M"));
    // parameter alpha_F for for generalized-alpha scheme
    scatratimeparams->set<double>("alpha_F",scatradyn.get<double>("ALPHA_F"));
    // parameter gamma for for generalized-alpha scheme
    scatratimeparams->set<double>("gamma",  scatradyn.get<double>("GAMMA"));

    //------------------------------------------------------------------
    // create instance of time integration class (call the constructor)
    //------------------------------------------------------------------
    scatra_ = rcp(new SCATRA::TimIntGenAlpha::TimIntGenAlpha(actdis, solver, scatratimeparams, output));
  }
  else
    dserror("Unknown time-integration scheme for scalar tranport problem");

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraBaseAlgorithm::~ScaTraBaseAlgorithm()
{
}


#endif  // #ifdef CCADISCRET
