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

#include "adapter_scatra_base_algorithm.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_solver.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_elch.H"
#include "../drt_lib/drt_globalproblem.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_scatra/scatra_timint_stat.H"
#include "../drt_scatra/scatra_timint_ost.H"
#include "../drt_scatra/scatra_timint_bdf2.H"
#include "../drt_scatra/scatra_timint_genalpha.H"
#include "../drt_scatra/scatra_timint_tg.H"
#include "../drt_scatra/scatra_resulttest.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraBaseAlgorithm::ScaTraBaseAlgorithm(
    const Teuchos::ParameterList& prbdyn,
    bool isale,
    const std::string disname,
    const Teuchos::ParameterList& solverparams,
    const bool reinitswitch
)
{
  // setup scalar transport algorithm (overriding some dynamic parameters
  // with values specified in given problem-dependent ParameterList prbdyn)

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->GetDis(disname);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  RCP<IO::DiscretizationWriter> output =
    Teuchos::rcp(new IO::DiscretizationWriter(actdis));
  if(!reinitswitch) output->WriteMesh(0,0.0); // don't write mesh for level set reinitialization time loops

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  const Teuchos::ParameterList& scatradyn =
    DRT::Problem::Instance()->ScalarTransportDynamicParams();

  // print out default parameters of scalar transport parameter lists
  if (actdis->Comm().MyPID()==0)
  {
    DRT::INPUT::PrintDefaultParameters(IO::cout, scatradyn);
    DRT::INPUT::PrintDefaultParameters(IO::cout, scatradyn.sublist("STABILIZATION"));
    DRT::INPUT::PrintDefaultParameters(IO::cout, scatradyn.sublist("NONLINEAR"));
    /*
    const Teuchos::ParameterList& solverparams =
        DRT::Problem::Instance()->ScalarTransportFluidSolverParams();
    DRT::INPUT::PrintDefaultParameters(IO::cout, solverparams);
    */
  }

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // TODO: TAW use of solverparams??? change input parameter to solver number instead of parameter list? -> no default paramter possible any more
  RCP<LINALG::Solver> solver =
    Teuchos::rcp(new LINALG::Solver(solverparams,
                           actdis->Comm(),
                           DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  // make a copy (inside an Teuchos::rcp) containing also all sublists
  RCP<Teuchos::ParameterList> scatratimeparams= Teuchos::rcp(new Teuchos::ParameterList(scatradyn));

  // -------------------------------------------------------------------
  // overrule certain parameters for coupled problems
  // -------------------------------------------------------------------
  // the default time step size
  scatratimeparams->set<double>   ("TIMESTEP"    ,prbdyn.get<double>("TIMESTEP"));
  // maximum simulation time
  scatratimeparams->set<double>   ("MAXTIME"     ,prbdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  scatratimeparams->set<int>      ("NUMSTEP"     ,prbdyn.get<int>("NUMSTEP"));
  // restart
  scatratimeparams->set           ("RESTARTEVRY" ,prbdyn.get<int>("RESTARTEVRY"));
  // solution output
  scatratimeparams->set           ("UPRES"       ,prbdyn.get<int>("UPRES"));

  // -------------------------------------------------------------------
  // overrule flag for form of convective term as well as stabilization
  // for solid-based scalar transport in FS3I-type problems, to allow
  // for simultaneously using convective formulation and stabilization
  // in fluid-based scalar transport, while conservative formulation
  // and no stabilization is mandatorily used in solid-based scalar
  // transport, for the time being
  // (assumed disname = "scatra2" for solid-based scalar transport)
  // -------------------------------------------------------------------
  if (disname == "scatra2")
  {
    scatratimeparams->set<string>("CONVFORM","conservative");
    scatratimeparams->sublist("STABILIZATION").set<string>("STABTYPE","no_stabilization");
    scatratimeparams->sublist("STABILIZATION").set<string>("DEFINITION_TAU","Zero");

    // some provisions not yet activated
    /*if (DRT::INPUT::IntegralValue<INPAR::SCATRA::StabType>(scatratimeparams->sublist("STABILIZATION"),"STABTYPE") == INPAR::SCATRA::stabtype_SUPG)
      scatratimeparams->sublist("STABILIZATION").set<string>("STABTYPE","USFEM");

    if (DRT::INPUT::IntegralValue<INPAR::SCATRA::TauType>(scatratimeparams->sublist("STABILIZATION"),"DEFINITION_TAU") == INPAR::SCATRA::tau_franca_valentin)
      scatratimeparams->sublist("STABILIZATION").set<string>("DEFINITION_TAU","Franca_Madureira_Valentin");
    else if (DRT::INPUT::IntegralValue<INPAR::SCATRA::TauType>(scatratimeparams->sublist("STABILIZATION"),"DEFINITION_TAU") == INPAR::SCATRA::tau_franca_valentin_wo_dt or
     DRT::INPUT::IntegralValue<INPAR::SCATRA::TauType>(scatratimeparams->sublist("STABILIZATION"),"DEFINITION_TAU") == INPAR::SCATRA::tau_exact_1d)
      scatratimeparams->sublist("STABILIZATION").set<string>("DEFINITION_TAU","Franca_Madureira_Valentin_wo_dt");*/
  }

  // -------------------------------------------------------------------
  // list for extra parameters
  // (put here everything that is not available in scatradyn or its sublists)
  // -------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> extraparams
    = Teuchos::rcp(new Teuchos::ParameterList());

  // ------------------------------pointer to the error file (for output)
  extraparams->set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

  // ----------------Eulerian or ALE formulation of transport equation(s)
  extraparams->set<bool>("isale",isale);

  // --------------sublist for combustion-specific gfunction parameters
  /* This sublist COMBUSTION DYNAMIC/GFUNCTION contains parameters for the gfunction field
   * which are only relevant for a combustion problem.                         07/08 henke */
  if (DRT::Problem::Instance()->ProblemType() == prb_combust)
  {
    extraparams->sublist("COMBUSTION GFUNCTION")=prbdyn.sublist("COMBUSTION GFUNCTION");

    if(reinitswitch==true)
    {
      extraparams->sublist("COMBUSTION PDE REINITIALIZATION")=prbdyn.sublist("COMBUSTION PDE REINITIALIZATION");
      extraparams->set<bool>("REINITSWITCH", true);
    }
  }

  // -------------------sublist for electrochemistry-specific parameters
  if (DRT::Problem::Instance()->ProblemType() == prb_elch)
  {
    // temperature of electrolyte solution
    extraparams->set<double>("TEMPERATURE",prbdyn.get<double>("TEMPERATURE"));

    // we provide all available electrochemistry-related parameters
    extraparams->sublist("ELCH CONTROL")=prbdyn;

    // create a 2nd solver for block-preconditioning if chosen from input
    if (DRT::INPUT::IntegralValue<int>(scatradyn,"BLOCKPRECOND"))
    {
      // set Inverse1 block (for primary variable), use Fluid Scatra Solver
      Teuchos::ParameterList& inv1 = solver->Params().sublist("Inverse1");
      inv1 = solver->Params();
      inv1.remove("SIMPLER",false);
      inv1.remove("Inverse1",false);
      // set Inverse2 block (for secondary variable), use ScalarTransportElectricPotential Solver
      // get the solver number used for SIMPLER SOLVER
      const int linsolvernumber_simpler = scatradyn.get<int>("SIMPLER_SOLVER");
      if (linsolvernumber_simpler == (-1))
        dserror("no SIMPLER_SOLVER number set for ELCH problem (solved with SIMPLER). Please set SIMPLER_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");
      solver->PutSolverParamsToSubParams("Inverse2", DRT::Problem::Instance()->SolverParams(linsolvernumber_simpler));
      // use CheapSIMPLE preconditioner (hardwired; change me for others)
      solver->Params().sublist("CheapSIMPLE Parameters").set("Prec Type","CheapSIMPLE");
      solver->Params().set("ELCH",true); // internal CheapSIMPLE modus for ML null space computation

      // print unused solver parameters to screen
      /*
      if (actdis->Comm().MyPID()==0)
      {
        const Teuchos::ParameterList& solverparams2
        = DRT::Problem::Instance()->SolverParams(linsolvernumber_simpler);
        DRT::INPUT::PrintDefaultParameters(IO::cout, solverparams2);
      }
      */
    }
  }

  // ------------------------------------get also fluid turbulence sublist
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  extraparams->sublist("TURBULENCE MODEL")=fdyn.sublist("TURBULENCE MODEL");
  extraparams->sublist("SUBGRID VISCOSITY")=fdyn.sublist("SUBGRID VISCOSITY");
  extraparams->sublist("MULTIFRACTAL SUBGRID SCALES")=fdyn.sublist("MULTIFRACTAL SUBGRID SCALES");
  extraparams->sublist("TURBULENT INFLOW")=fdyn.sublist("TURBULENT INFLOW");

  // ----------------------------- add some loma specific parameters
  // get also scatra stabilization sublist
  const Teuchos::ParameterList& lomadyn =
    DRT::Problem::Instance()->LOMAControlParams();
  extraparams->sublist("LOMA").set<bool>("update material",DRT::INPUT::IntegralValue<int>(lomadyn,"SGS_MATERIAL_UPDATE"));

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // respective time-integration (or stationary) scheme
  // -------------------------------------------------------------------
   INPAR::SCATRA::TimeIntegrationScheme timintscheme =
     DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(scatradyn,"TIMEINTEGR");

   if (reinitswitch==false)
   {
     switch(timintscheme)
     {
     case INPAR::SCATRA::timeint_stationary:
     {
       // create instance of time integration class (call the constructor)
       scatra_ = Teuchos::rcp(new SCATRA::TimIntStationary(actdis, solver, scatratimeparams, extraparams, output));
       break;
     }
     case INPAR::SCATRA::timeint_one_step_theta:
     {
       // create instance of time integration class (call the constructor)
       scatra_ = Teuchos::rcp(new SCATRA::TimIntOneStepTheta(actdis, solver, scatratimeparams, extraparams,output));
       break;
     }
     case INPAR::SCATRA::timeint_bdf2:
     {
       // create instance of time integration class (call the constructor)
       scatra_ = Teuchos::rcp(new SCATRA::TimIntBDF2(actdis, solver, scatratimeparams,extraparams, output));
       break;
     }
     case INPAR::SCATRA::timeint_gen_alpha:
     {
       // create instance of time integration class (call the constructor)
       scatra_ = Teuchos::rcp(new SCATRA::TimIntGenAlpha(actdis, solver, scatratimeparams,extraparams, output));
       break;
     }
     case INPAR::SCATRA::timeint_tg2:
     {
       // create instance of time integration class (call the constructor)
       scatra_ = Teuchos::rcp(new SCATRA::TimIntTaylorGalerkin(actdis, solver, scatratimeparams,extraparams, output));
       break;
     }
     case INPAR::SCATRA::timeint_tg3:
     {
       // create instance of time integration class (call the constructor)
       scatra_ = Teuchos::rcp(new SCATRA::TimIntTaylorGalerkin(actdis, solver, scatratimeparams,extraparams, output));
       break;
     }
     default:
       dserror("Unknown time-integration scheme for scalar transport problem");
       break;
     }// switch(timintscheme)
   }
   else{
     // -------------------------------------------------------------------
     // algorithm construction depending on
     // respective time-integration (or stationary) scheme
     // -------------------------------------------------------------------
     INPAR::SCATRA::TimeIntegrationScheme timintscheme_reinitialization = DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(prbdyn.sublist("COMBUSTION PDE REINITIALIZATION"),"REINIT_TIMEINTEGR");

     switch(timintscheme_reinitialization)
     {
     case INPAR::SCATRA::timeint_one_step_theta:
     {
       // create instance of time integration class (call the constructor)
       scatra_ = Teuchos::rcp(new SCATRA::TimIntOneStepTheta(actdis, solver, scatratimeparams, extraparams,output));
       break;
     }
     case INPAR::SCATRA::timeint_tg2:
     {
       // create instance of time integration class (call the constructor)
       scatra_ = Teuchos::rcp(new SCATRA::TimIntTaylorGalerkin(actdis, solver, scatratimeparams,extraparams, output));
       break;
     }
     default:
       dserror("Unknown time-integration scheme for reinitialization problem");
       break;
     }// switch(timintscheme_reinitialization)

   } // switch(reinitswitch)

  return;

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraBaseAlgorithm::~ScaTraBaseAlgorithm()
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntImpl& ADAPTER::ScaTraBaseAlgorithm::ScaTraField()
{
  return *scatra_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::ScaTraBaseAlgorithm::CreateScaTraFieldTest()
{
  return Teuchos::rcp(new SCATRA::ScaTraResultTest(*scatra_));
}


