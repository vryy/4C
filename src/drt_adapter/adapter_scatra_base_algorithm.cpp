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
#include "../drt_io/io.H"
#include "adapter_scatra_base_algorithm.H"
#include "../linalg/linalg_solver.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_elch.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_stat.H"
#include "../drt_scatra/scatra_timint_ost.H"
#include "../drt_scatra/scatra_timint_bdf2.H"
#include "../drt_scatra/scatra_timint_genalpha.H"
#include "../drt_scatra/scatra_timint_tg.H"             //schott 05/11
//#include "../drt_scatra/scatra_timint_reinitialization.H"       //schott 05/11
#include "../drt_scatra/scatra_resulttest.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraBaseAlgorithm::ScaTraBaseAlgorithm(
    const Teuchos::ParameterList& prbdyn,
    bool isale,
    const int disnum,
    const Teuchos::ParameterList& solverparams
)
{
  // setup scalar transport algorithm (overriding some dynamic parameters
  // with values specified in given problem-dependent ParameterList prbdyn)

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numscatra,disnum);

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

  // print out default parameters of scalar transport parameter lists
  if (actdis->Comm().MyPID()==0)
  {
    DRT::INPUT::PrintDefaultParameters(std::cout, scatradyn);
    DRT::INPUT::PrintDefaultParameters(std::cout, scatradyn.sublist("STABILIZATION"));
    DRT::INPUT::PrintDefaultParameters(std::cout, scatradyn.sublist("NONLINEAR"));
    /*
    const Teuchos::ParameterList& solverparams =
        DRT::Problem::Instance()->ScalarTransportFluidSolverParams();
    DRT::INPUT::PrintDefaultParameters(std::cout, solverparams);
    */
  }

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RCP<LINALG::Solver> solver =
    rcp(new LINALG::Solver(solverparams,   //DRT::Problem::Instance()->ScalarTransportFluidSolverParams(),
                           actdis->Comm(),
                           DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  // make a copy (inside an rcp) containing also all sublists
  RCP<ParameterList> scatratimeparams= rcp(new ParameterList(scatradyn));

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
  // list for extra parameters
  // (put here everything that is not available in scatradyn or its sublists)
  // -------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> extraparams
    = Teuchos::rcp(new Teuchos::ParameterList());

  // ----problem type (type of scalar transport problem we want to solve)
  extraparams->set<string>("problem type",DRT::Problem::Instance()->ProblemType());

  // ------------------------------pointer to the error file (for output)
  extraparams->set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

  // ----------------Eulerian or ALE formulation of transport equation(s)
  extraparams->set<bool>("isale",isale);

  // --------------sublist for combustion-specific gfunction parameters
  /* This sublist COMBUSTION DYNAMIC/GFUNCTION contains parameters for the gfunction field
   * which are only relevant for a combustion problem.                         07/08 henke */
  if (genprob.probtyp == prb_combust)
  {
    extraparams->sublist("COMBUSTION GFUNCTION")=prbdyn.sublist("COMBUSTION GFUNCTION");
  }

  // -------------------sublist for electrochemistry-specific parameters
  if (genprob.probtyp == prb_elch)
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
      solver->PutSolverParamsToSubParams("Inverse2", DRT::Problem::Instance()->ScalarTransportElectricPotentialSolverParams());
      // use CheapSIMPLE preconditioner (hardwired; change me for others)
      solver->Params().sublist("CheapSIMPLE Parameters").set("Prec Type","CheapSIMPLE");
      solver->Params().set("ELCH",true); // internal CheapSIMPLE modus for ML null space computation

      // print unused solver parameters to screen
      /*
      if (actdis->Comm().MyPID()==0)
      {
        const Teuchos::ParameterList& solverparams2
        = DRT::Problem::Instance()->ScalarTransportElectricPotentialSolverParams();
        DRT::INPUT::PrintDefaultParameters(std::cout, solverparams2);
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

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // respective time-integration (or stationary) scheme
  // -------------------------------------------------------------------
   INPAR::SCATRA::TimeIntegrationScheme timintscheme =
     DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(scatradyn,"TIMEINTEGR");

   switch(timintscheme)
   {
   case INPAR::SCATRA::timeint_stationary:
   {
     // create instance of time integration class (call the constructor)
     scatra_ = rcp(new SCATRA::TimIntStationary(actdis, solver, scatratimeparams, extraparams, output));
     break;
   }
   case INPAR::SCATRA::timeint_one_step_theta:
   {
     // create instance of time integration class (call the constructor)
     scatra_ = rcp(new SCATRA::TimIntOneStepTheta(actdis, solver, scatratimeparams, extraparams,output));
     break;
   }
   case INPAR::SCATRA::timeint_bdf2:
   {
     // create instance of time integration class (call the constructor)
     scatra_ = rcp(new SCATRA::TimIntBDF2(actdis, solver, scatratimeparams,extraparams, output));
     break;
   }
   case INPAR::SCATRA::timeint_gen_alpha:
   {
     // create instance of time integration class (call the constructor)
     scatra_ = rcp(new SCATRA::TimIntGenAlpha(actdis, solver, scatratimeparams,extraparams, output));
     break;
   }
   case INPAR::SCATRA::timeint_tg2: //schott 05/11
   {
     // create instance of time integration class (call the constructor)
	 scatra_ = rcp(new SCATRA::TimIntTaylorGalerkin(actdis, solver, scatratimeparams,extraparams, output));
     break;
   }
   case INPAR::SCATRA::timeint_tg2_LW: //schott 05/11
   {
     // create instance of time integration class (call the constructor)
	 scatra_ = rcp(new SCATRA::TimIntTaylorGalerkin(actdis, solver, scatratimeparams,extraparams, output));
     break;
   }
   case INPAR::SCATRA::timeint_tg3: //schott 05/11
   {
     // create instance of time integration class (call the constructor)
	 scatra_ = rcp(new SCATRA::TimIntTaylorGalerkin(actdis, solver, scatratimeparams,extraparams, output));
     break;
   }
   case INPAR::SCATRA::timeint_tg4_leapfrog: //schott 05/11
   {
     // create instance of time integration class (call the constructor)
     scatra_ = rcp(new SCATRA::TimIntTaylorGalerkin(actdis, solver, scatratimeparams,extraparams, output));
     break;
   }
   case INPAR::SCATRA::timeint_tg4_onestep: //schott 05/11
   {
     // create instance of time integration class (call the constructor)
     scatra_ = rcp(new SCATRA::TimIntTaylorGalerkin(actdis, solver, scatratimeparams,extraparams, output));
     break;
   }
   default:
     dserror("Unknown time-integration scheme for scalar transport problem");
   }// switch(timintscheme)


   if (genprob.probtyp == prb_combust)
   {

	   //TODO: Do this within the Combust Algorithm !!!

	   // -------------------------------------------------------------------
	   // create a solver
	   // -------------------------------------------------------------------
	   RCP<LINALG::Solver> solver_reinit =
	     rcp(new LINALG::Solver(solverparams, //DRT::Problem::Instance()->ScalarTransportSolverParams(),
	                            actdis->Comm(),
	                            DRT::Problem::Instance()->ErrorFile()->Handle()));
	   actdis->ComputeNullSpaceIfNecessary(solver_reinit->Params());


     extraparams->sublist("COMBUSTION GFUNCTION")=prbdyn.sublist("COMBUSTION GFUNCTION");
     extraparams->sublist("COMBUSTION PDE REINITIALIZATION")=prbdyn.sublist("COMBUSTION PDE REINITIALIZATION");
     extraparams->set<bool>("REINITSWITCH", true);
     // -------------------------------------------------------------------
     // set parameters in list required for all schemes
     // -------------------------------------------------------------------
     // make a copy (inside an rcp) containing also all sublists
     RCP<ParameterList> reinittimeparams= rcp(new ParameterList(scatradyn));


//     DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(*params,"REINIT_TIMEINTEGR")
     // -------------------------------------------------------------------
     // overrule certain parameters for coupled problems
     // -------------------------------------------------------------------
//     DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>
     // the default time step size
//     reinittimeparams->set<int>      ("TIMEINTEGR"  ,DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(prbdyn.sublist("COMBUSTION PDE REINITIALIZATION"),"REINIT_TIMEINTEGR"));
     // the default time step size
     reinittimeparams->set<double>   ("TIMESTEP"    ,prbdyn.sublist("COMBUSTION PDE REINITIALIZATION").get<double>("PSEUDOTIMESTEP_FACTOR"));
     // maximum simulation time
     reinittimeparams->set<double>   ("MAXTIME"     ,prbdyn.get<double>("MAXTIME"));
     // maximum number of timesteps
     reinittimeparams->set<int>      ("NUMSTEP"     ,prbdyn.sublist("COMBUSTION PDE REINITIALIZATION").get<int>("NUMPSEUDOSTEPS"));
     // restart
     reinittimeparams->set           ("RESTARTEVRY" ,prbdyn.get<int>("RESTARTEVRY"));
     // solution output
     reinittimeparams->set           ("UPRES"       ,prbdyn.get<int>("UPRES"));



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
       reinit_ = rcp(new SCATRA::TimIntOneStepTheta(actdis, solver_reinit, reinittimeparams, extraparams,output));
       break;
     }
     case INPAR::SCATRA::timeint_tg2:
     {
       // create instance of time integration class (call the constructor)
       reinit_ = rcp(new SCATRA::TimIntTaylorGalerkin(actdis, solver_reinit, reinittimeparams,extraparams, output));
       break;
     }
     default:
       dserror("Unknown time-integration scheme for reinitialization problem");
     }// switch(timintscheme_reinitialization)

   }

  return;

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntImpl& ADAPTER::ScaTraBaseAlgorithm::ScaTraField()
{
  return *scatra_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntImpl& ADAPTER::ScaTraBaseAlgorithm::ScaTraReinitField()
{
  return *reinit_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::ScaTraBaseAlgorithm::CreateScaTraFieldTest()
{
  return Teuchos::rcp(new SCATRA::ScaTraResultTest(*scatra_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraBaseAlgorithm::~ScaTraBaseAlgorithm()
{
}


#endif  // #ifdef CCADISCRET
