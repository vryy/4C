/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_base_algorithm.cpp

\brief Fluid Base Algorithm

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "adapter_fluid_base_algorithm.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_validparameters.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "adapter_fluid_impl.H"
#include "adapter_xfluid_impl.H"
#include "adapter_fluid_genalpha.H"
#include "adapter_fluid_combust.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidBaseAlgorithm::FluidBaseAlgorithm(const Teuchos::ParameterList& prbdyn, bool isale)
{
  SetupFluid(prbdyn, isale);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidBaseAlgorithm::~FluidBaseAlgorithm()
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidBaseAlgorithm::SetupFluid(const Teuchos::ParameterList& prbdyn, bool& isale)
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("ADAPTER::FluidBaseAlgorithm::SetupFluid");
  Teuchos::TimeMonitor monitor(*t);

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numff,0);

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
  //const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& probsize = DRT::Problem::Instance()->ProblemSizeParams();
  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();

  if (actdis->Comm().MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, fdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RCP<LINALG::Solver> solver =
    rcp(new LINALG::Solver(DRT::Problem::Instance()->FluidSolverParams(),
                           actdis->Comm(),
                           DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // -------------------------------------------------------------------
  // create a second solver for SIMPLER preconditioner if chosen from input
  // -------------------------------------------------------------------
  if (getIntegralValue<int>(fdyn,"SIMPLER"))
  {
    solver->PutSolverParamsToSubParams("SIMPLER",
                                       DRT::Problem::Instance()->FluidPressureSolverParams());
  }

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  RCP<ParameterList> fluidtimeparams = rcp(new ParameterList());

  fluidtimeparams->set<int>("Simple Preconditioner",Teuchos::getIntegralValue<int>(fdyn,"SIMPLER"));

  // -------------------------------------- number of degrees of freedom
  // number of degrees of freedom
  fluidtimeparams->set<int>              ("number of velocity degrees of freedom" ,probsize.get<int>("DIM"));

  // ---------------------------- low-Mach-number or incompressible flow
  fluidtimeparams->set<string>("low-Mach-number solver"   ,fdyn.get<string>("LOWMACH"));

  // ------------------------------------------------ basic scheme, i.e.
  // --------------------- solving nonlinear or linearised flow equation
  fluidtimeparams->set<int>("type of nonlinear solve" ,
                     Teuchos::getIntegralValue<int>(fdyn,"DYNAMICTYP"));

  // -------------------------------------------------- time integration
  // note: here, the values are taken out of the problem-dependent ParameterList prbdyn
  // (which also can be fluiddyn itself!)

  // the default time step size
  fluidtimeparams->set<double>           ("time step size"           ,prbdyn.get<double>("TIMESTEP"));
  // maximum simulation time
  fluidtimeparams->set<double>           ("total time"               ,prbdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  fluidtimeparams->set<int>              ("max number timesteps"     ,prbdyn.get<int>("NUMSTEP"));

  // ---------------------------------------------- nonlinear iteration
  // set linearisation scheme
  fluidtimeparams->set<string>           ("Linearisation"             ,fdyn.get<string>("NONLINITER"));
  // maximum number of nonlinear iteration steps
  fluidtimeparams->set<int>             ("max nonlin iter steps"     ,fdyn.get<int>("ITEMAX"));
  // stop nonlinear iteration when both incr-norms are below this bound
  fluidtimeparams->set<double>          ("tolerance for nonlin iter" ,fdyn.get<double>("CONVTOL"));
  // set convergence check
  fluidtimeparams->set<string>          ("CONVCHECK"  ,fdyn.get<string>("CONVCHECK"));
  // set adaptoive linear solver tolerance
  fluidtimeparams->set<bool>            ("ADAPTCONV",getIntegralValue<int>(fdyn,"ADAPTCONV")==1);
  fluidtimeparams->set<double>          ("ADAPTCONV_BETTER",fdyn.get<double>("ADAPTCONV_BETTER"));

  // ----------------------------------------------- restart and output
  // restart
  fluidtimeparams->set                  ("write restart every"       ,prbdyn.get<int>("RESTARTEVRY"));
  // solution output
  fluidtimeparams->set                  ("write solution every"      ,prbdyn.get<int>("UPRES"));
  // flag for writing stresses
  fluidtimeparams->set                  ("write stresses"            ,Teuchos::getIntegralValue<int>(ioflags,"FLUID_STRESS"));
  // ---------------------------------------------------- lift and drag
  fluidtimeparams->set<int>("liftdrag",Teuchos::getIntegralValue<int>(fdyn,"LIFTDRAG"));

  // -----------evaluate error for test flows with analytical solutions
  int init = Teuchos::getIntegralValue<int>(fdyn,"INITIALFIELD");
  fluidtimeparams->set                  ("eval err for analyt sol"   ,init);

  // ---------------------------- fine-scale subgrid viscosity approach
  fluidtimeparams->set<string>           ("fs subgrid viscosity"   ,fdyn.get<string>("FSSUGRVISC"));

  // -----------------------sublist containing stabilization parameters
  fluidtimeparams->sublist("STABILIZATION")=fdyn.sublist("STABILIZATION");

  // ------------------------------------------- Robin scheme parameters
  if (genprob.probtyp == prb_fsi)
  {
    INPUTPARAMS::FSIPartitionedCouplingMethod method =
      Teuchos::getIntegralValue<INPUTPARAMS::FSIPartitionedCouplingMethod>(prbdyn,"PARTITIONED");
    fluidtimeparams->set<bool>("fluidrobin",
        method==INPUTPARAMS::fsi_RobinNeumann or method==INPUTPARAMS::fsi_RobinRobin);
    fluidtimeparams->set<double>("alpharobinf",prbdyn.get<double>("ALPHA_F"));
  }

  // --------------------------sublist containing turbulence parameters
  {
    fluidtimeparams->sublist("TURBULENCE MODEL")=fdyn.sublist("TURBULENCE MODEL");

    fluidtimeparams->sublist("TURBULENCE MODEL").set<string>("statistics outfile",DRT::Problem::Instance()->OutputControlFile()->FileName());
  }

  // ----------------------------------------------- XFEM related stuff
  {
    const Teuchos::ParameterList& xdyn = DRT::Problem::Instance()->XFEMGeneralParams();
    fluidtimeparams->set<bool>("DLM_condensation", getIntegralValue<int>(xdyn,"DLM_CONDENSATION")==1 );
    fluidtimeparams->set<bool>("CONDEST", getIntegralValue<int>(xdyn,"CONDEST")==1 );
  }

  // --------------------------sublist for combustion-specific fluid parameters
  /* This sublist COMBUSTION DYNAMIC/FLUID contains parameters for the fluid field
   * which are only relevant for a combustion problem.                 07/08 henke */
  if (genprob.probtyp == prb_combust)
  {
    fluidtimeparams->sublist("COMBUSTION FLUID")=prbdyn.sublist("COMBUSTION FLUID");
  }

  // -------------------------------------------------------------------
  // additional parameters and algorithm call depending on respective
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  FLUID_TIMEINTTYPE iop = Teuchos::getIntegralValue<FLUID_TIMEINTTYPE>(fdyn,"TIMEINTEGR");

  // sanity checks and default flags
  if (genprob.probtyp == prb_fsi)
  {
    // in case of FSI calculations we do not want a stationary fluid solver
    if (iop == timeint_stationary)
      dserror("Stationary fluid solver not allowed for FSI.");

    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

    fluidtimeparams->set<bool>("interface second order", Teuchos::getIntegralValue<int>(fsidyn,"SECONDORDER"));
    fluidtimeparams->set<bool>("shape derivatives",
                               Teuchos::getIntegralValue<int>(fsidyn,"SHAPEDERIVATIVES"));

    int coupling = Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO");
    if (coupling == fsi_iter_monolithic or
        coupling == fsi_iter_monolithiclagrange or
        coupling == fsi_iter_monolithicstructuresplit)
    {
      // there are a couple of restrictions in monolithic FSI
      fluidtimeparams->set<bool>("do explicit predictor",false);
    }
  }
  // sanity checks and default flags
  if (genprob.probtyp == prb_fsi_xfem)
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    fluidtimeparams->set<bool>("interface second order", Teuchos::getIntegralValue<int>(fsidyn,"SECONDORDER"));

    int coupling = Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO");
    if (coupling == fsi_iter_monolithic or
        coupling == fsi_iter_monolithiclagrange or
        coupling == fsi_iter_monolithicstructuresplit)
    {
      // there are a couple of restrictions in monolithic FSI
      dserror("XFEM and monolithic FSI not tested!");
      fluidtimeparams->set<bool>("do explicit predictor",false);
    }
  }


  if(iop == timeint_stationary or
     iop == timeint_one_step_theta or
     iop == timeint_bdf2
    )
  {
    // -----------------------------------------------------------------
    // set additional parameters in list for OST/BDF2/stationary scheme
    // -----------------------------------------------------------------
    // type of time-integration (or stationary) scheme
    fluidtimeparams->set<FLUID_TIMEINTTYPE>("time int algo"            ,iop);
    // parameter theta for time-integration schemes
    fluidtimeparams->set<double>           ("theta"                    ,fdyn.get<double>("THETA"));
    // number of steps for potential start algorithm
    fluidtimeparams->set<int>              ("number of start steps"    ,fdyn.get<int>("NUMSTASTEPS"));
    // parameter theta for potential start algorithm
    fluidtimeparams->set<double>           ("start theta"              ,fdyn.get<double>("START_THETA"));
    // parameter for grid velocity interpolation
    fluidtimeparams->set<int>              ("order gridvel"            ,fdyn.get<int>("GRIDVEL"));

    fluidtimeparams->set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

    bool dirichletcond = true;
    if (genprob.probtyp == prb_fsi)
    {
      // FSI input parameters
      const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
      int coupling = Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO");
      if (coupling == fsi_iter_monolithic or
          coupling == fsi_iter_monolithiclagrange or
          coupling == fsi_iter_monolithicstructuresplit)
      {
        dirichletcond = false;
      }
      else
      {
        INPUTPARAMS::FSIPartitionedCouplingMethod method =
          Teuchos::getIntegralValue<INPUTPARAMS::FSIPartitionedCouplingMethod>(fsidyn,"PARTITIONED");
        if (method==INPUTPARAMS::fsi_RobinNeumann or
            method==INPUTPARAMS::fsi_RobinRobin)
          dirichletcond = false;
      }
    }

    //------------------------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor);
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom
    if (genprob.probtyp == prb_fsi_xfem)
    {
      // This is temporary until I found a solution how to model the non-exising interface mesh of a XFEM fluid
      // this is solely for the XFEM development and will go away (a.ger 04/08)
      RCP<DRT::Discretization> soliddis = null;
      soliddis = DRT::Problem::Instance()->Dis(genprob.numsf,0);

      if (!soliddis->Filled()) soliddis->FillComplete();

      fluid_ = rcp(new ADAPTER::XFluidImpl(actdis, soliddis, solver, fluidtimeparams, output, isale));
    }
    else if (genprob.probtyp == prb_combust)
    {
      cout << "Ich bin jetzt eine CombustFluidTimeInt" << &endl;
      //fluid_ = rcp(new ADAPTER::FluidImpl(actdis, solver, fluidtimeparams, output, isale));
      fluid_ = rcp(new ADAPTER::FluidCombust(actdis, solver, fluidtimeparams, output));
    }
    else
    {
      fluid_ = rcp(new ADAPTER::FluidImpl(actdis, solver, fluidtimeparams, output, isale, dirichletcond));
    }
  }
  else if (iop == timeint_gen_alpha)
  {
#if 1
    // -------------------------------------------------------------------
    // set additional parameters in list for generalized-alpha scheme
    // -------------------------------------------------------------------
    // parameter alpha_M for for generalized-alpha scheme
    fluidtimeparams->set<double>           ("alpha_M"                  ,fdyn.get<double>("ALPHA_M"));
    // parameter alpha_F for for generalized-alpha scheme
    fluidtimeparams->set<double>           ("alpha_F"                  ,fdyn.get<double>("ALPHA_F"));

    fluidtimeparams->set<double>           ("gamma"                    ,fdyn.get<double>("GAMMA"));
#else
    // -------------------------------------------------------------------
    // set additional parameters in list for generalized-alpha scheme
    // -------------------------------------------------------------------
    // parameter alpha_M for for generalized-alpha scheme
    fluidtimeparams->set<double>           ("alpha_M"                  ,1.-prbdyn.get<double>("ALPHA_M"));
    // parameter alpha_F for for generalized-alpha scheme
    fluidtimeparams->set<double>           ("alpha_F"                  ,1.-prbdyn.get<double>("ALPHA_F"));

    fluidtimeparams->set<double>           ("gamma"                    ,prbdyn.get<double>("GAMMA"));
#endif

    // create all vectors and variables associated with the time
    // integration (call the constructor);
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom
    fluid_ = rcp(new ADAPTER::FluidGenAlpha(actdis, solver, fluidtimeparams, output, isale));
  }
  else
  {
    dserror("Unknown time integration for fluid\n");
  }

  // set initial field by given function
  // we do this here, since we have direct access to all necessary parameters
  if(init>0)
  {
    int startfuncno = fdyn.get<int>("STARTFUNCNO");
    if (init!=2 and init!=3)
    {
      startfuncno=-1;
    }
    fluid_->SetInitialFlowField(init,startfuncno);
  }

  return;
}

#endif  // #ifdef CCADISCRET
