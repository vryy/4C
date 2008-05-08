#ifdef CCADISCRET

#include "adapter_condif_base_algorithm.H"
#include "adapter_condif_genalpha.H"
// further includes for ConDifBaseAlgorithm:
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>


#include "../drt_fluid/condifimplicitintegration.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ConDifBaseAlgorithm::ConDifBaseAlgorithm(const Teuchos::ParameterList& prbdyn)
{
  /// setup convection diffusion algorithm (overriding some fluid parameters with
  /// values specified in given problem-dependent ParameterList)

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numcdf,0);

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
  SOLVAR        *actsolv  = &solv[genprob.numcdf];

  //const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  //const Teuchos::ParameterList& probsize = DRT::Problem::Instance()->ProblemSizeParams();
  //const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();

  if (actdis->Comm().MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, fdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RCP<ParameterList> solveparams = rcp(new ParameterList());
  RCP<LINALG::Solver> solver =
    rcp(new LINALG::Solver(solveparams,actdis->Comm(),allfiles.out_err));
  solver->TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  RCP<ParameterList> condiftimeparams= rcp(new ParameterList());

  // -----------------------------------------------condif initial field
  condiftimeparams->set<int>              ("condif initial field"     ,Teuchos::getIntegralValue<int>(fdyn,"CD_INITIALFIELD"));
  condiftimeparams->set<int>              ("condif initial field func number",fdyn.get<int>("CD_INITFUNCNO"));
  
  // -----------------------------------------------------velocity field
  condiftimeparams->set<int>              ("condif velocity field"     ,Teuchos::getIntegralValue<int>(fdyn,"CD_VELOCITY"));
  condiftimeparams->set<int>              ("condif velocity function number",fdyn.get<int>("CD_VELFUNCNO"));

  // -------------------------------------------------- time integration
  // the default time step size
  condiftimeparams->set<double>           ("time step size"           ,fdyn.get<double>("TIMESTEP"));
  // maximum simulation time
  condiftimeparams->set<double>           ("total time"               ,fdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  condiftimeparams->set<int>              ("max number timesteps"     ,fdyn.get<int>("NUMSTEP"));

  // ----------------------------------------------- restart and output
  // restart
  condiftimeparams->set                  ("write restart every"       ,fdyn.get<int>("RESTARTEVRY"));
  // solution output
  condiftimeparams->set                  ("write solution every"      ,fdyn.get<int>("UPRES"));

  // ---------------------------------(fine-scale) subgrid diffusivity?
  condiftimeparams->set<string>           ("fs subgrid diffusivity"   ,fdyn.get<string>("FSSUGRVISC"));


  // -------------------------------------------------------------------
  // additional parameters and algorithm call depending on respective
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  FLUID_TIMEINTTYPE iop = Teuchos::getIntegralValue<FLUID_TIMEINTTYPE>(fdyn,"TIMEINTEGR");
  if(iop == timeint_stationary or
     iop == timeint_one_step_theta or
     iop == timeint_bdf2
  )
  {
    // -----------------------------------------------------------------
    // set additional parameters in list for OST/BDF2/stationary scheme
    // -----------------------------------------------------------------
    // type of time-integration (or stationary) scheme
    condiftimeparams->set<FLUID_TIMEINTTYPE>("time int algo",iop);
    // parameter theta for time-integration schemes
    condiftimeparams->set<double>           ("theta"                    ,fdyn.get<double>("THETA"));
    // number of steps for potential start algorithm
    condiftimeparams->set<int>              ("number of start steps"    ,fdyn.get<int>("NUMSTASTEPS"));
    // parameter theta for potential start algorithm
    condiftimeparams->set<double>           ("start theta"              ,fdyn.get<double>("START_THETA"));

    //------------------------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor)
    //------------------------------------------------------------------
    condif_ = rcp(new ADAPTER::ConDifImplicit(actdis, solver, condiftimeparams, output));

#if 0     
    // initial field from restart
    if (probtype.get<int>("RESTART"))
    {
      // read the restart information, set vectors and variables
      condifimplicit.ReadRestart(probtype.get<int>("RESTART"));
    }

#endif
  }
  else if (iop == timeint_gen_alpha)
  {
    dserror("no adapter for generalized alpha condif dynamic routine implemented.");
    // -------------------------------------------------------------------
    // set additional parameters in list for generalized-alpha scheme
    // -------------------------------------------------------------------
    // parameter alpha_M for for generalized-alpha scheme
    condiftimeparams->set<double>           ("alpha_M"                  ,fdyn.get<double>("ALPHA_M"));
    // parameter alpha_F for for generalized-alpha scheme
    condiftimeparams->set<double>           ("alpha_F"                  ,fdyn.get<double>("ALPHA_F"));

    //------------------------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor)
    //------------------------------------------------------------------
    condif_ = rcp(new ADAPTER::ConDifGenAlpha(actdis, solver, condiftimeparams, output));
    
  }
  else
  {
    dserror("Unknown solver type for convection-diffusion");
  }

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ConDifBaseAlgorithm::~ConDifBaseAlgorithm()
{
}


#endif  // #ifdef CCADISCRET
