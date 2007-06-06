
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fsi_dirichletneumann.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;

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
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;


void FSI::DirichletNeumannCoupling::Setup()
{
  SetupStructure();
  SetupFluid();
  SetupAle();
}

void FSI::DirichletNeumannCoupling::SetupStructure()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numsf,0);

  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  DiscretizationWriter output(actdis);
  output.WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR*         actsolv  = &solv[genprob.numsf];
  FSI_DYNAMIC *fsidyn     = alldyn[3].fsidyn;
  STRUCT_DYNAMIC* sdyn     = alldyn[genprob.numsf].sdyn;

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  LINALG::Solver solver(solveparams,actdis->Comm(),allfiles.out_err);
  solver.TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // -------------------------------------------------------------------
  // create a generalized alpha time integrator
  // -------------------------------------------------------------------
  ParameterList genalphaparams;
  StruGenAlpha::SetDefaults(genalphaparams);

  genalphaparams.set<bool>  ("damping",sdyn->damp);
  genalphaparams.set<double>("damping factor K",sdyn->k_damp);
  genalphaparams.set<double>("damping factor M",sdyn->m_damp);

  genalphaparams.set<double>("beta",sdyn->beta);
  genalphaparams.set<double>("gamma",sdyn->gamma);
  genalphaparams.set<double>("alpha m",sdyn->alpha_m);
  genalphaparams.set<double>("alpha f",sdyn->alpha_f);

  genalphaparams.set<double>("total time",0.0);
  genalphaparams.set<double>("delta time",fsidyn->dt);
  genalphaparams.set<int>   ("step",0);
  genalphaparams.set<int>   ("nstep",fsidyn->nstep);
  genalphaparams.set<int>   ("max iterations",sdyn->maxiter);
  genalphaparams.set<int>   ("num iterations",-1);
  genalphaparams.set<double>("tolerance displacements",sdyn->toldisp);

  genalphaparams.set<bool>  ("io structural disp",ioflags.struct_disp);
  genalphaparams.set<int>   ("io disp every nstep",sdyn->updevry_disp);
  genalphaparams.set<bool>  ("io structural stress",ioflags.struct_stress);
  genalphaparams.set<int>   ("io disp every nstep",sdyn->updevry_stress);
  genalphaparams.set<bool>  ("print to screen",true);
  genalphaparams.set<bool>  ("print to err",true);
  genalphaparams.set<FILE*> ("err file",allfiles.out_err);

  // takes values "full newton" , "modified newton" , "nonlinear cg"
  genalphaparams.set<string>("equilibrium iteration","full newton");

  structure_ = rcp(new Structure(genalphaparams,*actdis,solver,output));
}


void FSI::DirichletNeumannCoupling::SetupFluid()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numff,0);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  DiscretizationWriter output(actdis);
  output.WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR        *actsolv  = &solv[genprob.numff];

  FSI_DYNAMIC *fsidyn     = alldyn[3].fsidyn;
  FLUID_DYNAMIC *fdyn     = alldyn[genprob.numff].fdyn;
  fdyn->step              =   0;
  fdyn->acttime           = 0.0;

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  LINALG::Solver solver(solveparams,actdis->Comm(),allfiles.out_err);
  solver.TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // -------------------------------------------------------------------
  // create a fluid nonlinear time integrator
  // -------------------------------------------------------------------
  ParameterList fluidtimeparams;
  FluidImplicitTimeInt::SetDefaults(fluidtimeparams);

  // number of degrees of freedom
  fluidtimeparams.set<int>              ("number of velocity degrees of freedom" ,genprob.ndim);
  // the default time step size
  fluidtimeparams.set<double>           ("time step size"           ,fsidyn->dt);
  // max. sim. time
  fluidtimeparams.set<double>           ("total time"               ,fsidyn->maxtime);
  // parameter for time-integration
  fluidtimeparams.set<double>           ("theta"                    ,fdyn->theta);
  // which kind of time-integration
  fluidtimeparams.set<FLUID_TIMEINTTYPE>("time int algo"            ,fdyn->iop);
  // bound for the number of timesteps
  fluidtimeparams.set<int>              ("max number timesteps"     ,fsidyn->nstep);
  // number of steps with start algorithm
  fluidtimeparams.set<int>              ("number of start steps"    ,fdyn->nums);
  // parameter for start algo
  fluidtimeparams.set<double>           ("start theta"              ,fdyn->thetas);


  // ---------------------------------------------- nonlinear iteration
  // maximum number of nonlinear iteration steps
  fluidtimeparams.set<int>             ("max nonlin iter steps"     ,fdyn->itemax);
  // stop nonlinear iteration when both incr-norms are below this bound
  fluidtimeparams.set<double>          ("tolerance for nonlin iter" ,fdyn->ittol);

  // restart
  fluidtimeparams.set                  ("write restart every"       ,fdyn->uprestart);

  //--------------------------------------------------
  // evaluate error for test flows with analytical solutions
  fluidtimeparams.set                  ("eval err for analyt sol"   ,fdyn->init);


  //--------------------------------------------------
  // create all vectors and variables associated with the time
  // integration (call the constructor)
  // the only parameter from the list required here is the number of
  // velocity degrees of freedom
  fluid_ = rcp(new Fluid(actdis, solver, fluidtimeparams, output));
}


void FSI::DirichletNeumannCoupling::SetupAle()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numaf,0);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  DiscretizationWriter output(actdis);
  output.WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR        *actsolv  = &solv[genprob.numaf];

  FSI_DYNAMIC *fsidyn     = alldyn[3].fsidyn;
  ALE_DYNAMIC *adyn     = alldyn[genprob.numaf].adyn;
  adyn->step            =   0;
  adyn->time            = 0.0;

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  LINALG::Solver solver(solveparams,actdis->Comm(),allfiles.out_err);
  solver.TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  ParameterList params;
  params.set<int>("nstep", fsidyn->nstep);
  params.set<double>("maxtime", fsidyn->maxtime);
  params.set<double>("dt", fsidyn->dt);

  ale_ = rcp(new AleLinear(actdis, solver, params, output));
}

void FSI::DirichletNeumannCoupling::Timeloop()
{
}

#endif
#endif
