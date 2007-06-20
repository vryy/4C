/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include <ctime>
#include <cstdlib>
#include <iostream>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "stru_dyn_nln_drt.H"
#include "strugenalpha.H"
#include "../io/io_drt.H"
#include "drt_globalproblem.H"

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

/*----------------------------------------------------------------------*
  | structural nonlinear dynamics (gen-alpha)              m.gee 12/06  |
 *----------------------------------------------------------------------*/
void dyn_nlnstructural_drt()
{
  DSTraceHelper dst("dyn_nlnstructural_drt");

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

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR*         actsolv  = &solv[0];
  STRUCT_DYNAMIC* sdyn     = alldyn[0].sdyn;

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
  genalphaparams.set<double>("delta time",sdyn->dt);
  genalphaparams.set<int>   ("step",0);
  genalphaparams.set<int>   ("nstep",sdyn->nstep);
  genalphaparams.set<int>   ("max iterations",sdyn->maxiter);
  genalphaparams.set<int>   ("num iterations",-1);
  genalphaparams.set<double>("tolerance displacements",sdyn->toldisp);

  genalphaparams.set<bool>  ("io structural disp",ioflags.struct_disp);
  genalphaparams.set<int>   ("io disp every nstep",sdyn->updevry_disp);
  genalphaparams.set<bool>  ("io structural stress",ioflags.struct_stress);
  genalphaparams.set<int>   ("io stress every nstep",sdyn->updevry_stress);

  genalphaparams.set<int>   ("restart",genprob.restart);
  genalphaparams.set<int>   ("write restart every",sdyn->res_write_evry);

  genalphaparams.set<bool>  ("print to screen",true);
  genalphaparams.set<bool>  ("print to err",true);
  genalphaparams.set<FILE*> ("err file",allfiles.out_err);

  // takes values "full newton" , "modified newton" , "nonlinear cg"
  genalphaparams.set<string>("equilibrium iteration","full newton");

  // takes values "constant" consistent"  
  genalphaparams.set<string>("predictor","constant");

  StruGenAlpha timeintegrator(genalphaparams,*actdis,solver,output);

  if (genprob.restart)
  {
    timeintegrator.ReadRestart(genprob.restart);
  }
  
  // write mesh only if this is not a restart
  output.WriteMesh(genalphaparams.get<int>("step",0),
                   genalphaparams.get<double>("total time",0.0));

  // integrate
  timeintegrator.Integrate();

  return;
} // end of dyn_nlnstructural_drt()





#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
