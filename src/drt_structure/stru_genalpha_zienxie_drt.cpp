/*======================================================================*/
/*!
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include <ctime>
#include <cstdlib>
#include <iostream>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "stru_genalpha_zienxie_drt.H"
#include "../drt_timada/timeadaptivity.H"
#include "../drt_timada/ta_zienkiewiczxie.H"
#include "../io/io_drt.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern FIELD* field;

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern GENPROB genprob;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern FILES allfiles;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern IO_FLAGS ioflags;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern SOLVAR* solv;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA* alldyn;

/*======================================================================*/
/*!
  \brief Structural nonlinear dynamics with generalised-alpha and
         time step size adapivity with Zienkiewicz-Xie error
         indicator (GA/ZX)
  \author bborn
  \date 10/07
*/
void stru_genalpha_zienxie_drt() 
{
  DSTraceHelper dst("stru_genalpha_zienxie_drt");

  // ---------------------------------------------------------------------
  // set some pointers and variables
  const INT disnum = 0;
  SOLVAR* actsolv = &solv[disnum];
  STRUCT_DYNAMIC* sdyn = alldyn[disnum].sdyn;
  TIMADA_DYNAMIC* timada = &(sdyn->timada);

  // ---------------------------------------------------------------------
  // access the discretization
  // ---------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numsf,0);
  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();
  // processor ID
  int myrank = (*actdis).Comm().MyPID();

  // ---------------------------------------------------------------------
  // A word to the user
  if (myrank == 0)
  {
     printf("Adaptive Structural time integration with\n");
     printf("   Generalised-alpha (marching scheme)\n");
     printf("   Zienkiewicz-Xie (auxiliar scheme)\n");
  }
  
  // ---------------------------------------------------------------------
  // allocate adaptive time integrator
  ZienkiewiczXie adatimint
  (
     (double) timada->dt_max,
     (double) timada->dt_min,
     (double) timada->dt_scl_min,
     (double) timada->dt_scl_max,
     (double) timada->dt_scl_saf,
     (double) timada->err_tol,
     2
  );

  cout << adatimint << endl;


//   // -------------------------------------------------------------------
//   // context for output and restart
//   // -------------------------------------------------------------------
//   IO::DiscretizationWriter output(actdis);



//   // -------------------------------------------------------------------
//   // create a solver
//   // -------------------------------------------------------------------
//   RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
//   LINALG::Solver solver(solveparams,actdis->Comm(),allfiles.out_err);
//   solver.TranslateSolverParameters(*solveparams,actsolv);
//   actdis->ComputeNullSpaceIfNecessary(*solveparams);

//   // -------------------------------------------------------------------
//   // create a generalized alpha time integrator
//   // -------------------------------------------------------------------
//   switch(sdyn->Typ)
//   {
//     //==================================================================
//     // Generalized alpha time integration
//     //==================================================================
//     case STRUCT_DYNAMIC::gen_alfa :
//     {
//       ParameterList genalphaparams;
//       StruGenAlpha::SetDefaults(genalphaparams);

//       genalphaparams.set<bool>  ("damping",sdyn->damp);
//       genalphaparams.set<double>("damping factor K",sdyn->k_damp);
//       genalphaparams.set<double>("damping factor M",sdyn->m_damp);

//       genalphaparams.set<double>("beta",sdyn->beta);
//       genalphaparams.set<double>("gamma",sdyn->gamma);
//       genalphaparams.set<double>("alpha m",sdyn->alpha_m);
//       genalphaparams.set<double>("alpha f",sdyn->alpha_f);

//       genalphaparams.set<double>("total time",0.0);
//       genalphaparams.set<double>("delta time",sdyn->dt);
//       genalphaparams.set<int>   ("step",0);
//       genalphaparams.set<int>   ("nstep",sdyn->nstep);
//       genalphaparams.set<int>   ("max iterations",sdyn->maxiter);
//       genalphaparams.set<int>   ("num iterations",-1);
//       genalphaparams.set<double>("tolerance displacements",sdyn->toldisp);

//       genalphaparams.set<bool>  ("io structural disp",ioflags.struct_disp);
//       genalphaparams.set<int>   ("io disp every nstep",sdyn->updevry_disp);
//       genalphaparams.set<bool>  ("io structural stress",ioflags.struct_stress);
//       genalphaparams.set<int>   ("io stress every nstep",sdyn->updevry_stress);

//       genalphaparams.set<int>   ("restart",genprob.restart);
//       genalphaparams.set<int>   ("write restart every",sdyn->res_write_evry);

//       genalphaparams.set<bool>  ("print to screen",true);
//       genalphaparams.set<bool>  ("print to err",true);
//       genalphaparams.set<FILE*> ("err file",allfiles.out_err);

//       switch(sdyn->nlnSolvTyp)
//       {
//         case STRUCT_DYNAMIC::fullnewton:
//           genalphaparams.set<string>("equilibrium iteration","full newton");
//         break;
//         case STRUCT_DYNAMIC::modnewton:
//           genalphaparams.set<string>("equilibrium iteration","modified newton");
//         break;
//         case STRUCT_DYNAMIC::matfreenewton:
//           genalphaparams.set<string>("equilibrium iteration","matrixfree newton");
//         break;
//         case STRUCT_DYNAMIC::nlncg:
//           genalphaparams.set<string>("equilibrium iteration","nonlinear cg");
//         break;
//         case STRUCT_DYNAMIC::ptc:
//           genalphaparams.set<string>("equilibrium iteration","ptc");
//         break;
//         default:
//           genalphaparams.set<string>("equilibrium iteration","full newton");
//         break;
//       }

//       // takes values "constant" "consistent"
//       genalphaparams.set<string>("predictor","constant");

//       // create the time integrator
//       StruGenAlpha timeintegrator(genalphaparams,*actdis,solver,output);

//       // do restart if demanded from input file
//       // note that this changes time and step in genalphaparams
//       if (genprob.restart)
//         timeintegrator.ReadRestart(genprob.restart);

//       // write mesh always at beginning of calc or restart
//       {
//         int    step = genalphaparams.get<int>("step",0);
//         double time = genalphaparams.get<double>("total time",0.0);
//         output.WriteMesh(step,time);
//       }

//       // integrate in time and space
//       timeintegrator.Integrate();
//     }
//     break;
//     //==================================================================
//     // Generalized Energy Momentum Method
//     //==================================================================
//     case STRUCT_DYNAMIC::Gen_EMM :
//     {
//       dserror("Not yet impl.");
//     }
//     break;
//   } // end of switch(sdyn->Typ)

  return;
} // end of dyn_nlnstructural_drt()





#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
