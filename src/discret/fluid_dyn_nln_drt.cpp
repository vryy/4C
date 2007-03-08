/*!----------------------------------------------------------------------
\file
\brief Control routine for fluid time integration. Includes

     o Singele step one-step-theta time integration

     o Two step BDF2 Gear's methode with one-step-theta start step


     
<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
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

#include "fluid_dyn_nln_drt.H"


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
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | INT                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern INT            numcurve;
extern struct _CURVE *curve;



#include "fluid_dyn_nln_drt.H"


/*----------------------------------------------------------------------*
 * Time integration loop for fluid.
 * 
 *        o One-step-theta
 *        o BDF2
 *
 *----------------------------------------------------------------------*/

void dyn_fluid_drt()
{


 DSTraceHelper dst("dyn_fluid_drt");


  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  {
    vector<RefCountPtr<DRT::Discretization> >* fool =
              (vector<RefCountPtr<DRT::Discretization> >*)field[0].ccadis;
    actdis = (*fool)[0];
  }
  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();

  
  // -------------------------------------------------------------------
  // get a communicator and myrank
  // -------------------------------------------------------------------
  const Epetra_Comm& Comm = actdis->Comm();
  const int myrank  = Comm.MyPID();

  
  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR        *actsolv  = &solv[0];
  FLUID_DYNAMIC *fdyn     = alldyn[0].fdyn;

  
  //-----------------------------------------------------create a solver
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  LINALG::Solver solver(solveparams,actdis->Comm(),allfiles.out_err);
  solver.TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);


  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = actdis->DofRowMap();


  // -------------------------------------------------------------------
  // create empty system matrix --- stiffness and mass are assembled in
  // one system matrix!
  // -------------------------------------------------------------------
  RefCountPtr<Epetra_CrsMatrix> sys_mat
      = LINALG::CreateMatrix(*dofrowmap,81);


  // -------------------------------------------------------------------
  // create empty right hand side vector for linear system
  // -------------------------------------------------------------------
  RefCountPtr<Epetra_Vector> rhs  = LINALG::CreateVector(*dofrowmap,true);

  
  // -------------------------------------------------------------------
  // create empty vectors 
  // -------------------------------------------------------------------

  // accelerations at time n and n-1
  RefCountPtr<Epetra_Vector> accn  = LINALG::CreateVector(*dofrowmap,true);
  RefCountPtr<Epetra_Vector> accnm = LINALG::CreateVector(*dofrowmap,true);

  // velocities and pressures at time n+1, n and n-1
  RefCountPtr<Epetra_Vector> velnp = LINALG::CreateVector(*dofrowmap,true);
  RefCountPtr<Epetra_Vector> veln  = LINALG::CreateVector(*dofrowmap,true);
  RefCountPtr<Epetra_Vector> velnm = LINALG::CreateVector(*dofrowmap,true);
  
  // toggle vector indicating which dofs have Dirichlet BCs
  RefCountPtr<Epetra_Vector> dirichtoggle = LINALG::CreateVector(*dofrowmap,true);

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RefCountPtr<Epetra_Vector> zeros = LINALG::CreateVector(*dofrowmap,true);


  /*------------------------------------------- set initial step and time */
  fdyn->step    =   0;
  fdyn->acttime = 0.0;
  
  /*--------------------------------------- init all applied time curves -*/
  for (int actcurve=0; actcurve<numcurve; actcurve++)
  {
   /* the last three parameters are obsolete!!! */  
   dyn_init_curve(actcurve,fdyn->step,fdyn->dt,fdyn->maxtime);
  }

  /*------------------------------------------------------- printout head */
  if (myrank==0)
  {

  } /* end if (myrank==0) */
  
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  /*                               TIMELOOP                             */
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  {
   bool stop_timeloop=false;
   while (stop_timeloop==false)
   {
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*                NONLINEAR ITERATION (FLUID)                       */
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    bool stop_nonliniter=false;
    while (stop_nonliniter==false)
    {


	// check steady state, maxiter and maxtime
	stop_nonliniter=true;
    }
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*               END NONLINEAR ITERATION (FLUID)                    */
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
    /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/

    // check steady state, maxiter and maxtime
    stop_timeloop=true;
   }
  
  }
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  /*                            END TIMELOOP                            */
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  /*<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>*/
  
  
  //---------- this is the end. Beautiful friend. My only friend, The end.
  // thanks to RefCountPtr<> we do not need to delete anything here!
  
  return;

} // end of dyn_fluid_drt()

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
