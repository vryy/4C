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
  {
    vector<RefCountPtr<DRT::Discretization> >* fool =
              (vector<RefCountPtr<DRT::Discretization> >*)field[0].ccadis;
    actdis = (*fool)[0];
  }
  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // get a communicator and myrank and numproc
  // -------------------------------------------------------------------
  const Epetra_Comm& Comm = actdis->Comm();
  const int myrank  = Comm.MyPID();
  const int numproc = Comm.NumProc();
  
  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  FIELD*          actfield = &field[0];
  SOLVAR*         actsolv  = &solv[0];
  STRUCT_DYNAMIC* sdyn     = alldyn[0].sdyn;
  STRUCT_DYN_CALC dynvar;
  memset(&dynvar, 0, sizeof(STRUCT_DYN_CALC));
  double          acttime = 0.0;
  
  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = actdis->DofRowMap();

  // -------------------------------------------------------------------
  // create empty matrices
  // -------------------------------------------------------------------
  RefCountPtr<Epetra_CrsMatrix> stiff_mat = LINALG::CreateMatrix(*dofrowmap,81);
  RefCountPtr<Epetra_CrsMatrix> mass_mat  = LINALG::CreateMatrix(*dofrowmap,81);
  RefCountPtr<Epetra_CrsMatrix> damp_mat  = null;
  bool damping = false;
  if (sdyn->damp==1) 
  {
    damping = true;
    damp_mat = LINALG::CreateMatrix(*dofrowmap,81);
  }
  
  // -------------------------------------------------------------------
  // create empty vectors
  // -------------------------------------------------------------------
  // original load vector without curve
  //RefCountPtr<Epetra_Vector> rhs0 = LINALG::CreateVector(*dofrowmap,true);

  // load vectors at time t and t-dt
  RefCountPtr<Epetra_Vector> rhs1 = LINALG::CreateVector(*dofrowmap,true);
  //RefCountPtr<Epetra_Vector> rhs2 = LINALG::CreateVector(*dofrowmap,true);

  // interpolated load vector
  //RefCountPtr<Epetra_Vector> rhs3 = LINALG::CreateVector(*dofrowmap,true);

  // forces from dirichlet conditions
  //RefCountPtr<Epetra_Vector> dirich = LINALG::CreateVector(*dofrowmap,true);

  // solution at time t and t-dt
  RefCountPtr<Epetra_Vector> sol0 = LINALG::CreateVector(*dofrowmap,true);
  //RefCountPtr<Epetra_Vector> sol1 = LINALG::CreateVector(*dofrowmap,true);
  
  // incremental displacements
  RefCountPtr<Epetra_Vector> dispi = LINALG::CreateVector(*dofrowmap,true);
  
  // velocities
  //RefCountPtr<Epetra_Vector> vel = LINALG::CreateVector(*dofrowmap,true);
  
  // accelerations
  //RefCountPtr<Epetra_Vector> acc = LINALG::CreateVector(*dofrowmap,true);
  
  // internal forces
  //RefCountPtr<Epetra_Vector> intforce = LINALG::CreateVector(*dofrowmap,true);
  
  // internal forces at t, t-dt and interpolated  
  //RefCountPtr<Epetra_Vector> fie0 = LINALG::CreateVector(*dofrowmap,true);
  //RefCountPtr<Epetra_Vector> fie1 = LINALG::CreateVector(*dofrowmap,true);
  //RefCountPtr<Epetra_Vector> fie2 = LINALG::CreateVector(*dofrowmap,true);
  
  // working vectors
  //RefCountPtr<Epetra_Vector> work0 = LINALG::CreateVector(*dofrowmap,true);

  // -------------------------------------------------------------------
  // call elements to calculate stiffness and mass
  // -------------------------------------------------------------------
  {
    // create the parameters for the discretization
    ParameterList params;
    // action for elements
    params.set("action","calc_struct_nlnstiffmass");
    // choose what to assemble
    params.set("assemble matrix 1",true);
    params.set("assemble matrix 2",true);
    params.set("assemble vector 1",false);
    params.set("assemble vector 2",false);
    params.set("assemble vector 3",false);
    // other parameters that might be needed by the elements
    params.set("total time",acttime);
    params.set("delta time",sdyn->dt);
    // set vector values needed by elements
    actdis->ClearState();
    actdis->SetState("residual displacement",dispi);
    actdis->SetState("displacement",sol0);
    actdis->Evaluate(params,stiff_mat,mass_mat,null,null,null);
    actdis->ClearState();
  }

  // get dirichlet conditions from discretization and apply them to system of equations

  // complete stiffness and mass matrix
  LINALG::Complete(*stiff_mat);
  LINALG::Complete(*mass_mat);
  
  // build damping matrix if neccessary
  if (damping)
  {
    LINALG::Add(*stiff_mat,false,sdyn->k_damp,*damp_mat,0.0);
    LINALG::Add(*mass_mat ,false,sdyn->m_damp,*damp_mat,1.0);
    LINALG::Complete(*damp_mat);
  }

  /*------------------------------------------- set initial step and time */
  sdyn->step = -1;
  sdyn->time = 0.0;
  
  /*--------------------------------------- init all applied time curves -*/
  for (int actcurve=0; actcurve<numcurve; actcurve++)
    dyn_init_curve(actcurve,sdyn->nstep,sdyn->dt,sdyn->maxtime);
  
  /*------------------------------------------------------- printout head */
  if (myrank==0) dyn_nlnstruct_outhead(&dynvar,sdyn);
  
  /*----------------------------------------------------------------------*/
  /*                     START LOOP OVER ALL STEPS                        */
  /*----------------------------------------------------------------------*/
  timeloop:
  double t0 = ds_cputime();
  /*--------------------------------------------- increment step and time */
  sdyn->step++;
  sdyn->time += sdyn->dt;
  acttime = sdyn->time;
  /*-------------------------------------------------- set some constants */
  double constants[16];
  {
    const double beta   = sdyn->beta;
    const double gamma  = sdyn->gamma;
    const double alpham = sdyn->alpha_m;
    const double alphaf = sdyn->alpha_f;
    const double dt     = sdyn->dt;
    constants[6]  = 1.0-alphaf;
    constants[7]  = 1.0-alpham;
    constants[8]  = (1.0/beta)/(DSQR(dt));
    constants[9]  = -constants[8] * dt;
    constants[10] = 1.0 - 0.5/beta;
    constants[11] = (gamma/beta)/dt;
    constants[12] = 1.0 - gamma/beta;
    constants[13] = (1.0-(gamma/beta)/2.0)*dt;
    constants[0]  = constants[7]*constants[8];
    constants[1]  = constants[0]*dt;
    constants[2]  = (constants[7]/2.0)/beta - 1.0;
    constants[3]  = constants[6]*constants[11];
    constants[4]  = constants[3]*dt - 1.0;
    constants[5]  = ((gamma/beta)/2.0 - 1.0)*dt*constants[6];
    constants[14] = dt;
    constants[15] = beta;
    constants[16] = gamma;
  }
  
  
  /*------------------------- set incremental displacements dispi to zero */
  dispi->PutScalar(0.0);
  
  /*----------------------------------------------------------------------*/
  /*                     PREDICTOR                                        */
  /*----------------------------------------------------------------------*/
  
  //-------------------------------------------------- evaluate Neumann BCs
  {
    ParameterList params;
    // action for elements
    params.set("action","calc_struct_eleload");
    // choose what to assemble
    params.set("assemble matrix 1",false);
    params.set("assemble matrix 2",false);
    params.set("assemble vector 1",true);
    params.set("assemble vector 2",false);
    params.set("assemble vector 3",false);
    // other parameters needed by the elements
    params.set("total time",acttime);
    params.set("delta time",sdyn->dt);
    // set vector values needed by elements
    actdis->ClearState();
    actdis->SetState("displacement",sol0);
    actdis->EvaluateNeumann(params,*rhs1);
    actdis->ClearState();
  }




















  return;
} // end of dyn_nlnstructural_drt()





#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
