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
  if (!actdis->HaveDofs()) actdis->AssignDegreesOfFreedom();

  // -------------------------------------------------------------------
  // get a communicator and myrank and numproc
  // -------------------------------------------------------------------
  const Epetra_Comm& Comm = actdis->Comm();
  const int myrank = Comm.MyPID();
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
  const Epetra_Map* dofcolmap = actdis->DofColMap();

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
  RefCountPtr<Epetra_Vector> rhs0 = LINALG::CreateVector(*dofrowmap,true);

  // load vectors at time t and t-dt
  RefCountPtr<Epetra_Vector> rhs1 = LINALG::CreateVector(*dofrowmap,true);
  RefCountPtr<Epetra_Vector> rhs2 = LINALG::CreateVector(*dofrowmap,true);

  // interpolated load vector
  RefCountPtr<Epetra_Vector> rhs3 = LINALG::CreateVector(*dofrowmap,true);

  // forces from dirichlet conditions
  RefCountPtr<Epetra_Vector> dirich = LINALG::CreateVector(*dofrowmap,true);

  // solution at time t and t-dt
  RefCountPtr<Epetra_Vector> sol0 = LINALG::CreateVector(*dofrowmap,true);
  RefCountPtr<Epetra_Vector> sol1 = LINALG::CreateVector(*dofrowmap,true);
  
  // incremental displacements
  RefCountPtr<Epetra_Vector> dispi = LINALG::CreateVector(*dofrowmap,true);
  
  // velocities
  RefCountPtr<Epetra_Vector> vel = LINALG::CreateVector(*dofrowmap,true);
  
  // accelerations
  RefCountPtr<Epetra_Vector> acc = LINALG::CreateVector(*dofrowmap,true);
  
  // internal forces
  RefCountPtr<Epetra_Vector> intforce = LINALG::CreateVector(*dofrowmap,true);
  
  // internal forces at t, t-dt and interpolated  
  RefCountPtr<Epetra_Vector> fie0 = LINALG::CreateVector(*dofrowmap,true);
  RefCountPtr<Epetra_Vector> fie1 = LINALG::CreateVector(*dofrowmap,true);
  RefCountPtr<Epetra_Vector> fie2 = LINALG::CreateVector(*dofrowmap,true);
  
  // working vectors
  RefCountPtr<Epetra_Vector> work0 = LINALG::CreateVector(*dofrowmap,true);

  // -------------------------------------------------------------------
  // call elements to calculate stiffness and mass
  // -------------------------------------------------------------------
  {
    // create the parameters for the discretization
    ParameterList params;
    params.set("action","calc_struct_nlnstiffmass");
    params.set("total time",acttime);
    params.set("delta time",sdyn->dt);
    // set values needed by elements (have to be column map!)
    RefCountPtr<Epetra_Vector> cfie0 = LINALG::CreateVector(*dofcolmap,false);
    LINALG::Export(*fie0,*cfie0);
    RefCountPtr<Epetra_Vector> csol0 = LINALG::CreateVector(*dofcolmap,false);
    LINALG::Export(*sol0,*csol0);
    actdis->ClearState();
    actdis->SetState("residual",cfie0);
    actdis->SetState("displacement",csol0);
    actdis->Evaluate(params,stiff_mat,mass_mat,null,null,null);
  }



















  return;
} // end of dyn_nlnstructural_drt()





#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
