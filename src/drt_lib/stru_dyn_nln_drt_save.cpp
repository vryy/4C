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
  RefCountPtr<DRT::Discretization> actdis = 
    DRT::Problem::Instance()->Dis(genprob.numsf,0);

  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  DiscretizationWriter output(actdis);
  output.WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // get a communicator and myrank
  // -------------------------------------------------------------------
  const Epetra_Comm& Comm = actdis->Comm();
  const int myrank  = Comm.MyPID();

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR*         actsolv  = &solv[0];
  STRUCT_DYNAMIC* sdyn     = alldyn[0].sdyn;
  STRUCT_DYN_CALC dynvar;
  memset(&dynvar, 0, sizeof(STRUCT_DYN_CALC));
  double          acttime = 0.0;

  //-----------------------------------------------------create a solver
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  LINALG::Solver solver(solveparams,actdis->Comm(),allfiles.out_err);
  solver.TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

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

  // load vectors at time t-dt and t
  RefCountPtr<Epetra_Vector> rhs1 = LINALG::CreateVector(*dofrowmap,true);
  RefCountPtr<Epetra_Vector> rhs2 = LINALG::CreateVector(*dofrowmap,true);
  // interpolated load vector
  RefCountPtr<Epetra_Vector> rhs  = LINALG::CreateVector(*dofrowmap,true);

  // solution at time t-dt and t
  RefCountPtr<Epetra_Vector> sol0 = LINALG::CreateVector(*dofrowmap,true);
  RefCountPtr<Epetra_Vector> sol1 = LINALG::CreateVector(*dofrowmap,true);

  // incremental displacements
  RefCountPtr<Epetra_Vector> dx = LINALG::CreateVector(*dofrowmap,true);

  // residual incremental displacements
  RefCountPtr<Epetra_Vector> rdx = LINALG::CreateVector(*dofrowmap,true);

  // toggle vector indicating which dofs have Dirichlet BCs
  RefCountPtr<Epetra_Vector> dirichtoggle = LINALG::CreateVector(*dofrowmap,true);

  // velocities
  RefCountPtr<Epetra_Vector> vel = LINALG::CreateVector(*dofrowmap,true);

  // accelerations
  RefCountPtr<Epetra_Vector> acc = LINALG::CreateVector(*dofrowmap,true);

  // internal forces at t-dt, t
  RefCountPtr<Epetra_Vector> fie1 = LINALG::CreateVector(*dofrowmap,true);
  RefCountPtr<Epetra_Vector> fie2 = LINALG::CreateVector(*dofrowmap,true);
  // interpolated internal forces
  RefCountPtr<Epetra_Vector> fie  = LINALG::CreateVector(*dofrowmap,true);

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RefCountPtr<Epetra_Vector> zeros = LINALG::CreateVector(*dofrowmap,true);

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
    actdis->SetState("residual displacement",rdx);
    actdis->SetState("displacement",sol0);
    actdis->Evaluate(params,stiff_mat,mass_mat,null,null,null);
    actdis->ClearState();
  }

  // complete stiffness and mass matrix
  LINALG::Complete(*stiff_mat);
  LINALG::Complete(*mass_mat);
  const int maxentriesperrow = stiff_mat->MaxNumEntries();

  // build damping matrix if neccessary
  if (damping)
  {
    LINALG::Add(*stiff_mat,false,sdyn->k_damp,*damp_mat,0.0);
    LINALG::Add(*mass_mat ,false,sdyn->m_damp,*damp_mat,1.0);
    LINALG::Complete(*damp_mat);
  }
  stiff_mat = null;

  /*------------------------------------------- set initial step and time */
  sdyn->step = -1;
  sdyn->time = 0.0;


  /*------------------------------------------------------- printout head */
  if (myrank==0) dyn_nlnstruct_outhead(&dynvar,sdyn);

  //------------------------------------------------- output initial state
  output.NewStep(0,0.0);
  output.WriteVector("displacement", sol0);
  output.WriteVector("velocity", vel);
  output.WriteVector("acceleration", acc);

  /*----------------------------------------------------------------------*/
  /*                     START LOOP OVER ALL STEPS                        */
  /*----------------------------------------------------------------------*/
  timeloop:
  Epetra_Time timelooptime(Comm);
  const double t0 = timelooptime.WallTime();
  /*--------------------------------------------- increment step and time */
  sdyn->step++;
  sdyn->time += sdyn->dt;
  acttime = sdyn->time;
  /*-------------------------------------------------- set some constants */
  const double beta   = sdyn->beta;
  const double gamma  = sdyn->gamma;
  const double alpham = sdyn->alpha_m;
  const double alphaf = sdyn->alpha_f;
  const double dt     = sdyn->dt;

  /*------------------------- set incremental displacements dispi to zero */
  dx->PutScalar(0.0);
  //------------------------------------- set residual displacements to zero
  rdx->PutScalar(0.0);

  /*----------------------------------------------------------------------*/
  /*                     PREDICTOR                                        */
  /*----------------------------------------------------------------------*/

  //------------------------------------- evaluate Neumann and Dirichlet BCs
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
    // predicted dirichlet values
    // sol0 then also holds prescribed new dirichlet displacements
    actdis->EvaluateDirichlet(params,*sol1,*dirichtoggle);
    actdis->ClearState();
    actdis->SetState("displacement",sol0);
    // predicted rhs
    rhs2->PutScalar(0.0);
    actdis->EvaluateNeumann(params,*rhs2); cout << *rhs2;
    actdis->ClearState();
  }

  //------------------------------------ evaluate stiffness and internal forces
  {
    // zero out the stiffness matrix
    stiff_mat = LINALG::CreateMatrix(*dofrowmap,81);
    // zero out internal forces
    fie1->PutScalar(0.0);
    // create the parameters for the discretization
    ParameterList params;
    // action for elements
    params.set("action","calc_struct_nlnstiff");
    // choose what to assemble
    params.set("assemble matrix 1",true);
    params.set("assemble matrix 2",false);
    params.set("assemble vector 1",true);
    params.set("assemble vector 2",false);
    params.set("assemble vector 3",false);
    // other parameters that might be needed by the elements
    params.set("total time",acttime);
    params.set("delta time",sdyn->dt);
    // set vector values needed by elements
    actdis->ClearState();
    actdis->SetState("residual displacement",rdx);
    actdis->SetState("displacement",sol0);
    fie1->PutScalar(0.0);
    actdis->Evaluate(params,stiff_mat,null,fie1,null,null);
    actdis->ClearState();
  }

  //----------------------------------------------- interpolate external forces
  // rhs = (1-alphaf)*rhs2 + alphaf*rhs1
  rhs->Update((1.-alphaf),*rhs2,alphaf,*rhs1,0.0);

  //---------------- subtract internal forces from interpolated external forces
  // rhs = rhs - fie1;
  rhs->Update(-1.0,*fie1,1.0);

  //========================================================build effective RHS
  // rhs = (1-alphaf)*rhs1 + alphaf*rhs0 - fie1 (already done above)
  //       + M*(-a1*dx + a2*vel + a3*acc)
  //       + D*(-a4*dx + a5*vel + a6*acc)       (if present)
  //
  {
    const double a1 = (1.0-alpham) * (1./beta)/(dt*dt);
    const double a2 = ((1.0-alpham) * (1./beta)/(dt*dt))*dt;
    const double a3 = (1.0-alpham) / (2. * beta) - 1.0;
    const double a4 = (1.0-alphaf) * (gamma/(beta*dt));
    const double a5 = ((1.0-alphaf) * (gamma/(beta*dt)))*dt - 1.0;
    const double a6 = ((gamma/beta)/2.0 - 1.0)*dt*(1.0-alphaf);

    // build work = -a1*dx + a2*vel + a3*acc
    RefCountPtr<Epetra_Vector> work = LINALG::CreateVector(*dofrowmap,false);
    work->Update(-a1,*dx,a2,*vel,0.0);
    work->Update(a3,*acc,1.0);
    // rhs = rhs + mass_mat*work
    RefCountPtr<Epetra_Vector> tmp = LINALG::CreateVector(*dofrowmap,false);
    mass_mat->Multiply(false,*work,*tmp);
    rhs->Update(1.0,*tmp,1.0);
    if (damping)
    {
      // build work = -a4*dx + a5*vel + a6*acc (if present)
      work->Update(-a4,*dx,a5,*vel,0.0);
      work->Update(a6,*acc,1.0);
      // rhs = rhs + damp_mat*work2
      damp_mat->Multiply(false,*work,*tmp);
      rhs->Update(1.0,*tmp,1.0);
    }
  }

  //======================================================= build effective LHS
  // keff =   (1.0-alphaf)*K
  //        + (1.0-alpham)*(1.0/(beta*dt*dt))*M
  //        + (1.0-alphaf)*(gamma/(beta*dt))*D (if present)
  {
    const double a1 = (1.0-alphaf);
    const double a2 = (1.0-alpham)*(1.0/(beta*dt*dt));
    const double a3 = (1.0-alphaf)*(gamma/(beta*dt));
    LINALG::Add(*mass_mat,false,a2,*stiff_mat,a1);
    if (damping) LINALG::Add(*damp_mat,false,a3,*stiff_mat,1.0);
    LINALG::Complete(*stiff_mat);
  }

  //================ Apply dirichlet boundary conditions to system of equations
  // increment dx at Dirichlet BCs should be zero as we already have proper
  // dirichlet values in sol1
  // we use the vector dx to set zeros to dx and rhs respectively
  LINALG::ApplyDirichlettoSystem(stiff_mat,dx,rhs,zeros,dirichtoggle);

  //============================================== solve for dx = Keff^-1 * rhs
  solver.Solve(stiff_mat,dx,rhs,true,true);

  //====================================================== update displacements
  // new dirichlet values present on sol1
  sol1->Update(1.0,*dx,1.0);
  // make dx contain increment of dirichlet BCs
  dx->Update(1.0,*sol1,-1.0,*sol0,0.0);

  //------------------ start with residual discplacements zero as initial guess
  rdx->PutScalar(0.0);

  /*----------------------------------------------------------------------*/
  /*                     PERFORM EQUILLIBRIUM ITERATION                   */
  /*----------------------------------------------------------------------*/
  int itnum = 0;
  iterloop:

  //--------------- call elements to calculate stiffness and internal forces
  {
    // zero out the stiffness matrix
    stiff_mat = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow);
    // zero out internal forces
    fie2->PutScalar(0.0);
    // create the parameters for the discretization
    ParameterList params;
    // action for elements
    params.set("action","calc_struct_nlnstiff");
    // choose what to assemble
    params.set("assemble matrix 1",true);
    params.set("assemble matrix 2",false);
    params.set("assemble vector 1",true);
    params.set("assemble vector 2",false);
    params.set("assemble vector 3",false);
    // other parameters that might be needed by the elements
    params.set("total time",acttime);
    params.set("delta time",sdyn->dt);
    // set vector values needed by elements
    actdis->ClearState();
    actdis->SetState("residual displacement",rdx);
    actdis->SetState("displacement",sol1);
    fie2->PutScalar(0.0);
    actdis->Evaluate(params,stiff_mat,null,fie2,null,null);
    actdis->ClearState();
  }

  //----------------------------------------------- interpolate external forces
  // rhs = (1-alphaf)*rhs2 + alphaf*rhs1
  rhs->Update((1.-alphaf),*rhs2,alphaf,*rhs1,0.0);

  //----------------------------------------------- interpolate internal forces
  fie->Update((1.-alphaf),*fie2,alphaf,*fie1,0.0);

  //------------------------------ subtract internal forces from external force
  rhs->Update(-1.0,*fie,1.0);

  //---------------------------------------------- create effective load vector
  // rhs = (1-alphaf)*rhs1 + alphaf*rhs0 - {(1-alphaf)*fie2 + alphaf*fie1} (already done above)
  //       + M*(-a1*dx + a2*vel + a3*acc)
  //       + D*(-a4*dx + a5*vel + a6*acc)       (if present)
  //
  {
    const double a1 = (1.0-alpham) * (1./beta)/(dt*dt);
    const double a2 = ((1.0-alpham) * (1./beta)/(dt*dt))*dt;
    const double a3 = (1.0-alpham) / (2. * beta) - 1.0;
    const double a4 = (1.0-alphaf) * (gamma/(beta*dt));
    const double a5 = ((1.0-alphaf) * (gamma/(beta*dt)))*dt - 1.0;
    const double a6 = ((gamma/beta)/2.0 - 1.0)*dt*(1.0-alphaf);

    // build work = -a1*dx + a2*vel + a3*acc
    RefCountPtr<Epetra_Vector> work = LINALG::CreateVector(*dofrowmap,false);
    work->Update(-a1,*dx,a2,*vel,0.0);
    work->Update(a3,*acc,1.0);
    // rhs = rhs + mass_mat*work
    RefCountPtr<Epetra_Vector> tmp = LINALG::CreateVector(*dofrowmap,false);
    mass_mat->Multiply(false,*work,*tmp);
    rhs->Update(1.0,*tmp,1.0);
    if (damping)
    {
      // build work = -a4*dx + a5*vel + a6*acc (if present)
      work->Update(-a4,*dx,a5,*vel,0.0);
      work->Update(a6,*acc,1.0);
      // rhs = rhs + damp_mat*work
      damp_mat->Multiply(false,*work,*tmp);
      rhs->Update(1.0,*tmp,1.0);
    }
  }

  //------------------------------------------------------ create effective LHS
  // keff =   (1.0-alphaf)*K
  //        + (1.0-alpham)*(1.0/(beta*dt*dt))*M
  //        + (1.0-alphaf)*(gamma/(beta*dt))*D (if present)
  {
    const double a1 = (1.0-alphaf);
    const double a2 = (1.0-alpham)*(1.0/(beta*dt*dt));
    const double a3 = (1.0-alphaf)*(gamma/(beta*dt));
    LINALG::Add(*mass_mat ,false,a2,*stiff_mat,a1);
    if (damping) LINALG::Add(*damp_mat ,false,a3,*stiff_mat,1.0);
    LINALG::Complete(*stiff_mat);
  }

  //---------------- Apply dirichlet boundary conditions to system of equations
  // residual discplacements are supposed to be zero at boundary conditions
  LINALG::ApplyDirichlettoSystem(stiff_mat,rdx,rhs,zeros,dirichtoggle);

  //-------solve for residual displacements to correct incremental displacements
  rdx->PutScalar(0.0);
  solver.Solve(stiff_mat,rdx,rhs,true,false);

  //--------------- update incremental displacements by residual displacements
  dx->Update(1.0,*rdx,1.0);

  //------------------------------------------------ update total displacements
  sol1->Update(1.0,*rdx,1.0);

  //----------------------------------------------------- check for convergence
  {
    bool conv = false;
    double rdxnorm2;
    double dxnorm2;
    double rdxnorminf;
    rdx->Norm2(&rdxnorm2);
    rdx->NormInf(&rdxnorminf);
    dx->Norm2(&dxnorm2);
    if (myrank==0)
    {
      printf("  rdx %15.10E\n",rdxnorm2);
      fflush(stdout);
    }
    if (rdxnorm2 < sdyn->toldisp ||
        (rdxnorm2 < 1.0e-14 && rdxnorminf < 1.0e-12) ||
        dxnorm2 < 1.0e-14)
    {
      itnum++;
      conv = true;
    }
    else
    {
      itnum++;
      if (itnum==sdyn->maxiter) dserror("No convergence in Newton in max iterations steps");
      goto iterloop;
    }
  }
  /*----------------------------------------------------------------------*/
  /*                      END OF EQUILLIBRIUM ITERATION                   */
  /*----------------------------------------------------------------------*/
  //--------------------------------------------- do update for next time step
  {
    const double a1 = (1.0/beta)/(dt*dt);
    const double a2 = -(1.0/beta)/(dt);
    const double a3 = 1.0 - 0.5/beta;
    const double a4 = (gamma/beta)/dt;
    const double a5 = 1.0 - gamma/beta;
    const double a6 = (1.0-(gamma/beta)/2.0)*dt;

    // displacement increment is dx (already calculated in Newton iteration)
    // dx contains correct dirichlet BCs increment

    // copy load vector from rhs2 to rhs1
    rhs1->Update(1.0,*rhs2,0.0);

    // copy internal forces from fie2 to fie1
    fie1->Update(1.0,*fie2,0.0);

    // copy displacements from sol1 to sol0
    sol0->Update(1.0,*sol1,0.0);

    // temporary copy of acc
    RefCountPtr<Epetra_Vector> acc_old = LINALG::CreateVector(*dofrowmap,false);
    acc_old->Update(1.0,*acc,0.0);

    // new accelerations acc = a1*dx + a2*vel + a3*acc_old
    acc->Update(a1,*dx,a2,*vel,a3);

    // new velocities vel = a4*dx + a5*vel_old + a6*acc_old
    vel->Update(a4,*dx,a6,*acc_old,a5);

/*
    printf("displacements\n");
    for (int i=0; i<sol1->MyLength(); ++i)
      printf("%d     %20.15e\n",i,(*sol1)[i]);
    fflush(stdout);
    printf("velocities\n");
    for (int i=0; i<vel->MyLength(); ++i)
      printf("%d     %20.15e\n",i,(*vel)[i]);
    fflush(stdout);
    printf("accelerations\n");
    for (int i=0; i<acc->MyLength(); ++i)
      printf("%d     %20.15e\n",i,(*acc)[i]);
    fflush(stdout);
*/
  }

  //---------------------------------------check whether to write output or not
  int mod_disp   = sdyn->step % sdyn->updevry_disp;
  int mod_stress = sdyn->step % sdyn->updevry_stress;

  //--------------------------write displacements, velocities and accelerations
  if (!mod_disp && ioflags.struct_disp==1)
  {
    output.NewStep(sdyn->step+1, sdyn->time);
    output.WriteVector("displacement", sol0);
    output.WriteVector("velocity", vel);
    output.WriteVector("acceleration", acc);
  }

  //----------------------------------------------------- do stress calculation
  if (mod_stress==0 && ioflags.struct_stress==1)
  {
    // create the parameters for the discretization
    ParameterList params;
    // action for elements
    params.set("action","calc_struct_stress");
    // choose what to assemble
    params.set("assemble matrix 1",false);
    params.set("assemble matrix 2",false);
    params.set("assemble vector 1",false);
    params.set("assemble vector 2",false);
    params.set("assemble vector 3",false);
    // other parameters that might be needed by the elements
    params.set("total time",acttime);
    params.set("delta time",sdyn->dt);
    // set vector values needed by elements
    actdis->ClearState();
    rdx->PutScalar(0.0);
    actdis->SetState("residual displacement",rdx);
    actdis->SetState("displacement",sol1);
    actdis->Evaluate(params,null,null,null,null,null);
    actdis->ClearState();
  }

  //---------------- update any material history parameters or something else...
  {
    // create the parameters for the discretization
    ParameterList params;
    // action for elements
    params.set("action","calc_struct_update_istep");
    // choose what to assemble
    params.set("assemble matrix 1",false);
    params.set("assemble matrix 2",false);
    params.set("assemble vector 1",false);
    params.set("assemble vector 2",false);
    params.set("assemble vector 3",false);
    // other parameters that might be needed by the elements
    params.set("total time",acttime);
    params.set("delta time",sdyn->dt);
    actdis->Evaluate(params,null,null,null,null,null);
  }

  //--------------------------------------------------print output of time step
  if (Comm.MyPID()==0)
  {
    printf("STEP=%6d | NSTEP=%6d | TIME=%-14.8E | DT=%-14.8E | NUMITER=%3d\n",
           sdyn->step,sdyn->nstep,sdyn->time,dt,itnum);
    fprintf(allfiles.out_err,"STEP=%6d | NSTEP=%6d | TIME=%-14.8E | DT=%-14.8E | NUMITER=%3d\n",
            sdyn->step,sdyn->nstep,sdyn->time,dt,itnum);
    fflush(stdout);
    fflush(allfiles.out_err);

  }
  //--------------------------------------------- print timing to error file
  Comm.Barrier();
  const double t1 = timelooptime.WallTime();
  fprintf(allfiles.out_err,"TIME for step %d : %15.10e sec\n",sdyn->step,t1-t0);

  //========================================== check time and number of steps
  if (sdyn->step < sdyn->nstep-1 && sdyn->time <= sdyn->maxtime)
    goto timeloop;

  //------------- this is the end, no clean up needed thanxs to RefCountPtrs!
  return;
} // end of dyn_nlnstructural_drt()





#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
