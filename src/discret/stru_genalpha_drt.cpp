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

/*----------------------------------------------------------------------*
  | structural nonlinear dynamics (gen-alpha)              m.gee 12/06  |
 *----------------------------------------------------------------------*/
void stru_genalpha_drt()
{
  DSTraceHelper dst("stru_genalpha_drt");

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
  
  /*--------------------------------------- init all applied time curves -*/
  for (int actcurve=0; actcurve<numcurve; actcurve++)
    dyn_init_curve(actcurve,sdyn->nstep,sdyn->dt,sdyn->maxtime);

  // -------------------------------------------------------------------
  // create empty vectors
  // -------------------------------------------------------------------
  // a zero vector of full length
  RefCountPtr<Epetra_Vector> zeros = LINALG::CreateVector(*dofrowmap,true);
  // vector of full length; for each component 
  //                /  1   i-th DOF is supported, ie Dirichlet BC
  //    vector_i =  <
  //                \  0   i-th DOF is free
  RefCountPtr<Epetra_Vector> dirichtoggle = LINALG::CreateVector(*dofrowmap,true);
  // opposite of dirichtoggle vector, ie for each component
  //                /  0   i-th DOF is supported, ie Dirichlet BC
  //    vector_i =  <
  //                \  1   i-th DOF is free
  RefCountPtr<Epetra_Vector> invtoggle = LINALG::CreateVector(*dofrowmap,false);
  // displacements D_{n} at last time
  RefCountPtr<Epetra_Vector> dis = LINALG::CreateVector(*dofrowmap,true);
  // velocities V_{n} at last time
  RefCountPtr<Epetra_Vector> vel = LINALG::CreateVector(*dofrowmap,true);
  // accelerations A_{n} at last time
  RefCountPtr<Epetra_Vector> acc = LINALG::CreateVector(*dofrowmap,true);
  
  // displacements D_{n+1} at new time
  RefCountPtr<Epetra_Vector> disn = LINALG::CreateVector(*dofrowmap,false);
  // velocities V_{n+1} at new time 
  RefCountPtr<Epetra_Vector> veln = LINALG::CreateVector(*dofrowmap,false);
  // accelerations A_{n+1} at new time
  RefCountPtr<Epetra_Vector> accn = LINALG::CreateVector(*dofrowmap,false);

  // mid-displacements D_{n+1-alpha_f}
  RefCountPtr<Epetra_Vector> dism = LINALG::CreateVector(*dofrowmap,false);
  // mid-velocities V_{n+1-alpha_f}
  RefCountPtr<Epetra_Vector> velm = LINALG::CreateVector(*dofrowmap,false);
  // mid-accelerations A_{n+1-alpha_m}
  RefCountPtr<Epetra_Vector> accm = LINALG::CreateVector(*dofrowmap,false);

  // iterative displacement increments IncD_{n+1}
  // also known as residual displacements
  RefCountPtr<Epetra_Vector> disi = LINALG::CreateVector(*dofrowmap,true);

  // internal force vector F_int at different times
  RefCountPtr<Epetra_Vector> fint = LINALG::CreateVector(*dofrowmap,true);
  // external force vector F_ext at last times
  RefCountPtr<Epetra_Vector> fext = LINALG::CreateVector(*dofrowmap,true);
  // external mid-force vector F_{ext;n+1-alpha_f}
  RefCountPtr<Epetra_Vector> fextm = LINALG::CreateVector(*dofrowmap,false);
  // external force vector F_{n+1} at new time
  RefCountPtr<Epetra_Vector> fextn = LINALG::CreateVector(*dofrowmap,true);

  // dynamic force residual at mid-time R_{n+1-alpha}
  // also known at out-of-balance-force
  RefCountPtr<Epetra_Vector> fresm = LINALG::CreateVector(*dofrowmap,false);

  //-------------------------------------------- calculate external forces
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
    actdis->SetState("displacement",dis);
    // predicted dirichlet values
    // dis then also holds prescribed new dirichlet displacements
    actdis->EvaluateDirichlet(params,*dis,*dirichtoggle);
    actdis->ClearState();
    actdis->SetState("displacement",dis);
    // predicted rhs
    actdis->EvaluateNeumann(params,*fext);
    actdis->ClearState();
  }

  //----------------------- compute an inverse of the dirichtoggle vector
  invtoggle->PutScalar(1.0);
  invtoggle->Update(-1.0,*dirichtoggle,1.0);
  
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
    params.set("assemble vector 1",true);
    params.set("assemble vector 2",false);
    params.set("assemble vector 3",false);
    // other parameters that might be needed by the elements
    params.set("total time",acttime);
    params.set("delta time",sdyn->dt);
    // set vector values needed by elements
    actdis->ClearState();
    actdis->SetState("residual displacement",zeros);
    actdis->SetState("displacement",dis);
    //actdis->SetState("velocity",vel); // not used at the moment
    actdis->Evaluate(params,stiff_mat,mass_mat,fint,null,null);
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
  
  //--------------------------- calculate consistent initial accelerations
  {
    RefCountPtr<Epetra_Vector> rhs = LINALG::CreateVector(*dofrowmap,true);
    if (damping) damp_mat->Multiply(false,*vel,*rhs);
    rhs->Update(-1.0,*fint,1.0,*fext,-1.0);
    Epetra_Vector rhscopy(*rhs);
    rhs->Multiply(1.0,*invtoggle,rhscopy,0.0);
    solver.Solve(mass_mat,acc,rhs,true,true);
  }
  
  //----------------------------------- time integration scheme parameters
  const double beta   = sdyn->beta;
  const double gamma  = sdyn->gamma;
  const double alpham = sdyn->alpha_m;
  const double alphaf = sdyn->alpha_f;
  const double dt     = sdyn->dt;

  //==================================================== start of timeloop
  int istep = 0;
  while (acttime<=sdyn->maxtime && istep<=sdyn->nstep)
  { 
    //--------------------------------------------------- new time t_{n+1}
    acttime += dt;
    //--------------------------------------------------- predicting state
    // constant predictor displacement in domain
    disn->Update(1.0,*dis,0.0);      
    // apply new displacements at DBCs 
    // and get new external force vector
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
      actdis->SetState("displacement",disn);
      // predicted dirichlet values
      // disn then also holds prescribed new dirichlet displacements
      actdis->EvaluateDirichlet(params,*disn,*dirichtoggle);
      actdis->ClearState();
      actdis->SetState("displacement",disn);
      actdis->EvaluateNeumann(params,*fextn);
      actdis->ClearState();
    }
  
    // predicting velocity V_{n+1} (veln)
    // V_{n+1} := gamma/(beta*dt) * (D_{n+1} - D_n)
    //          + (beta-gamma)/beta * V_n
    //          + (2.*beta-gamma)/(2.*beta) * A_n
    veln->Update(1.0,*disn,-1.0,*dis,0.0);
    veln->Update((beta-gamma)/beta,*vel,
                (2.*beta-gamma)/(2.*beta),*acc,
                 gamma/(beta*dt));
    // predicting accelerations A_{n+1} (accn)
    // A_{n+1} := 1./(beta*dt*dt) * (D_{n+1} - D_n)
    //          - 1./(beta*dt) * V_n
    //          + (2.*beta-1.)/(2.*beta) * A_n
    accn->Update(1.0,*disn,-1.0,*dis,0.0);
    accn->Update(-1./(beta*dt),*vel,
                (2.*beta-1.)/(2.*beta),*acc,
                1./(beta*dt*dt));
  
    //------------------------------ compute interpolated dis, vel and acc
    // mid-displacements D_{n+1-alpha_f} (dism)
    //    D_{n+1-alpha_f} := (1.-alphaf) * D_{n+1} + alpha_f * D_{n}
    dism->Update(1.-alphaf,*disn,alphaf,*dis,0.0);
    // mid-velocities V_{n+1-alpha_f} (velm)
    //    V_{n+1-alpha_f} := (1.-alphaf) * V_{n+1} + alpha_f * V_{n}
    velm->Update(1.-alphaf,*veln,alphaf,*vel,0.0);
    // mid-accelerations A_{n+1-alpha_m} (accm)
    //    A_{n+1-alpha_m} := (1.-alpha_m) * A_{n+1} + alpha_m * A_{n}
    accm->Update(1.-alpham,*accn,alpham,*acc,0.0);
  
    //------------------------- ----- compute interpolated external forces
    // external mid-forces F_{ext;n+1-alpha_f} (fextm)
    //    F_{ext;n+1-alpha_f} := (1.-alphaf) * F_{ext;n+1} + alpha_f * F_{ext;n}
    fextm->Update(1.-alphaf,*fextn,alphaf,*fext,0.0);
  
    //------------- eval fint at interpolated state, eval stiffness matrix
    {
      // zero out stiffness
      stiff_mat = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow);
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
      actdis->SetState("residual displacement",disi);
      actdis->SetState("displacement",dism);
      //actdis->SetState("velocity",velm); // not used at the moment
      fint->PutScalar(0.0);  // initialise internal force vector
      actdis->Evaluate(params,stiff_mat,null,fint,null,null);
      actdis->ClearState();
      // finalize the stiffness matrix
      LINALG::Complete(*stiff_mat);
    }
    
    //-------------------------------------------- compute residual forces
    // Res = M . A_{n+1-alpha_m} 
    //     + C . V_{n+1-alpha_f} 
    //     + F_int(D_{n+1-alpha_f})
    //     - F_{ext;n+1-alpha_f}
    // add mid-inertial force
    mass_mat->Multiply(false,*accm,*fresm);
    // add mid-viscous damping force 
    if (damping)
    {
    	RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,false);
    	damp_mat->Multiply(false,*velm,*fviscm);
    	fresm->Update(1.0,*fviscm,1.0);
    }
    // add static mid-balance
    fresm->Update(1.0,*fint,-1.0,*fextm,1.0);
    // blank residual at DOFs on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm);
      fresm->Multiply(1.0,*invtoggle,fresmcopy,0.0);
    }

    //------------------------------------------------ build residual norm
    double norm;
    fresm->Norm2(&norm);
    cout << "Norm bevor Newton " << norm << endl;
  
    //=================================================== equilibrium loop
    int numiter=0;
    while (norm>sdyn->toldisp && numiter<=sdyn->maxiter)
    {
      //------------------------------------------- effective rhs is fresm
      //---------------------------------------------- build effective lhs
      RefCountPtr<Epetra_CrsMatrix> eff_mat = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow);
      LINALG::Add(*stiff_mat,false,1.-alphaf,*eff_mat,0.0);
      if (damping)
      {
        LINALG::Add(*damp_mat,false,(1.-alphaf)*gamma/(beta*dt),*eff_mat,1.0);
      }
      LINALG::Add(*mass_mat,false,(1.-alpham)/(beta*dt*dt),*eff_mat,1.0);
      LINALG::Complete(*eff_mat);
    
      //----------------------- apply dirichlet BCs to system of equations
      fresm->Scale(-1.0);  // delete this by building fresm with other sign
      disi->PutScalar(0.0);  // Useful? depends on solver and more
      LINALG::ApplyDirichlettoSystem(eff_mat,disi,fresm,zeros,dirichtoggle);
    
      //--------------------------------------------------- solve for disi
      // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
      if (numiter==0)
      {
        solver.Solve(eff_mat,disi,fresm,true,true);
      }
      else
      {
        solver.Solve(eff_mat,disi,fresm,true,false);
      }

      //---------------------------------- update mid configuration values
      // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
      dism->Update(1.-alphaf,*disi,1.0);
      // V_{n+1-alpha_f} := V_{n+1-alpha_f} 
      //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
      velm->Update((1.-alphaf)*gamma/(beta*dt),*disi,1.0);
      // A_{n+1-alpha_m} := A_{n+1-alpha_m} 
      //                  + (1-alpha_m)/beta/dt^2*IncD_{n+1}
      accm->Update((1.-alpham)/(beta*dt*dt),*disi,1.0);
    
      //---------------------------- compute internal forces and stiffness
      {
        // zero out stiffness
        stiff_mat = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow);
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
        actdis->SetState("residual displacement",disi);
        actdis->SetState("displacement",dism);
        //actdis->SetState("velocity",velm); // not used at the moment
        fint->PutScalar(0.0);  // initialise internal force vector
        actdis->Evaluate(params,stiff_mat,null,fint,null,null);
        actdis->ClearState();
        // finalize the stiffness matrix
        LINALG::Complete(*stiff_mat);
      }
    
      //------------------------------------------ compute residual forces
      // Res = M . A_{n+1-alpha_m} 
      //     + C . V_{n+1-alpha_f} 
      //     + F_int(D_{n+1-alpha_f})
      //     - F_{ext;n+1-alpha_f}
      // add interia mid-forces
      mass_mat->Multiply(false,*accm,*fresm);
      // add viscous mid-forces
      if (damping)
      {
        RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,false);
        damp_mat->Multiply(false,*velm,*fviscm);
        fresm->Update(1.0,*fviscm,1.0);
      }
      // add static mid-balance
      fresm->Update(1.0,*fint,-1.0,*fextm,1.0);
      // blank residual DOFs with are on Dirichlet BC
      {
        Epetra_Vector fresmcopy(*fresm);
        fresm->Multiply(1.0,*invtoggle,fresmcopy,0.0);
      }
    
      //---------------------------------------------- build residual norm
      double disinorm;
      disi->Norm2(&disinorm);
    
      fresm->Norm2(&norm);
      cout << "numiter " << numiter 
           << " norm " << norm
           << " disinorm " << disinorm << endl;;
    
    
      //--------------------------------- increment equilibrium loop index
      ++numiter;
    } // end equilibrium loop

    //---------------------------- determine new end-quantities and update
    // new displacements at t_{n+1} -> t_n
    //    D_{n} := D_{n+1} = 1./(1.-alphaf) * D_{n+1-alpha_f} 
    //                     - alphaf/(1.-alphaf) * D_n
    dis->Update(1./(1.-alphaf),*dism,-alphaf/(1.-alphaf));
    // new velocities at t_{n+1} -> t_n
    //    V_{n} := V_{n+1} = 1./(1.-alphaf) * V_{n+1-alpha_f} 
    //                     - alphaf/(1.-alphaf) * V_n
    vel->Update(1./(1.-alphaf),*velm,-alphaf/(1.-alphaf));
    // new accelerations at t_{n+1} -> t_n
    //    A_{n} := A_{n+1} = 1./(1.-alpham) * A_{n+1-alpha_m} 
    //                     - alpham/(1.-alpham) * A_n
    acc->Update(1./(1.-alpham),*accm,-alpham/(1.-alpham));
    // update new external force
    //    F_{ext;n} := F_{ext;n+1}
    fext->Update(1.0,*fextn,0.0);

    //------------------------------------------------ increment time step
    ++istep;
  }  // end time loop

  //----------------------------- this is the end my lonely friend the end
  return;
} // end of dyn_nlnstructural_drt()





#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
