/*!----------------------------------------------------------------------
\file strugenalpha.cpp
\brief Generalized Alpha time integration for structural problems

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "strugenalpha.H"
#include "iostream"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 03/07|
 *----------------------------------------------------------------------*/
StruGenAlpha::StruGenAlpha(ParameterList& params,
                           DRT::Discretization& dis,
                           LINALG::Solver& solver,
                           IO::DiscretizationWriter& output) :
params_(params),
discret_(dis),
solver_(solver),
output_(output),
myrank_(discret_.Comm().MyPID()),
maxentriesperrow_(81)
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time    = params_.get<double>("total time"      ,0.0);
  double dt      = params_.get<double>("delta time"      ,0.01);
  bool damping   = params_.get<bool>  ("damping"         ,false);
  double kdamp   = params_.get<double>("damping factor K",0.0);
  double mdamp   = params_.get<double>("damping factor M",0.0);
  int step       = params_.get<int>   ("step"            ,0);
  bool outerr    = params_.get<bool>  ("print to err"    ,false);
  FILE* errfile  = params_.get<FILE*> ("err file"        ,NULL);
  if (!errfile) outerr = false;

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  // -------------------------------------------------------------------
  if (!discret_.Filled()) discret_.FillComplete();
  const Epetra_Map* dofrowmap = discret_.DofRowMap();

  // -------------------------------------------------------------------
  // create empty matrices
  // -------------------------------------------------------------------
  stiff_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,81,true,false));
  mass_  = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,81,true,true));
  if (damping) damp_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,81,true,true));

  // -------------------------------------------------------------------
  // create empty vectors
  // -------------------------------------------------------------------
  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*dofrowmap,true);
  // vector of full length; for each component
  //                /  1   i-th DOF is supported, ie Dirichlet BC
  //    vector_i =  <
  //                \  0   i-th DOF is free
  dirichtoggle_ = LINALG::CreateVector(*dofrowmap,true);
  // opposite of dirichtoggle vector, ie for each component
  //                /  0   i-th DOF is supported, ie Dirichlet BC
  //    vector_i =  <
  //                \  1   i-th DOF is free
  invtoggle_ = LINALG::CreateVector(*dofrowmap,false);

  // displacements D_{n} at last time
  dis_ = LINALG::CreateVector(*dofrowmap,true);
  // velocities V_{n} at last time
  vel_ = LINALG::CreateVector(*dofrowmap,true);
  // accelerations A_{n} at last time
  acc_ = LINALG::CreateVector(*dofrowmap,true);

  // displacements D_{n+1} at new time
  disn_ = LINALG::CreateVector(*dofrowmap,true);
  // velocities V_{n+1} at new time
  veln_ = LINALG::CreateVector(*dofrowmap,true);
  // accelerations A_{n+1} at new time
  accn_ = LINALG::CreateVector(*dofrowmap,true);

  // mid-displacements D_{n+1-alpha_f}
  dism_ = LINALG::CreateVector(*dofrowmap,true);
  // mid-velocities V_{n+1-alpha_f}
  velm_ = LINALG::CreateVector(*dofrowmap,true);
  // mid-accelerations A_{n+1-alpha_m}
  accm_ = LINALG::CreateVector(*dofrowmap,true);

  // iterative displacement increments IncD_{n+1}
  // also known as residual displacements
  disi_ = LINALG::CreateVector(*dofrowmap,true);

  // internal force vector F_int at different times
  fint_ = LINALG::CreateVector(*dofrowmap,true);
  // inertial force vector F_inert at different times
  finert_ = LINALG::CreateVector(*dofrowmap,true);
  // viscous force vector F_visc at different times
  fvisc_ = LINALG::CreateVector(*dofrowmap,true);
  // external force vector F_ext at last times
  fext_ = LINALG::CreateVector(*dofrowmap,true);
  // external mid-force vector F_{ext;n+1-alpha_f}
  fextm_ = LINALG::CreateVector(*dofrowmap,true);
  // external force vector F_{n+1} at new time
  fextn_ = LINALG::CreateVector(*dofrowmap,true);

  // dynamic force residual at mid-time R_{n+1-alpha}
  // also known as out-of-balance-force
  fresm_ = LINALG::CreateVector(*dofrowmap,false);

  //-------------------------------------------- calculate external forces
  {
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_eleload");
    // choose what to assemble
    p.set("assemble matrix 1",false);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",false);
    p.set("assemble vector 3",false);
    // other parameters needed by the elements
    p.set("total time",time);
    p.set("delta time",dt);

    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("displacement",dis_);
    // predicted dirichlet values
    // dis then also holds prescribed new dirichlet displacements
    discret_.EvaluateDirichlet(p,dis_,null,null,dirichtoggle_);
    discret_.ClearState();
    discret_.SetState("displacement",dis_);
    // predicted rhs
    discret_.EvaluateNeumann(p,*fext_);
    discret_.ClearState();
  }

  //----------------------- compute an inverse of the dirichtoggle vector
  invtoggle_->PutScalar(1.0);
  invtoggle_->Update(-1.0,*dirichtoggle_,1.0);

  // -------------------------------------------------------------------
  // call elements to calculate stiffness and mass
  // -------------------------------------------------------------------
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_nlnstiffmass");
    // choose what to assemble
    p.set("assemble matrix 1",true);
    p.set("assemble matrix 2",true);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",false);
    p.set("assemble vector 3",false);
    // other parameters that might be needed by the elements
    p.set("total time",time);
    p.set("delta time",dt);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("residual displacement",zeros_);
    discret_.SetState("displacement",dis_);
    //discret_.SetState("velocity",vel_); // not used at the moment
    discret_.Evaluate(p,stiff_,mass_,fint_,null,null);
    discret_.ClearState();

    //initialize Constraint Manager
    ConstrMan_=rcp(new ConstrManager(discret_, dis_));

    // Check for surface stress conditions due to interfacial phenomena
    vector<DRT::Condition*> surfstresscond(0);
    discret_.GetCondition("SurfaceStress",surfstresscond);
    if (surfstresscond.size())
    {
      // Determine local number of associated surface elements
      int localnumsurf=0;
      for (unsigned i=0;i<surfstresscond.size();++i)
      {
        localnumsurf+=surfstresscond[i]->Geometry().size();
      }

      // Create a map of the associated surface elements
      // This is needed because the SurfStressManager creates
      // several vectors containing history variables of the
      // corresponding surface elements. These vectors need to be
      // distributed over the processors in the same way as the
      // surface elements are.

      int localsurf[localnumsurf];
      int count=0;
      for (unsigned i=0;i<surfstresscond.size();++i)
      {
        for (map<int,RefCountPtr<DRT::Element> >::iterator iter=surfstresscond[i]->Geometry().begin();
             iter!=surfstresscond[i]->Geometry().end();++iter)
        {
          localsurf[count]=iter->first;
          count++;
        }
      }
      Epetra_Map surfmap(-1, localnumsurf, localsurf, 0, discret_.Comm());

      int numsurf;
      discret_.Comm().SumAll(&localnumsurf, &numsurf, 1);
      surf_stress_man_=rcp(new DRT::SurfStressManager(discret_, numsurf, surfmap));
    }
  }

  // close mass matrix
  mass_->Complete();

  // build damping matrix if desired
  if (damping)
  {
    stiff_->Complete();
    damp_->Add(*stiff_,false,kdamp,0.0);
    damp_->Add(*mass_,false,mdamp,1.0);
    damp_->Complete();
  }

  //--------------------------- calculate consistent initial accelerations
  {
    RefCountPtr<Epetra_Vector> rhs = LINALG::CreateVector(*dofrowmap,true);
    if (damping) damp_->Multiply(false,*vel_,*rhs);
    rhs->Update(-1.0,*fint_,1.0,*fext_,-1.0);
    Epetra_Vector rhscopy(*rhs);
    rhs->Multiply(1.0,*invtoggle_,rhscopy,0.0);
    solver.Solve(mass_->EpetraMatrix(),acc_,rhs,true,true);
  }

  //------------------------------------------------------ time step index
  step = 0;
  params_.set<int>("step",step);
  return;
} // StruGenAlpha::StruGenAlpha

/*----------------------------------------------------------------------*
 |  do constant predictor step (public)                      mwgee 03/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::ConstantPredictor()
{

  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time        = params_.get<double>("total time"     ,0.0);
  double dt          = params_.get<double>("delta time"     ,0.01);
  int    step        = params_.get<int>   ("step"           ,0);
  bool   damping     = params_.get<bool>  ("damping"        ,false);
  double alphaf      = params_.get<double>("alpha f"        ,0.459);
  bool   printscreen = params_.get<bool>  ("print to screen",false);
  string convcheck   = params_.get<string>("convcheck"      ,"AbsRes_Or_AbsDis");

  // store norms of old displacements and maximum of norms of
  // internal, external and inertial forces if a relative convergence
  // check is desired
  if (step != 0 and (convcheck != "AbsRes_And_AbsDis" or convcheck != "AbsRes_Or_AbsDis"))
  {
    CalcRefNorms();
  }

  // increment time and step
  double timen = time + dt;
  //int istep = step + 1;

  //--------------------------------------------------- predicting state
  // constant predictor : displacement in domain
  disn_->Update(1.0,*dis_,0.0);

  // apply new displacements at DBCs
  // and get new external force vector
  {
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_eleload");
    // choose what to assemble
    p.set("assemble matrix 1",false);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",false);
    p.set("assemble vector 3",false);
    // other parameters needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("displacement",disn_);
    // predicted dirichlet values
    // disn then also holds prescribed new dirichlet displacements
    discret_.EvaluateDirichlet(p,disn_,null,null,dirichtoggle_);
    discret_.ClearState();
    discret_.SetState("displacement",disn_);
    fextn_->PutScalar(0.0);  // initialize external force vector (load vect)
    discret_.EvaluateNeumann(p,*fextn_);
    discret_.ClearState();
  }

  // constant predictor
  veln_->Update(1.0,*vel_,0.0);
  accn_->Update(1.0,*acc_,0.0);

  //------------------------------ compute interpolated dis, vel and acc
  // constant predictor
  // mid-displacements D_{n+1-alpha_f} (dism)
  //    D_{n+1-alpha_f} := (1.-alphaf) * D_{n+1} + alpha_f * D_{n}
  dism_->Update(1.-alphaf,*disn_,alphaf,*dis_,0.0);
  velm_->Update(1.0,*vel_,0.0);
  accm_->Update(1.0,*acc_,0.0);

  //------------------------------- compute interpolated external forces
  // external mid-forces F_{ext;n+1-alpha_f} (fextm)
  //    F_{ext;n+1-alpha_f} := (1.-alphaf) * F_{ext;n+1}
  //                         + alpha_f * F_{ext;n}
  fextm_->Update(1.-alphaf,*fextn_,alphaf,*fext_,0.0);

  //------------- eval fint at interpolated state, eval stiffness matrix
  {
    // zero out stiffness
    stiff_->Zero();
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_nlnstiff");
    // choose what to assemble
    p.set("assemble matrix 1",true);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",false);
    p.set("assemble vector 3",false);
    // other parameters that might be needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("residual displacement",disi_);
    discret_.SetState("displacement",dism_);
    //discret_.SetState("velocity",velm_); // not used at the moment
    fint_->PutScalar(0.0);  // initialise internal force vector
    discret_.Evaluate(p,stiff_,null,fint_,null,null);
    discret_.ClearState();

    if (surf_stress_man_!=null)
    {
      p.set("surfstr_man", surf_stress_man_);
      surf_stress_man_->EvaluateSurfStress(p,dism_,fint_,stiff_);
    }

    // do NOT finalize the stiffness matrix, add mass and damping to it later
  }

  //-------------------------------------------- compute residual forces
  // Res = M . A_{n+1-alpha_m}
  //     + C . V_{n+1-alpha_f}
  //     + F_int(D_{n+1-alpha_f})
  //     - F_{ext;n+1-alpha_f}
  // add mid-inertial force
  mass_->Multiply(false,*accm_,*finert_);
  fresm_->Update(1.0,*finert_,0.0);
  // add mid-viscous damping force
  if (damping)
  {
    //RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,true);
    damp_->Multiply(false,*velm_,*fvisc_);
    fresm_->Update(1.0,*fvisc_,1.0);
  }

  // add static mid-balance
  fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);

  // blank residual at DOFs on Dirichlet BC
  {
    Epetra_Vector fresmcopy(*fresm_);
    fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
  }

  //------------------------------------------------ build residual norm
  double fresmnorm = 1.0;

  // store norms of displacements and maximum of norms of internal,
  // external and inertial forces if a relative convergence check
  // is desired and we are in the first time step
  if (step == 0 and (convcheck != "AbsRes_And_AbsDis" or convcheck != "AbsRes_Or_AbsDis"))
  {
    CalcRefNorms();
  }

  if (printscreen)
    fresm_->Norm2(&fresmnorm);
  if (!myrank_ and printscreen)
  {
    PrintPredictor(convcheck, fresmnorm);
  }

  return;
} // StruGenAlpha::ConstantPredictor()

/*----------------------------------------------------------------------*
 |  do matrix free constant predictor step (public)             lw 10/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::MatrixFreeConstantPredictor()
{
  // currently no preconditioner can be used in combination with this
  // matrix-free method
  ParameterList& solvparams = solver_.Params();
  string solvertype = solvparams.get("solver", "aztec");
  if (solvertype == "aztec")
  {
    ParameterList& azlist = solvparams.sublist("Aztec Parameters");
    int azprectype = azlist.get("AZ_precond", AZ_none);
    if (azprectype != AZ_none)
      dserror("Matrix-free Newton does not support preconditioning");
  }

  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time        = params_.get<double>("total time"     ,0.0);
  double dt          = params_.get<double>("delta time"     ,0.01);
  int    step        = params_.get<int>   ("step"           ,0);
  bool   damping     = params_.get<bool>  ("damping"        ,false);
  double alphaf      = params_.get<double>("alpha f"        ,0.459);
  bool   printscreen = params_.get<bool>  ("print to screen",false);
  string convcheck   = params_.get<string>("convcheck"      ,"AbsRes_Or_AbsDis");

  // store norms of old displacements and maximum of norms of
  // internal, external and inertial forces if a relative convergence
  // check is desired
  if (step != 0 and (convcheck != "AbsRes_And_AbsDis" or convcheck != "AbsRes_Or_AbsDis"))
  {
    CalcRefNorms();
  }

  // increment time and step
  double timen = time + dt;  // t_{n+1}
  //int istep = step + 1;  // n+1

  //--------------------------------------------------- predicting state
  // constant predictor : displacement in domain
  disn_->Update(1.0,*dis_,0.0);

  // apply new displacements at DBCs
  // and get new external force vector
  {
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_eleload");
    // choose what to assemble
    p.set("assemble matrix 1",false);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",false);
    p.set("assemble vector 3",false);
    // other parameters needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("displacement",disn_);
    // predicted dirichlet values
    // disn then also holds prescribed new dirichlet displacements
    discret_.EvaluateDirichlet(p,disn_,null,null,dirichtoggle_);
    discret_.ClearState();
    discret_.SetState("displacement",disn_);
    fextn_->PutScalar(0.0);  // initialize external force vector (load vect)
    discret_.EvaluateNeumann(p,*fextn_);
    discret_.ClearState();
  }

  // constant predictor
  veln_->Update(1.0,*vel_,0.0);
  accn_->Update(1.0,*acc_,0.0);

  //------------------------------ compute interpolated dis, vel and acc
  // constant predictor
  // mid-displacements D_{n+1-alpha_f} (dism)
  //    D_{n+1-alpha_f} := (1.-alphaf) * D_{n+1} + alpha_f * D_{n}
  dism_->Update(1.-alphaf,*disn_,alphaf,*dis_,0.0);
  velm_->Update(1.0,*vel_,0.0);
  accm_->Update(1.0,*acc_,0.0);

  //------------------------------- compute interpolated external forces
  // external mid-forces F_{ext;n+1-alpha_f} (fextm)
  //    F_{ext;n+1-alpha_f} := (1.-alphaf) * F_{ext;n+1}
  //                         + alpha_f * F_{ext;n}
  fextm_->Update(1.-alphaf,*fextn_,alphaf,*fext_,0.0);

  //------------------------------------ eval fint at interpolated state
  {
    stiff_->Zero();

    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_nlnstiff");
    // choose what to assemble
//     p.set("assemble matrix 1",false);
    p.set("assemble matrix 1",true);                           // test only!
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",false);
    p.set("assemble vector 3",false);
    // other parameters that might be needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("residual displacement",disi_);
    discret_.SetState("displacement",dism_);
    //discret_.SetState("velocity",velm_); // not used at the moment
    fint_->PutScalar(0.0);  // initialise internal force vector
//     discret_.Evaluate(p,null,null,fint_,null,null);
    discret_.Evaluate(p,stiff_,null,fint_,null,null);          // test only!
    discret_.ClearState();
    stiff_->Complete();                                 // test only!
  }

  //-------------------------------------------- compute residual forces
  // Res = M . A_{n+1-alpha_m}
  //     + C . V_{n+1-alpha_f}
  //     + F_int(D_{n+1-alpha_f})
  //     - F_{ext;n+1-alpha_f}
  // add mid-inertial force
  mass_->Multiply(false,*accm_,*finert_);
  fresm_->Update(1.0,*finert_,0.0);
  // add mid-viscous damping force
  if (damping)
  {
      //RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,true);
      damp_->Multiply(false,*velm_,*fvisc_);
      fresm_->Update(1.0,*fvisc_,1.0);
  }

  // add static mid-balance
  fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);

  // blank residual at DOFs on Dirichlet BC
  {
    Epetra_Vector fresmcopy(*fresm_);
    fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
  }

  //------------------------------------------------ build residual norm
  double fresmnorm = 1.0;

  // store norms of displacements and maximum of norms of internal,
  // external and inertial forces if a relative convergence check
  // is desired and we are in the first time step
  if (step == 0 and (convcheck != "AbsRes_And_AbsDis" or convcheck != "AbsRes_Or_AbsDis"))
  {
    CalcRefNorms();
  }

  if (printscreen)
    fresm_->Norm2(&fresmnorm);
  if (!myrank_ and printscreen)
  {
    PrintPredictor(convcheck, fresmnorm);
  }

  return;
} // StruGenAlpha::MatrixFreeConstantPredictor()


/*----------------------------------------------------------------------*
 |  do consistent predictor step (public)                    mwgee 07/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::ConsistentPredictor()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time        = params_.get<double>("total time"     ,0.0);
  double dt          = params_.get<double>("delta time"     ,0.01);
  int    step        = params_.get<int>   ("step"           ,0);
  bool   damping     = params_.get<bool>  ("damping"        ,false);
  double alphaf      = params_.get<double>("alpha f"        ,0.459);
  double alpham      = params_.get<double>("alpha m"        ,0.378);
  double beta        = params_.get<double>("beta"           ,0.292);
  double gamma       = params_.get<double>("gamma"          ,0.581);
  bool   printscreen = params_.get<bool>  ("print to screen",false);
  string convcheck   = params_.get<string>("convcheck"      ,"AbsRes_Or_AbsDis");

  // store norms of old displacements and maximum of norms of
  // internal, external and inertial forces if a relative convergence
  // check is desired
  if (step != 0 and (convcheck != "AbsRes_And_AbsDis" or convcheck != "AbsRes_Or_AbsDis"))
  {
    CalcRefNorms();
  }

  // increment time and step
  double timen = time + dt;  // t_{n+1}
  //int istep = step + 1;  // n+1

  //--------------------------------------------------- predicting state
  // constant predictor : displacement in domain
  disn_->Update(1.0,*dis_,0.0);

  // apply new displacements at DBCs
  // and get new external force vector
  {
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_eleload");
    // choose what to assemble
    p.set("assemble matrix 1",false);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",false);
    p.set("assemble vector 3",false);
    // other parameters needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("displacement",disn_);
    // predicted dirichlet values
    // disn then also holds prescribed new dirichlet displacements
    discret_.EvaluateDirichlet(p,disn_,null,null,dirichtoggle_);
    discret_.ClearState();
    discret_.SetState("displacement",disn_);
    fextn_->PutScalar(0.0);  // initialize external force vector (load vect)
    discret_.EvaluateNeumann(p,*fextn_);
    discret_.ClearState();
  }

  //cout << *disn_ << endl;

  // consistent predictor
  // predicting velocity V_{n+1} (veln)
  // V_{n+1} := gamma/(beta*dt) * (D_{n+1} - D_n)
  //          + (beta-gamma)/beta * V_n
  //          + (2.*beta-gamma)/(2.*beta) * A_n
  veln_->Update(1.0,*disn_,-1.0,*dis_,0.0);
  veln_->Update((beta-gamma)/beta,*vel_,
                (2.*beta-gamma)*dt/(2.*beta),*acc_,gamma/(beta*dt));


#ifdef STRUGENALPHA_STRONGDBC
  // apply new velocities at DBCs
  {
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_eleload");
    // choose what to assemble
    p.set("assemble matrix 1",false);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",false);
    p.set("assemble vector 3",false);
    // other parameters needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    //p.set("time derivative degree",1);  // we want velocities
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("velocity",veln_);
    // predicted dirichlet values
    // veln_ then also holds prescribed new Dirichlet velocities
    discret_.EvaluateDirichlet(p,null,veln_,null,dirichtoggle_);
    discret_.ClearState();
  }
#endif

  //cout << *veln_ << endl;

  // predicting accelerations A_{n+1} (accn)
  // A_{n+1} := 1./(beta*dt*dt) * (D_{n+1} - D_n)
  //          - 1./(beta*dt) * V_n
  //          + (2.*beta-1.)/(2.*beta) * A_n
  accn_->Update(1.0,*disn_,-1.0,*dis_,0.0);
  accn_->Update(-1./(beta*dt),*vel_,
                (2.*beta-1.)/(2.*beta),*acc_,1./(beta*dt*dt));

#ifdef STRUGENALPHA_STRONGDBC
  // apply new accelerations at DBCs
  {
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_eleload");
    // choose what to assemble
    p.set("assemble matrix 1",false);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",false);
    p.set("assemble vector 3",false);
    // other parameters needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    //p.set("time derivative degree",2);  // we want accelerations
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("acceleration",accn_);
    // predicted dirichlet values
    // accn_ then also holds prescribed new Dirichlet accelerations
    discret_.EvaluateDirichlet(p,null,null,accn_,dirichtoggle_);
    discret_.ClearState();
  }
#endif

  //cout << *accn_ << endl;


  //------------------------------ compute interpolated dis, vel and acc
  // consistent predictor
  // mid-displacements D_{n+1-alpha_f} (dism)
  //    D_{n+1-alpha_f} := (1.-alphaf) * D_{n+1} + alpha_f * D_{n}
  dism_->Update(1.-alphaf,*disn_,alphaf,*dis_,0.0);
  // mid-velocities V_{n+1-alpha_f} (velm)
  //    V_{n+1-alpha_f} := (1.-alphaf) * V_{n+1} + alpha_f * V_{n}
  velm_->Update(1.-alphaf,*veln_,alphaf,*vel_,0.0);
  // mid-accelerations A_{n+1-alpha_m} (accm)
  //    A_{n+1-alpha_m} := (1.-alpha_m) * A_{n+1} + alpha_m * A_{n}
  accm_->Update(1.-alpham,*accn_,alpham,*acc_,0.0);

  //------------------------------- compute interpolated external forces
  // external mid-forces F_{ext;n+1-alpha_f} (fextm)
  //    F_{ext;n+1-alpha_f} := (1.-alphaf) * F_{ext;n+1}
  //                         + alpha_f * F_{ext;n}
  fextm_->Update(1.-alphaf,*fextn_,alphaf,*fext_,0.0);

  //------------- eval fint at interpolated state, eval stiffness matrix
  {
    // zero out stiffness
    stiff_->Zero();
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_nlnstiff");
    // choose what to assemble
    p.set("assemble matrix 1",true);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",false);
    p.set("assemble vector 3",false);
    // other parameters that might be needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("residual displacement",disi_);
    discret_.SetState("displacement",dism_);
    //discret_.SetState("velocity",velm_); // not used at the moment
    fint_->PutScalar(0.0);  // initialise internal force vector
    discret_.Evaluate(p,stiff_,null,fint_,null,null);
    discret_.ClearState();

    if (surf_stress_man_!=null)
    {
      p.set("surfstr_man", surf_stress_man_);
      surf_stress_man_->EvaluateSurfStress(p,dism_,fint_,stiff_);
    }

    // do NOT finalize the stiffness matrix, add mass and damping to it later
  }

  //-------------------------------------------- compute residual forces
  // Res = M . A_{n+1-alpha_m}
  //     + C . V_{n+1-alpha_f}
  //     + F_int(D_{n+1-alpha_f})
  //     - F_{ext;n+1-alpha_f}
  // add mid-inertial force
  mass_->Multiply(false,*accm_,*finert_);
  fresm_->Update(1.0,*finert_,0.0);
  // add mid-viscous damping force
  if (damping)
  {
    //RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,true);
    damp_->Multiply(false,*velm_,*fvisc_);
    fresm_->Update(1.0,*fvisc_,1.0);
  }

  // add static mid-balance
  fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);

  // blank residual at DOFs on Dirichlet BC
  {
    Epetra_Vector fresmcopy(*fresm_);
    fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
  }

  //------------------------------------------------ build residual norm
  double fresmnorm = 1.0;

  // store norms of displacements and maximum of norms of internal,
  // external and inertial forces if a relative convergence check
  // is desired and we are in the first time step
  if (step == 0 and (convcheck != "AbsRes_And_AbsDis" or convcheck != "AbsRes_Or_AbsDis"))
  {
    CalcRefNorms();
  }

  if (printscreen)
    fresm_->Norm2(&fresmnorm);
  if (!myrank_ and printscreen)
  {
    PrintPredictor(convcheck, fresmnorm);
  }

  return;
} // StruGenAlpha::ConsistentPredictor()


/*----------------------------------------------------------------------*
 |  setup equilibrium with additional external forces        u.kue 02/08|
 *----------------------------------------------------------------------*/
void StruGenAlpha::ApplyExternalForce(const LINALG::MapExtractor& extractor,
                                      Teuchos::RCP<Epetra_Vector> iforce)
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time        = params_.get<double>("total time"     ,0.0);
  double dt          = params_.get<double>("delta time"     ,0.01);
  int    step        = params_.get<int>   ("step"           ,0);
  bool   damping     = params_.get<bool>  ("damping"        ,false);
  double alphaf      = params_.get<double>("alpha f"        ,0.459);
  //double alpham      = params_.get<double>("alpha m"        ,0.378);
  //double beta        = params_.get<double>("beta"           ,0.292);
  //double gamma       = params_.get<double>("gamma"          ,0.581);
  bool   printscreen = params_.get<bool>  ("print to screen",false);
  string convcheck   = params_.get<string>("convcheck"      ,"AbsRes_Or_AbsDis");

  // increment time and step
  double timen = time + dt;

  // get new external force vector
  {
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_eleload");
    // other parameters needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("displacement",disn_);
    fextn_->PutScalar(0.0);  // initialize external force vector (load vect)
    discret_.EvaluateNeumann(p,*fextn_);
    discret_.ClearState();
  }

  // Add iforce to fextn_
  // there might already be (body) forces at the interface nodes
  extractor.AddCondVector(iforce,fextn_);

  //------------------------------- compute interpolated external forces
  // external mid-forces F_{ext;n+1-alpha_f} (fextm)
  //    F_{ext;n+1-alpha_f} := (1.-alphaf) * F_{ext;n+1}
  //                         + alpha_f * F_{ext;n}
  fextm_->Update(1.-alphaf,*fextn_,alphaf,*fext_,0.0);

  //------------- eval fint at interpolated state, eval stiffness matrix
  {
    // zero out stiffness
    stiff_->Zero();
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_nlnstiff");
    // choose what to assemble
    p.set("assemble matrix 1",true);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",false);
    p.set("assemble vector 3",false);
    // other parameters that might be needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("residual displacement",disi_);
    discret_.SetState("displacement",dism_);
    //discret_.SetState("velocity",velm_); // not used at the moment
    fint_->PutScalar(0.0);  // initialise internal force vector
    discret_.Evaluate(p,stiff_,null,fint_,null,null);
    discret_.ClearState();

    if (surf_stress_man_!=null)
    {
      p.set("surfstr_man", surf_stress_man_);
      surf_stress_man_->EvaluateSurfStress(p,dism_,fint_,stiff_);
    }

    // do NOT finalize the stiffness matrix, add mass and damping to it later
  }

  //-------------------------------------------- compute residual forces
  // Res = M . A_{n+1-alpha_m}
  //     + C . V_{n+1-alpha_f}
  //     + F_int(D_{n+1-alpha_f})
  //     - F_{ext;n+1-alpha_f}
  // add mid-inertial force
  mass_->Multiply(false,*accm_,*finert_);
  fresm_->Update(1.0,*finert_,0.0);
  // add mid-viscous damping force
  if (damping)
  {
    //RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,true);
    damp_->Multiply(false,*velm_,*fvisc_);
    fresm_->Update(1.0,*fvisc_,1.0);
  }

  // add static mid-balance
  fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);

  // blank residual at DOFs on Dirichlet BC
  {
    Epetra_Vector fresmcopy(*fresm_);
    fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
  }

  //------------------------------------------------ build residual norm
  double fresmnorm = 1.0;

  // store norms of displacements and maximum of norms of internal,
  // external and inertial forces if a relative convergence check
  // is desired and we are in the first time step
  if (step == 0 and (convcheck != "AbsRes_And_AbsDis" or convcheck != "AbsRes_Or_AbsDis"))
  {
    CalcRefNorms();
  }

  if (printscreen)
    fresm_->Norm2(&fresmnorm);
  if (!myrank_ and printscreen)
  {
    PrintPredictor(convcheck, fresmnorm);
  }
}


/*----------------------------------------------------------------------*
 |  build linear system matrix and rhs (public)              u.kue 03/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::Evaluate(Teuchos::RCP<const Epetra_Vector> disp)
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time      = params_.get<double>("total time"             ,0.0);
  double dt        = params_.get<double>("delta time"             ,0.01);
  double timen     = time + dt;
  bool   damping   = params_.get<bool>  ("damping"                ,false);
  double beta      = params_.get<double>("beta"                   ,0.292);
  double gamma     = params_.get<double>("gamma"                  ,0.581);
  double alpham    = params_.get<double>("alpha m"                ,0.378);
  double alphaf    = params_.get<double>("alpha f"                ,0.459);

  // On the first call in a time step we have to have
  // disp==Teuchos::null. Then we just finished one of our predictors,
  // than contains the element loop, so we can fast forward and finish
  // up the linear system.

  if (disp!=Teuchos::null)
  {
    // set the new solution we just got
    disi_->Update(1.0,*disp,0.0);

    //---------------------------------- update mid configuration values
    // displacements
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
    dism_->Update(1.-alphaf,*disi_,1.0);
    // velocities
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // V_{n+1-alpha_f} := V_{n+1-alpha_f}
    //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
    velm_->Update((1.-alphaf)*gamma/(beta*dt),*disi_,1.0);
#else
    // incremental (required for constant predictor)
    velm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
    velm_->Update((beta-(1.0-alphaf)*gamma)/beta,*vel_,
                  (1.0-alphaf)*(2.*beta-gamma)*dt/(2.*beta),*acc_,
                  gamma/(beta*dt));
#endif
    // accelerations
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // A_{n+1-alpha_m} := A_{n+1-alpha_m}
    //                  + (1-alpha_m)/beta/dt^2*IncD_{n+1}
    accm_->Update((1.-alpham)/(beta*dt*dt),*disi_,1.0);
#else
    // incremental (required for constant predictor)
    accm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
    accm_->Update(-(1.-alpham)/(beta*dt),*vel_,
                  (2.*beta-1.+alpham)/(2.*beta),*acc_,
                  (1.-alpham)/((1.-alphaf)*beta*dt*dt));
#endif

    //---------------------------- compute internal forces and stiffness
    {
      // zero out stiffness
      stiff_->Zero();
      // create the parameters for the discretization
      ParameterList p;
      // action for elements
      p.set("action","calc_struct_nlnstiff");

      // other parameters that might be needed by the elements
      p.set("total time",timen);
      p.set("delta time",dt);

      // set vector values needed by elements
      discret_.ClearState();
      discret_.SetState("residual displacement",disi_);
      discret_.SetState("displacement",dism_);

      //discret_.SetState("velocity",velm_); // not used at the moment
      fint_->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,stiff_,fint_);
      discret_.ClearState();
      // do NOT finalize the stiffness matrix to add masses to it later
    }

    //------------------------------------------ compute residual forces
    // Res = M . A_{n+1-alpha_m}
    //     + C . V_{n+1-alpha_f}
    //     + F_int(D_{n+1-alpha_f})
    //     - F_{ext;n+1-alpha_f}
    // add inertia mid-forces
    mass_->Multiply(false,*accm_,*finert_);
    fresm_->Update(1.0,*finert_,0.0);
    // add viscous mid-forces
    if (damping)
    {
      //RCP<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,false);
      damp_->Multiply(false,*velm_,*fvisc_);
      fresm_->Update(1.0,*fvisc_,1.0);
    }
    // add static mid-balance
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
    // blank residual DOFs that are on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
    }
  }

  //------------------------------------------- effective rhs is fresm
  //---------------------------------------------- build effective lhs
  // (using matrix stiff_ as effective matrix)
  stiff_->Add(*mass_,false,(1.-alpham)/(beta*dt*dt),1.-alphaf);
  if (damping)
  {
    stiff_->Add(*damp_,false,(1.-alphaf)*gamma/(beta*dt),1.0);
  }
  stiff_->Complete();

  //----------------------- apply dirichlet BCs to system of equations
  disi_->PutScalar(0.0);  // Useful? depends on solver and more
  LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

} // StruGenAlpha::Evaluate()

/*----------------------------------------------------------------------*
 |  do Newton iteration (public)                             mwgee 03/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::FullNewton()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time      = params_.get<double>("total time"             ,0.0);
  double dt        = params_.get<double>("delta time"             ,0.01);
  double timen     = time + dt;
  int    maxiter   = params_.get<int>   ("max iterations"         ,10);
  bool   damping   = params_.get<bool>  ("damping"                ,false);
  double beta      = params_.get<double>("beta"                   ,0.292);
  double gamma     = params_.get<double>("gamma"                  ,0.581);
  double alpham    = params_.get<double>("alpha m"                ,0.378);
  double alphaf    = params_.get<double>("alpha f"                ,0.459);
  string convcheck = params_.get<string>("convcheck"              ,"AbsRes_Or_AbsDis");
  double toldisp   = params_.get<double>("tolerance displacements",1.0e-07);
  double tolres    = params_.get<double>("tolerance residual"     ,1.0e-07);
  bool printscreen = params_.get<bool>  ("print to screen",true);
  bool printerr    = params_.get<bool>  ("print to err",false);
  FILE* errfile    = params_.get<FILE*> ("err file",NULL);
  if (!errfile) printerr = false;

  // check whether we have a stiffness matrix, that is not filled yet
  // and mass and damping are present
  if (stiff_->Filled()) dserror("stiffness matrix may not be filled here");
  if (!mass_->Filled()) dserror("mass matrix must be filled here");
  if (damping)
    if (!damp_->Filled()) dserror("damping matrix must be filled here");

  //=================================================== equilibrium loop
  int numiter=0;
  double fresmnorm = 1.0e6;
  double disinorm = 1.0e6;
  fresm_->Norm2(&fresmnorm);
  Epetra_Time timer(discret_.Comm());
  timer.ResetStartTime();
  bool print_unconv = true;

  while (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) and numiter<=maxiter)
  {
    //------------------------------------------- effective rhs is fresm
    //---------------------------------------------- build effective lhs
    // (using matrix stiff_ as effective matrix)
    stiff_->Add(*mass_,false,(1.-alpham)/(beta*dt*dt),1.-alphaf);
    if (damping)
    {
      stiff_->Add(*damp_,false,(1.-alphaf)*gamma/(beta*dt),1.0);
    }
    stiff_->Complete();

    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

    //--------------------------------------------------- solve for disi
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (!numiter)
      solver_.Solve(stiff_->EpetraMatrix(),disi_,fresm_,true,true);
    else
      solver_.Solve(stiff_->EpetraMatrix(),disi_,fresm_,true,false);

    //---------------------------------- update mid configuration values
    // displacements
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
    dism_->Update(1.-alphaf,*disi_,1.0);
    // velocities
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // V_{n+1-alpha_f} := V_{n+1-alpha_f}
    //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
    velm_->Update((1.-alphaf)*gamma/(beta*dt),*disi_,1.0);
#else
    // incremental (required for constant predictor)
    velm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
    velm_->Update((beta-(1.0-alphaf)*gamma)/beta,*vel_,
                  (1.0-alphaf)*(2.*beta-gamma)*dt/(2.*beta),*acc_,
                  gamma/(beta*dt));
#endif
    // accelerations
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // A_{n+1-alpha_m} := A_{n+1-alpha_m}
    //                  + (1-alpha_m)/beta/dt^2*IncD_{n+1}
    accm_->Update((1.-alpham)/(beta*dt*dt),*disi_,1.0);
#else
    // incremental (required for constant predictor)
    accm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
    accm_->Update(-(1.-alpham)/(beta*dt),*vel_,
                  (2.*beta-1.+alpham)/(2.*beta),*acc_,
                  (1.-alpham)/((1.-alphaf)*beta*dt*dt));
#endif

    //---------------------------- compute internal forces and stiffness
    {
      // zero out stiffness
      stiff_->Zero();
      // create the parameters for the discretization
      ParameterList p;
      // action for elements
      p.set("action","calc_struct_nlnstiff");
      // choose what to assemble
      p.set("assemble matrix 1",true);
      p.set("assemble matrix 2",false);
      p.set("assemble vector 1",true);
      p.set("assemble vector 2",false);
      p.set("assemble vector 3",false);
      // other parameters that might be needed by the elements
      p.set("total time",timen);
      p.set("delta time",dt);
      // set vector values needed by elements
      discret_.ClearState();
      discret_.SetState("residual displacement",disi_);
      discret_.SetState("displacement",dism_);
      //discret_.SetState("velocity",velm_); // not used at the moment
      fint_->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,stiff_,null,fint_,null,null);
      discret_.ClearState();

      if (surf_stress_man_!=null)
      {
        p.set("surfstr_man", surf_stress_man_);
        surf_stress_man_->EvaluateSurfStress(p,dism_,fint_,stiff_);
      }

      if (ConstrMan_->HaveConstraint())
      {
          ConstrMan_->StiffnessAndInternalForces(time+dt,disn_,fint_,stiff_);
      }
      // do NOT finalize the stiffness matrix to add masses to it later
    }

    //------------------------------------------ compute residual forces
    // Res = M . A_{n+1-alpha_m}
    //     + C . V_{n+1-alpha_f}
    //     + F_int(D_{n+1-alpha_f})
    //     - F_{ext;n+1-alpha_f}
    // add inertia mid-forces
    mass_->Multiply(false,*accm_,*finert_);
    fresm_->Update(1.0,*finert_,0.0);
    // add viscous mid-forces
    if (damping)
    {
      //RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,false);
      damp_->Multiply(false,*velm_,*fvisc_);
      fresm_->Update(1.0,*fvisc_,1.0);
    }
    // add static mid-balance
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
    // blank residual DOFs that are on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
    }

    //---------------------------------------------- build residual norm
    disi_->Norm2(&disinorm);

    fresm_->Norm2(&fresmnorm);

    // a short message
    if (!myrank_ and (printscreen or printerr))
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }

    //--------------------------------- increment equilibrium loop index
    ++numiter;

  }
  //=================================================================== end equilibrium loop
  print_unconv = false;

  //-------------------------------- test whether max iterations was hit
  if (numiter>=maxiter)
  {
     dserror("Newton unconverged in %d iterations",numiter);
  }
  else
  {
    if (ConstrMan_->HaveMonitor())
    {
      ConstrMan_->ComputeMonitorValues(dism_);
    }
    if (!myrank_ and printscreen)
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }
  }

  params_.set<int>("num iterations",numiter);

  return;
} // StruGenAlpha::FullNewton()

/*----------------------------------------------------------------------*
 |  do Uzawa iteration (public)                                         |
 *----------------------------------------------------------------------*/
void StruGenAlpha::NonLinearUzawaFullNewton(int predictor)
{
    int  maxiterUzawa  = params_.get<int>   ("uzawa maxiter"          ,50);
    double Uzawa_param = params_.get<double>("uzawa parameter"        ,1);
    double tolconstr   = params_.get<double>("tolerance constraint"   ,1.0e-07);
    double alphaf      = params_.get<double>("alpha f"                ,0.459);
    double time        = params_.get<double>("total time"             ,0.0);
    double dt          = params_.get<double>("delta time"             ,0.01);


    FullNewton();
    //--------------------update end configuration
    disn_->Update(1./(1.-alphaf),*dism_,-alphaf/(1.-alphaf));
    //--------------------compute constraint error
    ConstrMan_->ComputeError(time+dt,disn_);
    double constrnorm=ConstrMan_->GetErrorNorm();
    cout<<"Constraint error for Newton solution: "<<constrnorm<<endl;
    int numiter_uzawa=0;
    while (constrnorm>tolconstr and numiter_uzawa <= maxiterUzawa)
    {
        // Lagrange multiplier is increased by Uzawa_param*ConstrErr
        ConstrMan_->UpdateLagrMult(Uzawa_param);
        // Keep new Lagrange multiplier fixed and solve for new displacements
        if      (predictor==1) ConstantPredictor();
        else if (predictor==2) ConsistentPredictor();
        ConstrMan_->StiffnessAndInternalForces(time+dt,disn_,fint_,stiff_);
        FullNewton();
        //--------------------update end configuration
        disn_->Update(1./(1.-alphaf),*dism_,-alphaf/(1.-alphaf));
        //--------------------compute constraint error
        ConstrMan_->ComputeError(time+dt,disn_);
        constrnorm=ConstrMan_->GetErrorNorm();
        cout<<"Constraint error for computed displacement: "<<constrnorm<<endl;
        numiter_uzawa++;
    }
    params_.set<int>("num iterations",numiter_uzawa+1);
}

/*----------------------------------------------------------------------*
 |                                                              tk 11/07|
 | Newton iteration respecting constraints.                             |
 | Uzawa algorithm is used to deal with lagrange multipliers            |
 *----------------------------------------------------------------------*/
void StruGenAlpha::FullNewtonLinearUzawa()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time      = params_.get<double>("total time"             ,0.0);
  double dt        = params_.get<double>("delta time"             ,0.01);
  double timen     = time + dt;
  int    maxiter   = params_.get<int>   ("max iterations"         ,10);
  bool   damping   = params_.get<bool>  ("damping"                ,false);
  double beta      = params_.get<double>("beta"                   ,0.292);
  double gamma     = params_.get<double>("gamma"                  ,0.581);
  double alpham    = params_.get<double>("alpha m"                ,0.378);
  double alphaf    = params_.get<double>("alpha f"                ,0.459);
  string convcheck = params_.get<string>("convcheck"              ,"AbsRes_Or_AbsDis");
  double toldisp   = params_.get<double>("tolerance displacements",1.0e-07);
  double tolres    = params_.get<double>("tolerance residual"     ,1.0e-07);
  double tolconstr    = params_.get<double>("tolerance constraint"     ,1.0e-07);
  bool printscreen = params_.get<bool>  ("print to screen",true);
  bool printerr    = params_.get<bool>  ("print to err",false);
  int  maxiterUzawa   = params_.get<int>   ("uzawa maxiter"         ,50);
  double Uzawa_param=params_.get<double>("uzawa parameter",1);
  FILE* errfile    = params_.get<FILE*> ("err file",NULL);
  if (!errfile) printerr = false;
  const Epetra_Map* dofrowmap = discret_.DofRowMap();

  // check whether we have a stiffness matrix, that is not filled yet
  // and mass and damping are present
  if (stiff_->Filled()) dserror("stiffness matrix may not be filled here");
  if (!mass_->Filled()) dserror("mass matrix must be filled here");
  if (damping)
  if (!damp_->Filled()) dserror("damping matrix must be filled here");

  //=================================================== equilibrium loop
  int numiter=0;
  double fresmnorm = 1.0e6;
  double disinorm = 1.0e6;
  fresm_->Norm2(&fresmnorm);

  double constrnorm=ConstrMan_->GetErrorNorm();
  Epetra_Time timer(discret_.Comm());
  timer.ResetStartTime();
  bool print_unconv = true;  
  while (!Converged(convcheck, disinorm, fresmnorm, constrnorm, toldisp, tolres, tolconstr) and numiter<=maxiter){
    //------------------------------------------- effective rhs is fresm
    //---------------------------------------------- build effective lhs
    // (using matrix stiff_ as effective matrix)
    stiff_->Add(*mass_,false,(1.-alpham)/(beta*dt*dt),1.-alphaf);
    if (damping)
    {
      stiff_->Add(*damp_,false,(1.-alphaf)*gamma/(beta*dt),1.0);
    }
    stiff_->Complete();
    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

    //===================================================uzawa loop
    // For every iteration step an uzawa algorithm is used to solve the linear system.
    //Preparation of uzawa method to solve the linear system.
    double norm_uzawa;
    double norm_uzawa_alt;
    double quotient;
    double norm_constr_uzawa;
    int numiter_uzawa=0;
    //counter used for adaptivity
    int count_paramadapt=1;

    ConstrMan_->ScaleLagrIncr(0.0);
    RCP<Epetra_Vector> constrVecWeight=LINALG::CreateVector(*dofrowmap,true);
    
    RCP<Epetra_Map> domainmap= ConstrMan_->GetConstraintMap();
    RCP<Epetra_Vector> dotprod= rcp(new Epetra_Vector(*domainmap));

    // Compute residual of the uzawa algorithm

    ConstrMan_->ComputeConstrTimesLagrIncr(constrVecWeight);
    RCP<Epetra_Vector> fresmcopy=rcp(new Epetra_Vector(*fresm_));
    fresmcopy->Update(1.0,*constrVecWeight,1.0);
    Epetra_Vector uzawa_res(*fresmcopy);
    (*stiff_).Multiply(false,*disi_,uzawa_res);
    uzawa_res.Update(1.0,*fresmcopy,-1.0);
    // blank residual DOFs which are on Dirichlet BC
    {
        Epetra_Vector rescopy(uzawa_res);
        uzawa_res.Multiply(1.0,*invtoggle_,rescopy,0.0);
    }
    uzawa_res.Norm2(&norm_uzawa);
    Epetra_Vector constr_res(*domainmap);
    ConstrMan_->ComputeConstrTimesDisi(*disi_,dotprod);
    constr_res.Update(1.0,*dotprod,1.0,*(ConstrMan_->GetError()),0.0);
    constr_res.Norm2(&norm_constr_uzawa);
    quotient =1;
    //Solve one iteration step with augmented lagrange
    //Since we calculate displacement norm as well, at least one step has to be taken
    while (((norm_uzawa > tolres/10 or norm_constr_uzawa>tolconstr/10)
            and numiter_uzawa < maxiterUzawa) or numiter_uzawa<1)
    {
        LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresmcopy,zeros_,dirichtoggle_);
        // solve for disi
        // Solve K . IncD = -R  ===>  IncD_{n+1}
        if (numiter_uzawa==0 and numiter==0)
        {
            solver_.Solve(stiff_->EpetraMatrix(),disi_,fresmcopy,true,true);
        }
        else
        {
            solver_.Solve(stiff_->EpetraMatrix(),disi_,fresmcopy,true,false);
        }

        //compute Lagrange multiplier increment
        ConstrMan_->ComputeConstrTimesDisi(*disi_,dotprod);
        ConstrMan_->UpdateLagrIncr(Uzawa_param, *dotprod);

        //Compute residual of the uzawa algorithm
        ConstrMan_->ComputeConstrTimesLagrIncr(constrVecWeight);
        fresmcopy->Update(1.0,*constrVecWeight,1.0,*fresm_,0.0);
        Epetra_Vector uzawa_res(*fresmcopy);
        (*stiff_).Multiply(false,*disi_,uzawa_res);
        uzawa_res.Update(1.0,*fresmcopy,-1.0);
        // blank residual DOFs which are on Dirichlet BC
        {
            Epetra_Vector rescopy(uzawa_res);
            uzawa_res.Multiply(1.0,*invtoggle_,rescopy,0.0);
        }
        norm_uzawa_alt=norm_uzawa;
        uzawa_res.Norm2(&norm_uzawa);
        Epetra_Vector constr_res(*domainmap);

        constr_res.Update(1.0,*dotprod,1.0,*(ConstrMan_->GetError()),0.0);
        constr_res.Norm2(& norm_constr_uzawa);

        //-------------Adapt Uzawa parameter--------------
        // For a constant parameter the quotient of two successive residual norms
        // stays constant during the computation. So this quotient seems to be a good
        // measure for the parameter choice

        // Adaptivity only takes place every second step. Otherwise the quotient is not significant.
        if (count_paramadapt>=2)
        {
            double quotient_neu=norm_uzawa/norm_uzawa_alt;
            // In case of divergence the parameter must be choosen too high
            if (quotient_neu>1)
            {
                Uzawa_param=Uzawa_param/2.;
                count_paramadapt=0;
                quotient=1;
            }
            else
            {
                // In case the newly computed quotient is better than the one obtained from the
                // previous parameter, the parameter is increased by a factor (1+quotient_neu)
                if (quotient>quotient_neu)
                {
                    Uzawa_param=Uzawa_param*(1+quotient_neu);
                    quotient=quotient_neu;
                    count_paramadapt=0;
                }
                // In case the newly computed quotient is worse than the one obtained from the
                // previous parameter, the parameter is decreased by a factor 1/(1+quotient_neu)
                else
                {
                    Uzawa_param=Uzawa_param/(1+quotient_neu);
                    quotient=quotient_neu;
                    count_paramadapt=0;
                }
            }
        }
        count_paramadapt++;
        numiter_uzawa++;
    } //Uzawa loop

    if (!myrank_)
    {
        cout<<"Uzawa steps "<<numiter_uzawa<<endl;
    }

    //update lagrange multiplier
    ConstrMan_->UpdateLagrMult();

    //---------------------------------- update mid and end configuration values
    // displacements
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
    dism_->Update(1.-alphaf,*disi_,1.0);
    disn_->Update(1.,*disi_,1.0);
    // velocities
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // V_{n+1-alpha_f} := V_{n+1-alpha_f}
    //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
    velm_->Update((1.-alphaf)*gamma/(beta*dt),*disi_,1.0);
#else
    // incremental (required for constant predictor)
    velm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
    velm_->Update((beta-(1.0-alphaf)*gamma)/beta,*vel_,
                  (1.0-alphaf)*(2.*beta-gamma)*dt/(2.*beta),*acc_,
                  gamma/(beta*dt));
#endif
    // accelerations
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // A_{n+1-alpha_m} := A_{n+1-alpha_m}
    //                  + (1-alpha_m)/beta/dt^2*IncD_{n+1}
    accm_->Update((1.-alpham)/(beta*dt*dt),*disi_,1.0);
#else
    // incremental (required for constant predictor)
    accm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
    accm_->Update(-(1.-alpham)/(beta*dt),*vel_,
                  (2.*beta-1.+alpham)/(2.*beta),*acc_,
                  (1.-alpham)/((1.-alphaf)*beta*dt*dt));
#endif

    //---------------------------- compute internal forces and stiffness
    {
      // zero out stiffness
      stiff_->Zero();
      // create the parameters for the discretization
      ParameterList p;
      // action for elements
      p.set("action","calc_struct_nlnstiff");
      // choose what to assemble
      p.set("assemble matrix 1",true);
      p.set("assemble matrix 2",false);
      p.set("assemble vector 1",true);
      p.set("assemble vector 2",false);
      p.set("assemble vector 3",false);
      // other parameters that might be needed by the elements
      p.set("total time",timen);
      p.set("delta time",dt);
      // set vector values needed by elements
      discret_.ClearState();
      discret_.SetState("residual displacement",disi_);
      discret_.SetState("displacement",dism_);
      //discret_.SetState("velocity",velm_); // not used at the moment
      fint_->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,stiff_,null,fint_,null,null);
      discret_.ClearState();

      ConstrMan_->StiffnessAndInternalForces(timen,disn_,fint_,stiff_);
      constrnorm=ConstrMan_->GetErrorNorm();
      // do NOT finalize the stiffness matrix to add masses to it later
    }

    //------------------------------------------ compute residual forces
    // Res = M . A_{n+1-alpha_m}
    //     + C . V_{n+1-alpha_f}
    //     + F_int(D_{n+1-alpha_f})
    //     - F_{ext;n+1-alpha_f}
    // add inertia mid-forces
    mass_->Multiply(false,*accm_,*finert_);
    fresm_->Update(1.0,*finert_,0.0);
    // add viscous mid-forces
    if (damping)
    {
      //RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,false);
      damp_->Multiply(false,*velm_,*fvisc_);
      fresm_->Update(1.0,*fvisc_,1.0);
    }
    // add static mid-balance
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
    // blank residual DOFs that are on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
    }

    //---------------------------------------------- build residual norm
    disi_->Norm2(&disinorm);

    fresm_->Norm2(&fresmnorm);

    // a short message
    if (!myrank_ and (printscreen or printerr))
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck,constrnorm,Uzawa_param);
    }

    //--------------------------------- increment equilibrium loop index
    ++numiter;


  }
  //=================================================================== end equilibrium loop
  print_unconv = false;

  //-------------------------------- test whether max iterations was hit
  if (numiter>=maxiter)
  {
     dserror("Newton unconverged in %d iterations",numiter);
  }
  else
  {
    if (ConstrMan_->HaveMonitor())
    {
      ConstrMan_->ComputeMonitorValues(dism_);
    }
    if (!myrank_ and printscreen)
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck,constrnorm, Uzawa_param);
    }
  }

  params_.set<int>("num iterations",numiter);
  params_.set<double>("uzawa parameter",Uzawa_param);

  return;
} // StruGenAlpha::FullNewtonUzawa()



/*----------------------------------------------------------------------*
 |  do modified Newton iteration (public)                    mwgee 03/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::ModifiedNewton()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time      = params_.get<double>("total time"             ,0.0);
  double dt        = params_.get<double>("delta time"             ,0.01);
  double timen     = time + dt;
  int    maxiter   = params_.get<int>   ("max iterations"         ,10);
  bool   damping   = params_.get<bool>  ("damping"                ,false);
  double beta      = params_.get<double>("beta"                   ,0.292);
  double gamma     = params_.get<double>("gamma"                  ,0.581);
  double alpham    = params_.get<double>("alpha m"                ,0.378);
  double alphaf    = params_.get<double>("alpha f"                ,0.459);
  string convcheck = params_.get<string>("convcheck"              ,"AbsRes_Or_AbsDis");
  double toldisp   = params_.get<double>("tolerance displacements",1.0e-07);
  double tolres    = params_.get<double>("tolerance residual"     ,1.0e-07);
  bool printscreen = params_.get<bool>  ("print to screen",true);
  bool printerr    = params_.get<bool>  ("print to err",false);
  FILE* errfile    = params_.get<FILE*> ("err file",NULL);
  if (!errfile) printerr = false;
  //const Epetra_Map* dofrowmap = discret_.DofRowMap();

  // check whether we have a stiffness matrix, that is not filled yet
  // and mass and damping are present
  if (stiff_->Filled()) dserror("stiffness matrix may not be filled here");
  if (!mass_->Filled()) dserror("mass matrix must be filled here");
  if (damping)
    if (!damp_->Filled()) dserror("damping matrix must be filled here");

  int numiter=0;
  //---------------------------------------------- build effective lhs
  // (using matrix stiff_ as effective matrix)
  stiff_->Add(*mass_,false,(1.-alpham)/(beta*dt*dt),1.-alphaf);
  if (damping)
    stiff_->Add(*damp_,false,(1.-alphaf)*gamma/(beta*dt),1.0);
  stiff_->Complete();
  LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

  //=================================================== equilibrium loop
  double fresmnorm = 1.0e6;
  double disinorm = 1.0e6;
  fresm_->Norm2(&fresmnorm);
  Epetra_Time timer(discret_.Comm());
  timer.ResetStartTime();
  bool print_unconv = true;

  while (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres)  and numiter<=maxiter)
  {
    //------------------------------------------- effective rhs is fresm
    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more

    //--------------------------------------------------- solve for disi
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (!numiter)
      solver_.Solve(stiff_->EpetraMatrix(),disi_,fresm_,true,true);
    else
      solver_.Solve(stiff_->EpetraMatrix(),disi_,fresm_,false,false);

    //---------------------------------- update mid configuration values
    // displacements
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
    dism_->Update(1.-alphaf,*disi_,1.0);
    // velocities
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // V_{n+1-alpha_f} := V_{n+1-alpha_f}
    //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
    velm_->Update((1.-alphaf)*gamma/(beta*dt),*disi_,1.0);
#else
    // incremental (required for constant predictor)
    velm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
    velm_->Update((beta-(1.0-alphaf)*gamma)/beta,*vel_,
                  (1.0-alphaf)*(2.*beta-gamma)*dt/(2.*beta),*acc_,
                  gamma/(beta*dt));
#endif
    // accelerations
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // A_{n+1-alpha_m} := A_{n+1-alpha_m}
    //                  + (1-alpha_m)/beta/dt^2*IncD_{n+1}
    accm_->Update((1.-alpham)/(beta*dt*dt),*disi_,1.0);
#else
    // incremental (required for constant predictor)
    accm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
    accm_->Update(-(1.-alpham)/(beta*dt),*vel_,
                  (2.*beta-1.+alpham)/(2.*beta),*acc_,
                  (1.-alpham)/((1.-alphaf)*beta*dt*dt));
#endif

    //----------------------------------------- compute internal forces
    {
      // create the parameters for the discretization
      ParameterList p;
      // action for elements
      p.set("action","calc_struct_nlnstiff");
      // choose what to assemble
      p.set("assemble matrix 1",false);
      p.set("assemble matrix 2",false);
      p.set("assemble vector 1",true);
      p.set("assemble vector 2",false);
      p.set("assemble vector 3",false);
      // other parameters that might be needed by the elements
      p.set("total time",timen);
      p.set("delta time",dt);
      // set vector values needed by elements
      discret_.ClearState();
      discret_.SetState("residual displacement",disi_);
      discret_.SetState("displacement",dism_);
      //discret_.SetState("velocity",velm_); // not used at the moment
      fint_->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,null,null,fint_,null,null);
      discret_.ClearState();
    }

    //------------------------------------------ compute residual forces
    // Res = M . A_{n+1-alpha_m}
    //     + C . V_{n+1-alpha_f}
    //     + F_int(D_{n+1-alpha_f})
    //     - F_{ext;n+1-alpha_f}
    // add inertia mid-forces
    mass_->Multiply(false,*accm_,*finert_);
    fresm_->Update(1.0,*finert_,0.0);
    // add viscous mid-forces
    if (damping)
    {
      //RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,false);
      damp_->Multiply(false,*velm_,*fvisc_);
      fresm_->Update(1.0,*fvisc_,1.0);
    }
    // add static mid-balance
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
    // blank residual DOFs with are on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
    }

    //---------------------------------------------- build residual norm
    disi_->Norm2(&disinorm);

    fresm_->Norm2(&fresmnorm);
    // a short message
    if (!myrank_ and (printscreen or printerr))
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }

    //--------------------------------- increment equilibrium loop index
    ++numiter;
  }
  //============================================= end equilibrium loop
  print_unconv = false;

  //-------------------------------- test whether max iterations was hit
  if (numiter>=maxiter)
  {
     dserror("Newton unconverged in %d iterations",numiter);
  }
  else
  {
    if (ConstrMan_->HaveMonitor())
    {
      ConstrMan_->ComputeMonitorValues(dism_);
    }
    if (!myrank_ and printscreen)
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }
  }

  params_.set<int>("num iterations",numiter);

  return;
} // StruGenAlpha::ModifiedNewton()


/*----------------------------------------------------------------------*
 |  do matrix free Newton iteration (public)                    lw 10/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::MatrixFreeNewton()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time      = params_.get<double>("total time"             ,0.0);
  double dt        = params_.get<double>("delta time"             ,0.01);
  double timen     = time + dt;
  int    maxiter   = params_.get<int>   ("max iterations"         ,10);
  bool   damping   = params_.get<bool>  ("damping"                ,false);
  double beta      = params_.get<double>("beta"                   ,0.292);
  double gamma     = params_.get<double>("gamma"                  ,0.581);
  double alpham    = params_.get<double>("alpha m"                ,0.378);
  double alphaf    = params_.get<double>("alpha f"                ,0.459);
  string convcheck = params_.get<string>("convcheck"              ,"AbsRes_Or_AbsDis");
  double toldisp   = params_.get<double>("tolerance displacements",1.0e-07);
  double tolres    = params_.get<double>("tolerance residual"     ,1.0e-07);
  bool printscreen = params_.get<bool>  ("print to screen",true);
  bool printerr    = params_.get<bool>  ("print to err",false);
  FILE* errfile    = params_.get<FILE*> ("err file",NULL);
  if (!errfile) printerr = false;
  //const Epetra_Map* dofrowmap = discret_.DofRowMap();

//   int numiter=0;
//   //---------------------------------------------- build effective lhs
//   LINALG::ApplyDirichlettoSystem(disi_,fresm_,zeros_,dirichtoggle_);

//   //=================================================== equilibrium loop
//   double fresmnorm = 1.0e6;
//   double disinorm = 1.0e6;
//   fresm_->Norm2(&fresmnorm);
//   Epetra_Time timer(discret_.Comm());
//   timer.ResetStartTime();
//   bool print_unconv = true;

  // check whether we have a stiffness matrix, that is not filled yet
  // and mass and damping are present          -------------------> test only !!!!!!!!!!!!!!
  if (stiff_->Filled()) dserror("stiffness matrix may not be filled here");
  if (!mass_->Filled()) dserror("mass matrix must be filled here");
  if (damping)
    if (!damp_->Filled()) dserror("damping matrix must be filled here");

  //=================================================== equilibrium loop
  int numiter=0;
  double fresmnorm = 1.0e6;
  double disinorm = 1.0e6;
  fresm_->Norm2(&fresmnorm);
  Epetra_Time timer(discret_.Comm());
  timer.ResetStartTime();
  bool print_unconv = true;

  while (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres)  and numiter<=maxiter)
  {

    // (using matrix stiff_ as effective matrix) ----------> test only!!!!!!!
    stiff_->Add(*mass_,false,(1.-alpham)/(beta*dt*dt),1.-alphaf);
    if (damping)
    {
      stiff_->Add(*damp_,false,(1.-alphaf)*gamma/(beta*dt),1.0);
    }
    stiff_->Complete();

    //------------------------------------------- effective rhs is fresm
    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

    //----------------------------------------- build MatrixFreeOperator
    RCP<Epetra_Operator> mfop = rcp(new LINALG::MatrixFreeOperator(*this, stiff_));

    //--------------------------------------------------- solve for disi
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (!numiter)
      solver_.Solve(mfop,disi_,fresm_,true,true);
    else
      solver_.Solve(mfop,disi_,fresm_,false,false);

    //---------------------------------- update mid configuration values
    // displacements
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
    dism_->Update(1.-alphaf,*disi_,1.0);
    // velocities
//#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // V_{n+1-alpha_f} := V_{n+1-alpha_f}
    //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
//    velm_->Update((1.-alphaf)*gamma/(beta*dt),*disi_,1.0);
//#else
    // incremental (required for constant predictor)
    velm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
    velm_->Update((beta-(1.0-alphaf)*gamma)/beta,*vel_,
                  (1.0-alphaf)*(2.*beta-gamma)*dt/(2.*beta),*acc_,
                  gamma/(beta*dt));
//#endif
    // accelerations
//#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // A_{n+1-alpha_m} := A_{n+1-alpha_m}
    //                  + (1-alpha_m)/beta/dt^2*IncD_{n+1}
//    accm_->Update((1.-alpham)/(beta*dt*dt),*disi_,1.0);
//#else
    // incremental (required for constant predictor)
    accm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
    accm_->Update(-(1.-alpham)/(beta*dt),*vel_,
                  (2.*beta-1.+alpham)/(2.*beta),*acc_,
                  (1.-alpham)/((1.-alphaf)*beta*dt*dt));
//#endif

    //----------------------------------------- compute internal forces
    {
      // create the parameters for the discretization
      ParameterList p;
      // action for elements
      p.set("action","calc_struct_nlnstiff");
      // choose what to assemble
      p.set("assemble matrix 1",false);
      p.set("assemble matrix 2",false);
      p.set("assemble vector 1",true);
      p.set("assemble vector 2",false);
      p.set("assemble vector 3",false);
      // other parameters that might be needed by the elements
      p.set("total time",timen);
      p.set("delta time",dt);
      // set vector values needed by elements
      discret_.ClearState();
      discret_.SetState("residual displacement",disi_);
      discret_.SetState("displacement",dism_);
      //discret_.SetState("velocity",velm_); // not used at the moment
      fint_->PutScalar(0.0);  // initialise internal force vector
      //discret_.Evaluate(p,null,null,fint_,null,null);
      discret_.Evaluate(p,stiff_,null,fint_,null,null);
      discret_.ClearState();

//       if surface stress etc should be incorporated, be sure that
//       they are also incorporated in computeFmatrixfree!!!!!!!!
//
//       if (surf_stress_man_!=null)
//       {
//         p.set("surfstr_man", surf_stress_man_);
//         surf_stress_man_->EvaluateSurfStress(p,dism_,fint_,stiff_);
//       }

//       if (ConstrMan_->HaveConstraint())
//       {
//     	  ConstrMan_->StiffnessAndInternalForces(time+dt,disn_,fint_,stiff_);
//       }
      // do NOT finalize the stiffness matrix to add masses to it later
    }

    //------------------------------------------ compute residual forces
    // Res = M . A_{n+1-alpha_m}
    //     + C . V_{n+1-alpha_f}
    //     + F_int(D_{n+1-alpha_f})
    //     - F_{ext;n+1-alpha_f}
    // add inertia mid-forces
    mass_->Multiply(false,*accm_,*finert_);
    fresm_->Update(1.0,*finert_,0.0);
    // add viscous mid-forces
    if (damping)
    {
      //RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,false);
      damp_->Multiply(false,*velm_,*fvisc_);
      fresm_->Update(1.0,*fvisc_,1.0);
    }
    // add static mid-balance
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
    // blank residual DOFs with are on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
    }

    //---------------------------------------------- build residual norm
    disi_->Norm2(&disinorm);

    fresm_->Norm2(&fresmnorm);
    // a short message
    if (!myrank_ and (printscreen or printerr))
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }

    //--------------------------------- increment equilibrium loop index
    ++numiter;
  }
  //============================================= end equilibrium loop
  print_unconv = false;

  //-------------------------------- test whether max iterations was hit
  if (numiter>=maxiter)
  {
     dserror("Newton unconverged in %d iterations",numiter);
  }
  else
  {
     if (!myrank_ and printscreen)
     {
       PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                   fresmnorm,disinorm,convcheck);
     }
  }

  params_.set<int>("num iterations",numiter);

  return;
} // StruGenAlpha::MatrixFreeNewton()


/*----------------------------------------------------------------------*
 |  do nonlinear cg iteration (public)                       mwgee 03/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::NonlinearCG()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  //double time      = params_.get<double>("total time"             ,0.0);
  double dt        = params_.get<double>("delta time"             ,0.01);
  int    maxiter   = params_.get<int>   ("max iterations"         ,10);
  bool   damping   = params_.get<bool>  ("damping"                ,false);
  double beta      = params_.get<double>("beta"                   ,0.292);
  double gamma     = params_.get<double>("gamma"                  ,0.581);
  double alpham    = params_.get<double>("alpha m"                ,0.378);
  double alphaf    = params_.get<double>("alpha f"                ,0.459);
  double toldisp   = params_.get<double>("tolerance displacements",1.0e-07);
  double tolres    = params_.get<double>("tolerance residual"     ,1.0e-07);
  //bool printscreen = params_.get<bool>  ("print to screen",true);
  bool printerr    = params_.get<bool>  ("print to err",false);
  FILE* errfile    = params_.get<FILE*> ("err file",NULL);
  if (!errfile) printerr = false;

  // get ml stuff from input file
  bool hasml = solver_.Params().isSublist("ML Parameters");
  if (!hasml) dserror("Solver does not have ML configured from input file");
  ParameterList& linearmllist = solver_.Params().sublist("ML Parameters");
  int outlevel = linearmllist.get<int>("output",0);
  int maxlevel = linearmllist.get<int>("max levels",3);
  int maxcsize = linearmllist.get<int>("coarse: max size",1);

  // check whether we have a stiffness matrix, that is not filled yet
  // and mass and damping are present
  if (stiff_==null)     dserror("stiffness matrix = null");
  if (stiff_->Filled()) dserror("stiffness matrix may not be filled here");
  if (!mass_->Filled()) dserror("mass matrix must be filled here");
  if (damping)
    if (!damp_->Filled()) dserror("damping matrix must be filled here");

  //---------------------------------------------- set initial guess
  disi_->PutScalar(0.0);

  // create the nox parameter list if it does not exist
  ParameterList& noxparams = params_.sublist("nox parameters");
  RefCountPtr<Teuchos::ParameterList> rcpparams = rcp(&noxparams);
  rcpparams.release();

  ParameterList& printParams = noxparams.sublist("Printing");
  printParams.set("MyPID",myrank_);
  printParams.set("Output Precision", 9);
  printParams.set("Output Processor", 0);
  if (outlevel)
    printParams.set("Output Information",
                    NOX::Utils::OuterIteration +
                    //NOX::Utils::OuterIterationStatusTest +
                    //NOX::Utils::InnerIteration +
                    //NOX::Utils::Parameters +
                    //NOX::Utils::Details +
                    NOX::Utils::Warning
                    );
  else
    printParams.set("Output Information",0);

  // Set the nonlinear solver method as line search
  noxparams.set("Nonlinear Solver","Line Search Based");

  // get sublist for type of linesearch
  ParameterList& searchParams = noxparams.sublist("Line Search");
  searchParams.set("Method","NonlinearCG");

  // Sublist for direction
  ParameterList& dirParams = noxparams.sublist("Direction");
  dirParams.set("Method", "NonlinearCG");

  ParameterList& nlcgParams = dirParams.sublist("Nonlinear CG");
  nlcgParams.set("Precondition","On");
  nlcgParams.set("Orthogonalize","Polak-Ribiere");
  //nlcgParams.set("Orthogonalize", "Fletcher-Reeves");
  nlcgParams.set("Restart Frequency", 25);

  ParameterList& lsParams = nlcgParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");
  lsParams.set("Max Iterations", 100);
  lsParams.set("Tolerance", 1e-11);
  lsParams.set("Output Frequency", 10);
  //lsParams.set("Preconditioning", "None");
  //lsParams.set("Preconditioner","None");
  lsParams.set("Preconditioning", "User Supplied Preconditioner");
  lsParams.set("Preconditioner","User Defined");

  // create the nonlinear ml parameter list
  ParameterList& mlparams = params_.sublist("ml parameters");
  mlparams.set("nlnML output",                                      outlevel   ); // ML-output-level (0-10)
  mlparams.set("nlnML max levels",                                  maxlevel   ); // max. # levels (minimum = 2 !)
  mlparams.set("nlnML coarse: max size",                            maxcsize   ); // the size ML stops generating coarser levels
  mlparams.set("nlnML is linear preconditioner",                    false      );
  mlparams.set("nlnML is matrixfree",                               false      );
  mlparams.set("nlnML apply constraints",                           false      );
  mlparams.set("nlnML Jacobian fix diagonal",                       false      );
  mlparams.set("nlnML finite difference fine level",                false      );
  mlparams.set("nlnML finite difference alpha",                     1.0e-08    );
  mlparams.set("nlnML finite difference beta",                      1.0e-07    );
  mlparams.set("nlnML finite difference centered",                  false      );

  mlparams.set("nlnML absolute residual tolerance",                 tolres     );
  mlparams.set("nlnML max cycles",                                  maxiter    );
  mlparams.set("nlnML adaptive recompute",                          0.0        ); // recompute if residual is larger then this value
  mlparams.set("nlnML offset recompute",                            0          ); // every offset this preconditioner is recomputed
  mlparams.set("nlnML additional adaptive nullspace",               0          ); // compute adaptive nullspace (additional kernel vectors)
  mlparams.set("nlnML PDE equations",                               6          ); // dof per node
  mlparams.set("nlnML null space: dimension",                       6          ); // dimension of nullspace
  mlparams.set("nlnML spatial dimension",                           3          );
  mlparams.set("nlnML coarse: type",                                "Uncoupled"); // Uncoupled METIS VBMETIS
  mlparams.set("nlnML nodes per aggregate",                         9          ); // # nodes per agg for coarsening METIS and VBMETIS

  mlparams.set("nlnML use nlncg on fine level",                     true); // use nlnCG or mod. Newton's method
  mlparams.set("nlnML use nlncg on medium level",                   true);
  mlparams.set("nlnML use nlncg on coarsest level",                 true);

  mlparams.set("nlnML max iterations newton-krylov fine level",     100); // # iterations of lin. CG in mod. Newton's method
  mlparams.set("nlnML max iterations newton-krylov medium level" ,  50);
  mlparams.set("nlnML max iterations newton-krylov coarsest level", 150);

  mlparams.set("nlnML linear smoother type fine level",             "MLS"); // MLS SGS BSGS Jacobi MLS Bcheby AmesosKLU
  mlparams.set("nlnML linear smoother type medium level",           "MLS");
  mlparams.set("nlnML linear smoother type coarsest level",         "AmesosKLU");
  mlparams.set("nlnML linear smoother sweeps fine level",           24);
  mlparams.set("nlnML linear smoother sweeps medium level",         24);
  mlparams.set("nlnML linear smoother sweeps coarsest level",       1);

  mlparams.set("nlnML nonlinear presmoothing sweeps fine level",    0);
  mlparams.set("nlnML nonlinear presmoothing sweeps medium level",  0);
  mlparams.set("nlnML nonlinear smoothing sweeps coarse level",     15);
  mlparams.set("nlnML nonlinear postsmoothing sweeps medium level", 3);
  mlparams.set("nlnML nonlinear postsmoothing sweeps fine level",   10);

  // create the fine level interface if it does not exist
  if (fineinterface_==null)
  {
    int printlevel = mlparams.get("nlnML output",6);
    fineinterface_ = rcp(new NoxInterface(*this,printlevel));
  }

  // create the nonlinear ml preconditioner if it does not exist
  if (prec_==null)
    prec_ = rcp(new NLNML::NLNML_Preconditioner(fineinterface_,mlparams,discret_.Comm()));
  // tell preconditioner to recompute from scratch
  prec_->setinit(false);

#if 0 // use the nonlinear preconditioner as a solver without the outer nox loop
  {
    Epetra_Time timer(discret_.Comm());
    double t0 = timer.ElapsedTime();
    //prec_->solve_variant();
    prec_->solve();
    double t1 = timer.ElapsedTime();

    // get status and print output message
    if (outlevel and myrank_==0)
    {
      printf("NOX/ML :============solve time incl. setup : %15.4f sec\n",t1-t0);
      double appltime = fineinterface_->getsumtime();
      printf("NOX/ML :===========of which time in ccarat : %15.4f sec\n",appltime);
      cout << "NOX/ML :======number calls to computeF in this solve : "
           << fineinterface_->getnumcallscomputeF() << "\n\n\n";
      fflush(stdout);
    }
    fineinterface_->resetsumtime();
    fineinterface_->setnumcallscomputeF(0);

  //---------------------------------- update mid configuration values
  // displacements
  // incremental (disi is now D_{n+1}-D_{n})
  dism_->Update((1.-alphaf),*disi_,1.0);

  // velocities
#ifdef STRUGENALPHA_INCRUPDT
  // incremental (required for constant predictor)
  velm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
  velm_->Update((beta-(1.0-alphaf)*gamma)/beta,*vel_,
                (1.0-alphaf)*(2.*beta-gamma)*dt/(2.*beta),*acc_,
                gamma/(beta*dt));
#else
  dserror("Please verify that incremental update is not needed here!");
#endif

  // accelerations
#ifdef STRUGENALPHA_INCRUPDT
  // incremental (required for constant predictor)
  accm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
  accm_->Update(-(1.-alpham)/(beta*dt),*vel_,
                (2.*beta-1.+alpham)/(2.*beta),*acc_,
                (1.-alpham)/((1.-alphaf)*beta*dt*dt));
#else
  dserror("Please verify that incremental update is not needed here!");
#endif

    return;
  }
#endif


  // create a matrix free operator if it does not exist
  if (matfreeoperator_==null)
    matfreeoperator_ = rcp(new NOX::Epetra::MatrixFree(printParams,fineinterface_,*disi_,false));

  // create a linear system if it does not exist
  if (rcpazlinsys_==null)
  {
    NOX::Epetra::Vector initialGuess(*disi_);
    const RefCountPtr<NOX::Epetra::Interface::Jacobian>       ijac  = matfreeoperator_;
    const RefCountPtr<Epetra_Operator>                        jac   = matfreeoperator_;
    const RefCountPtr<NOX::Epetra::Interface::Preconditioner> iprec = prec_;
    const RefCountPtr<Epetra_Operator>                        prec  = prec_;
    rcpazlinsys_ = rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,lsParams,
                                                            ijac,
                                                            jac,
                                                            iprec,
                                                            prec,
                                                            initialGuess));
  }

  // create a group if it does not exist
  NOX::Epetra::Vector initialGuess(*disi_);
  RefCountPtr<NOX::Epetra::Group> rcpgrp = rcp(new NOX::Epetra::Group(printParams,
                                                                      fineinterface_,
                                                                      initialGuess,
                                                                      rcpazlinsys_));

  // create a convergence test if it does not exist
  if (combo_==null)
  {
    RefCountPtr<NOX::StatusTest::NormF> absresid =
      rcp( new NOX::StatusTest::NormF(tolres));
    RefCountPtr<NOX::StatusTest::NormUpdate> nupdate =
      rcp(new NOX::StatusTest::NormUpdate(toldisp));
    RefCountPtr<NOX::StatusTest::Combo> converged =
      rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
    converged->addStatusTest(absresid);
    converged->addStatusTest(nupdate);
    RefCountPtr<NOX::StatusTest::FiniteValue> fv =
      rcp(new NOX::StatusTest::FiniteValue());
    int maxcycle = mlparams.get("nlnML max cycles",200);
    RefCountPtr<NOX::StatusTest::MaxIters> maxiters =
      rcp(new NOX::StatusTest::MaxIters(maxcycle));
    combo_ = rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo_->addStatusTest(maxiters);
    combo_->addStatusTest(converged);
    combo_->addStatusTest(fv);
  }

  // create nox solver manager if it does not exist
  RCP<NOX::Solver::Generic> noxsolver = NOX::Solver::buildSolver(rcpgrp,combo_,rcpparams);
  prec_->SetNoxSolver(noxsolver);

  // solve nonlinear problem
  Epetra_Time timer(discret_.Comm());
  double t0 = timer.ElapsedTime();
  NOX::StatusTest::StatusType status = noxsolver->solve();
  double t1 = timer.ElapsedTime();

  bool converged = true;
  if (status != NOX::StatusTest::Converged) converged = false;

  // get status and print output message
  if (outlevel and myrank_==0)
  {
    printf("NOX/ML :============solve time incl. setup : %15.4f sec\n",t1-t0);
    double appltime = fineinterface_->getsumtime();
    printf("NOX/ML :===========of which time in ccarat : %15.4f sec\n",appltime);
    cout << "NOX/ML :======number calls to computeF in this solve : "
         << fineinterface_->getnumcallscomputeF() << "\n\n\n";
    if (!converged)
      cout << "***WRN***: NOX not converged!\n";
    fflush(stdout);
  }
  fineinterface_->resetsumtime();
  fineinterface_->setnumcallscomputeF(0);

  //---------------------------------- update mid configuration values
  // displacements
  // incremental (disi is now D_{n+1}-D_{n})
  dism_->Update((1.-alphaf),*disi_,1.0);

  // velocities
#ifdef STRUGENALPHA_INCRUPDT
  // incremental (required for constant predictor)
  velm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
  velm_->Update((beta-(1.0-alphaf)*gamma)/beta,*vel_,
                (1.0-alphaf)*(2.*beta-gamma)*dt/(2.*beta),*acc_,
                gamma/(beta*dt));
#else
  dserror("Please verify that incremental update is not needed here!");
#endif

  // accelerations

  // incremental (required for constant predictor)
#ifdef STRUGENALPHA_INCRUPDT
  accm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
  accm_->Update(-(1.-alpham)/(beta*dt),*vel_,
                (2.*beta-1.+alpham)/(2.*beta),*acc_,
                (1.-alpham)/((1.-alphaf)*beta*dt*dt));
#else
  dserror("Please verify that incremental update is not needed here!");
#endif

  return;
} // StruGenAlpha::NonlinearCG()

/*----------------------------------------------------------------------*
 |  do Newton iteration (public)                             mwgee 03/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::PTC()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time      = params_.get<double>("total time"             ,0.0);
  double dt        = params_.get<double>("delta time"             ,0.01);
  double timen     = time + dt;
  int    maxiter   = params_.get<int>   ("max iterations"         ,10);
  bool   damping   = params_.get<bool>  ("damping"                ,false);
  double beta      = params_.get<double>("beta"                   ,0.292);
  double gamma     = params_.get<double>("gamma"                  ,0.581);
  double alpham    = params_.get<double>("alpha m"                ,0.378);
  double alphaf    = params_.get<double>("alpha f"                ,0.459);
  string convcheck = params_.get<string>("convcheck"              ,"AbsRes_Or_AbsDis");
  double toldisp   = params_.get<double>("tolerance displacements",1.0e-07);
  double tolres    = params_.get<double>("tolerance residual"     ,1.0e-07);
  bool printscreen = params_.get<bool>  ("print to screen",true);
  bool printerr    = params_.get<bool>  ("print to err",false);
  FILE* errfile    = params_.get<FILE*> ("err file",NULL);
  if (!errfile) printerr = false;

  // check whether we have a stiffness matrix, that is not filled yet
  // and mass and damping are present
  if (stiff_->Filled()) dserror("stiffness matrix may not be filled here");
  if (!mass_->Filled()) dserror("mass matrix must be filled here");
  if (damping)
    if (!damp_->Filled()) dserror("damping matrix must be filled here");

  // hard wired ptc parameters
  double ptcdt = 1.0e-03;
  double nc;
  fresm_->NormInf(&nc);
  double dti = 1/ptcdt;
  double dti0 = dti;
  RCP<Epetra_Vector> x0 = rcp(new Epetra_Vector(*disi_));

  //=================================================== equilibrium loop
  int numiter=0;
  double fresmnorm = 1.0e6;
  double disinorm = 1.0e6;
  fresm_->Norm2(&fresmnorm);
  Epetra_Time timer(discret_.Comm());
  timer.ResetStartTime();
  bool print_unconv = true;

  while (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) and numiter<=maxiter)
  {
    double dtim = dti0;
    dti0 = dti;
    RCP<Epetra_Vector> xm = rcp(new Epetra_Vector(*x0));
    x0->Update(1.0,*disi_,0.0);
    //------------------------------------------- effective rhs is fresm
    //---------------------------------------------- build effective lhs
    // (using matrix stiff_ as effective matrix)
    stiff_->Add(*mass_,false,(1.-alpham)/(beta*dt*dt),1.-alphaf);
    if (damping)
      stiff_->Add(*damp_,false,(1.-alphaf)*gamma/(beta*dt),1.0);
    stiff_->Complete();

    //------------------------------- do ptc modification to effective LHS
    {
      RCP<Epetra_Vector> tmp = LINALG::CreateVector(stiff_->RowMap(),false);
      tmp->PutScalar(dti);
      RCP<Epetra_Vector> diag = LINALG::CreateVector(stiff_->RowMap(),false);
      stiff_->ExtractDiagonalCopy(*diag);
      diag->Update(1.0,*tmp,1.0);
      stiff_->ReplaceDiagonalValues(*diag);
    }

    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

    //--------------------------------------------------- solve for disi
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (!numiter)
      solver_.Solve(stiff_->EpetraMatrix(),disi_,fresm_,true,true);
    else
      solver_.Solve(stiff_->EpetraMatrix(),disi_,fresm_,true,false);

    //---------------------------------- update mid configuration values
    // displacements
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
    dism_->Update(1.-alphaf,*disi_,1.0);
    // velocities
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // V_{n+1-alpha_f} := V_{n+1-alpha_f}
    //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
    velm_->Update((1.-alphaf)*gamma/(beta*dt),*disi_,1.0);
#else
    // incremental (required for constant predictor)
    velm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
    velm_->Update((beta-(1.0-alphaf)*gamma)/beta,*vel_,
                  (1.0-alphaf)*(2.*beta-gamma)*dt/(2.*beta),*acc_,
                  gamma/(beta*dt));
#endif
    // accelerations
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // A_{n+1-alpha_m} := A_{n+1-alpha_m}
    //                  + (1-alpha_m)/beta/dt^2*IncD_{n+1}
    accm_->Update((1.-alpham)/(beta*dt*dt),*disi_,1.0);
#else
    // incremental (required for constant predictor)
    accm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
    accm_->Update(-(1.-alpham)/(beta*dt),*vel_,
                  (2.*beta-1.+alpham)/(2.*beta),*acc_,
                  (1.-alpham)/((1.-alphaf)*beta*dt*dt));
#endif

    //---------------------------- compute internal forces and stiffness
    {
      // zero out stiffness
      stiff_->Zero();
      // create the parameters for the discretization
      ParameterList p;
      // action for elements
      p.set("action","calc_struct_nlnstiff");
      // choose what to assemble
      p.set("assemble matrix 1",true);
      p.set("assemble matrix 2",false);
      p.set("assemble vector 1",true);
      p.set("assemble vector 2",false);
      p.set("assemble vector 3",false);
      // other parameters that might be needed by the elements
      p.set("total time",timen);
      p.set("delta time",dt);
      // set vector values needed by elements
      discret_.ClearState();
      discret_.SetState("residual displacement",disi_);
      discret_.SetState("displacement",dism_);
      //discret_.SetState("velocity",velm_); // not used at the moment
      fint_->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,stiff_,null,fint_,null,null);
      discret_.ClearState();
      // do NOT finalize the stiffness matrix to add masses to it later
    }

    //------------------------------------------ compute residual forces
    // Res = M . A_{n+1-alpha_m}
    //     + C . V_{n+1-alpha_f}
    //     + F_int(D_{n+1-alpha_f})
    //     - F_{ext;n+1-alpha_f}
    // add inertia mid-forces
    mass_->Multiply(false,*accm_,*finert_);
    fresm_->Update(1.0,*finert_,0.0);
    // add viscous mid-forces
    if (damping)
    {
      //RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,false);
      damp_->Multiply(false,*velm_,*fvisc_);
      fresm_->Update(1.0,*fvisc_,1.0);
    }
    // add static mid-balance
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
    // blank residual DOFs that are on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
    }

    // compute inf norm of residual
    double np;
    fresm_->NormInf(&np);

    //---------------------------------------------- build residual norm
    disi_->Norm2(&disinorm);

    fresm_->Norm2(&fresmnorm);
    // a short message
    if (!myrank_ and (printscreen or printerr))
    {
      PrintPTC(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck,dti);
    }

    //------------------------------------ PTC update of artificial time
#if 1
    // SER step size control
    dti *= (np/nc);
    dti = max(dti,0.0);
    nc = np;
#else
    {
      // TTI step size control
      double ttau=0.75;
      RCP<Epetra_Vector> d1 = LINALG::CreateVector(stiff_->RowMap(),false);
      d1->Update(1.0,*disi_,-1.0,*x0,0.0);
      d1->Scale(dti0);
      RCP<Epetra_Vector> d0 = LINALG::CreateVector(stiff_->RowMap(),false);
      d0->Update(1.0,*x0,-1.0,*xm,0.0);
      d0->Scale(dtim);
      double dt0 = 1/dti0;
      double dtm = 1/dtim;
      RCP<Epetra_Vector> xpp = LINALG::CreateVector(stiff_->RowMap(),false);
      xpp->Update(2.0/(dt0+dtm),*d1,-2.0/(dt0+dtm),*d0,0.0);
      RCP<Epetra_Vector> xtt = LINALG::CreateVector(stiff_->RowMap(),false);
      for (int i=0; i<xtt->MyLength(); ++i) (*xtt)[i] = abs((*xpp)[i])/(1.0+abs((*disi_)[i]));
      double ett;
      xtt->MaxValue(&ett);
      ett = ett / (2.*ttau);
      dti = sqrt(ett);
      nc = np;
    }
#endif

    //--------------------------------- increment equilibrium loop index
    ++numiter;

  }
  //============================================= end equilibrium loop
  print_unconv = false;

  //-------------------------------- test whether max iterations was hit
  if (numiter>=maxiter)
  {
     dserror("PTC unconverged in %d iterations",numiter);
  }
  else
  {
     if (!myrank_ and printscreen)
     {
       PrintPTC(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                   fresmnorm,disinorm,convcheck,dti);
     }
  }

  params_.set<int>("num iterations",numiter);

  return;
} // StruGenAlpha::PTC()


/*----------------------------------------------------------------------*
 |  compute residual to given state x (public)               mwgee 10/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::computeF(const Epetra_Vector& x, Epetra_Vector& F)
{
    //==================== compute residual for given displcaement increment x
    double time    = params_.get<double>("total time",0.0);
    double dt      = params_.get<double>("delta time",0.01);
    double timen   = time + dt;
    double beta    = params_.get<double>("beta"      ,0.292);
    double gamma   = params_.get<double>("gamma"     ,0.581);
    double alpham  = params_.get<double>("alpha m"   ,0.378);
    double alphaf  = params_.get<double>("alpha f"   ,0.459);
    bool   damping = params_.get<bool>  ("damping"   ,false);
    const Epetra_Map* dofrowmap = discret_.DofRowMap();

    // cast away constness of x
    fresm_->PutScalar(0.0);
    Epetra_Vector& dx = const_cast<Epetra_Vector&>(x);
    RefCountPtr<Epetra_Vector> disi = rcp(&dx);
    disi.release();


    // blank increment on Dirichlet BC
    LINALG::ApplyDirichlettoSystem(disi,fresm_,zeros_,dirichtoggle_);

    //---------------------------------- update mid configuration values
    RefCountPtr<Epetra_Vector> dism = LINALG::CreateVector(*dofrowmap,false);
    RefCountPtr<Epetra_Vector> velm = LINALG::CreateVector(*dofrowmap,false);
    RefCountPtr<Epetra_Vector> accm = LINALG::CreateVector(*dofrowmap,false);
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
    dism->Update(1.-alphaf,*disi,1.0,*dism_,0.0);

    // velocities
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // V_{n+1-alpha_f} := V_{n+1-alpha_f}
    //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
    velm->Update((1.-alphaf)*gamma/(beta*dt),*disi,1.0,*velm_,0.0);
#else
    // incremental (required for constant predictor)
    velm->Update(1.0,*dism,-1.0,*dis_,0.0);
    velm->Update((beta-(1.0-alphaf)*gamma)/beta,*vel_,
                  (1.0-alphaf)*(2.*beta-gamma)*dt/(2.*beta),*acc_,
                  gamma/(beta*dt));
#endif

    // accelerations
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // A_{n+1-alpha_m} := A_{n+1-alpha_m}
    //                  + (1-alpha_m)/beta/dt^2*IncD_{n+1}
    accm->Update((1.-alpham)/(beta*dt*dt),*disi,1.0,*accm_,0.0);
#else
    // incremental (required for constant predictor)
    accm->Update(1.0,*dism,-1.0,*dis_,0.0);
    accm->Update(-(1.-alpham)/(beta*dt),*vel_,
                  (2.*beta-1.+alpham)/(2.*beta),*acc_,
                  (1.-alpham)/((1.-alphaf)*beta*dt*dt));
#endif

    //----------------------------------------- compute internal forces
    {
      // create the parameters for the discretization
      ParameterList p;
      // action for elements
      p.set("action","calc_struct_nlnstiff");
      // choose what to assemble
      p.set("assemble matrix 1",false);
      p.set("assemble matrix 2",false);
      p.set("assemble vector 1",true);
      p.set("assemble vector 2",false);
      p.set("assemble vector 3",false);
      // other parameters that might be needed by the elements
      p.set("total time",timen);
      p.set("delta time",dt);
      // set vector values needed by elements
      discret_.ClearState();
      discret_.SetState("residual displacement",disi);
      discret_.SetState("displacement",dism);
      //discret_.SetState("velocity",velm); // not used at the moment
      fint_->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,null,null,fint_,null,null);
      discret_.ClearState();
    }

    //------------------------------------------ compute residual forces
    // Res = M . A_{n+1-alpha_m}
    //     + C . V_{n+1-alpha_f}
    //     + F_int(D_{n+1-alpha_f})
    //     - F_{ext;n+1-alpha_f}
    // add inertia mid-forces
    mass_->Multiply(false,*accm,*finert_);
    fresm_->Update(1.0,*finert_,0.0);
    // add viscous mid-forces
    if (damping)
    {
      //RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,false);
      damp_->Multiply(false,*velm,*fvisc_);
      fresm_->Update(1.0,*fvisc_,1.0);
    }
    // add static mid-balance
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
    // blank residual at DOFs that are on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
    }

    // copy residual forces to nox vector
    F.Update(1.0,*fresm_,0.0);
    // copy disi to time integrator
    disi_->Update(1.0,*disi,0.0);
    //========================================================================

  return;
} // void StruGenAlpha::computeF(const Epetra_Vector& x, Epetra_Vector& F)


/*----------------------------------------------------------------------*
 |  compute residual to given state dism_+x (public)            lw 10/07|
 |  in this routine, none of the state vectors is changed               |
 |  -> necessary in case of matrix-free methods where the residual is   |
 |  evaluated also for pertubations                                     |
 *----------------------------------------------------------------------*/
void StruGenAlpha::computeFmatrixfree(const Epetra_Vector& x, Epetra_Vector& F)
{
    //==================== compute residual for given displcaement increment x
    double time    = params_.get<double>("total time",0.0);
    double dt      = params_.get<double>("delta time",0.01);
    double beta    = params_.get<double>("beta"      ,0.292);
    double gamma   = params_.get<double>("gamma"     ,0.581);
    double alpham  = params_.get<double>("alpha m"   ,0.378);
    double alphaf  = params_.get<double>("alpha f"   ,0.459);
    bool   damping = params_.get<bool>  ("damping"   ,false);
    const Epetra_Map* dofrowmap = discret_.DofRowMap();

    RCP<Epetra_Vector> fresm = rcp(new Epetra_Vector(*dofrowmap, true));
    RCP<Epetra_Vector> fint = rcp(new Epetra_Vector(*dofrowmap, true));

    // cast away constness of x
    Epetra_Vector& dx = const_cast<Epetra_Vector&>(x);
    RefCountPtr<Epetra_Vector> disi = rcp(&dx);
    disi.release();


    // blank increment on Dirichlet BC
    LINALG::ApplyDirichlettoSystem(disi,fresm,zeros_,dirichtoggle_);

    //---------------------------------- update mid configuration values
    RefCountPtr<Epetra_Vector> dism = LINALG::CreateVector(*dofrowmap,false);
    RefCountPtr<Epetra_Vector> velm = LINALG::CreateVector(*dofrowmap,false);
    RefCountPtr<Epetra_Vector> accm = LINALG::CreateVector(*dofrowmap,false);
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
    dism->Update(1.-alphaf,*disi,1.0,*dism_,0.0);


    // velocities
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // V_{n+1-alpha_f} := V_{n+1-alpha_f}
    //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
    velm->Update((1.-alphaf)*gamma/(beta*dt),*disi,1.0,*velm_,0.0);
#else
    // incremental (required for constant predictor)
    velm->Update(1.0,*dism,-1.0,*dis_,0.0);
    velm->Update((beta-(1.0-alphaf)*gamma)/beta,*vel_,
                  (1.0-alphaf)*(2.*beta-gamma)*dt/(2.*beta),*acc_,
                  gamma/(beta*dt));
#endif

    // accelerations
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // A_{n+1-alpha_m} := A_{n+1-alpha_m}
    //                  + (1-alpha_m)/beta/dt^2*IncD_{n+1}
    accm->Update((1.-alpham)/(beta*dt*dt),*disi,1.0,*accm_,0.0);
#else
    // incremental (required for constant predictor)
    accm->Update(1.0,*dism,-1.0,*dis_,0.0);
    accm->Update(-(1.-alpham)/(beta*dt),*vel_,
                  (2.*beta-1.+alpham)/(2.*beta),*acc_,
                  (1.-alpham)/((1.-alphaf)*beta*dt*dt));
#endif

    //----------------------------------------- compute internal forces
    {
      // create the parameters for the discretization
      ParameterList p;
      // action for elements
      p.set("action","calc_struct_nlnstiff");
      // choose what to assemble
      p.set("assemble matrix 1",false);
      p.set("assemble matrix 2",false);
      p.set("assemble vector 1",true);
      p.set("assemble vector 2",false);
      p.set("assemble vector 3",false);
      // other parameters that might be needed by the elements
      p.set("total time",time);
      p.set("delta time",dt);
      // set vector values needed by elements
      discret_.ClearState();
      discret_.SetState("residual displacement",disi);
      discret_.SetState("displacement",dism);
      //discret_.SetState("velocity",velm); // not used at the moment
      fint->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,null,null,fint,null,null);
      discret_.ClearState();
    }

    //------------------------------------------ compute residual forces
    // Res = M . A_{n+1-alpha_m}
    //     + C . V_{n+1-alpha_f}
    //     + F_int(D_{n+1-alpha_f})
    //     - F_{ext;n+1-alpha_f}
    // add inertia mid-forces
    mass_->Multiply(false,*accm,*fresm);
    // add viscous mid-forces
    if (damping)
    {
      RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,false);
      damp_->Multiply(false,*velm,*fviscm);
      fresm->Update(1.0,*fviscm,1.0);
    }
    // add static mid-balance
    fresm->Update(-1.0,*fint,1.0,*fextm_,-1.0);
    // blank residual at DOFs that are on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm);
      fresm->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
    }

    // copy residual forces to nox vector
    F.Update(1.0,*fresm,0.0);
    // copy disi to time integrator
    // disi_->Update(1.0,*disi,0.0);
    //========================================================================

  return;
} // void StruGenAlpha::computeFmatrixfree(const Epetra_Vector& x, Epetra_Vector& F)


/*----------------------------------------------------------------------*
 |  compute Jacobian to given state x (public)               mwgee 10/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::computeJacobian(const Epetra_Vector& x)
{
    double time    = params_.get<double>("total time",0.0);
    double dt      = params_.get<double>("delta time",0.01);
    double timen   = time + dt;
    double beta    = params_.get<double>("beta"      ,0.292);
    double gamma   = params_.get<double>("gamma"     ,0.581);
    double alpham  = params_.get<double>("alpha m"   ,0.378);
    double alphaf  = params_.get<double>("alpha f"   ,0.459);
    bool   damping = params_.get<bool>  ("damping"   ,false);
    const Epetra_Map* dofrowmap = discret_.DofRowMap();

    // cast away constness of x
    Epetra_Vector& dx = const_cast<Epetra_Vector&>(x);
    RefCountPtr<Epetra_Vector> disi = rcp(&dx);
    disi.release();

    //---------------------------------- update mid configuration values
    RefCountPtr<Epetra_Vector> dism = LINALG::CreateVector(*dofrowmap,false);
    RefCountPtr<Epetra_Vector> velm = LINALG::CreateVector(*dofrowmap,false);
    RefCountPtr<Epetra_Vector> accm = LINALG::CreateVector(*dofrowmap,false);
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
    dism->Update(1.-alphaf,*disi,1.0,*dism_,0.0);
    // V_{n+1-alpha_f} := V_{n+1-alpha_f}
    //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
    velm->Update((1.-alphaf)*gamma/(beta*dt),*disi,1.0,*velm_,0.0);
    // A_{n+1-alpha_m} := A_{n+1-alpha_m}
    //                  + (1-alpha_m)/beta/dt^2*IncD_{n+1}
    accm->Update((1.-alpham)/(beta*dt*dt),*disi,1.0,*accm_,0.0);

    //------------------------------------------------ compute stiffness
    {
      // zero out stiffness
      stiff_->Zero();
      // create the parameters for the discretization
      ParameterList p;
      // action for elements
      p.set("action","calc_struct_nlnstiff");
      // choose what to assemble
      p.set("assemble matrix 1",true);
      p.set("assemble matrix 2",false);
      p.set("assemble vector 1",false);
      p.set("assemble vector 2",false);
      p.set("assemble vector 3",false);
      // other parameters that might be needed by the elements
      p.set("total time",timen);
      p.set("delta time",dt);
      // set vector values needed by elements
      discret_.ClearState();
      discret_.SetState("residual displacement",disi);
      discret_.SetState("displacement",dism);
      fint_->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,stiff_,null,null,null,null);
      discret_.ClearState();
      // do NOT finalize the stiffness matrix to add masses to it later
    }

    //---------------------------------------------- build effective lhs
    // (using matrix stiff_ as effective matrix)
    stiff_->Add(*mass_,false,(1.-alpham)/(beta*dt*dt),1.-alphaf);
    if (damping)
      stiff_->Add(*damp_,false,(1.-alphaf)*gamma/(beta*dt),1.0);
    stiff_->Complete();
    LINALG::ApplyDirichlettoSystem(stiff_,
                                   disi,
                                   fresm_,
                                   zeros_,
                                   dirichtoggle_);
  return;
} // StruGenAlpha::computeJacobian

/*----------------------------------------------------------------------*
 |  do update and output (public)                            mwgee 03/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::UpdateandOutput()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time          = params_.get<double>("total time"             ,0.0);
  double dt            = params_.get<double>("delta time"             ,0.01);
  double timen         = time + dt;  // t_{n+1}
  int    step          = params_.get<int>   ("step"                   ,0);
  int    istep         = step + 1;  // n+1
  int    nstep         = params_.get<int>   ("nstep"                  ,5);
  int    numiter       = params_.get<int>   ("num iterations"         ,-1);

  double alpham        = params_.get<double>("alpha m"                ,0.378);
  double alphaf        = params_.get<double>("alpha f"                ,0.459);

  bool   iodisp        = params_.get<bool>  ("io structural disp"     ,true);
  int    updevrydisp   = params_.get<int>   ("io disp every nstep"    ,10);
  bool   iostress      = params_.get<bool>  ("io structural stress"   ,false);
  int    updevrystress = params_.get<int>   ("io stress every nstep"  ,10);

  int    writeresevry  = params_.get<int>   ("write restart every"    ,0);

  bool   printscreen   = params_.get<bool>  ("print to screen"        ,true);
  bool   printerr      = params_.get<bool>  ("print to err"           ,true);
  FILE*  errfile       = params_.get<FILE*> ("err file"               ,NULL);
  if (!errfile) printerr = false;

  //----------------------------------------------- update time and step
  params_.set<double>("total time", timen);
  params_.set<int>("step", istep);

  //---------------------------- determine new end-quantities and update
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1} = 1./(1.-alphaf) * D_{n+1-alpha_f}
  //                     - alphaf/(1.-alphaf) * D_n
  dis_->Update(1./(1.-alphaf),*dism_,-alphaf/(1.-alphaf));
  // new velocities at t_{n+1} -> t_n
  //    V_{n} := V_{n+1} = 1./(1.-alphaf) * V_{n+1-alpha_f}
  //                     - alphaf/(1.-alphaf) * V_n
  vel_->Update(1./(1.-alphaf),*velm_,-alphaf/(1.-alphaf));
  // new accelerations at t_{n+1} -> t_n
  //    A_{n} := A_{n+1} = 1./(1.-alpham) * A_{n+1-alpha_m}
  //                     - alpham/(1.-alpham) * A_n
  acc_->Update(1./(1.-alpham),*accm_,-alpham/(1.-alpham));
  // update new external force
  //    F_{ext;n} := F_{ext;n+1}
  fext_->Update(1.0,*fextn_,0.0);

  //----- update anything that needs to be updated at the element level
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_update_istep");
    // other parameters that might be needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    discret_.Evaluate(p,null,null,null,null,null);
  }

  bool isdatawritten = false;

  //------------------------------------------------- write restart step
  if (writeresevry and istep%writeresevry==0)
  {
    output_.WriteMesh(istep,timen);
    output_.NewStep(istep, timen);
    output_.WriteVector("displacement",dis_);
    output_.WriteVector("velocity",vel_);
    output_.WriteVector("acceleration",acc_);
    output_.WriteVector("fexternal",fext_);
    isdatawritten = true;

    if (discret_.Comm().MyPID()==0 and printscreen)
    {
      cout << "====== Restart written in step " << istep << endl;
      fflush(stdout);
    }
    if (errfile and printerr)
    {
      fprintf(errfile,"====== Restart written in step %d\n",istep);
      fflush(errfile);
    }
  }

  //----------------------------------------------------- output results
  if (iodisp and updevrydisp and istep%updevrydisp==0 and !isdatawritten)
  {
    output_.NewStep(istep, timen);
    output_.WriteVector("displacement",dis_);
    output_.WriteVector("velocity",vel_);
    output_.WriteVector("acceleration",acc_);
    output_.WriteElementData();
    isdatawritten = true;
  }

  //------------------------------------- do stress calculation and output
  if (updevrystress and !(istep%updevrystress) and iostress)
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_stress");
    // other parameters that might be needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    Teuchos::RCP<std::vector<char> > stress = Teuchos::rcp(new std::vector<char>());
    p.set("stress", stress);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("residual displacement",zeros_);
    discret_.SetState("displacement",dis_);
    discret_.Evaluate(p,null,null,null,null,null);
    discret_.ClearState();
    if (!isdatawritten) output_.NewStep(istep, timen);
    isdatawritten = true;
    output_.WriteVector("gauss_stresses_xyz",*stress,*discret_.ElementColMap());
  }

  //---------------------------------------------------------- print out
  if (!myrank_)
  {
    if (printscreen)
    {
      printf("step %6d | nstep %6d | time %-14.8E | dt %-14.8E | numiter %3d\n",
             istep,nstep,timen,dt,numiter);
      printf("----------------------------------------------------------------------------------\n");
      fflush(stdout);
    }
    if (printerr)
    {
      fprintf(errfile,"step %6d | nstep %6d | time %-14.8E | dt %-14.8E | numiter %3d\n",
              istep,nstep,timen,dt,numiter);
      fprintf(errfile,"----------------------------------------------------------------------------------\n");
      fflush(errfile);
    }
  }

  return;
} // StruGenAlpha::UpdateandOutput()


/*----------------------------------------------------------------------*
 |  determine new state at t_{n+1} (public)                bborn 11/07  |
 *----------------------------------------------------------------------*/
void StruGenAlpha::ExtrapolateEndState()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double alpham        = params_.get<double>("alpha m"                ,0.378);
  double alphaf        = params_.get<double>("alpha f"                ,0.459);

  //---------------------------- determine new end-quantities and update
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1} = 1./(1.-alphaf) * D_{n+1-alpha_f}
  //                     - alphaf/(1.-alphaf) * D_n
  disn_->Update(1./(1.-alphaf),*dism_,-alphaf/(1.-alphaf),*dis_,0.0);
  // new velocities at t_{n+1} -> t_n
  //    V_{n} := V_{n+1} = 1./(1.-alphaf) * V_{n+1-alpha_f}
  //                     - alphaf/(1.-alphaf) * V_n
  veln_->Update(1./(1.-alphaf),*velm_,-alphaf/(1.-alphaf),*vel_,0.0);
  // new accelerations at t_{n+1} -> t_n
  //    A_{n} := A_{n+1} = 1./(1.-alpham) * A_{n+1-alpha_m}
  //                     - alpham/(1.-alpham) * A_n
  accn_->Update(1./(1.-alpham),*accm_,-alpham/(1.-alpham),*acc_,0.0);
  // update new external force
  //    F_{ext;n} := F_{ext;n+1}
  //fext_->Update(1.0,*fextn_,0.0);

  return;
} // StruGenAlpha::ExtrapolateEndState


/*----------------------------------------------------------------------*
 |  integrate in time          (static/public)               mwgee 03/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::Integrate()
{
  int    step  = params_.get<int>   ("step" ,0);
  int    nstep = params_.get<int>   ("nstep",5);

  // can have values "full newton" , "modified newton" , "nonlinear cg"
  string equil = params_.get<string>("equilibrium iteration","full newton");

  // can have values takes values "constant" consistent"
  string pred  = params_.get<string>("predictor","constant");
  int predictor=-1;
  if      (pred=="constant")   predictor = 1;
  else if (pred=="consistent") predictor = 2;
  else dserror("Unknown type of predictor");

  //in case a constraint is defined, use defined algorithm
  if (ConstrMan_->HaveConstraint())
  {
	  string algo = params_.get<string>("uzawa algorithm","newtonlinuzawa");
	  ConstrMan_->ScaleLagrMult(0.0);
	  for (int i=step; i<nstep; ++i)
      {
        if      (predictor==1) ConstantPredictor();
        else if (predictor==2) ConsistentPredictor();
        //Does predicted displacement satisfy constraint?
        double time          = params_.get<double>("total time"             ,0.0);
        double dt            = params_.get<double>("delta time"             ,0.01);
        // what algorithm is used?
        // - "newtonlinuzawa": 			Potential is linearized wrt displacements and Lagrange multipliers
        //					   			Linear problem is solved with Uzawa algorithm
        // - "augmentedlagrange":		Potential is linearized wrt displacements keeping Lagrange multiplier fixed
        //								Until convergence Lagrange multiplier increased by Uzawa_param*(Vol_err)
        if (algo=="newtonlinuzawa")
        {
        	ConstrMan_->ScaleLagrMult(0.0);
        	ConstrMan_->StiffnessAndInternalForces(time+dt,disn_,fint_,stiff_);
        	FullNewtonLinearUzawa();
        }
        else if (algo=="augmentedlagrange")
        {
        	ConstrMan_->StiffnessAndInternalForces(time+dt,disn_,fint_,stiff_);
        	NonLinearUzawaFullNewton(predictor);
        }
        else dserror("Unknown type of algorithm to deal with constraints");
        UpdateandOutput();
      }
  }
  else if (equil=="full newton")
  {
    for (int i=step; i<nstep; ++i)
    {
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      FullNewton();
      UpdateandOutput();
    }
  }
  else if (equil=="modified newton")
  {
    for (int i=step; i<nstep; ++i)
    {
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      ModifiedNewton();
      UpdateandOutput();
    }
  }
  else if (equil=="matrixfree newton")
  {
    for (int i=step; i<nstep; ++i)
    {
      //dserror("Matrix-free Newton not yet implemented");
      MatrixFreeConstantPredictor();
      MatrixFreeNewton();
      UpdateandOutput();
    }
  }
  else if (equil=="nonlinear cg")
  {
    for (int i=step; i<nstep; ++i)
    {
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      NonlinearCG();
      UpdateandOutput();
    }
  }
  else if (equil=="ptc")
  {
    for (int i=step; i<nstep; ++i)
    {
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      PTC();
      UpdateandOutput();
    }
  }
  else dserror("Unknown type of equilibrium iteration");

  return;
} // void StruGenAlpha::Integrate()


/*----------------------------------------------------------------------*
 |  integrate one step in time (static/public)               bborn 10/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::IntegrateStep()
{
   //int    step = params_.get<int>   ("step" ,0);

  // can have values "full newton" , "modified newton" , "nonlinear cg"
  string equil = params_.get<string>("equilibrium iteration","full newton");

  // can have values takes values "constant" consistent"
  string pred  = params_.get<string>("predictor","constant");
  int predictor=-1;
  if      (pred=="constant")   predictor = 1;
  else if (pred=="consistent") predictor = 2;
  else dserror("Unknown type of predictor");

  if (equil=="full newton")
  {
    if      (predictor==1) ConstantPredictor();
    else if (predictor==2) ConsistentPredictor();
    FullNewton();
  }
  else if (equil=="modified newton")
  {
    if      (predictor==1) ConstantPredictor();
    else if (predictor==2) ConsistentPredictor();
    ModifiedNewton();
  }
  else if (equil=="nonlinear cg")
  {
    if      (predictor==1) ConstantPredictor();
    else if (predictor==2) ConsistentPredictor();
    NonlinearCG();
  }
  else if (equil=="ptc")
  {
    if      (predictor==1) ConstantPredictor();
    else if (predictor==2) ConsistentPredictor();
    PTC();
  }
  else if (equil=="matrixfree newton")
  {
    //dserror("Matfree Newton not yet impl.");
    MatrixFreeConstantPredictor();
    MatrixFreeNewton();
  }
  else dserror("Unknown type of equilibrium iteration");

  return;
} // void StruGenAlpha::IntegrateStep()


/*----------------------------------------------------------------------*
 |  set default parameter list (static/public)               mwgee 03/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::SetDefaults(ParameterList& params)
{
  params.set<bool>  ("print to screen"        ,true);
  params.set<bool>  ("print to err"           ,false);
  params.set<FILE*> ("err file"               ,NULL);
  params.set<bool>  ("damping"                ,false);
  params.set<double>("damping factor K"       ,0.00001);
  params.set<double>("damping factor M"       ,0.00001);
  params.set<double>("beta"                   ,0.292);
  params.set<double>("gamma"                  ,0.581);
  params.set<double>("alpha m"                ,0.378);
  params.set<double>("alpha f"                ,0.459);
  params.set<double>("total time"             ,0.0);
  params.set<double>("delta time"             ,0.01);
  params.set<int>   ("step"                   ,0);
  params.set<int>   ("nstep"                  ,5);
  params.set<int>   ("max iterations"         ,10);
  params.set<int>   ("num iterations"         ,-1);
  params.set<double>("tolerance displacements",1.0e-07);
  params.set<bool>  ("io structural disp"     ,false);
  params.set<int>   ("io disp every nstep"    ,10);
  params.set<bool>  ("io structural stress"   ,false);
  params.set<int>   ("io disp every nstep"    ,10);
  params.set<int>   ("restart"                ,0);
  params.set<int>   ("write restart every"    ,0);
  params.set<bool>  ("contact"                ,false);
  // takes values "constant" consistent"
  params.set<string>("predictor"              ,"constant");
  // takes values "full newton" , "modified newton" , "matrixfree newton", "nonlinear cg" "ptc"
  params.set<string>("equilibrium iteration"  ,"full newton");
  return;
}

/*----------------------------------------------------------------------*
 |  read restart (public)                                    mwgee 06/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::ReadRestart(int step)
{
  RefCountPtr<DRT::Discretization> rcpdiscret = rcp(&discret_);
  rcpdiscret.release();
  IO::DiscretizationReader reader(rcpdiscret,step);
  double time  = reader.ReadDouble("time");
  int    rstep = reader.ReadInt("step");
  if (rstep != step) dserror("Time step on file not equal to given step");

  reader.ReadVector(dis_, "displacement");
  reader.ReadVector(vel_, "velocity");
  reader.ReadVector(acc_, "acceleration");
  reader.ReadVector(fext_,"fexternal");
  reader.ReadMesh(step);

  // iverride current time and step with values from file
  params_.set<double>("total time",time);
  params_.set<int>   ("step",rstep);

  return;
}

/*----------------------------------------------------------------------*
 |  get time integration scheme (TIS) parameters           bborn 10/07  |
 *----------------------------------------------------------------------*/
void StruGenAlpha::GetTISPara(double& beta,
                              double& gamma,
                              double& alpham,
                              double& alphaf)
{
  beta = params_.get<double>("beta", 0.292);
  gamma = params_.get<double>("gamma", 0.581);
  alpham = params_.get<double>("alpha m", 0.378);
  alphaf = params_.get<double>("alpha f", 0.459);
  return;
}

/*----------------------------------------------------------------------*
 |  set time step size                                     bborn 10/07  |
 *----------------------------------------------------------------------*/
void StruGenAlpha::SetTimeStepSize(const double& timstpsiz)
{
  params_.set<double>("delta time", timstpsiz);
  return;
}

/*----------------------------------------------------------------------*
 |  set time t_{n+1}                                       bborn 11/07  |
 *----------------------------------------------------------------------*/
void StruGenAlpha::SetTime(const double& tim)
{
  params_.set<double>("total time", tim);
  return;
}

/*----------------------------------------------------------------------*
 |  set step n+1                                           bborn 11/07  |
 *----------------------------------------------------------------------*/
void StruGenAlpha::SetTimeStep(const int& stp)
{
  params_.set<int>("step", stp);
  return;
}

/*----------------------------------------------------------------------*
 |  check convergence of Newton iteration (public)              lw 12/07|
 *----------------------------------------------------------------------*/
bool StruGenAlpha::Converged(const string type, const double disinorm,
                             const double resnorm, const double toldisp,
                             const double tolres)
{
  if (type == "AbsRes_Or_AbsDis")
  {
    return (disinorm<toldisp or resnorm<tolres);
  }
  else if (type == "AbsRes_And_AbsDis")
  {
    return (disinorm<toldisp and resnorm<tolres);
  }
  else if (type == "RelRes_Or_AbsDis")
  {
    if (ref_fnorm_ == 0.) ref_fnorm_ = 1.0;
    return (disinorm<toldisp or (resnorm/ref_fnorm_)<tolres);
  }
  else if (type == "RelRes_And_AbsDis")
  {
    if (ref_fnorm_ == 0.) ref_fnorm_ = 1.0;
    return (disinorm<toldisp and (resnorm/ref_fnorm_)<tolres);
  }
  else if (type == "RelRes_Or_RelDis")
  {
    if (ref_fnorm_ == 0.) ref_fnorm_ = 1.0;
    if (ref_disnorm_ == 0.) ref_disnorm_ = 1.0;
    return ((disinorm/ref_disnorm_)<toldisp or (resnorm/ref_fnorm_)<tolres);
  }
  else if (type == "RelRes_And_RelDis")
  {
    if (ref_fnorm_ == 0.) ref_fnorm_ = 1.0;
    if (ref_disnorm_ == 0.) ref_disnorm_ = 1.0;
    return ((disinorm/ref_disnorm_)<toldisp and (resnorm/ref_fnorm_)<tolres);
  }
  else
  {
    dserror("Requested convergence check not (yet) implemented");
    return true;
  }
}

/*----------------------------------------------------------------------*
 |  check convergence of Newton iteration (public)              tk 01/08|
 |  take the constraints into account as well                     |
 *----------------------------------------------------------------------*/
bool StruGenAlpha::Converged(const string type, const double disinorm,
        const double resnorm, const double constrnorm,
        const double toldisp, const double tolres,
        const double tolconstr)
{
	return (Converged(type,disinorm, resnorm, toldisp,tolres) and (constrnorm<tolconstr));
}

/*----------------------------------------------------------------------*
 |  calculate reference norms for relative convergence checks   lw 12/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::CalcRefNorms()
{
  // The reference norms are used to scale the calculated iterative
  // displacement norm and/or the residual force norm. For this
  // purpose we only need the right order of magnitude, so we don't
  // mind evaluating the corresponding norms at possibly different
  // points within the timestep (end point, generalized midpoint).

  dis_->Norm2(&ref_disnorm_);

  bool damping   = params_.get<bool>  ("damping",false);

  double fintnorm, fextnorm, finertnorm;
  fint_->Norm2(&fintnorm);  // norm of the internal forces
  fextm_->Norm2(&fextnorm);  // norm of the external forces
  finert_->Norm2(&finertnorm);  // norm of the inertial forces

  if (damping)
  {
    double fviscnorm;
    fvisc_->Norm2(&fviscnorm);
    ref_fnorm_=max(fviscnorm, max(finertnorm, max(fintnorm, fextnorm)));
  }
  else
  {
    ref_fnorm_=max(finertnorm, max(fintnorm, fextnorm));
  }
}

/*----------------------------------------------------------------------*
 |  print to screen and/or error file                           lw 12/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::PrintNewton(bool printscreen, bool printerr, bool print_unconv,
                               FILE* errfile, Epetra_Time timer, int numiter,
                               int maxiter, double fresmnorm, double disinorm,
                               string convcheck)
{
  bool relres        = (convcheck == "RelRes_And_AbsDis" or convcheck == "RelRes_Or_AbsDis");
  bool relres_reldis = (convcheck == "RelRes_And_RelDis" or convcheck == "RelRes_Or_RelDis");

  if (relres)
  {
    fresmnorm /= ref_fnorm_;
  }
  if (relres_reldis)
  {
    fresmnorm /= ref_fnorm_;
    disinorm  /= ref_disnorm_;
  }

  if (print_unconv)
  {
    if (printscreen)
    {
      if (relres)
      {
        printf("numiter %2d scaled res-norm %10.5e absolute dis-norm %20.15E\n",numiter+1, fresmnorm, disinorm);
        fflush(stdout);
      }
      else if (relres_reldis)
      {
        printf("numiter %2d scaled res-norm %10.5e scaled dis-norm %20.15E\n",numiter+1, fresmnorm, disinorm);
        fflush(stdout);
      }
      else
        {
        printf("numiter %2d absolute res-norm %10.5e absolute dis-norm %20.15E\n",numiter+1, fresmnorm, disinorm);
        fflush(stdout);
      }
    }
    if (printerr)
    {
      if (relres)
      {
        fprintf(errfile, "numiter %2d scaled res-norm %10.5e absolute dis-norm %20.15E\n",numiter+1, fresmnorm, disinorm);
        fflush(errfile);
      }
      else if (relres_reldis)
      {
        fprintf(errfile, "numiter %2d scaled res-norm %10.5e scaled dis-norm %20.15E\n",numiter+1, fresmnorm, disinorm);
        fflush(errfile);
      }
      else
        {
        fprintf(errfile, "numiter %2d absolute res-norm %10.5e absolute dis-norm %20.15E\n",numiter+1, fresmnorm, disinorm);
        fflush(errfile);
      }
    }
  }
  else
  {
    if (ConstrMan_->HaveMonitor())
    {
      ConstrMan_->PrintMonitorValues();
    }
    double timepernlnsolve = timer.ElapsedTime();

    if (relres)
    {
      printf("Newton iteration converged: numiter %d scaled res-norm %e absolute dis-norm %e time %10.5f\n",
             numiter,fresmnorm,disinorm,timepernlnsolve);
      fflush(stdout);
    }
    else if (relres_reldis)
    {
      printf("Newton iteration converged: numiter %d scaled res-norm %e scaled dis-norm %e time %10.5f\n",
             numiter,fresmnorm,disinorm,timepernlnsolve);
      fflush(stdout);
    }
    else
    {
      printf("Newton iteration converged: numiter %d absolute res-norm %e absolute dis-norm %e time %10.5f\n",
             numiter,fresmnorm,disinorm,timepernlnsolve);
      fflush(stdout);
    }
  }
}

/*------------------------------------------------------------------------------*
 |  print to screen and/or error file considering constraints           tk 01/08|
 *------------------------------------------------------------------------------*/
void StruGenAlpha::PrintNewton(bool printscreen, bool printerr, bool print_unconv,
                               FILE* errfile, Epetra_Time timer, int numiter,
                               int maxiter, double fresmnorm, double disinorm,
                               string convcheck, double constrnorm, double UzawaPara)
{
  bool relres        = (convcheck == "RelRes_And_AbsDis" or convcheck == "RelRes_Or_AbsDis");
  bool relres_reldis = (convcheck == "RelRes_And_RelDis" or convcheck == "RelRes_Or_RelDis");


  if (relres)
  {
    fresmnorm /= ref_fnorm_;
  }
  if (relres_reldis)
  {
    fresmnorm /= ref_fnorm_;
    disinorm  /= ref_disnorm_;
  }

  if (print_unconv)
  {
    if (printscreen)
    {

      if (relres)
      {
        printf("numiter %2d scaled res-norm %10.5e absolute dis-norm %20.15E absolute constr-norm %10.5e current Uzawa parameter %10.5e\n",
                numiter+1, fresmnorm, disinorm, constrnorm, UzawaPara);
        fflush(stdout);
      }
      else if (relres_reldis)
      {
        printf("numiter %2d scaled res-norm %10.5e scaled dis-norm %20.15E absolute constr_norm %10.5e current Uzawa parameter %10.5e\n",
                numiter+1, fresmnorm, disinorm, constrnorm, UzawaPara);
        fflush(stdout);
      }
      else
        {
        printf("numiter %2d absolute res-norm %10.5e absolute dis-norm %20.15E absolute constr_norm %10.5e current Uzawa parameter %10.5e\n",
                numiter+1, fresmnorm, disinorm, constrnorm, UzawaPara);
        fflush(stdout);
      }
    }
    if (printerr)
    {
      if (relres)
      {
        fprintf(errfile, "numiter %2d scaled res-norm %10.5e absolute dis-norm %20.15E absolute constr_norm %10.5e current Uzawa parameter %10.5e\n",
                numiter+1, fresmnorm, disinorm, constrnorm, UzawaPara);
        fflush(errfile);
      }
      else if (relres_reldis)
      {
        fprintf(errfile, "numiter %2d scaled res-norm %10.5e scaled dis-norm %20.15E absolute constr_norm %10.5e current Uzawa parameter %10.5e\n",
                numiter+1, fresmnorm, disinorm, constrnorm, UzawaPara);
        fflush(errfile);
      }
      else
        {
        fprintf(errfile, "numiter %2d absolute res-norm %10.5e absolute dis-norm %20.15E absolute constr_norm %10.5e current Uzawa parameter %10.5e\n",
                numiter+1, fresmnorm, disinorm, constrnorm, UzawaPara);
        fflush(errfile);
      }
    }
  }
  else
  {
    if (ConstrMan_->HaveMonitor())
    {
      ConstrMan_->PrintMonitorValues();
    }
    double timepernlnsolve = timer.ElapsedTime();

    if (relres)
    {
      printf("Newton iteration converged: numiter %d scaled res-norm %e absolute dis-norm %e absolute constr_norm %e time %10.5f\n",
             numiter,fresmnorm,disinorm, constrnorm,timepernlnsolve);
      fflush(stdout);
    }
    else if (relres_reldis)
    {
      printf("Newton iteration converged: numiter %d scaled res-norm %e scaled dis-norm %e absolute constr_norm %e time %10.5f\n",
             numiter,fresmnorm,disinorm, constrnorm,timepernlnsolve);
      fflush(stdout);
    }
    else
    {
      printf("Newton iteration converged: numiter %d absolute res-norm %e absolute dis-norm %e absolute constr_norm %e time %10.5f\n",
             numiter,fresmnorm,disinorm, constrnorm,timepernlnsolve);
      fflush(stdout);
    }
  }
}

/*----------------------------------------------------------------------*
 |  print to screen and/or error file                          gee 01/08|
 *----------------------------------------------------------------------*/
void StruGenAlpha::PrintPTC(bool printscreen, bool printerr, bool print_unconv,
                               FILE* errfile, Epetra_Time timer, int numiter,
                               int maxiter, double fresmnorm, double disinorm,
                               string convcheck, double dti)
{
  bool relres        = (convcheck == "RelRes_And_AbsDis" or convcheck == "RelRes_Or_AbsDis");
  bool relres_reldis = (convcheck == "RelRes_And_RelDis" or convcheck == "RelRes_Or_RelDis");

  if (relres)
  {
    fresmnorm /= ref_fnorm_;
  }
  if (relres_reldis)
  {
    fresmnorm /= ref_fnorm_;
    disinorm  /= ref_disnorm_;
  }

  if (print_unconv)
  {
    if (printscreen)
    {
      if (relres)
      {
        printf("numiter %2d scaled res-norm %10.5e absolute dis-norm %20.15E dti %20.15E\n",numiter+1, fresmnorm, disinorm,dti);
        fflush(stdout);
      }
      else if (relres_reldis)
      {
        printf("numiter %2d scaled res-norm %10.5e scaled dis-norm %20.15E dti %20.15E\n",numiter+1, fresmnorm, disinorm,dti);
        fflush(stdout);
      }
      else
        {
        printf("numiter %2d absolute res-norm %10.5e absolute dis-norm %20.15E dti %20.15E\n",numiter+1, fresmnorm, disinorm,dti);
        fflush(stdout);
      }
    }
    if (printerr)
    {
      if (relres)
      {
        fprintf(errfile, "numiter %2d scaled res-norm %10.5e absolute dis-norm %20.15E dti %20.15E\n",numiter+1, fresmnorm, disinorm,dti);
        fflush(errfile);
      }
      else if (relres_reldis)
      {
        fprintf(errfile, "numiter %2d scaled res-norm %10.5e scaled dis-norm %20.15E dti %20.15E\n",numiter+1, fresmnorm, disinorm,dti);
        fflush(errfile);
      }
      else
        {
        fprintf(errfile, "numiter %2d absolute res-norm %10.5e absolute dis-norm %20.15E dti %20.15E\n",numiter+1, fresmnorm, disinorm,dti);
        fflush(errfile);
      }
    }
  }
  else
  {

	double timepernlnsolve = timer.ElapsedTime();

    if (relres)
    {
      printf("Psitc iteration converged: numiter %d scaled res-norm %e absolute dis-norm %e time %10.5f\n",
             numiter,fresmnorm,disinorm,timepernlnsolve);
      fflush(stdout);
    }
    else if (relres_reldis)
    {
      printf("Psitc iteration converged: numiter %d scaled res-norm %e scaled dis-norm %e time %10.5f\n",
             numiter,fresmnorm,disinorm,timepernlnsolve);
      fflush(stdout);
    }
    else
    {
      printf("Psitc iteration converged: numiter %d absolute res-norm %e absolute dis-norm %e time %10.5f\n",
             numiter,fresmnorm,disinorm,timepernlnsolve);
      fflush(stdout);
    }
  }
} // end of StruGenAlpha::PrintPTC


/*----------------------------------------------------------------------*
 |  print to screen                                             lw 12/07|
 *----------------------------------------------------------------------*/
void StruGenAlpha::PrintPredictor(string convcheck, double fresmnorm)
{
  if (convcheck != "AbsRes_And_AbsDis" and convcheck != "AbsRes_Or_AbsDis")
  {
    fresmnorm /= ref_fnorm_;
    cout << "Predictor scaled res-norm " << fresmnorm << endl;
  }
  else
  {
    cout << "Predictor absolute res-norm " << fresmnorm << endl;
  }
  fflush(stdout);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> StruGenAlpha::LinearRelaxationSolve(Teuchos::RCP<Epetra_Vector> relax)
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time      = params_.get<double>("total time",0.0);
  double dt        = params_.get<double>("delta time"             ,0.01);
  double timen     = time + dt;  // t_{n+1}
  const Epetra_Map* dofrowmap = discret_.DofRowMap();
  bool   damping   = params_.get<bool>  ("damping"                ,false);
  double beta      = params_.get<double>("beta"                   ,0.292);
  double gamma     = params_.get<double>("gamma"                  ,0.581);
  double alpham    = params_.get<double>("alpha m"                ,0.378);
  double alphaf    = params_.get<double>("alpha f"                ,0.459);

  // we start from zero
  fextm_->Update(1.-alphaf,*relax,0.0);

  // This (re)creates the stiffness matrix at the current
  // configuration.
  //------------- eval fint at interpolated state, eval stiffness matrix
  {
    // zero out stiffness
    stiff_->Zero();
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_nlnstiff");
    // other parameters that might be needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("residual displacement",disi_);
    discret_.SetState("displacement",dism_);
    //discret_.SetState("velocity",velm_); // not used at the moment
    fint_->PutScalar(0.0);  // initialise internal force vector
    discret_.Evaluate(p,stiff_,fint_);
    discret_.ClearState();
    // do NOT finalize the stiffness matrix, add mass and damping to it later
  }

  //-------------------------------------------- compute residual forces
  // Res = M . A_{n+1-alpha_m}
  //     + C . V_{n+1-alpha_f}
  //     + F_int(D_{n+1-alpha_f})
  //     - F_{ext;n+1-alpha_f}
  // add mid-inertial force
  mass_->Multiply(false,*accm_,*fresm_);
  // add mid-viscous damping force
  if (damping)
  {
      RCP<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,true);
      damp_->Multiply(false,*velm_,*fviscm);
      fresm_->Update(1.0,*fviscm,1.0);
  }

  // add static mid-balance
  fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);

  // blank residual at DOFs on Dirichlet BC
  Epetra_Vector fresmcopy(*fresm_);
  fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);

  //------------------------------------------- effective rhs is fresm
  //---------------------------------------------- build effective lhs
  // (using matrix stiff_ as effective matrix)
  stiff_->Add(*mass_,false,(1.-alpham)/(beta*dt*dt),1.-alphaf);
  if (damping)
    stiff_->Add(*damp_,false,(1.-alphaf)*gamma/(beta*dt),1.0);
  stiff_->Complete();

  //----------------------- apply dirichlet BCs to system of equations
  disi_->PutScalar(0.0);  // Useful? depends on solver and more
  LINALG::ApplyDirichlettoSystem(stiff_,disi_,fextm_,zeros_,dirichtoggle_);

  //--------------------------------------------------- solve for disi
  // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
  solver_.Solve(stiff_->EpetraMatrix(),disi_,fextm_,true,true);

  return disi_;
}


#endif  // #ifdef CCADISCRET
