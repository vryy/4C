/*----------------------------------------------------------------------*/
/*!
\file strutimint_ab2.cpp
\brief Structural time integration with Adams-Bashforth 2nd order

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "strutimint_ab2.H"

/*----------------------------------------------------------------------*/
/* Constructor */
StruTimIntAB2::StruTimIntAB2
(
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& sdynparams,
  const Teuchos::ParameterList& xparams,
//  const Teuchos::ParameterList& ab2params,
  DRT::Discretization& actis,
  LINALG::Solver& solver,
  IO::DiscretizationWriter& output
)
: StruTimIntExpl
  (
    ioparams,
    sdynparams,
    xparams,
    actis,
    solver,
    output
  ),
  diso_(Teuchos::null),
  velo_(Teuchos::null),
  acco_(Teuchos::null),
  fextn_(Teuchos::null),
  fintn_(Teuchos::null),
  fviscn_(Teuchos::null),
  frimpn_(Teuchos::null)
{
  // info to user : AB2 --- your federal highway "Warschauer Allee"
  if (myrank_ == 0)
  {
    std::cout << "with Adams-Bashforth 2nd order"
              << std::endl;
  }

  // resize state vector
  state_->Resize(-1, 0, true);
  // associate conveniance pointers
  diso_ = state_->Dis(-1);
  velo_ = state_->Vel(-1);
  acco_ = state_->Acc(-1);

  // allocate force vectors
  fextn_ = LINALG::CreateVector(*dofrowmap_, true);
  fintn_ = LINALG::CreateVector(*dofrowmap_, true);
  fviscn_ = LINALG::CreateVector(*dofrowmap_, true);
  frimpn_ = LINALG::CreateVector(*dofrowmap_, true);

  // let it rain
  return;
}

/*----------------------------------------------------------------------*/
/* Integrate step */
void StruTimIntAB2::IntegrateStep()
{
  const double dt = dt_;  // \f$\Delta t_{n}\f$
  const double dto = dt_;  // \f$\Delta t_{n-1}\f$

  // new displacements \f$D_{n+}\f$
  disn_->Update(1.0, *dis_, 0.0);
  disn_->Update((2.0*dt*dto+dt*dt)/(2.0*dto), *vel_,
                -(dt*dt)/(2.0*dto), *velo_,
                1.0);

  // new velocities \f$V_{n+1}\f$
  veln_->Update(1.0, *vel_, 0.0);
  veln_->Update((2.0*dt*dto+dt*dt)/(2.0*dto), *acc_,
                -(dt*dt)/(2.0*dto), *acco_,
                1.0);

  // apply Dirichlet BCs
  ApplyDirichletBC(timen_, disn_, veln_, Teuchos::null);

  // build new external forces
  fextn_->PutScalar(0.0);
  ApplyForceExternal(timen_, disn_, fextn_);

  // initialise internal forces
  fintn_->PutScalar(0.0);

  // initialise stiffness matrix to zero
  stiff_->Zero();

  // ordinary internal force and stiffness
  {
    // displacement increment in step
    Epetra_Vector disinc = Epetra_Vector(*disn_);
    disinc.Update(-1.0, *dis_, 1.0);
    // internal force
    ApplyForceStiffInternal(timen_, dt,
                            disn_, zeros_, //Teuchos::rcp(&disinc,false),
                            fintn_, stiff_);
  }

  // viscous forces due Rayleigh damping
  if (damping_)
  {
    damp_->Multiply(false, *veln_, *fviscn_);
  }

  // determine time derivative of linear momentum vector,
  // ie \f$\dot{P} = M \dot{V}_{n=1}\f$
  frimpn_->Update(1.0, *fextn_, -1.0, *fintn_, 0.0);
  double some, more, less;
  frimpn_->Norm2(&some);
  fextn_->Norm2(&more);
  fintn_->Norm2(&less);
  cout << some << " " << more << " " << less << endl;
  if (damping_)
  {
    frimpn_->Update(-1.0, *fviscn_, 1.0);
  }

  // obtain new accelerations \f$A_{n+1}\f$
  {
    dsassert(mass_->Filled(), "Mass matrix has to be completed");
    // blank linear momentum zero on DOFs subjected to DBCs
    Epetra_Vector rhscopy = Epetra_Vector(*frimpn_);
    frimpn_->Multiply(1.0, *invtoggle_, rhscopy, 0.0);
    // get accelerations 
    accn_->PutScalar(0.0);
    // refactor==false: This is not necessary, because we always
    // use the same constant mass matrix, which was firstly factorised
    // in StruTimInt::DetermineMassDampConsistAccel
    solver_.Solve(mass_->EpetraMatrix(), accn_, frimpn_, false, true);
  }

  // apply Dirichlet BCs on accelerations
  ApplyDirichletBC(timen_, Teuchos::null, Teuchos::null, accn_);

  // wassup?
  return;
}

/*----------------------------------------------------------------------*/
/* Update step */
void StruTimIntAB2::UpdateStep()
{
  // update all old state at t_{n-1} etc
  // important for step size adaptivity
  state_->UpdateStep();
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1}
  dis_->Update(1.0, *disn_, 0.0);
  // new velocities at t_{n+1} -> t_n
  //    V_{n} := V_{n+1}
  vel_->Update(1.0, *veln_, 0.0);
  // new accelerations at t_{n+1} -> t_n
  //    A_{n} := A_{n+1}
  acc_->Update(1.0, *accn_, 0.0);

  // update anything that needs to be updated at the element level
  {
    // create the parameters for the discretization
    ParameterList p;
    // other parameters that might be needed by the elements
    p.set("total time", timen_);
    p.set("delta time", dt_);
    // action for elements
    p.set("action", "calc_struct_update_istep");    
    // go to elements
    discret_.Evaluate(p, null, null, null, null, null);
  }

  // bye
  return;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
