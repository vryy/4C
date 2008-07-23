/*----------------------------------------------------------------------*/
/*!
\file strtimint_ab2.cpp
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
#include "strtimint_ab2.H"


/*----------------------------------------------------------------------*/
/* Constructor */
STR::StruTimIntAB2::StruTimIntAB2
(
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& sdynparams,
  const Teuchos::ParameterList& xparams,
  //const Teuchos::ParameterList& ab2params,
  Teuchos::RCP<DRT::Discretization> actdis,
  Teuchos::RCP<LINALG::Solver> solver,
  Teuchos::RCP<IO::DiscretizationWriter> output
)
: StruTimIntExpl
  (
    ioparams,
    sdynparams,
    xparams,
    actdis,
    solver,
    output
  ),
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

  // resize of multi-step quantities
  ResizeMStep();

  // allocate force vectors
  fextn_ = LINALG::CreateVector(*dofrowmap_, true);
  fintn_ = LINALG::CreateVector(*dofrowmap_, true);
  fviscn_ = LINALG::CreateVector(*dofrowmap_, true);
  frimpn_ = LINALG::CreateVector(*dofrowmap_, true);

  // let it rain
  return;
}

/*----------------------------------------------------------------------*/
/* Resizing of multi-step quantities */
void STR::StruTimIntAB2::ResizeMStep()
{
  // resize time and stepsize fields
  time_->Resize(-1, 0, (*time_)[0]);
  dt_->Resize(-1, 0, (*dt_)[0]);

  // resize state vectors, AB2 is a 2-step method, thus we need two
  // past steps at t_{n} and t_{n-1}
  dis_->Resize(-1, 0, dofrowmap_, true);
  vel_->Resize(-1, 0, dofrowmap_, true);
  acc_->Resize(-1, 0, dofrowmap_, true);
}

/*----------------------------------------------------------------------*/
/* Integrate step */
void STR::StruTimIntAB2::IntegrateStep()
{
  const double dt = (*dt_)[0];  // \f$\Delta t_{n}\f$
  const double dto = (*dt_)[-1];  // \f$\Delta t_{n-1}\f$

  // new displacements \f$D_{n+}\f$
  disn_->Update(1.0, *(*dis_)(0), 0.0);
  disn_->Update((2.0*dt*dto+dt*dt)/(2.0*dto), *(*vel_)(0),
                -(dt*dt)/(2.0*dto), *(*vel_)(-1),
                1.0);

  // new velocities \f$V_{n+1}\f$
  veln_->Update(1.0, *(*vel_)(0), 0.0);
  veln_->Update((2.0*dt*dto+dt*dt)/(2.0*dto), *(*acc_)(0),
                -(dt*dt)/(2.0*dto), *(*acc_)(-1),
                1.0);

  // apply Dirichlet BCs
  ApplyDirichletBC(timen_, disn_, veln_, Teuchos::null);

  // build new external forces
  fextn_->PutScalar(0.0);
  ApplyForceExternal(timen_, disn_, veln_, fextn_);

  // initialise internal forces
  fintn_->PutScalar(0.0);

  // initialise stiffness matrix to zero
  stiff_->Zero();

  // ordinary internal force and stiffness
  {
    // displacement increment in step
    Epetra_Vector disinc = Epetra_Vector(*disn_);
    disinc.Update(-1.0, *(*dis_)(0), 1.0);
    // internal force
    ApplyForceInternal(timen_, dt,
                       disn_, Teuchos::rcp(&disinc,false), veln_,
                       fintn_);
  }

  // viscous forces due Rayleigh damping
  if (damping_ == damp_rayleigh)
  {
    damp_->Multiply(false, *veln_, *fviscn_);
  }

  // determine time derivative of linear momentum vector,
  // ie \f$\dot{P} = M \dot{V}_{n=1}\f$
  frimpn_->Update(1.0, *fextn_, -1.0, *fintn_, 0.0);
//   double some, more, less;
//   frimpn_->Norm2(&some);
//   fextn_->Norm2(&more);
//   fintn_->Norm2(&less);
//   cout << some << " " << more << " " << less << endl;
  if (damping_ == damp_rayleigh)
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
    solver_->Solve(mass_->EpetraMatrix(), accn_, frimpn_, false, true);
  }

  // apply Dirichlet BCs on accelerations
  ApplyDirichletBC(timen_, Teuchos::null, Teuchos::null, accn_);

  // wassup?
  return;
}

/*----------------------------------------------------------------------*/
/* Update step */
void STR::StruTimIntAB2::UpdateStep()
{
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1}, D_{n-1} := D_{n}
  dis_->UpdateSteps(*disn_);
  // new velocities at t_{n+1} -> t_n
  //    V_{n} := V_{n+1}, V_{n-1} := V_{n}
  vel_->UpdateSteps(*veln_);
  // new accelerations at t_{n+1} -> t_n
  //    A_{n} := A_{n+1}, A_{n-1} := A_{n}
  acc_->UpdateSteps(*accn_);

  // update anything that needs to be updated at the element level
  {
    // create the parameters for the discretization
    ParameterList p;
    // other parameters that might be needed by the elements
    p.set("total time", timen_);
    p.set("delta time", (*dt_)[0]);
    // action for elements
    p.set("action", "calc_struct_update_istep");    
    // go to elements
    discret_->Evaluate(p, null, null, null, null, null);
  }

  // bye
  return;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
