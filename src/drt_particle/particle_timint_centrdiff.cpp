/*----------------------------------------------------------------------*/
/*!
\file particle_timint_centrdiff.cpp
\brief Particle time integration with central difference scheme 2nd order (explicit),
       also known as Velocity-Verlet algorithm

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-152537
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_centrdiff.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*/
/* Constructor */
PARTICLE::TimIntCentrDiff::TimIntCentrDiff(
    const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& sdynparams,
    const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<LINALG::Solver> contactsolver,
    Teuchos::RCP<IO::DiscretizationWriter> output
  ) : STR::TimIntCentrDiff
  (
    ioparams,
    sdynparams,
    xparams,
    actdis,
    solver,
    contactsolver,
    output,
    false // DetermineMassDampConsistAccel() must not be called
  )
{
  // allocate vectors
  radius_  = LINALG::CreateVector(*discret_->NodeRowMap(), true);
  radiusn_  = LINALG::CreateVector(*discret_->NodeRowMap(), true);
  density_ = LINALG::CreateVector(*discret_->NodeRowMap(), true);
  densityn_  = LINALG::CreateVector(*discret_->NodeRowMap(), true);

  // initialize displacement field with nodal positions
  for(int n=0; n<discret_->NumMyRowNodes(); n++)
  {
    DRT::Node* actnode = discret_->lRowNode(n);
    // get the first gid of a node and convert it into a LID
    int gid = discret_->Dof(actnode, 0);
    int lid = discret_->DofRowMap()->LID(gid);
    for (int dim=0; dim<3; dim++)
    {
      (*(*dis_)(0))[lid+dim] = actnode->X()[dim];
    }

    // for testing reasons
    (*radiusn_)[n] = 1.3;
    (*densityn_)[n] = 2.4;

    double mass = (*densityn_)[n] * 4.0/3.0 * M_PI * pow((*radiusn_)[n], 3.0);
    for (int dim=2; dim<3; dim++)
    {
      double force = 2.3;
      (*(*acc_)(0))[lid+dim] = force / mass;
    }
    // end: for testing reasons

  }

  return;
}


/*----------------------------------------------------------------------*/
/* Integrate step */
void PARTICLE::TimIntCentrDiff::IntegrateStep()
{
  // time this step
  timer_->ResetStartTime();

  const double dt = (*dt_)[0];   // \f$\Delta t_{n}\f$
  const double dthalf = dt/2.0;  // \f$\Delta t_{n+1/2}\f$

  // new velocities \f$V_{n+1/2}\f$
  veln_->Update(1.0, *(*vel_)(0), 0.0);
  veln_->Update(dthalf, *(*acc_)(0), 1.0);

  // new displacements \f$D_{n+1}\f$
  disn_->Update(1.0, *(*dis_)(0), 0.0);
  disn_->Update(dt, *veln_, 1.0);

  // apply Dirichlet BCs
//  ApplyDirichletBC(timen_, disn_, veln_, Teuchos::null, false);

  // build new external forces
  fextn_->PutScalar(0.0);
//  ApplyForceExternal(timen_, disn_, veln_, fextn_);

//   TIMING
//  double dtcpu = timer_->WallTime();
//
//  // initialise internal forces
//  fintn_->PutScalar(0.0);
//
//  // ordinary internal force and stiffness
//  {
//    // displacement increment in step
//    Epetra_Vector disinc = Epetra_Vector(*disn_);
//    disinc.Update(-1.0, *(*dis_)(0), 1.0);
//    // internal force
//    ApplyForceInternal(timen_, dt,
//                       disn_, Teuchos::rcp(&disinc,false), veln_,
//                       fintn_);
//  }
//
//   TIMING
//  if (!myrank_) cout << "\nT_internal: " << timer_->WallTime() -dtcpu << endl;
//
//   viscous forces due Rayleigh damping
//  if (damping_ == INPAR::STR::damp_rayleigh)
//  {
//    damp_->Multiply(false, *veln_, *fviscn_);
//  }
//
//   TIMING
//  dtcpu = timer_->WallTime();



  // TIMING
  //if (!myrank_) cout << "T_contact:  " << timer_->WallTime() - dtcpu  << endl;

//   determine time derivative of linear momentum vector,
//   ie \f$\dot{P} = M \dot{V}_{n=1}\f$
//  frimpn_->Update(1.0, *fextn_, -1.0, *fintn_, 0.0);

//  if (damping_ == INPAR::STR::damp_rayleigh)
//  {
//    frimpn_->Update(-1.0, *fviscn_, 1.0);
//  }

  // TIMING
  //dtcpu = timer_->WallTime();

  // obtain new accelerations \f$A_{n+1}\f$
//  {
//    dsassert(mass_->Filled(), "Mass matrix has to be completed");
//    // blank linear momentum zero on DOFs subjected to DBCs
//    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), frimpn_);
//    // get accelerations
//    accn_->PutScalar(0.0);
//
//    // in case of no lumping or if mass matrix is a BlockSparseMatrix, use solver
//    if (lumpmass_==false || Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(mass_)==Teuchos::null)
//    {
//      // linear solver call
//      // refactor==false: This is not necessary, because we always
//      // use the same constant mass matrix, which was firstly factorised
//      // in TimInt::DetermineMassDampConsistAccel
//      solver_->Solve(mass_->EpetraOperator(), accn_, frimpn_, false, true);
//    }
//
//    // direct inversion based on lumped mass matrix
//    else
//    {
//      RCP<LINALG::SparseMatrix> massmatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(mass_);
//      RCP<Epetra_Vector> diagonal = LINALG::CreateVector(*dofrowmap_, true);
//      int error = massmatrix->ExtractDiagonalCopy(*diagonal);
//      if (error!=0) dserror("ERROR: ExtractDiagonalCopy went wrong");
//      accn_->ReciprocalMultiply(1.0,*diagonal,*frimpn_,0.0);
//    }
//  }

  // obtain new accelerations \f$A_{n+1}\f$
  for(int n=0; n<discret_->NumMyRowNodes(); n++)
  {
    DRT::Node* actnode = discret_->lRowNode(n);

    double mass = (*densityn_)[n] * 4.0/3.0 * M_PI * pow((*radiusn_)[n], 3.0);

    // get the first gid of a node and convert it into a LID
    int gid = discret_->Dof(actnode, 0);
    int lid = discret_->DofRowMap()->LID(gid);
    // Vorsicht: hier nur z momentan
    for (int i=2; i<3; i++)
    {
      double force = 2.3;
      double zahl=1.0; //1.0+std::sin(timen_*2*M_PI/7.0);
      (*accn_)[lid+i] = zahl*force / mass;
    }
  }


  // TIMING
  //if (!myrank_) cout << "T_linsolve: " << timer_->WallTime() - dtcpu << endl;

  // apply Dirichlet BCs on accelerations
//  ApplyDirichletBC(timen_, Teuchos::null, Teuchos::null, accn_, false);

  // update of end-velocities \f$V_{n+1}\f$
  veln_->Update(dthalf, *accn_, 1.0);

  return;
}


/*----------------------------------------------------------------------*/
/* State vectors are updated according to the new distribution of particles */
void PARTICLE::TimIntCentrDiff::UpdateStatesAfterParticleTransfer()
{
  Teuchos::RCP<Epetra_Vector> old;

  if (disn_ != Teuchos::null)
  {
    old = disn_;
    disn_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
    LINALG::Export(*old, *disn_);
  }

  if (veln_ != Teuchos::null)
  {
    old = veln_;
    veln_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
    LINALG::Export(*old, *veln_);
  }

  if (accn_ != Teuchos::null)
  {
    old = accn_;
    accn_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
    LINALG::Export(*old, *accn_);
  }

  if ((*dis_)(0) != Teuchos::null)
  {
    const Teuchos::RCP<Epetra_Vector> oldvec = Teuchos::rcp(new Epetra_Vector(*(*dis_)(0)));
    dis_->ReplaceMaps(discret_->DofRowMap());
    LINALG::Export(*oldvec, *(*dis_)(0));
  }

  if ((*vel_)(0) != Teuchos::null)
  {
    const Teuchos::RCP<Epetra_Vector> oldvec = Teuchos::rcp(new Epetra_Vector(*(*vel_)(0)));
    vel_->ReplaceMaps(discret_->DofRowMap());
    LINALG::Export(*oldvec, *(*vel_)(0));
  }

  if ((*acc_)(0) != Teuchos::null)
  {
    const Teuchos::RCP<Epetra_Vector> oldvec = Teuchos::rcp(new Epetra_Vector(*(*acc_)(0)));
    acc_->ReplaceMaps(discret_->DofRowMap());
    LINALG::Export(*oldvec, *(*acc_)(0));
  }

  if (radiusn_ != Teuchos::null)
  {
    old = radiusn_;
    radiusn_ = LINALG::CreateVector(*discret_->NodeRowMap(),true);
    LINALG::Export(*old, *radiusn_);
  }

  if (densityn_ != Teuchos::null)
  {
    old = densityn_;
    densityn_ = LINALG::CreateVector(*discret_->NodeRowMap(),true);
    LINALG::Export(*old, *densityn_);
  }

  return;
}


/*----------------------------------------------------------------------*/
