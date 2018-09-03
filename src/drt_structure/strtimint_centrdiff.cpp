/*----------------------------------------------------------------------*/
/*!
\file strtimint_centrdiff.cpp
\brief Structural time integration with central difference scheme 2nd order (explicit)
\level 1

<pre>
\maintainer Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint_centrdiff.H"
#include "../drt_mortar/mortar_manager_base.H"
#include "../drt_mortar/mortar_strategy_base.H"
#include "../drt_contact/meshtying_contact_bridge.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"

/*----------------------------------------------------------------------*/
/* Constructor */
STR::TimIntCentrDiff::TimIntCentrDiff(const Teuchos::ParameterList& timeparams,
    const Teuchos::ParameterList& ioparams, const Teuchos::ParameterList& sdynparams,
    const Teuchos::ParameterList& xparams, Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<LINALG::Solver> contactsolver,
    Teuchos::RCP<IO::DiscretizationWriter> output)
    : TimIntExpl(timeparams, ioparams, sdynparams, xparams, actdis, solver, contactsolver, output),
      fextn_(Teuchos::null),
      fintn_(Teuchos::null),
      fviscn_(Teuchos::null),
      fcmtn_(Teuchos::null),
      frimpn_(Teuchos::null)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g.
  // redistribution of elements. Only then call the setup to this class. This will call the setup to
  // all classes in the inheritance hierarchy. This way, this class may also override a method that
  // is called during Setup() in a base class.
  return;
}

/*----------------------------------------------------------------------------------------------*
 * Initialize this class                                                            rauch 09/16 |
 *----------------------------------------------------------------------------------------------*/
void STR::TimIntCentrDiff::Init(const Teuchos::ParameterList& timeparams,
    const Teuchos::ParameterList& sdynparams, const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<LINALG::Solver> solver)
{
  // call Init() in base class
  STR::TimIntExpl::Init(timeparams, sdynparams, xparams, actdis, solver);

  // info to user
  if (myrank_ == 0)
  {
    std::cout << "with central differences" << std::endl;
  }

  // let it rain
  return;
}

/*----------------------------------------------------------------------------------------------*
 * Setup this class                                                                 rauch 09/16 |
 *----------------------------------------------------------------------------------------------*/
void STR::TimIntCentrDiff::Setup()
{
  // call Setup() in base class
  STR::TimIntExpl::Setup();

  // determine mass, damping and initial accelerations
  DetermineMassDampConsistAccel();

  // resize of multi-step quantities
  ResizeMStep();

  // allocate force vectors
  fextn_ = LINALG::CreateVector(*DofRowMapView(), true);
  fintn_ = LINALG::CreateVector(*DofRowMapView(), true);
  fviscn_ = LINALG::CreateVector(*DofRowMapView(), true);
  fcmtn_ = LINALG::CreateVector(*DofRowMapView(), true);
  frimpn_ = LINALG::CreateVector(*DofRowMapView(), true);


  return;
}

/*----------------------------------------------------------------------*/
/* Resizing of multi-step quantities */
void STR::TimIntCentrDiff::ResizeMStep()
{
  // nothing to do, because CentrDiff is a 1-step method
  return;
}

/*----------------------------------------------------------------------*/
/* Integrate step */
int STR::TimIntCentrDiff::IntegrateStep()
{
  // things to be done before integrating
  PreSolve();

  // time this step
  timer_->ResetStartTime();

  const double dt = (*dt_)[0];     // \f$\Delta t_{n}\f$
  const double dthalf = dt / 2.0;  // \f$\Delta t_{n+1/2}\f$

  // new velocities \f$V_{n+1/2}\f$
  veln_->Update(1.0, *(*vel_)(0), 0.0);
  veln_->Update(dthalf, *(*acc_)(0), 1.0);

  // new displacements \f$D_{n+1}\f$
  disn_->Update(1.0, *(*dis_)(0), 0.0);
  disn_->Update(dt, *veln_, 1.0);

  // *********** time measurement ***********
  double dtcpu = timer_->WallTime();
  // *********** time measurement ***********

  // apply Dirichlet BCs
  ApplyDirichletBC(timen_, disn_, veln_, Teuchos::null, false);

  // initialise stiffness matrix to zero
  stiff_->Zero();

  // build new external forces
  fextn_->PutScalar(0.0);
  ApplyForceExternal(timen_, disn_, veln_, fextn_);

  // TIMING
  // double dtcpu = timer_->WallTime();

  // initialise internal forces
  fintn_->PutScalar(0.0);

  // ordinary internal force and stiffness
  {
    // displacement increment in step
    Epetra_Vector disinc = Epetra_Vector(*disn_);
    disinc.Update(-1.0, *(*dis_)(0), 1.0);
    // internal force
    ApplyForceInternal(timen_, dt, disn_, Teuchos::rcp(&disinc, false), veln_, fintn_);
  }

  // *********** time measurement ***********
  dtele_ = timer_->WallTime() - dtcpu;
  // *********** time measurement ***********

  // viscous forces due Rayleigh damping
  if (damping_ == INPAR::STR::damp_rayleigh)
  {
    damp_->Multiply(false, *veln_, *fviscn_);
  }

  // *********** time measurement ***********
  dtcpu = timer_->WallTime();
  // *********** time measurement ***********

  // contact or meshtying forces
  if (HaveContactMeshtying())
  {
    fcmtn_->PutScalar(0.0);

    if (cmtbridge_->HaveMeshtying())
      cmtbridge_->MtManager()->GetStrategy().ApplyForceStiffCmt(
          disn_, stiff_, fcmtn_, stepn_, 0, false);
    if (cmtbridge_->HaveContact())
      cmtbridge_->ContactManager()->GetStrategy().ApplyForceStiffCmt(
          disn_, stiff_, fcmtn_, stepn_, 0, false);
  }

  // *********** time measurement ***********
  dtcmt_ = timer_->WallTime() - dtcpu;
  // *********** time measurement ***********

  // determine time derivative of linear momentum vector,
  // ie \f$\dot{P} = M \dot{V}_{n=1}\f$
  frimpn_->Update(1.0, *fextn_, -1.0, *fintn_, 0.0);

  if (damping_ == INPAR::STR::damp_rayleigh)
  {
    frimpn_->Update(-1.0, *fviscn_, 1.0);
  }

  if (HaveContactMeshtying())
  {
    frimpn_->Update(1.0, *fcmtn_, 1.0);
  }

  // *********** time measurement ***********
  dtcpu = timer_->WallTime();
  // *********** time measurement ***********

  // obtain new accelerations \f$A_{n+1}\f$
  {
    dsassert(mass_->Filled(), "Mass matrix has to be completed");
    // blank linear momentum zero on DOFs subjected to DBCs
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), frimpn_);
    // get accelerations
    accn_->PutScalar(0.0);

    // in case of no lumping or if mass matrix is a BlockSparseMatrix, use solver
    if (lumpmass_ == false ||
        Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(mass_) == Teuchos::null)
    {
      // linear solver call
      // refactor==false: This is not necessary, because we always
      // use the same constant mass matrix, which was firstly factorised
      // in TimInt::DetermineMassDampConsistAccel
      solver_->Solve(mass_->EpetraOperator(), accn_, frimpn_, false, true);
    }

    // direct inversion based on lumped mass matrix
    else
    {
      Teuchos::RCP<LINALG::SparseMatrix> massmatrix =
          Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(mass_);
      Teuchos::RCP<Epetra_Vector> diagonal = LINALG::CreateVector(*DofRowMapView(), true);
      int error = massmatrix->ExtractDiagonalCopy(*diagonal);
      if (error != 0) dserror("ERROR: ExtractDiagonalCopy went wrong");
      accn_->ReciprocalMultiply(1.0, *diagonal, *frimpn_, 0.0);
    }
  }

  // apply Dirichlet BCs on accelerations
  ApplyDirichletBC(timen_, Teuchos::null, Teuchos::null, accn_, false);

  // *********** time measurement ***********
  dtsolve_ = timer_->WallTime() - dtcpu;
  // *********** time measurement ***********

  // update of end-velocities \f$V_{n+1}\f$
  veln_->Update(dthalf, *accn_, 1.0);

  // things to be done after integrating
  PostSolve();

  // wassup?
  return 0;
}

/*----------------------------------------------------------------------*/
/* Update step */
void STR::TimIntCentrDiff::UpdateStepState()
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

  // update contact and meshtying
  UpdateStepContactMeshtying();

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/* update after time step after output on element level*/
// update anything that needs to be updated at the element level
void STR::TimIntCentrDiff::UpdateStepElement()
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // other parameters that might be needed by the elements
  p.set("total time", timen_);
  p.set("delta time", (*dt_)[0]);
  // action for elements
  p.set("action", "calc_struct_update_istep");
  // go to elements
  discret_->Evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
}

/*----------------------------------------------------------------------*/
/* read restart forces */
void STR::TimIntCentrDiff::ReadRestartForce()
{
  dserror("No restart ability for central differences time integrator!");
  return;
}

/*----------------------------------------------------------------------*/
/* write internal and external forces for restart */
void STR::TimIntCentrDiff::WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output)
{
  return;
}
