/*----------------------------------------------------------------------*/
/*! \file
\brief Structural time integration with forward Euler (explicit)
\level 2
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_structure_timint_expleuler.hpp"

#include "4C_contact_meshtying_contact_bridge.hpp"
#include "4C_io.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_mortar_manager_base.hpp"
#include "4C_mortar_strategy_base.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/* Constructor */
STR::TimIntExplEuler::TimIntExplEuler(const Teuchos::ParameterList& timeparams,
    const Teuchos::ParameterList& ioparams, const Teuchos::ParameterList& sdynparams,
    const Teuchos::ParameterList& xparams, Teuchos::RCP<Discret::Discretization> actdis,
    Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Core::LinAlg::Solver> contactsolver,
    Teuchos::RCP<Core::IO::DiscretizationWriter> output)
    : TimIntExpl(timeparams, ioparams, sdynparams, xparams, actdis, solver, contactsolver, output),
      modexpleuler_(Core::UTILS::IntegralValue<int>(sdynparams, "MODIFIEDEXPLEULER") == 1),
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
void STR::TimIntExplEuler::Init(const Teuchos::ParameterList& timeparams,
    const Teuchos::ParameterList& sdynparams, const Teuchos::ParameterList& xparams,
    Teuchos::RCP<Discret::Discretization> actdis, Teuchos::RCP<Core::LinAlg::Solver> solver)
{
  // call Init() in base class
  STR::TimIntExpl::Init(timeparams, sdynparams, xparams, actdis, solver);


  // info to user
  if (myrank_ == 0)
  {
    std::cout << "with " << (modexpleuler_ ? "modified" : "standard") << " forward Euler"
              << std::endl
              << "lumping activated: " << (lumpmass_ ? "true" : "false") << std::endl
              << std::endl;
  }

  // let it rain
  return;
}

/*----------------------------------------------------------------------------------------------*
 * Setup this class                                                                 rauch 09/16 |
 *----------------------------------------------------------------------------------------------*/
void STR::TimIntExplEuler::Setup()
{
  // call Setup() in base class
  STR::TimIntExpl::Setup();

  // determine mass, damping and initial accelerations
  determine_mass_damp_consist_accel();

  // resize of multi-step quantities
  ResizeMStep();

  // allocate force vectors
  fextn_ = Core::LinAlg::CreateVector(*dof_row_map_view(), true);
  fintn_ = Core::LinAlg::CreateVector(*dof_row_map_view(), true);
  fviscn_ = Core::LinAlg::CreateVector(*dof_row_map_view(), true);
  fcmtn_ = Core::LinAlg::CreateVector(*dof_row_map_view(), true);
  frimpn_ = Core::LinAlg::CreateVector(*dof_row_map_view(), true);

  return;
}

/*----------------------------------------------------------------------*/
/* Resizing of multi-step quantities */
void STR::TimIntExplEuler::ResizeMStep()
{
  // nothing to do, because ExplEuler is a 1-step method
  return;
}

/*----------------------------------------------------------------------*/
/* Integrate step */
int STR::TimIntExplEuler::IntegrateStep()
{
  // things to be done before integrating
  PreSolve();

  // time this step
  timer_->reset();

  const double dt = (*dt_)[0];  // \f$\Delta t_{n}\f$

  // new velocities \f$V_{n+1}\f$
  veln_->Update(1.0, *(*vel_)(0), 0.0);
  veln_->Update(dt, *(*acc_)(0), 1.0);

  // new displacements \f$D_{n+1}\f$, modified expl Euler uses veln_ for
  // updating disn_
  disn_->Update(1.0, *(*dis_)(0), 0.0);
  if (modexpleuler_ == true)
    disn_->Update(dt, *(*veln_)(0), 1.0);
  else
    disn_->Update(dt, *(*vel_)(0), 1.0);

  // *********** time measurement ***********
  double dtcpu = timer_->wallTime();
  // *********** time measurement ***********

  // apply Dirichlet BCs
  apply_dirichlet_bc(timen_, disn_, veln_, Teuchos::null, false);

  // initialise stiffness matrix to zero
  stiff_->Zero();

  // build new external forces
  fextn_->PutScalar(0.0);
  apply_force_external(timen_, disn_, veln_, fextn_);

  // additional external forces are added (e.g. interface forces)
  fextn_->Update(1.0, *fifc_, 1.0);

  // TIMING
  // double dtcpu = timer_->wallTime();

  // initialise internal forces
  fintn_->PutScalar(0.0);

  // ordinary internal force and stiffness
  {
    // displacement increment in step
    Epetra_Vector disinc = Epetra_Vector(*disn_);
    disinc.Update(-1.0, *(*dis_)(0), 1.0);
    // internal force
    apply_force_internal(timen_, dt, disn_, Teuchos::rcp(&disinc, false), veln_, fintn_);
  }

  // *********** time measurement ***********
  dtele_ = timer_->wallTime() - dtcpu;
  // *********** time measurement ***********

  // viscous forces due Rayleigh damping
  if (damping_ == Inpar::STR::damp_rayleigh)
  {
    damp_->Multiply(false, *veln_, *fviscn_);
  }

  // *********** time measurement ***********
  dtcpu = timer_->wallTime();
  // *********** time measurement ***********

  // contact or meshtying forces
  if (have_contact_meshtying())
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
  dtcmt_ = timer_->wallTime() - dtcpu;
  // *********** time measurement ***********

  // determine time derivative of linear momentum vector,
  // ie \f$\dot{P} = M \dot{V}_{n=1}\f$
  frimpn_->Update(1.0, *fextn_, -1.0, *fintn_, 0.0);

  if (damping_ == Inpar::STR::damp_rayleigh)
  {
    frimpn_->Update(-1.0, *fviscn_, 1.0);
  }

  if (have_contact_meshtying())
  {
    frimpn_->Update(1.0, *fcmtn_, 1.0);
  }

  // *********** time measurement ***********
  dtcpu = timer_->wallTime();
  // *********** time measurement ***********

  // obtain new accelerations \f$A_{n+1}\f$
  {
    FOUR_C_ASSERT(mass_->Filled(), "Mass matrix has to be completed");
    // blank linear momentum zero on DOFs subjected to DBCs
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), frimpn_);
    // get accelerations
    accn_->PutScalar(0.0);

    // in case of no lumping or if mass matrix is a BlockSparseMatrix, use solver
    if (lumpmass_ == false ||
        Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(mass_) == Teuchos::null)
    {
      // refactor==false: This is not necessary, because we always
      // use the same constant mass matrix, which was firstly factorised
      // in TimInt::determine_mass_damp_consist_accel
      Core::LinAlg::SolverParams solver_params;
      solver_params.reset = true;
      solver_->Solve(mass_->EpetraOperator(), accn_, frimpn_, solver_params);
    }
    // direct inversion based on lumped mass matrix
    else
    {
      // extract the diagonal values of the mass matrix
      Teuchos::RCP<Epetra_Vector> diag = Core::LinAlg::CreateVector(
          (Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(mass_))->RowMap(), false);
      (Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(mass_))->ExtractDiagonalCopy(*diag);
      // A_{n+1} = M^{-1} . ( -fint + fext )
      accn_->ReciprocalMultiply(1.0, *diag, *frimpn_, 0.0);
    }
  }

  // apply Dirichlet BCs on accelerations
  apply_dirichlet_bc(timen_, Teuchos::null, Teuchos::null, accn_, false);

  // *********** time measurement ***********
  dtsolve_ = timer_->wallTime() - dtcpu;
  // *********** time measurement ***********

  // wassup?
  return 0;
}

/*----------------------------------------------------------------------*/
/* Update step */
void STR::TimIntExplEuler::UpdateStepState()
{
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1}
  dis_->UpdateSteps(*disn_);
  // new velocities at t_{n+1} -> t_n
  //    V_{n} := V_{n+1}
  vel_->UpdateSteps(*veln_);
  // new accelerations at t_{n+1} -> t_n
  //    A_{n} := A_{n+1}
  acc_->UpdateSteps(*accn_);

  // update contact and meshtying
  update_step_contact_meshtying();

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/* update after time step after output on element level*/
// update anything that needs to be updated at the element level
void STR::TimIntExplEuler::UpdateStepElement()
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
void STR::TimIntExplEuler::ReadRestartForce() { return; }

/*----------------------------------------------------------------------*/
/* write internal and external forces for restart */
void STR::TimIntExplEuler::WriteRestartForce(Teuchos::RCP<Core::IO::DiscretizationWriter> output)
{
  return;
}

FOUR_C_NAMESPACE_CLOSE
