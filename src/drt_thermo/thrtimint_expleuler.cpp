/*----------------------------------------------------------------------*/
/*!
\brief Thermal time integration with forward Euler
\level 3
\maintainer Christoph Meier
*/


/*----------------------------------------------------------------------*
 | headers                                                   dano 01/12 |
 *----------------------------------------------------------------------*/
#include "thrtimint_expleuler.H"
#include "thermo_ele_action.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"


/*----------------------------------------------------------------------*
 | constructor                                               dano 01/12 |
 *----------------------------------------------------------------------*/
THR::TimIntExplEuler::TimIntExplEuler(const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& tdynparams, const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<IO::DiscretizationWriter> output)
    : TimIntExpl(ioparams, tdynparams, xparams, actdis, solver, output),
      fextn_(Teuchos::null),
      fintn_(Teuchos::null)
{
  // info to user
  if (myrank_ == 0)
  {
    std::cout << "with forward Euler" << std::endl
              << "lumping activated: " << (lumpcapa_ ? "true" : "false") << std::endl
              << std::endl;
  }

  // determine capacity
  DetermineCapaConsistTempRate();

  // allocate force vectors
  fextn_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
  fintn_ = LINALG::CreateVector(*discret_->DofRowMap(), true);

  // let it rain
  return;

}  // TimIntExplEuler()


/*----------------------------------------------------------------------*
 | integrate step                                            dano 01/12 |
 *----------------------------------------------------------------------*/
void THR::TimIntExplEuler::IntegrateStep()
{
  const double dt = (*dt_)[0];  // \f$\Delta t_{n}\f$

  // new displacements \f$D_{n+}\f$
  // T_{n+1} = T_n + dt * r_n
  tempn_->Update(1.0, *(*temp_)(0), 0.0);
  tempn_->Update(dt, *(*rate_)(0), 1.0);

  // apply Dirichlet BCs
  ApplyDirichletBC(timen_, tempn_, raten_, false);

  // build new external forces
  fextn_->PutScalar(0.0);
  ApplyForceExternal(timen_, tempn_, fextn_);

  // interface forces to external forces
  fextn_->Update(1.0, *fifc_, 1.0);

  // TIMING
  // double dtcpu = timer_->WallTime();

  // initialise internal forces
  fintn_->PutScalar(0.0);

  // ordinary internal force and conductivity matrix
  {
    // temperature increment in step
    Teuchos::RCP<Epetra_Vector> tempinc = Teuchos::rcp(new Epetra_Vector(*tempn_));
    tempinc->Update(-1.0, *(*temp_)(0), 1.0);
    // create an empty parameter list for the discretisation
    Teuchos::ParameterList p;
    // internal force
    ApplyForceInternal(p, timen_, dt, tempn_, tempinc, fintn_);
  }

  // determine time derivative of capacity vector, ie \f$\dot{P} = C . \dot{T}_{n=1}\f$
  Teuchos::RCP<Epetra_Vector> frimpn = LINALG::CreateVector(*discret_->DofRowMap(), true);
  frimpn->Update(1.0, *fextn_, -1.0, *fintn_, 0.0);

  // obtain new temperature rates \f$R_{n+1}\f$
  {
    dsassert(tang_->Filled(), "capacity matrix has to be completed");
    // get accelerations
    raten_->PutScalar(0.0);
  }

  if ((lumpcapa_ == false) or
      (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(tang_) == Teuchos::null))
  {
    // refactor==false: This is not necessary, because we always
    // use the same constant capacity matrix, which was firstly factorised
    // in TimInt::DetermineCapaConsistTempRate
    solver_->Solve(tang_->EpetraOperator(), raten_, frimpn, false, true);
  }
  // direct inversion based on lumped capacity matrix
  else
  {
    // extract the diagonal values of the mass matrix
    Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(
        (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(tang_))->RowMap(), false);
    (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(tang_))->ExtractDiagonalCopy(*diag);
    // R_{n+1} = C^{-1} . ( -fint + fext )
    raten_->ReciprocalMultiply(1.0, *diag, *frimpn, 0.0);
  }

  // apply Dirichlet BCs on temperature rates
  ApplyDirichletBC(timen_, Teuchos::null, raten_, false);

  // wassup?
  return;

}  // IntegrateStep()


/*----------------------------------------------------------------------*
 | update step                                               dano 01/12 |
 *----------------------------------------------------------------------*/
void THR::TimIntExplEuler::UpdateStepState()
{
  // new temperatures at t_{n+1} -> t_n
  // T_n := T_{n+1}
  temp_->UpdateSteps(*tempn_);
  // new temperature rates at t_{n+1} -> t_n
  // R_n := R_{n+1}
  rate_->UpdateSteps(*raten_);

  // bye
  return;
}  // UpdateStepState()


/*----------------------------------------------------------------------*
 | update after time step after output on element level      dano 01/12 |
 | update anything that needs to be updated at the element level        |
 *----------------------------------------------------------------------*/
void THR::TimIntExplEuler::UpdateStepElement()
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // other parameters that might be needed by the elements
  p.set("total time", timen_);
  p.set("delta time", (*dt_)[0]);
  // action for elements
  // --> be careful: this action does nothing
  p.set<int>("action", THR::calc_thermo_update_istep);
  // go to elements and do nothing
  discret_->Evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

}  // UpdateStepElement()


/*----------------------------------------------------------------------*
 | read restart forces                                       dano 01/12 |
 *----------------------------------------------------------------------*/
void THR::TimIntExplEuler::ReadRestartForce()
{
  // do nothing
  return;

}  // ReadRestartForce


/*----------------------------------------------------------------------*
 | read restart forces                                       dano 07/13 |
 *----------------------------------------------------------------------*/
void THR::TimIntExplEuler::WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output)
{
  // do nothing
  return;

}  // WriteRestartForce()


/*----------------------------------------------------------------------*/
