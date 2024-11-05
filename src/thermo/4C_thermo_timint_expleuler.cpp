// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_thermo_timint_expleuler.hpp"

#include "4C_io.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_thermo_ele_action.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | constructor                                               dano 01/12 |
 *----------------------------------------------------------------------*/
Thermo::TimIntExplEuler::TimIntExplEuler(const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& tdynparams, const Teuchos::ParameterList& xparams,
    std::shared_ptr<Core::FE::Discretization> actdis, std::shared_ptr<Core::LinAlg::Solver> solver,
    std::shared_ptr<Core::IO::DiscretizationWriter> output)
    : TimIntExpl(ioparams, tdynparams, xparams, actdis, solver, output),
      fextn_(nullptr),
      fintn_(nullptr)
{
  // info to user
  if (myrank_ == 0)
  {
    std::cout << "with forward Euler" << std::endl
              << "lumping activated: " << (lumpcapa_ ? "true" : "false") << std::endl
              << std::endl;
  }

  // determine capacity
  determine_capa_consist_temp_rate();

  // allocate force vectors
  fextn_ = Core::LinAlg::create_vector(*discret_->dof_row_map(), true);
  fintn_ = Core::LinAlg::create_vector(*discret_->dof_row_map(), true);

  // let it rain
  return;

}  // TimIntExplEuler()


/*----------------------------------------------------------------------*
 | integrate step                                            dano 01/12 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntExplEuler::integrate_step()
{
  const double dt = (*dt_)[0];  // \f$\Delta t_{n}\f$

  // new displacements \f$D_{n+}\f$
  // T_{n+1} = T_n + dt * r_n
  tempn_->Update(1.0, *(*temp_)(0), 0.0);
  tempn_->Update(dt, *(*rate_)(0), 1.0);

  // apply Dirichlet BCs
  apply_dirichlet_bc(timen_, tempn_, raten_, false);

  // build new external forces
  fextn_->PutScalar(0.0);
  apply_force_external(timen_, tempn_, *fextn_);

  // interface forces to external forces
  fextn_->Update(1.0, *fifc_, 1.0);

  // TIMING
  // double dtcpu = timer_->WallTime();

  // initialise internal forces
  fintn_->PutScalar(0.0);

  // ordinary internal force and conductivity matrix
  {
    // temperature increment in step
    std::shared_ptr<Core::LinAlg::Vector<double>> tempinc =
        std::make_shared<Core::LinAlg::Vector<double>>(*tempn_);
    tempinc->Update(-1.0, *(*temp_)(0), 1.0);
    // create an empty parameter list for the discretisation
    Teuchos::ParameterList p;
    // internal force
    apply_force_internal(p, timen_, dt, tempn_, tempinc, fintn_);
  }

  // determine time derivative of capacity vector, ie \f$\dot{P} = C . \dot{T}_{n=1}\f$
  std::shared_ptr<Core::LinAlg::Vector<double>> frimpn =
      Core::LinAlg::create_vector(*discret_->dof_row_map(), true);
  frimpn->Update(1.0, *fextn_, -1.0, *fintn_, 0.0);

  // obtain new temperature rates \f$R_{n+1}\f$
  {
    FOUR_C_ASSERT(tang_->filled(), "capacity matrix has to be completed");
    // get accelerations
    raten_->PutScalar(0.0);
  }

  if ((lumpcapa_ == false) or
      (std::dynamic_pointer_cast<Core::LinAlg::SparseMatrix>(tang_) == nullptr))
  {
    // refactor==false: This is not necessary, because we always
    // use the same constant capacity matrix, which was firstly factorised
    // in TimInt::determine_capa_consist_temp_rate
    Core::LinAlg::SolverParams solver_params;
    solver_params.reset = true;
    solver_->solve(tang_->epetra_operator(), raten_, frimpn, solver_params);
  }
  // direct inversion based on lumped capacity matrix
  else
  {
    // extract the diagonal values of the mass matrix
    std::shared_ptr<Core::LinAlg::Vector<double>> diag = Core::LinAlg::create_vector(
        (std::dynamic_pointer_cast<Core::LinAlg::SparseMatrix>(tang_))->row_map());
    (std::dynamic_pointer_cast<Core::LinAlg::SparseMatrix>(tang_))->extract_diagonal_copy(*diag);
    // R_{n+1} = C^{-1} . ( -fint + fext )
    raten_->ReciprocalMultiply(1.0, *diag, *frimpn, 0.0);
  }

  // apply Dirichlet BCs on temperature rates
  apply_dirichlet_bc(timen_, nullptr, raten_, false);

  // wassup?
  return;

}  // IntegrateStep()


/*----------------------------------------------------------------------*
 | update step                                               dano 01/12 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntExplEuler::update_step_state()
{
  // new temperatures at t_{n+1} -> t_n
  // T_n := T_{n+1}
  temp_->update_steps(*tempn_);
  // new temperature rates at t_{n+1} -> t_n
  // R_n := R_{n+1}
  rate_->update_steps(*raten_);

  // bye
  return;
}  // update_step_state()


/*----------------------------------------------------------------------*
 | update after time step after output on element level      dano 01/12 |
 | update anything that needs to be updated at the element level        |
 *----------------------------------------------------------------------*/
void Thermo::TimIntExplEuler::update_step_element()
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // other parameters that might be needed by the elements
  p.set("total time", timen_);
  p.set("delta time", (*dt_)[0]);
  // action for elements
  // --> be careful: this action does nothing
  p.set<Thermo::Action>("action", Thermo::calc_thermo_update_istep);
  // go to elements and do nothing
  discret_->evaluate(p, nullptr, nullptr, nullptr, nullptr, nullptr);

}  // update_step_element()


/*----------------------------------------------------------------------*
 | read restart forces                                       dano 01/12 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntExplEuler::read_restart_force()
{
  // do nothing
  return;

}  // ReadRestartForce


/*----------------------------------------------------------------------*
 | read restart forces                                       dano 07/13 |
 *----------------------------------------------------------------------*/
void Thermo::TimIntExplEuler::write_restart_force(
    std::shared_ptr<Core::IO::DiscretizationWriter> output)
{
  // do nothing
  return;

}  // WriteRestartForce()


/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
