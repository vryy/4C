// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_porofluidmultiphase_meshtying_strategy_std.hpp"

#include "4C_linear_solver_method_linalg.hpp"
#include "4C_porofluidmultiphase_utils.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                (public) kremheller 04/18 |
 *----------------------------------------------------------------------*/
POROFLUIDMULTIPHASE::MeshtyingStrategyStd::MeshtyingStrategyStd(
    POROFLUIDMULTIPHASE::TimIntImpl* porofluidmultitimint, const Teuchos::ParameterList& probparams,
    const Teuchos::ParameterList& poroparams)
    : MeshtyingStrategyBase(porofluidmultitimint, probparams, poroparams)
{
  return;
}



/*----------------------------------------------------------------------*
 | prepare time loop                                   kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyStd::prepare_time_loop() { return; }

/*----------------------------------------------------------------------*
 | setup the variables to do a new time step  (public) kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyStd::prepare_time_step() { return; }

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                     kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyStd::update() { return; }

/*--------------------------------------------------------------------------*
 | initialize the linear solver                            kremheller 07/20 |
 *--------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyStd::initialize_linear_solver(
    std::shared_ptr<Core::LinAlg::Solver> solver)
{
  porofluidmultitimint_->discretization()->compute_null_space_if_necessary(solver->params());
}

/*--------------------------------------------------------------------------*
 | solve linear system of equations                        kremheller 04/18 |
 *--------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyStd::linear_solve(
    std::shared_ptr<Core::LinAlg::Solver> solver,
    std::shared_ptr<Core::LinAlg::SparseOperator> sysmat,
    std::shared_ptr<Core::LinAlg::Vector<double>> increment,
    std::shared_ptr<Core::LinAlg::Vector<double>> residual,
    Core::LinAlg::SolverParams& solver_params)
{
  solver_params.refactor = true;
  solver_params.reset = true;
  solver->solve(sysmat->epetra_operator(), increment, residual, solver_params);

  return;
}

/*----------------------------------------------------------------------*
 | Calculate problem specific norm                     kremheller 03/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyStd::calculate_norms(std::vector<double>& preresnorm,
    std::vector<double>& incprenorm, std::vector<double>& prenorm,
    const std::shared_ptr<const Core::LinAlg::Vector<double>> increment)
{
  preresnorm.resize(1);
  incprenorm.resize(1);
  prenorm.resize(1);

  preresnorm[0] = Utils::calculate_vector_norm(vectornormfres_, *porofluidmultitimint_->rhs());
  incprenorm[0] = Utils::calculate_vector_norm(vectornorminc_, *increment);
  prenorm[0] = Utils::calculate_vector_norm(vectornorminc_, *porofluidmultitimint_->phinp());

  return;
}

/*----------------------------------------------------------------------*
 | create result test for this field                   kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyStd::create_field_test() { return; }

/*----------------------------------------------------------------------*
 |  read restart data                                  kremheller 04/18 |
 -----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyStd::read_restart(const int step) { return; }

/*----------------------------------------------------------------------*
 | output of solution vector to BINIO                  kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyStd::output() { return; }

/*----------------------------------------------------------------------*
 | evaluate matrix and rhs                             kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyStd::evaluate() { return; }

/*----------------------------------------------------------------------*
 | extract and update                                  kremheller 04/18 |
 *----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>>
POROFLUIDMULTIPHASE::MeshtyingStrategyStd::extract_and_update_iter(
    const std::shared_ptr<const Core::LinAlg::Vector<double>> inc)
{
  return inc;
}

/*----------------------------------------------------------------------*
 | extract and update                                  kremheller 04/18 |
 *----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>>
POROFLUIDMULTIPHASE::MeshtyingStrategyStd::combined_increment(
    const std::shared_ptr<const Core::LinAlg::Vector<double>> inc) const
{
  return inc;
}

/*----------------------------------------------------------------------*
 | check initial fields                                kremheller 06/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyStd::check_initial_fields(
    std::shared_ptr<const Core::LinAlg::Vector<double>> vec_cont) const
{
  return;
}

/*-------------------------------------------------------------------------*
 | set element pairs that are close                       kremheller 03/19 |
 *------------------------------------------------------------------------ */
void POROFLUIDMULTIPHASE::MeshtyingStrategyStd::set_nearby_ele_pairs(
    const std::map<int, std::set<int>>* nearbyelepairs)
{
  return;
}

/*-------------------------------------------------------------------------*
 | setup the strategy                                     kremheller 03/19 |
 *------------------------------------------------------------------------ */
void POROFLUIDMULTIPHASE::MeshtyingStrategyStd::setup() { return; }

/*----------------------------------------------------------------------*
 | apply mesh movement                                 kremheller 06/18 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::MeshtyingStrategyStd::apply_mesh_movement() const { return; }

FOUR_C_NAMESPACE_CLOSE
