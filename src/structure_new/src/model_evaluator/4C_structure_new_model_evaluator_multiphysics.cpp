// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_structure_new_model_evaluator_multiphysics.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Solid::ModelEvaluator::Multiphysics::Multiphysics() : active_mt_(mt_none)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Multiphysics::init(
    const std::shared_ptr<Solid::ModelEvaluator::Data>& eval_data_ptr,
    const std::shared_ptr<Solid::TimeInt::BaseDataGlobalState>& gstate_ptr,
    const std::shared_ptr<Solid::TimeInt::BaseDataIO>& gio_ptr,
    const std::shared_ptr<Solid::Integrator>& int_ptr,
    const std::shared_ptr<const Solid::TimeInt::Base>& timint_ptr, const int& dof_offset)
{
  Solid::ModelEvaluator::Generic::init(
      eval_data_ptr, gstate_ptr, gio_ptr, int_ptr, timint_ptr, dof_offset);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Multiphysics::setup()
{
  check_init();
  issetup_ = true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::ModelEvaluator::Multiphysics::reset(const Core::LinAlg::Vector<double>& x)
{
  check_active_model_type();

  get_model_evaluator_from_map(active_mt_)->reset(x);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Solid::ModelEvaluator::Multiphysics::assemble_force(
    Core::LinAlg::Vector<double>& f, const double& timefac_np) const
{
  check_active_model_type();

  get_model_evaluator_from_map(active_mt_)->assemble_force(f, timefac_np);
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Solid::ModelEvaluator::Multiphysics::assemble_jacobian(
    Core::LinAlg::SparseOperator& jac, const double& timefac_np) const
{
  check_active_model_type();

  get_model_evaluator_from_map(active_mt_)->assemble_jacobian(jac, timefac_np);
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Solid::ModelEvaluator::Multiphysics::evaluate_force()
{
  check_active_model_type();

  get_model_evaluator_from_map(active_mt_)->evaluate_force();
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Solid::ModelEvaluator::Multiphysics::evaluate_stiff()
{
  check_active_model_type();

  get_model_evaluator_from_map(active_mt_)->evaluate_stiff();
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Solid::ModelEvaluator::Multiphysics::evaluate_force_stiff()
{
  check_active_model_type();

  get_model_evaluator_from_map(active_mt_)->evaluate_force_stiff();
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::ModelEvaluator::Multiphysics::update_step_state(const double& timefac_n)
{
  check_active_model_type();

  get_model_evaluator_from_map(active_mt_)->update_step_state(timefac_n);
}

FOUR_C_NAMESPACE_CLOSE
