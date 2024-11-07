// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beamcontact_str_model_evaluator_beaminteraction_old.hpp"

#include "4C_beamcontact_beam3contact_manager.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_io.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Solid::ModelEvaluator::BeamInteractionOld::BeamInteractionOld()
    : disnp_ptr_(nullptr),
      stiff_beaminteract_ptr_(nullptr),
      f_beaminteract_np_ptr_(nullptr),
      beamcman_(nullptr)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::ModelEvaluator::BeamInteractionOld::setup()
{
  if (not is_init()) FOUR_C_THROW("init() has not been called, yet!");

  // setup the pointers for displacement and stiffness
  disnp_ptr_ = global_state().get_dis_np();
  stiff_beaminteract_ptr_ = std::make_shared<Core::LinAlg::SparseMatrix>(
      *global_state().dof_row_map_view(), 81, true, true);
  f_beaminteract_np_ptr_ =
      std::make_shared<Core::LinAlg::Vector<double>>(*global_state().dof_row_map(), true);

  // create beam contact manager
  beamcman_ = std::make_shared<CONTACT::Beam3cmanager>(*discret_ptr(), 0.0);

  // gmsh output at beginning of simulation
#ifdef GMSHTIMESTEPS
  beamcman_->GmshOutput(*disnp_ptr_, 0, 0, true);
#endif

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::ModelEvaluator::BeamInteractionOld::reset(const Core::LinAlg::Vector<double>& x)
{
  check_init_setup();

  // update the structural displacement vector
  disnp_ptr_ = global_state().get_dis_np();

  // Zero out force and stiffness contributions
  f_beaminteract_np_ptr_->PutScalar(0.0);
  stiff_beaminteract_ptr_->zero();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Solid::ModelEvaluator::BeamInteractionOld::evaluate_force()
{
  check_init_setup();

  // ToDo replace the parameter list by a beam interaction model data container
  Teuchos::ParameterList beamcontactparams;
  beamcontactparams.set("iter", eval_data().get_nln_iter());
  beamcontactparams.set("dt", eval_data().get_delta_time());
  beamcontactparams.set("numstep", eval_data().get_step_np());

  beamcman_->evaluate(*stiff_beaminteract_ptr_, *f_beaminteract_np_ptr_, *disnp_ptr_,
      beamcontactparams, true, eval_data().get_total_time());

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Solid::ModelEvaluator::BeamInteractionOld::evaluate_stiff()
{
  check_init_setup();


  // ToDo replace the parameter list by a beam interaction model data container
  Teuchos::ParameterList beamcontactparams;
  beamcontactparams.set("iter", eval_data().get_nln_iter());
  beamcontactparams.set("dt", eval_data().get_delta_time());
  beamcontactparams.set("numstep", eval_data().get_step_np());

  beamcman_->evaluate(*stiff_beaminteract_ptr_, *f_beaminteract_np_ptr_, *disnp_ptr_,
      beamcontactparams, true, eval_data().get_total_time());

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Solid::ModelEvaluator::BeamInteractionOld::evaluate_force_stiff()
{
  check_init_setup();

  // ToDo replace the parameter list by a beam interaction model data container
  Teuchos::ParameterList beamcontactparams;
  beamcontactparams.set("iter", eval_data().get_nln_iter());
  beamcontactparams.set("dt", eval_data().get_delta_time());
  beamcontactparams.set("numstep", eval_data().get_step_np());

  beamcman_->evaluate(*stiff_beaminteract_ptr_, *f_beaminteract_np_ptr_, *disnp_ptr_,
      beamcontactparams, true, eval_data().get_total_time());

  // visualization of current Newton step
#ifdef GMSHNEWTONSTEPS
  beamcman_->GmshOutput(*disnp_ptr_, EvalData().get_step_np(), EvalData().GetNlnIter());
  beamcman_->ConsoleOutput();
#endif

  // update constraint norm
  beamcman_->update_constr_norm();  // ToDo

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Solid::ModelEvaluator::BeamInteractionOld::assemble_force(
    Core::LinAlg::Vector<double>& f, const double& timefac_np) const
{
  // Todo take care of the minus sign in front of timefac_np
  Core::LinAlg::assemble_my_vector(1.0, f, -timefac_np, *f_beaminteract_np_ptr_);
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Solid::ModelEvaluator::BeamInteractionOld::assemble_jacobian(
    Core::LinAlg::SparseOperator& jac, const double& timefac_np) const
{
  std::shared_ptr<Core::LinAlg::SparseMatrix> jac_dd_ptr = global_state().extract_displ_block(jac);
  jac_dd_ptr->add(*stiff_beaminteract_ptr_, false, timefac_np, 1.0);
  // no need to keep it
  stiff_beaminteract_ptr_->zero();

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::ModelEvaluator::BeamInteractionOld::write_restart(
    Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  beamcman_->write_restart(iowriter);  // ToDo

  // since the global output_step_state() routine is not called, if the
  // restart is written, we have to do it here manually.
  output_step_state(iowriter);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::ModelEvaluator::BeamInteractionOld::read_restart(
    Core::IO::DiscretizationReader& ioreader)
{
  beamcman_->read_restart(ioreader);  // ToDo
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::ModelEvaluator::BeamInteractionOld::run_post_compute_x(
    const Core::LinAlg::Vector<double>& xold, const Core::LinAlg::Vector<double>& dir,
    const Core::LinAlg::Vector<double>& xnew)
{
  // empty ToDo
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::ModelEvaluator::BeamInteractionOld::update_step_state(const double& timefac_n)
{
  beamcman_->update(*disnp_ptr_, eval_data().get_step_np(), eval_data().get_nln_iter());

  // add the old time factor scaled contributions to the residual
  std::shared_ptr<Core::LinAlg::Vector<double>>& fstructold_ptr =
      global_state().get_fstructure_old();

  // Todo take care of the minus sign in front of timefac_np
  fstructold_ptr->Update(-timefac_n, *f_beaminteract_np_ptr_, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::ModelEvaluator::BeamInteractionOld::update_step_element()
{
  // empty
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::ModelEvaluator::BeamInteractionOld::determine_stress_strain()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::ModelEvaluator::BeamInteractionOld::determine_energy()
{
  // ToDo
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::ModelEvaluator::BeamInteractionOld::determine_optional_quantity()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::ModelEvaluator::BeamInteractionOld::output_step_state(
    Core::IO::DiscretizationWriter& iowriter) const
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::ModelEvaluator::BeamInteractionOld::reset_step_state() { return; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::shared_ptr<const Epetra_Map>
Solid::ModelEvaluator::BeamInteractionOld::get_block_dof_row_map_ptr() const
{
  check_init_setup();
  return global_state().dof_row_map();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>>
Solid::ModelEvaluator::BeamInteractionOld::get_current_solution_ptr() const
{
  // there are no model specific solution entries
  return nullptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>>
Solid::ModelEvaluator::BeamInteractionOld::get_last_time_step_solution_ptr() const
{
  // there are no model specific solution entries
  return nullptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::ModelEvaluator::BeamInteractionOld::post_output()
{
  check_init_setup();
  // empty

  return;
}  // post_output()

FOUR_C_NAMESPACE_CLOSE
