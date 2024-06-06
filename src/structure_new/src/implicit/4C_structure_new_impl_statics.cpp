/*-----------------------------------------------------------*/
/*! \file

\brief Static (time) integrator.


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_impl_statics.hpp"

#include "4C_io.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_structure_new_dbc.hpp"
#include "4C_structure_new_model_evaluator.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_model_evaluator_structure.hpp"
#include "4C_structure_new_predict_generic.hpp"
#include "4C_structure_new_timint_implicit.hpp"

#include <Epetra_Vector.h>
#include <NOX_Epetra_Vector.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::IMPLICIT::Statics::Statics()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Statics::Setup()
{
  check_init();

  // Call the Setup() of the abstract base class first.
  Generic::Setup();

  // check for valid parameter combinations:
  if (eval_data().get_damping_type() != Inpar::STR::damp_none)
    FOUR_C_THROW("ERROR: Damping not provided for statics time integration!");

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Statics::post_setup()
{
  check_init_setup();
  // DO NOTHING
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Statics::set_state(const Epetra_Vector& x)
{
  check_init_setup();
  if (is_predictor_state()) return;

  Teuchos::RCP<Epetra_Vector> disnp_ptr = global_state().extract_displ_entries(x);
  global_state().get_dis_np()->Scale(1.0, *disnp_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::Statics::apply_force(const Epetra_Vector& x, Epetra_Vector& f)
{
  check_init_setup();
  reset_eval_params();
  return model_eval().apply_force(x, f, 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::Statics::apply_stiff(const Epetra_Vector& x, Core::LinAlg::SparseOperator& jac)
{
  check_init_setup();
  reset_eval_params();
  bool ok = model_eval().apply_stiff(x, jac, 1.0);
  jac.Complete();
  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::Statics::apply_force_stiff(
    const Epetra_Vector& x, Epetra_Vector& f, Core::LinAlg::SparseOperator& jac)
{
  check_init_setup();
  reset_eval_params();
  bool ok = model_eval().apply_force_stiff(x, f, jac, 1.0);
  jac.Complete();
  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::Statics::assemble_force(
    Epetra_Vector& f, const std::vector<Inpar::STR::ModelType>* without_these_models) const
{
  check_init_setup();
  return model_eval().assemble_force(1.0, f, without_these_models);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Statics::write_restart(
    Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  check_init_setup();

  // create empty dynamic forces
  auto finertialn = Core::LinAlg::CreateVector(*global_state().dof_row_map_view(), true);
  auto fviscon = Core::LinAlg::CreateVector(*global_state().dof_row_map_view(), true);

  // write dynamic forces, so that it can be used later on for restart dynamics analysis
  iowriter.WriteVector("finert", finertialn);
  iowriter.WriteVector("fvisco", fviscon);

  model_eval().write_restart(iowriter, forced_writerestart);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Statics::read_restart(Core::IO::DiscretizationReader& ioreader)
{
  check_init_setup();
  model_eval().read_restart(ioreader);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::IMPLICIT::Statics::calc_ref_norm_force(
    const enum ::NOX::Abstract::Vector::NormType& type) const
{
  check_init_setup();

  const Teuchos::RCP<Epetra_Vector> fintnp =
      Teuchos::rcp_const_cast<Epetra_Vector>(global_state().get_fint_np());
  const Teuchos::RCP<Epetra_Vector> fextnp =
      Teuchos::rcp_const_cast<Epetra_Vector>(global_state().get_fext_np());
  const Teuchos::RCP<Epetra_Vector> freactnp =
      Teuchos::rcp_const_cast<Epetra_Vector>(global_state().get_freact_np());

  // switch from Epetra_Vector to ::NOX::Epetra::Vector (view but read-only)
  Teuchos::RCP<const ::NOX::Epetra::Vector> fintnp_nox_ptr =
      Teuchos::rcp(new ::NOX::Epetra::Vector(fintnp, ::NOX::Epetra::Vector::CreateView));
  Teuchos::RCP<const ::NOX::Epetra::Vector> fextnp_nox_ptr =
      Teuchos::rcp(new ::NOX::Epetra::Vector(fextnp, ::NOX::Epetra::Vector::CreateView));
  Teuchos::RCP<const ::NOX::Epetra::Vector> freactnp_nox_ptr =
      Teuchos::rcp(new ::NOX::Epetra::Vector(freactnp, ::NOX::Epetra::Vector::CreateView));

  // norm of the internal forces
  double fintnorm = fintnp_nox_ptr->norm(type);

  // norm of the external forces
  double fextnorm = fextnp_nox_ptr->norm(type);

  // norm of reaction forces
  double freactnorm = freactnp_nox_ptr->norm(type);

  // return characteristic norm
  return std::max(fintnorm, std::max(fextnorm, freactnorm));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::IMPLICIT::Statics::get_int_param() const { return 0.0; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Statics::pre_update()
{
  check_init_setup();
  const STR::TimeInt::Implicit* impl_ptr = dynamic_cast<const STR::TimeInt::Implicit*>(&tim_int());
  if (impl_ptr == nullptr) return;

  // get the time step size
  const double dt = (*global_state().get_delta_time())[0];

  const Inpar::STR::PredEnum& pred_type = impl_ptr->Predictor().get_type();
  Teuchos::RCP<Epetra_Vector>& accnp_ptr = global_state().get_acc_np();
  Teuchos::RCP<Epetra_Vector>& velnp_ptr = global_state().get_vel_np();

  switch (pred_type)
  {
    // case: constant acceleration
    case Inpar::STR::pred_constacc:
    {
      // read-only access
      Teuchos::RCP<const Epetra_Vector> veln_ptr = global_state().get_vel_n();
      // update the pseudo acceleration (statics!)
      accnp_ptr->Update(1.0 / dt, *velnp_ptr, -1.0 / dt, *veln_ptr, 0.0);

      [[fallthrough]];
    }
    // case: constant acceleration OR constant velocity
    case Inpar::STR::pred_constvel:
    {
      // read-only access
      Teuchos::RCP<const Epetra_Vector> disn_ptr = global_state().get_dis_n();
      Teuchos::RCP<const Epetra_Vector> disnp_ptr = global_state().get_dis_np();
      // update the pseudo velocity (statics!)
      velnp_ptr->Update(1.0 / dt, *disnp_ptr, -1.0 / dt, *disn_ptr, 0.0);
      // ATTENTION: Break for both cases!
      break;
    }
    default:
      /* do nothing */
      break;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Statics::update_step_state()
{
  check_init_setup();
  // update model specific variables
  model_eval().update_step_state(0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Statics::update_step_element()
{
  check_init_setup();
  model_eval().update_step_element();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Statics::predict_const_dis_consist_vel_acc(
    Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const
{
  check_init_setup();
  // constant predictor : displacement in domain
  disnp.Update(1.0, *global_state().get_dis_n(), 0.0);
  // new end-point velocities, these stay zero in static calculation
  velnp.PutScalar(0.0);
  // new end-point accelerations, these stay zero in static calculation
  accnp.PutScalar(0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::Statics::predict_const_vel_consist_acc(
    Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const
{
  check_init_setup();
  // If there is not enough history information, return a fail status.
  if (global_state().get_step_n() == 0) return false;

  // Displacement increment over last time step
  Teuchos::RCP<Epetra_Vector> disp_inc =
      Teuchos::rcp(new Epetra_Vector(*global_state().dof_row_map_view(), true));
  disp_inc->Update((*global_state().get_delta_time())[0], *global_state().get_vel_n(), 0.);
  // apply the dbc on the auxiliary vector
  tim_int().get_dbc().apply_dirichlet_to_vector(disp_inc);
  // update the solution variables
  disnp.Update(1.0, *global_state().get_dis_n(), 0.0);
  disnp.Update(1.0, *disp_inc, 1.0);
  velnp.Update(1.0, *global_state().get_vel_n(), 0.0);
  accnp.Update(1.0, *global_state().get_acc_n(), 0.0);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::Statics::predict_const_acc(
    Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const
{
  check_init_setup();
  // If there is not enough history information try a different predictor with
  // less requirements.
  if (global_state().get_step_n() < 2) return predict_const_vel_consist_acc(disnp, velnp, accnp);

  // Displacement increment over last time step
  Teuchos::RCP<Epetra_Vector> disp_inc =
      Teuchos::rcp(new Epetra_Vector(*global_state().dof_row_map_view(), true));
  const double& dt = (*global_state().get_delta_time())[0];
  disp_inc->Update(dt, *global_state().get_vel_n(), 0.);
  disp_inc->Update(0.5 * dt * dt, *global_state().get_acc_n(), 1.0);
  // apply the dbc on the auxiliary vector
  tim_int().get_dbc().apply_dirichlet_to_vector(disp_inc);
  // update the solution variables
  disnp.Update(1.0, *global_state().get_dis_n(), 0.0);
  disnp.Update(1., *disp_inc, 1.);
  velnp.Update(1.0, *global_state().get_vel_n(), 0.0);
  accnp.Update(1.0, *global_state().get_acc_n(), 0.0);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Statics::reset_eval_params()
{
  // call base class
  STR::IMPLICIT::Generic::reset_eval_params();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::IMPLICIT::Statics::get_model_value(const Epetra_Vector& x)
{
  Teuchos::RCP<const Epetra_Vector> disnp_ptr = global_state().extract_displ_entries(x);
  const Epetra_Vector& disnp = *disnp_ptr;

  set_state(disnp);

  eval_data().clear_values_for_all_energy_types();
  STR::MODELEVALUATOR::Structure& str_model =
      dynamic_cast<STR::MODELEVALUATOR::Structure&>(evaluator(Inpar::STR::model_structure));

  str_model.determine_strain_energy(disnp, true);
  const double int_energy_np = eval_data().get_energy_data(STR::internal_energy);
  double ext_energy_np = 0.0;
  global_state().get_fext_np()->Dot(disnp, &ext_energy_np);
  const double total = int_energy_np - ext_energy_np;

  std::ostream& os = Core::IO::cout.os(Core::IO::debug);
  os << __LINE__ << __PRETTY_FUNCTION__ << "\n";
  os << "internal/strain energy       = " << int_energy_np << "\n"
     << "external energy              = " << ext_energy_np << "\n";
  os << std::string(80, '-') << "\n";
  os << "Total                     = " << total << "\n";
  os << std::string(80, '-') << "\n";


  return total;
}

FOUR_C_NAMESPACE_CLOSE
