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
  if (EvalData().GetDampingType() != INPAR::STR::damp_none)
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
  if (IsPredictorState()) return;

  Teuchos::RCP<Epetra_Vector> disnp_ptr = global_state().ExtractDisplEntries(x);
  global_state().GetDisNp()->Scale(1.0, *disnp_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::Statics::ApplyForce(const Epetra_Vector& x, Epetra_Vector& f)
{
  check_init_setup();
  reset_eval_params();
  return ModelEval().ApplyForce(x, f, 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::Statics::ApplyStiff(const Epetra_Vector& x, CORE::LINALG::SparseOperator& jac)
{
  check_init_setup();
  reset_eval_params();
  bool ok = ModelEval().ApplyStiff(x, jac, 1.0);
  jac.Complete();
  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::Statics::ApplyForceStiff(
    const Epetra_Vector& x, Epetra_Vector& f, CORE::LINALG::SparseOperator& jac)
{
  check_init_setup();
  reset_eval_params();
  bool ok = ModelEval().ApplyForceStiff(x, f, jac, 1.0);
  jac.Complete();
  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::Statics::assemble_force(
    Epetra_Vector& f, const std::vector<INPAR::STR::ModelType>* without_these_models) const
{
  check_init_setup();
  return ModelEval().assemble_force(1.0, f, without_these_models);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Statics::write_restart(
    IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  check_init_setup();

  // create empty dynamic forces
  auto finertialn = CORE::LINALG::CreateVector(*global_state().DofRowMapView(), true);
  auto fviscon = CORE::LINALG::CreateVector(*global_state().DofRowMapView(), true);

  // write dynamic forces, so that it can be used later on for restart dynamics analysis
  iowriter.WriteVector("finert", finertialn);
  iowriter.WriteVector("fvisco", fviscon);

  ModelEval().write_restart(iowriter, forced_writerestart);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Statics::read_restart(IO::DiscretizationReader& ioreader)
{
  check_init_setup();
  ModelEval().read_restart(ioreader);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::IMPLICIT::Statics::CalcRefNormForce(
    const enum ::NOX::Abstract::Vector::NormType& type) const
{
  check_init_setup();

  const Teuchos::RCP<Epetra_Vector> fintnp =
      Teuchos::rcp_const_cast<Epetra_Vector>(global_state().GetFintNp());
  const Teuchos::RCP<Epetra_Vector> fextnp =
      Teuchos::rcp_const_cast<Epetra_Vector>(global_state().GetFextNp());
  const Teuchos::RCP<Epetra_Vector> freactnp =
      Teuchos::rcp_const_cast<Epetra_Vector>(global_state().GetFreactNp());

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
double STR::IMPLICIT::Statics::GetIntParam() const { return 0.0; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Statics::PreUpdate()
{
  check_init_setup();
  const STR::TIMINT::Implicit* impl_ptr = dynamic_cast<const STR::TIMINT::Implicit*>(&tim_int());
  if (impl_ptr == nullptr) return;

  // get the time step size
  const double dt = (*global_state().GetDeltaTime())[0];

  const INPAR::STR::PredEnum& pred_type = impl_ptr->Predictor().GetType();
  Teuchos::RCP<Epetra_Vector>& accnp_ptr = global_state().GetAccNp();
  Teuchos::RCP<Epetra_Vector>& velnp_ptr = global_state().GetVelNp();

  switch (pred_type)
  {
    // case: constant acceleration
    case INPAR::STR::pred_constacc:
    {
      // read-only access
      Teuchos::RCP<const Epetra_Vector> veln_ptr = global_state().GetVelN();
      // update the pseudo acceleration (statics!)
      accnp_ptr->Update(1.0 / dt, *velnp_ptr, -1.0 / dt, *veln_ptr, 0.0);

      [[fallthrough]];
    }
    // case: constant acceleration OR constant velocity
    case INPAR::STR::pred_constvel:
    {
      // read-only access
      Teuchos::RCP<const Epetra_Vector> disn_ptr = global_state().GetDisN();
      Teuchos::RCP<const Epetra_Vector> disnp_ptr = global_state().GetDisNp();
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
void STR::IMPLICIT::Statics::UpdateStepState()
{
  check_init_setup();
  // update model specific variables
  ModelEval().UpdateStepState(0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Statics::UpdateStepElement()
{
  check_init_setup();
  ModelEval().UpdateStepElement();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Statics::predict_const_dis_consist_vel_acc(
    Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const
{
  check_init_setup();
  // constant predictor : displacement in domain
  disnp.Update(1.0, *global_state().GetDisN(), 0.0);
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
  if (global_state().GetStepN() == 0) return false;

  // Displacement increment over last time step
  Teuchos::RCP<Epetra_Vector> disp_inc =
      Teuchos::rcp(new Epetra_Vector(*global_state().DofRowMapView(), true));
  disp_inc->Update((*global_state().GetDeltaTime())[0], *global_state().GetVelN(), 0.);
  // apply the dbc on the auxiliary vector
  tim_int().GetDBC().apply_dirichlet_to_vector(disp_inc);
  // update the solution variables
  disnp.Update(1.0, *global_state().GetDisN(), 0.0);
  disnp.Update(1.0, *disp_inc, 1.0);
  velnp.Update(1.0, *global_state().GetVelN(), 0.0);
  accnp.Update(1.0, *global_state().GetAccN(), 0.0);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::Statics::PredictConstAcc(
    Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const
{
  check_init_setup();
  // If there is not enough history information try a different predictor with
  // less requirements.
  if (global_state().GetStepN() < 2) return predict_const_vel_consist_acc(disnp, velnp, accnp);

  // Displacement increment over last time step
  Teuchos::RCP<Epetra_Vector> disp_inc =
      Teuchos::rcp(new Epetra_Vector(*global_state().DofRowMapView(), true));
  const double& dt = (*global_state().GetDeltaTime())[0];
  disp_inc->Update(dt, *global_state().GetVelN(), 0.);
  disp_inc->Update(0.5 * dt * dt, *global_state().GetAccN(), 1.0);
  // apply the dbc on the auxiliary vector
  tim_int().GetDBC().apply_dirichlet_to_vector(disp_inc);
  // update the solution variables
  disnp.Update(1.0, *global_state().GetDisN(), 0.0);
  disnp.Update(1., *disp_inc, 1.);
  velnp.Update(1.0, *global_state().GetVelN(), 0.0);
  accnp.Update(1.0, *global_state().GetAccN(), 0.0);

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
double STR::IMPLICIT::Statics::GetModelValue(const Epetra_Vector& x)
{
  Teuchos::RCP<const Epetra_Vector> disnp_ptr = global_state().ExtractDisplEntries(x);
  const Epetra_Vector& disnp = *disnp_ptr;

  set_state(disnp);

  EvalData().clear_values_for_all_energy_types();
  STR::MODELEVALUATOR::Structure& str_model =
      dynamic_cast<STR::MODELEVALUATOR::Structure&>(Evaluator(INPAR::STR::model_structure));

  str_model.determine_strain_energy(disnp, true);
  const double int_energy_np = EvalData().GetEnergyData(STR::internal_energy);
  double ext_energy_np = 0.0;
  global_state().GetFextNp()->Dot(disnp, &ext_energy_np);
  const double total = int_energy_np - ext_energy_np;

  std::ostream& os = IO::cout.os(IO::debug);
  os << __LINE__ << __PRETTY_FUNCTION__ << "\n";
  os << "internal/strain energy       = " << int_energy_np << "\n"
     << "external energy              = " << ext_energy_np << "\n";
  os << std::string(80, '-') << "\n";
  os << "Total                     = " << total << "\n";
  os << std::string(80, '-') << "\n";


  return total;
}

FOUR_C_NAMESPACE_CLOSE
