/*----------------------------------------------------------------------------*/
/*! \file
\brief different strategies for the update/correction of the regularization
parameter cn

\level 3

*/
/*----------------------------------------------------------------------------*/


#include "4C_contact_aug_penalty_update.hpp"

#include "4C_contact_aug_lagrange_multiplier_function.hpp"
#include "4C_contact_aug_potential.hpp"
#include "4C_contact_aug_steepest_ascent_strategy.hpp"
#include "4C_contact_paramsinterface.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_structure_new_model_evaluator_contact.hpp"
#include "4C_utils_epetra_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::Aug::PenaltyUpdate* CONTACT::Aug::PenaltyUpdate::Create(
    const Teuchos::ParameterList& sa_params)
{
  const Inpar::CONTACT::PenaltyUpdate update_type =
      Teuchos::getIntegralValue<Inpar::CONTACT::PenaltyUpdate>(sa_params, "PENALTY_UPDATE");

  return Create(update_type);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::Aug::PenaltyUpdate* CONTACT::Aug::PenaltyUpdate::Create(
    const Inpar::CONTACT::PenaltyUpdate update_type, const PenaltyUpdate* pu_src)
{
  switch (update_type)
  {
    case Inpar::CONTACT::PenaltyUpdate::sufficient_lin_reduction:
    {
      // call copy constructor
      if (pu_src) return new PenaltyUpdateSufficientLinReduction(*pu_src);
      return new PenaltyUpdateSufficientLinReduction();
    }
    case Inpar::CONTACT::PenaltyUpdate::sufficient_angle:
    {
      // call copy constructor
      if (pu_src) return new PenaltyUpdateSufficientAngle(*pu_src);
      return new PenaltyUpdateSufficientAngle();
    }
    case Inpar::CONTACT::PenaltyUpdate::none:
    {
      // call copy constructor
      if (pu_src) return new PenaltyUpdateEmpty(*pu_src);
      return new PenaltyUpdateEmpty();
    }
    case Inpar::CONTACT::PenaltyUpdate::vague:
    {
      FOUR_C_THROW(
          "You specified no PENALTY_UPDATE routine. Fix you input file!\n"
          "enum = %d | \"%s\"",
          update_type, Inpar::CONTACT::PenaltyUpdate2String(update_type).c_str());
      exit(EXIT_FAILURE);
    }
    default:
    {
      FOUR_C_THROW("Unknown penalty update type! (enum = %d | \"%s\")", update_type,
          Inpar::CONTACT::PenaltyUpdate2String(update_type).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::PenaltyUpdate::init(
    CONTACT::Aug::Strategy* const strategy, CONTACT::Aug::DataContainer* const data)
{
  strategy_ptr_ = strategy;
  data_ptr_ = data;
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::PenaltyUpdate::throw_if_not_initialized() const
{
  if (not isinit_) FOUR_C_THROW("Call init() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::PenaltyUpdate::set_state(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold, const Epetra_Vector& dir)
{
  throw_if_not_initialized();

  Core::IO::cout(Core::IO::debug) << std::string(40, '*') << Core::IO::endl;
  Core::IO::cout(Core::IO::debug) << __LINE__ << " -- " << CONTACT_FUNC_NAME << Core::IO::endl;
  Core::IO::cout(Core::IO::debug) << std::string(40, '*') << Core::IO::endl;

  state_.set(xold, dir, data());

  double dir_nrm2 = 0.0;
  dir.Norm2(&dir_nrm2);
  dir_norm2_ = dir_nrm2;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::PenaltyUpdate::post_update()
{
  PrintUpdate(Core::IO::cout.os(Core::IO::standard));
  reset();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::PenaltyUpdate::reset()
{
  state_.reset();
  status_ = Status::unevaluated;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::PenaltyUpdate::PrintInfo(std::ostream& os) const
{
  os << "\nCONTACT::Aug::PenaltyUpdate\n";
  os << "Type = " << Inpar::CONTACT::PenaltyUpdate2String(Type()) << "\n";
  os << "isinit_ = " << (isinit_ ? "TRUE" : "FALSE") << "\n";
  os << "strategy_ptr_ = " << strategy_ptr_ << "\n";
  os << "data_ptr_ = " << data_ptr_ << "\n";
  os << "dir_norm2_ = " << dir_norm2_ << "\n";
  os << "ratio_ = " << ratio_ << "\n";
  state_.Print(os);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::PenaltyUpdate::State::set(
    const Epetra_Vector& xold, const Epetra_Vector& dir, const CONTACT::Aug::DataContainer& data)
{
  xold_ = Teuchos::rcp(new Epetra_Vector(xold));
  full_direction_ = Teuchos::rcp(new Epetra_Vector(dir));
  wgap_ = Teuchos::rcp(new Epetra_Vector(*data.WGapPtr()));
  tributary_area_active_ = Teuchos::rcp(new Epetra_Vector(*data.KappaVecPtr()));
  tributary_area_inactive_ = Teuchos::rcp(new Epetra_Vector(*data.AVecPtr()));

  CONTACT::Aug::Potential& pot = pu_.data().Potential();
  pot.Compute();
  gn_gn_ = pot.get(POTENTIAL::Type::infeasibility_measure, POTENTIAL::SetType::active);
  pot.ComputeLin(dir);
  gn_dgn_ = pot.GetLin(POTENTIAL::Type::infeasibility_measure, POTENTIAL::SetType::active,
      POTENTIAL::LinTerm::wrt_d);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::PenaltyUpdate::State::reset()
{
  xold_ = Teuchos::null;
  full_direction_ = Teuchos::null;
  wgap_ = Teuchos::null;
  tributary_area_active_ = Teuchos::null;
  tributary_area_inactive_ = Teuchos::null;

  gn_gn_ = 0.0;
  gn_dgn_ = 0.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& CONTACT::Aug::PenaltyUpdate::State::GetDirection() const
{
  if (full_direction_.is_null()) FOUR_C_THROW("Set the state first!");

  return *full_direction_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::PenaltyUpdate::State::Print(std::ostream& os) const
{
  os << "--- CONTACT::Aug::PenaltyUpdate::State object\n";
  os << "    <RCP> fulldirection_:" << full_direction_ << ")\n";
  os << "    <RCP> xold_: " << xold_ << "\n";
  os << "    <RCP> wgap_: " << wgap_ << "\n";
  os << "    <RCP> tributary_area_active_: " << tributary_area_active_ << "\n";
  os << "    <RCP> tributary_area_inactive_: " << tributary_area_inactive_ << "\n";
  os << "    <double> gn_gn_: " << gn_gn_ << "\n";
  os << "    <double> gn_dgn_: " << gn_dgn_ << "\n" << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& CONTACT::Aug::PenaltyUpdate::State::get_previously_accepted_state() const
{
  if (xold_.is_null()) FOUR_C_THROW("Set the state first!");

  return *xold_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& CONTACT::Aug::PenaltyUpdate::State::GetWGap() const
{
  if (wgap_.is_null()) FOUR_C_THROW("Set the state first!");

  return *wgap_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& CONTACT::Aug::PenaltyUpdate::State::get_active_tributary_area() const
{
  if (tributary_area_active_.is_null()) FOUR_C_THROW("Set the state first!");

  return *tributary_area_active_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& CONTACT::Aug::PenaltyUpdate::State::get_inactive_tributary_area() const
{
  if (tributary_area_inactive_.is_null()) FOUR_C_THROW("Set the state first!");

  return *tributary_area_inactive_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::PenaltyUpdate::Update(const CONTACT::ParamsInterface& cparams)
{
  pre_update();

  status_ = execute(cparams);

  post_update();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::PenaltyUpdate::Decrease(const CONTACT::ParamsInterface& cparams)
{
  status_ = execute_decrease(cparams);

  post_decrease();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum CONTACT::Aug::PenaltyUpdate::Status CONTACT::Aug::PenaltyUpdate::execute_decrease(
    const CONTACT::ParamsInterface& cparams)
{
  Core::IO::cout(Core::IO::standard) << std::string(80, '!') << "\n"
                                     << "WARNING: The current PenaltyUpdate strategy does "
                                        "not support a decrease of the regularization parameter!\n"
                                     << std::string(80, '!') << Core::IO::endl;

  return Status::unevaluated;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::PenaltyUpdate::post_decrease()
{
  PrintUpdate(Core::IO::cout.os(Core::IO::standard));
  reset();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::Aug::PenaltyUpdate::get_d_gap_n(
    const Epetra_Vector& dincr_slma) const
{
  Teuchos::RCP<Epetra_Vector> dgapn_ptr =
      Teuchos::rcp(new Epetra_Vector(data().d_lm_nw_gap_lin_matrix_ptr()->RangeMap(), true));

  CATCH_EPETRA_ERROR(data().d_lm_nw_gap_lin_matrix_ptr()->Multiply(false, dincr_slma, *dgapn_ptr));

  return dgapn_ptr.getConst();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::Aug::PenaltyUpdate::get_inconsistent_d_gap_n(
    const Epetra_Vector& dincr_slma) const
{
  Teuchos::RCP<Epetra_Vector> dgapn_ptr =
      Teuchos::rcp(new Epetra_Vector(data().BMatrixPtr()->RangeMap(), true));

  CATCH_EPETRA_ERROR(data().BMatrixPtr()->Multiply(false, dincr_slma, *dgapn_ptr));

  return dgapn_ptr.getConst();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::Aug::PenaltyUpdate::get_problem_rhs(
    const CONTACT::ParamsInterface& cparams,
    const std::vector<Inpar::STR::ModelType>* without_these_models) const
{
  const STR::MODELEVALUATOR::Generic& model = cparams.get_model_evaluator();
  const STR::MODELEVALUATOR::Contact& cmodel =
      dynamic_cast<const STR::MODELEVALUATOR::Contact&>(model);

  return cmodel.assemble_force_of_models(without_these_models, true).getConst();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Core::LinAlg::SparseMatrix>
CONTACT::Aug::PenaltyUpdate::get_structural_stiffness_matrix(
    const CONTACT::ParamsInterface& cparams) const
{
  const STR::MODELEVALUATOR::Generic& model = cparams.get_model_evaluator();
  const STR::MODELEVALUATOR::Contact& cmodel =
      dynamic_cast<const STR::MODELEVALUATOR::Contact&>(model);

  // access the full stiffness matrix
  Teuchos::RCP<const Core::LinAlg::SparseMatrix> full_stiff_ptr =
      cmodel.get_jacobian_block(STR::MatBlockType::displ_displ);
  return full_stiff_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::PenaltyUpdate::PrintUpdate(std::ostream& os) const
{
  double cnmax = 0.0;
  data().Cn().MaxValue(&cnmax);

  os << "\n"
     << std::string(80, '=') << "\n"
     << "Update of the regularization parameter\n"
     << std::string(80, '=') << "\n";
  os << "Type   = " << Inpar::CONTACT::PenaltyUpdate2String(Type()) << "\n"
     << "Increase Parameter = " << data().SaData().get_penalty_correction_parameter() << "\n"
     << "Decrease Parameter = " << data().SaData().get_penalty_decrease_correction_parameter()
     << "\n"
     << "Status = " << status2_string(status_) << "\n"
     << "Ratio (cN_new / cn_old) = " << std::setw(10) << std::setprecision(4) << std::scientific
     << ratio() << "\n"
     << "New cN = " << std::setw(10) << std::setprecision(4) << std::scientific << cnmax << "\n"
     << std::string(80, '=') << "\n"
     << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Inpar::CONTACT::PenaltyUpdate CONTACT::Aug::PenaltyUpdateEmpty::Type() const
{
  return Inpar::CONTACT::PenaltyUpdate::none;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Inpar::CONTACT::PenaltyUpdate CONTACT::Aug::PenaltyUpdateSufficientLinReduction::Type() const
{
  return Inpar::CONTACT::PenaltyUpdate::sufficient_lin_reduction;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::PenaltyUpdateSufficientLinReduction::set_state(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold, const Epetra_Vector& dir)
{
  PenaltyUpdate::set_state(cparams, xold, dir);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum CONTACT::Aug::PenaltyUpdate::Status CONTACT::Aug::PenaltyUpdateSufficientLinReduction::execute(
    const CONTACT::ParamsInterface& cparams)
{
  const State& state = get_state();

  if (std::abs(state.gn_gn_) < 1.0e-30)
  {
    ratio() = 1.0;
    return Status::unchanged;
  }

  const double beta_theta_v = beta_theta();
  ratio() = 1.0 / (1.0 - beta_theta_v) * (1.0 + state.gn_dgn_ / state.gn_gn_);

  if (ratio() > 1.0) data().Cn().Scale(ratio());

  return (ratio() > 1.0 ? Status::increased : Status::unchanged);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum CONTACT::Aug::PenaltyUpdate::Status
CONTACT::Aug::PenaltyUpdateSufficientLinReduction::execute_decrease(
    const CONTACT::ParamsInterface& cparams)
{
  const State& state = get_state();
  const double beta_theta_decrease_v = beta_theta_decrease();

  if (std::abs(state.gn_gn_) < 1.0e-30 or beta_theta_decrease_v >= 1.0)
  {
    ratio() = 1.0;
    return Status::unchanged;
  }

  ratio() = 1.0 / (1.0 - beta_theta_decrease_v) * (1.0 + state.gn_dgn_ / state.gn_gn_);

  if (ratio() > 0.0 and ratio() < 1.0)
  {
    data().Cn().Scale(ratio());
    return Status::decreased;
  }

  return Status::unchanged;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::PenaltyUpdateSufficientLinReduction::beta_theta() const
{
  return data().SaData().get_penalty_correction_parameter();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::PenaltyUpdateSufficientLinReduction::beta_theta_decrease() const
{
  return data().SaData().get_penalty_decrease_correction_parameter();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Inpar::CONTACT::PenaltyUpdate CONTACT::Aug::PenaltyUpdateSufficientAngle::Type() const
{
  return Inpar::CONTACT::PenaltyUpdate::sufficient_angle;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::PenaltyUpdateSufficientAngle::set_state(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold, const Epetra_Vector& dir)
{
  PenaltyUpdate::set_state(cparams, xold, dir);

  if (std::abs(get_state().gn_gn_) < 1.0e-30)
  {
    ratio() = 1.0;
    return;
  }

  // split the previously accepted state into its different parts
  Teuchos::RCP<Epetra_Vector> displ_slma, zn_active, zn_inactive;
  strategy().SplitStateVector(
      get_state().get_previously_accepted_state(), displ_slma, zn_active, zn_inactive);

  // split the direction into its different parts
  Teuchos::RCP<Epetra_Vector> dincr_slma, znincr_active, znincr_inactive;
  strategy().SplitStateVector(
      get_state().GetDirection(), dincr_slma, znincr_active, znincr_inactive);

  Epetra_Vector tmp(dincr_slma->Map(), true);
  data().DGLmLinMatrix().Multiply(false, *dincr_slma, tmp);
  double d_ddglm_d = 0.0;
  tmp.Dot(*dincr_slma, &d_ddglm_d);

  // directional derivative of the structural gradient
  const STR::MODELEVALUATOR::Contact& cmodel =
      dynamic_cast<const STR::MODELEVALUATOR::Contact&>(cparams.get_model_evaluator());

  const std::vector<Inpar::STR::ModelType> without_contact(1, cmodel.Type());
  Teuchos::RCP<Epetra_Vector> str_gradient =
      cmodel.assemble_force_of_models(&without_contact, true);

  Epetra_Vector dincr_str(str_gradient->Map());
  Core::LinAlg::ExtractMyVector(dir, dincr_str);

  double dstr_grad = 0.0;
  CATCH_EPETRA_ERROR(str_gradient->Dot(dincr_str, &dstr_grad));

  // directional derivative of the active constraint gradient
  Teuchos::RCP<const Epetra_Vector> dgapn_ptr = get_inconsistent_d_gap_n(*dincr_slma);
  Epetra_Vector dgapn_active(*data().g_active_n_dof_row_map_ptr());
  Core::LinAlg::ExtractMyVector(*dgapn_ptr, dgapn_active);

  double dgapn_zn = 0.0;
  CATCH_EPETRA_ERROR(dgapn_active.Dot(*zn_active, &dgapn_zn));

  Epetra_Vector sc_dgapn(dgapn_active);
  MultiplyElementwise(
      get_state().get_active_tributary_area(), data().GActiveNodeRowMap(), sc_dgapn, true);

  double sc_dgapn_dgapn = 0.0;
  CATCH_EPETRA_ERROR(sc_dgapn.Dot(dgapn_active, &sc_dgapn_dgapn));

  Epetra_Vector awgapn(*get_state().wgap_);
  MultiplyElementwise(
      get_state().get_active_tributary_area(), data().GActiveNodeRowMap(), awgapn, true);

  double awgapn_nrm2 = 0.0;
  awgapn.Norm2(&awgapn_nrm2);

  double dgapn_nrm2 = 0.0;
  dgapn_active.Norm2(&dgapn_nrm2);

  double cn_old = 0.0;
  data().CnPtr()->MaxValue(&cn_old);

  const double gamma_phi = beta_angle();

  const double cn_new =
      (cn_old * gamma_phi * dgapn_nrm2 * awgapn_nrm2 + d_ddglm_d - dstr_grad + dgapn_zn) /
      sc_dgapn_dgapn;

  ratio() = cn_new / cn_old;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum CONTACT::Aug::PenaltyUpdate::Status CONTACT::Aug::PenaltyUpdateSufficientAngle::execute(
    const CONTACT::ParamsInterface& cparams)
{
  if (ratio() > 1.0)
  {
    data().Cn().Scale(ratio());
    return Status::increased;
  }
  return Status::unchanged;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::PenaltyUpdateSufficientAngle::beta_angle() const
{
  return data().SaData().get_penalty_correction_parameter();
}

FOUR_C_NAMESPACE_CLOSE
