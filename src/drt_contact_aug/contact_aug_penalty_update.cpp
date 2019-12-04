/*----------------------------------------------------------------------------*/
/*! \file
\brief different strategies for the update/correction of the regularization
parameter cn

\level 3

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------------*/


#include "contact_aug_penalty_update.H"
#include "contact_aug_potential.H"
#include "contact_aug_steepest_ascent_strategy.H"
#include "contact_aug_lagrange_multiplier_function.H"
#include "../drt_contact/contact_paramsinterface.H"

#include "../drt_structure_new/str_model_evaluator_contact.H"
#include "../drt_inpar/inpar_structure.H"

#include "../drt_inpar/inpar_contact.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../drt_io/io_pstream.H"
#include "../linalg/linalg_multiply.H"
#include "../drt_lib/epetra_utils.H"

#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::PenaltyUpdate* CONTACT::AUG::PenaltyUpdate::Create(
    const Teuchos::ParameterList& sa_params)
{
  const INPAR::CONTACT::PenaltyUpdate update_type =
      Teuchos::getIntegralValue<INPAR::CONTACT::PenaltyUpdate>(sa_params, "PENALTY_UPDATE");

  return Create(update_type);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::PenaltyUpdate* CONTACT::AUG::PenaltyUpdate::Create(
    const INPAR::CONTACT::PenaltyUpdate update_type, const PenaltyUpdate* pu_src)
{
  switch (update_type)
  {
    case INPAR::CONTACT::PenaltyUpdate::sufficient_lin_reduction:
    {
      // call copy constructor
      if (pu_src) return new PenaltyUpdate_SufficientLinReduction(*pu_src);
      return new PenaltyUpdate_SufficientLinReduction();
    }
    case INPAR::CONTACT::PenaltyUpdate::sufficient_angle:
    {
      // call copy constructor
      if (pu_src) return new PenaltyUpdate_SufficientAngle(*pu_src);
      return new PenaltyUpdate_SufficientAngle();
    }
    case INPAR::CONTACT::PenaltyUpdate::none:
    {
      // call copy constructor
      if (pu_src) return new PenaltyUpdate_Empty(*pu_src);
      return new PenaltyUpdate_Empty();
    }
    case INPAR::CONTACT::PenaltyUpdate::vague:
    {
      dserror(
          "You specified no PENALTY_UPDATE routine. Fix you input file!\n"
          "enum = %d | \"%s\"",
          update_type, INPAR::CONTACT::PenaltyUpdate2String(update_type).c_str());
      exit(EXIT_FAILURE);
    }
    default:
    {
      dserror("Unknown penalty update type! (enum = %d | \"%s\")", update_type,
          INPAR::CONTACT::PenaltyUpdate2String(update_type).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::Init(
    CONTACT::AUG::Strategy* const strategy, CONTACT::AUG::DataContainer* const data)
{
  strategy_ptr_ = strategy;
  data_ptr_ = data;
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::ThrowIfNotInitialized() const
{
  if (not isinit_) dserror("Call Init() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::SetState(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold, const Epetra_Vector& dir)
{
  ThrowIfNotInitialized();

  IO::cout(IO::debug) << std::string(40, '*') << IO::endl;
  IO::cout(IO::debug) << __LINE__ << " -- " << CONTACT_FUNC_NAME << IO::endl;
  IO::cout(IO::debug) << std::string(40, '*') << IO::endl;

  state_.Set(xold, dir, Data());

  double dir_nrm2 = 0.0;
  dir.Norm2(&dir_nrm2);
  dir_norm2_ = dir_nrm2;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::PostUpdate()
{
  PrintUpdate(IO::cout.os(IO::standard));
  Reset();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::Reset()
{
  state_.Reset();
  status_ = Status::unevaluated;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::PrintInfo(std::ostream& os) const
{
  os << "\nCONTACT::AUG::PenaltyUpdate\n";
  os << "Type = " << INPAR::CONTACT::PenaltyUpdate2String(Type()) << "\n";
  os << "isinit_ = " << (isinit_ ? "TRUE" : "FALSE") << "\n";
  os << "strategy_ptr_ = " << strategy_ptr_ << "\n";
  os << "data_ptr_ = " << data_ptr_ << "\n";
  os << "dir_norm2_ = " << dir_norm2_ << "\n";
  os << "ratio_ = " << ratio_ << "\n";
  state_.Print(os);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::State::Set(
    const Epetra_Vector& xold, const Epetra_Vector& dir, const CONTACT::AUG::DataContainer& data)
{
  xold_ = Teuchos::rcp(new Epetra_Vector(xold));
  full_direction_ = Teuchos::rcp(new Epetra_Vector(dir));
  wgap_ = Teuchos::rcp(new Epetra_Vector(*data.WGapPtr()));
  tributary_area_active_ = Teuchos::rcp(new Epetra_Vector(*data.KappaVecPtr()));
  tributary_area_inactive_ = Teuchos::rcp(new Epetra_Vector(*data.AVecPtr()));

  CONTACT::AUG::Potential& pot = pu_.Data().Potential();
  pot.Compute();
  gn_gn_ = pot.Get(POTENTIAL::Type::infeasibility_measure, POTENTIAL::SetType::active);
  pot.ComputeLin(dir);
  gn_dgn_ = pot.GetLin(POTENTIAL::Type::infeasibility_measure, POTENTIAL::SetType::active,
      POTENTIAL::LinTerm::wrt_d);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::State::Reset()
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
const Epetra_Vector& CONTACT::AUG::PenaltyUpdate::State::GetDirection() const
{
  if (full_direction_.is_null()) dserror("Set the state first!");

  return *full_direction_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::State::Print(std::ostream& os) const
{
  os << "--- CONTACT::AUG::PenaltyUpdate::State object\n";
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
const Epetra_Vector& CONTACT::AUG::PenaltyUpdate::State::GetPreviouslyAcceptedState() const
{
  if (xold_.is_null()) dserror("Set the state first!");

  return *xold_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& CONTACT::AUG::PenaltyUpdate::State::GetWGap() const
{
  if (wgap_.is_null()) dserror("Set the state first!");

  return *wgap_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& CONTACT::AUG::PenaltyUpdate::State::GetActiveTributaryArea() const
{
  if (tributary_area_active_.is_null()) dserror("Set the state first!");

  return *tributary_area_active_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& CONTACT::AUG::PenaltyUpdate::State::GetInactiveTributaryArea() const
{
  if (tributary_area_inactive_.is_null()) dserror("Set the state first!");

  return *tributary_area_inactive_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::Update(const CONTACT::ParamsInterface& cparams)
{
  PreUpdate();

  status_ = Execute(cparams);

  PostUpdate();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::Decrease(const CONTACT::ParamsInterface& cparams)
{
  status_ = ExecuteDecrease(cparams);

  PostDecrease();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum CONTACT::AUG::PenaltyUpdate::Status CONTACT::AUG::PenaltyUpdate::ExecuteDecrease(
    const CONTACT::ParamsInterface& cparams)
{
  IO::cout(IO::standard) << std::string(80, '!') << "\n"
                         << "WARNING: The current PenaltyUpdate strategy does "
                            "not support a decrease of the regularization parameter!\n"
                         << std::string(80, '!') << IO::endl;

  return Status::unevaluated;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::PostDecrease()
{
  PrintUpdate(IO::cout.os(IO::standard));
  Reset();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::AUG::PenaltyUpdate::Get_DGapN(
    const Epetra_Vector& dincr_slma) const
{
  Teuchos::RCP<Epetra_Vector> dgapn_ptr =
      Teuchos::rcp(new Epetra_Vector(Data().DLmNWGapLinMatrixPtr()->RangeMap(), true));

  CATCH_EPETRA_ERROR(Data().DLmNWGapLinMatrixPtr()->Multiply(false, dincr_slma, *dgapn_ptr));

  return dgapn_ptr.getConst();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::AUG::PenaltyUpdate::Get_inconsistent_DGapN(
    const Epetra_Vector& dincr_slma) const
{
  Teuchos::RCP<Epetra_Vector> dgapn_ptr =
      Teuchos::rcp(new Epetra_Vector(Data().BMatrixPtr()->RangeMap(), true));

  CATCH_EPETRA_ERROR(Data().BMatrixPtr()->Multiply(false, dincr_slma, *dgapn_ptr));

  return dgapn_ptr.getConst();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::AUG::PenaltyUpdate::GetProblemRhs(
    const CONTACT::ParamsInterface& cparams,
    const std::vector<INPAR::STR::ModelType>* without_these_models) const
{
  const STR::MODELEVALUATOR::Generic& model = cparams.GetModelEvaluator();
  const STR::MODELEVALUATOR::Contact& cmodel =
      dynamic_cast<const STR::MODELEVALUATOR::Contact&>(model);

  return cmodel.AssembleForceOfModels(without_these_models, true).getConst();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::SparseMatrix> CONTACT::AUG::PenaltyUpdate::GetStructuralStiffnessMatrix(
    const CONTACT::ParamsInterface& cparams) const
{
  const STR::MODELEVALUATOR::Generic& model = cparams.GetModelEvaluator();
  const STR::MODELEVALUATOR::Contact& cmodel =
      dynamic_cast<const STR::MODELEVALUATOR::Contact&>(model);

  // access the full stiffness matrix
  Teuchos::RCP<const LINALG::SparseMatrix> full_stiff_ptr =
      cmodel.GetJacobianBlock(DRT::UTILS::block_displ_displ);
  return full_stiff_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate::PrintUpdate(std::ostream& os) const
{
  double cnmax = 0.0;
  Data().Cn().MaxValue(&cnmax);

  os << "\n"
     << std::string(80, '=') << "\n"
     << "Update of the regularization parameter\n"
     << std::string(80, '=') << "\n";
  os << "Type   = " << INPAR::CONTACT::PenaltyUpdate2String(Type()) << "\n"
     << "Increase Parameter = " << Data().SaData().GetPenaltyCorrectionParameter() << "\n"
     << "Decrease Parameter = " << Data().SaData().GetPenaltyDecreaseCorrectionParameter() << "\n"
     << "Status = " << Status2String(status_) << "\n"
     << "Ratio (cN_new / cn_old) = " << std::setw(10) << std::setprecision(4) << std::scientific
     << Ratio() << "\n"
     << "New cN = " << std::setw(10) << std::setprecision(4) << std::scientific << cnmax << "\n"
     << std::string(80, '=') << "\n"
     << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::CONTACT::PenaltyUpdate CONTACT::AUG::PenaltyUpdate_Empty::Type() const
{
  return INPAR::CONTACT::PenaltyUpdate::none;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::CONTACT::PenaltyUpdate CONTACT::AUG::PenaltyUpdate_SufficientLinReduction::Type() const
{
  return INPAR::CONTACT::PenaltyUpdate::sufficient_lin_reduction;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate_SufficientLinReduction::SetState(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold, const Epetra_Vector& dir)
{
  PenaltyUpdate::SetState(cparams, xold, dir);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum CONTACT::AUG::PenaltyUpdate::Status
CONTACT::AUG::PenaltyUpdate_SufficientLinReduction::Execute(const CONTACT::ParamsInterface& cparams)
{
  const State& state = GetState();

  if (std::abs(state.gn_gn_) < 1.0e-30)
  {
    Ratio() = 1.0;
    return Status::unchanged;
  }

  const double beta_theta = BetaTheta();
  Ratio() = 1.0 / (1.0 - beta_theta) * (1.0 + state.gn_dgn_ / state.gn_gn_);

  if (Ratio() > 1.0) Data().Cn().Scale(Ratio());

  return (Ratio() > 1.0 ? Status::increased : Status::unchanged);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum CONTACT::AUG::PenaltyUpdate::Status
CONTACT::AUG::PenaltyUpdate_SufficientLinReduction::ExecuteDecrease(
    const CONTACT::ParamsInterface& cparams)
{
  const State& state = GetState();
  const double beta_theta_decrease = BetaThetaDecrease();

  if (std::abs(state.gn_gn_) < 1.0e-30 or beta_theta_decrease >= 1.0)
  {
    Ratio() = 1.0;
    return Status::unchanged;
  }

  Ratio() = 1.0 / (1.0 - beta_theta_decrease) * (1.0 + state.gn_dgn_ / state.gn_gn_);

  if (Ratio() > 0.0 and Ratio() < 1.0)
  {
    Data().Cn().Scale(Ratio());
    return Status::decreased;
  }

  return Status::unchanged;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::PenaltyUpdate_SufficientLinReduction::BetaTheta() const
{
  return Data().SaData().GetPenaltyCorrectionParameter();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::PenaltyUpdate_SufficientLinReduction::BetaThetaDecrease() const
{
  return Data().SaData().GetPenaltyDecreaseCorrectionParameter();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::CONTACT::PenaltyUpdate CONTACT::AUG::PenaltyUpdate_SufficientAngle::Type() const
{
  return INPAR::CONTACT::PenaltyUpdate::sufficient_angle;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::PenaltyUpdate_SufficientAngle::SetState(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold, const Epetra_Vector& dir)
{
  PenaltyUpdate::SetState(cparams, xold, dir);

  if (std::abs(GetState().gn_gn_) < 1.0e-30)
  {
    Ratio() = 1.0;
    return;
  }

  // split the previously accepted state into its different parts
  Teuchos::RCP<Epetra_Vector> displ_slma, zn_active, zn_inactive;
  Strategy().SplitStateVector(
      GetState().GetPreviouslyAcceptedState(), displ_slma, zn_active, zn_inactive);

  // split the direction into its different parts
  Teuchos::RCP<Epetra_Vector> dincr_slma, znincr_active, znincr_inactive;
  Strategy().SplitStateVector(
      GetState().GetDirection(), dincr_slma, znincr_active, znincr_inactive);

  Epetra_Vector tmp(dincr_slma->Map(), true);
  Data().DGLmLinMatrix().Multiply(false, *dincr_slma, tmp);
  double d_ddglm_d = 0.0;
  tmp.Dot(*dincr_slma, &d_ddglm_d);

  // directional derivative of the structural gradient
  const STR::MODELEVALUATOR::Contact& cmodel =
      dynamic_cast<const STR::MODELEVALUATOR::Contact&>(cparams.GetModelEvaluator());

  const std::vector<INPAR::STR::ModelType> without_contact(1, cmodel.Type());
  Teuchos::RCP<Epetra_Vector> str_gradient = cmodel.AssembleForceOfModels(&without_contact, true);

  Epetra_Vector dincr_str(str_gradient->Map());
  LINALG::ExtractMyVector(dir, dincr_str);

  double dstr_grad = 0.0;
  CATCH_EPETRA_ERROR(str_gradient->Dot(dincr_str, &dstr_grad));

  // directional derivative of the active constraint gradient
  Teuchos::RCP<const Epetra_Vector> dgapn_ptr = Get_inconsistent_DGapN(*dincr_slma);
  Epetra_Vector dgapn_active(*Data().GActiveNDofRowMapPtr());
  LINALG::ExtractMyVector(*dgapn_ptr, dgapn_active);

  double dgapn_zn = 0.0;
  CATCH_EPETRA_ERROR(dgapn_active.Dot(*zn_active, &dgapn_zn));

  Epetra_Vector sc_dgapn(dgapn_active);
  MultiplyElementwise(
      GetState().GetActiveTributaryArea(), Data().GActiveNodeRowMap(), sc_dgapn, true);

  double sc_dgapn_dgapn = 0.0;
  CATCH_EPETRA_ERROR(sc_dgapn.Dot(dgapn_active, &sc_dgapn_dgapn));

  Epetra_Vector awgapn(*GetState().wgap_);
  MultiplyElementwise(
      GetState().GetActiveTributaryArea(), Data().GActiveNodeRowMap(), awgapn, true);

  double awgapn_nrm2 = 0.0;
  awgapn.Norm2(&awgapn_nrm2);

  double dgapn_nrm2 = 0.0;
  dgapn_active.Norm2(&dgapn_nrm2);

  double cn_old = 0.0;
  Data().CnPtr()->MaxValue(&cn_old);

  const double gamma_phi = BetaAngle();

  const double cn_new =
      (cn_old * gamma_phi * dgapn_nrm2 * awgapn_nrm2 + d_ddglm_d - dstr_grad + dgapn_zn) /
      sc_dgapn_dgapn;

  Ratio() = cn_new / cn_old;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum CONTACT::AUG::PenaltyUpdate::Status CONTACT::AUG::PenaltyUpdate_SufficientAngle::Execute(
    const CONTACT::ParamsInterface& cparams)
{
  if (Ratio() > 1.0)
  {
    Data().Cn().Scale(Ratio());
    return Status::increased;
  }
  return Status::unchanged;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::PenaltyUpdate_SufficientAngle::BetaAngle() const
{
  return Data().SaData().GetPenaltyCorrectionParameter();
}
