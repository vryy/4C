/*---------------------------------------------------------------------*/
/*!
\file contact_aug_combo_strategy_switching.cpp

\brief This file contains the implementation of the switching strategy
       for the AUG::ComboStrategy.

\level 3

\maintainer Michael Hiermeier

\date Mar 24, 2017

*/
/*---------------------------------------------------------------------*/


#include "contact_aug_combo_strategy.H"
#include "contact_augmented_strategy.H"
#include "contact_augmented_interface.H"

#include "../drt_contact/contact_paramsinterface.H"
#include "../drt_contact/contact_strategy_factory.H"

#include "../drt_io/io_pstream.H"

#include "../linalg/linalg_utils.H"

#include "../drt_structure_new/str_model_evaluator_contact.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::AUG::ComboStrategy::Switching> CONTACT::AUG::ComboStrategy::Switching::Create(
    ComboStrategy& combo)
{
  const Teuchos::ParameterList& p_combo = combo.Params().sublist("AUGMENTED").sublist("COMBO");

  const enum INPAR::CONTACT::SwitchingStrategy switch_type =
      DRT::INPUT::IntegralValue<enum INPAR::CONTACT::SwitchingStrategy>(
          p_combo, "SWITCHING_STRATEGY");

  switch (switch_type)
  {
    case INPAR::CONTACT::switch_preasymptotic:
      return Teuchos::rcp(new PreAsymptoticSwitching(combo, p_combo));
    default:
      dserror("Unknown switching strategy! (switch_type = %d)", switch_type);
      exit(EXIT_FAILURE);
  }

  dserror("Impossible to reach this point!");
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::ComboStrategy::Switching::Switching(
    ComboStrategy& combo, const Teuchos::ParameterList& p_combo)
    : combo_(combo), strat_types_(0)
{
  GetStrategyTypes(combo_.strategies_, strat_types_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::Switching::GetStrategyTypes(
    const plain_strategy_set& strategies, plain_strattype_set& strat_types) const
{
  for (plain_strategy_set::const_iterator cit = strategies.begin(); cit != strategies.end(); ++cit)
  {
    const CONTACT::CoAbstractStrategy& s = (**cit);

    switch (s.Type())
    {
      case INPAR::CONTACT::solution_augmented:
      case INPAR::CONTACT::solution_steepest_ascent:
      case INPAR::CONTACT::solution_steepest_ascent_sp:
      case INPAR::CONTACT::solution_std_lagrange:
        strat_types.push_back(s.Type());
        break;
      default:
        dserror(
            "The strategy is of a non-supported type! ( type = "
            "%s | %d )",
            INPAR::CONTACT::SolvingStrategy2String(s.Type()).c_str(), s.Type());
        exit(EXIT_FAILURE);
    }
  }

  if (strategies.size() != strat_types.size()) dserror("Size mismatch! Something went wrong.");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
unsigned CONTACT::AUG::ComboStrategy::Switching::Id(
    enum INPAR::CONTACT::SolvingStrategy sol_type) const
{
  return FindId(sol_type);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
unsigned CONTACT::AUG::ComboStrategy::Switching::FindId(
    INPAR::CONTACT::SolvingStrategy sol_type) const
{
  unsigned id = 0;
  for (plain_strattype_set::const_iterator cit = strat_types_.begin(); cit != strat_types_.end();
       ++cit)
  {
    if (*cit == sol_type) return id;
    ++id;
  }

  dserror("Couldn't find the given SolvingStrategy! (sol_type = %s | %d)",
      INPAR::CONTACT::SolvingStrategy2String(sol_type).c_str(), sol_type);
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::ComboStrategy::PreAsymptoticSwitching::PreAsymptoticSwitching(
    ComboStrategy& combo, const Teuchos::ParameterList& p_combo)
    : Switching(combo, p_combo),
      preasymptotic_id_(0),
      asymptotic_id_(0),
      is_asymptotic_(false),
      maxabsawgap_(*this)
{
  if (combo_.strategies_.size() > 2)
    dserror(
        "This basic switching strategy supports maximal a number of "
        "two strategies. Feel free to add a new switching strategy, if you "
        "need more.");

  const enum INPAR::CONTACT::SolvingStrategy preasymptotic =
      DRT::INPUT::IntegralValue<enum INPAR::CONTACT::SolvingStrategy>(p_combo, "STRATEGY_0");
  preasymptotic_id_ = FindId(preasymptotic);

  const enum INPAR::CONTACT::SolvingStrategy asymptotic =
      DRT::INPUT::IntegralValue<enum INPAR::CONTACT::SolvingStrategy>(p_combo, "STRATEGY_1");
  asymptotic_id_ = FindId(asymptotic);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
unsigned CONTACT::AUG::ComboStrategy::PreAsymptoticSwitching::Id() const
{
  if (is_asymptotic_) return asymptotic_id_;

  return preasymptotic_id_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::PreAsymptoticSwitching::Update(
    CONTACT::ParamsInterface& cparams, std::ostream& os)
{
  PrintUpdateHead(os);

  const bool is_predict = cparams.IsPredictor();
  const bool check_pen = CheckPenetration(os);
  bool is_asymptotic =
      (not is_predict and check_pen and (CheckResidual(cparams, os) or CheckCnBound(os)));

  // if the status is the same as before, do nothing
  if (is_asymptotic == is_asymptotic_) return;

  // --------------------------------------------------------------------------
  /* switch back to the pre-asymptotic phase:
   *
   * Only possible, if the penetration criterion is not fulfilled or the
   * asymptotic penetration is two times larger as the pre-asymptotic
   * penetration.
   *
   * We switch also back to the pre-asymptotic phase at the beginning of a new
   * time/load step ( prediction phase ). */
  // --------------------------------------------------------------------------
  if (not is_asymptotic and
      (is_predict or
          (not check_pen or maxabsawgap_.asymptotic_ > 2.0 * maxabsawgap_.pre_asymptotic_)))
  {
    os << "Switching back to the pre-asymptotic phase since";
    if (is_predict)
      os << " a new time/load step starts... \n";
    else if (not check_pen)
      os << " the penetration bound criterion is hurt... \n";
    else
      os << " the asymptotic constraint violation exceeds the pre-asymptotic"
            " by at least a factor of two...\n";

    is_asymptotic_ = false;
    STRATEGY::Factory::PrintStrategyBanner(strat_types_[Id()]);
  }
  // --------------------------------------------------------------------------
  // switch to the asymptotic phase
  // --------------------------------------------------------------------------
  else if (is_asymptotic)
  {
    is_asymptotic_ = true;
    STRATEGY::Factory::PrintStrategyBanner(strat_types_[Id()]);
  }

#ifdef DEBUG_COMBO_STRATEGY
  std::cout << __LINE__ << " -- " << __PRETTY_FUNCTION__ << std::endl;
  std::cout << "is_asymptotic_ = " << (is_asymptotic_ ? "TRUE" : "FALSE") << "\n";
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::PreAsymptoticSwitching::CheckResidual(
    CONTACT::ParamsInterface& cparams, std::ostream& os)
{
  Teuchos::RCP<Epetra_Vector> force_no_dbc_ptr = GetStructuralForceWithoutDbcDofs(cparams);

  Teuchos::RCP<Epetra_Vector> str_slmaforce = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> constr_slmaforce = Teuchos::null;
  GetActiveSlMaForces(*force_no_dbc_ptr, str_slmaforce, constr_slmaforce);

  const bool angle_check =
      CheckAngleBetweenStrForceAndContactForce(*str_slmaforce, *constr_slmaforce, os);
  const bool res_check = CheckContactResidualNorm(*str_slmaforce, *constr_slmaforce, os);

  return (res_check and angle_check);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::PreAsymptoticSwitching::CheckContactResidualNorm(
    const Epetra_Vector& str_slmaforce, const Epetra_Vector& constr_slmaforce,
    std::ostream& os) const
{
  double str_nrm = 0.0;
  double res_nrm = 0.0;

#ifdef DEBUG_COMBO_STRATEGY
  force.Norm2(&res_nrm);
  std::cout << __LINE__ << " -- " << __PRETTY_FUNCTION__ << std::endl;
  std::cout << "res_nrm = " << res_nrm << "\n";
#endif

  Epetra_Vector res_slma_active(str_slmaforce.Map(), false);
  EPETRA_CHK_ERR(res_slma_active.Update(1.0, str_slmaforce, -1.0, constr_slmaforce, 0.0));

  res_slma_active.Norm2(&res_nrm);

  str_slmaforce.Norm2(&str_nrm);
  const bool check_res = (res_nrm <= 1.0e-3 * str_nrm);

  os << std::setfill('.') << std::setw(80)
     << "# Checking residual between str. force and contact force "
     << (check_res ? "SUCCEEDED" : "FAILED") << "\n"
     << std::setfill(' ') << "  rel. sl/ma residual = " << std::setw(14) << std::setprecision(4)
     << std::scientific << res_nrm << " <= " << 1.0e-3 * str_nrm << "\n";

  return check_res;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::PreAsymptoticSwitching::CheckAngleBetweenStrForceAndContactForce(
    const Epetra_Vector& str_slmaforce, const Epetra_Vector& constr_slmaforce,
    std::ostream& os) const
{
  double constr_slmaforce_nrm = 0.0;
  constr_slmaforce.Norm2(&constr_slmaforce_nrm);

  double str_slmaforce_nrm = 0.0;
  str_slmaforce.Norm2(&str_slmaforce_nrm);

  const double nrm_prod = constr_slmaforce_nrm * str_slmaforce_nrm;

  double inner_prod = 0.0;
  constr_slmaforce.Dot(str_slmaforce, &inner_prod);

  const double angle_tol = 1.0e-6;

  const bool check_angle = (nrm_prod - inner_prod <= nrm_prod * angle_tol);

  static const double conv_to_deg = 180 / (std::atan(1.0) * 4.0);
  static const double angle_bound = acos(1.0 - angle_tol) * conv_to_deg;

  os << std::setfill('.') << std::setw(80)
     << "# Checking angle between str. force and contact force "
     << (check_angle ? "SUCCEEDED" : "FAILED") << "\n"
     << std::setfill(' ') << "  angle: 0.0 <= " << std::setw(14) << std::setprecision(4)
     << std::scientific;
  if (nrm_prod > 0.0)
  {
    const double cosine = inner_prod / nrm_prod;
    if (std::abs(cosine) <= 1.0)
      os << acos(cosine) * conv_to_deg;
    else
      os << 0.0;
  }
  else
    os << "undef.";
  os << " <= " << angle_bound << " [deg]\n";

  return check_angle;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::PreAsymptoticSwitching::CheckCnBound(std::ostream& os) const
{
  static const double cn_bound = 1.0e+15;
  double max_cn = 0.0;
  combo_.data_.Cn().MaxValue(&max_cn);

  const bool check_bound = (max_cn > cn_bound);

  os << std::setfill('.') << std::setw(80)
     << "# Checking upper bound of the regularization parameter "
     << (check_bound ? "SUCCEEDED" : "FAILED") << "\n"
     << std::setfill(' ') << "  current cn [" << std::setw(14) << std::setprecision(4)
     << std::scientific << max_cn << "] > bound [" << std::setw(14) << std::setprecision(4)
     << std::scientific << cn_bound << "]\n";

  return check_bound;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::PreAsymptoticSwitching::GetActiveSlMaForces(
    const Epetra_Vector& str_force, Teuchos::RCP<Epetra_Vector>& str_slmaforce,
    Teuchos::RCP<Epetra_Vector>& constr_slmaforce) const
{
  Epetra_Vector& slforce = *combo_.no_dbc_.slForce_;
  Epetra_Vector& maforce = *combo_.no_dbc_.maForce_;

  slforce.PutScalar(0.0);
  maforce.PutScalar(0.0);

  LINALG::ExtractMyVector(*combo_.data_.SlForceLmPtr(), slforce);
  LINALG::ExtractMyVector(*combo_.data_.MaForceLmPtr(), maforce);

  Teuchos::RCP<Epetra_Map> gSlActiveForceMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> gMaActiveForceMap = Teuchos::null;
  GetGlobalSlMaActiveForceMaps(slforce, maforce, gSlActiveForceMap, gMaActiveForceMap);
  Teuchos::RCP<Epetra_Map> gSlMaActiveForceMap =
      LINALG::MergeMap(*gSlActiveForceMap, *gMaActiveForceMap);

  Epetra_Vector& slmaforce = *combo_.no_dbc_.slMaForce_;
  slmaforce.PutScalar(0.0);

  LINALG::AssembleMyVector(1.0, slmaforce, 1.0, slforce);
  LINALG::AssembleMyVector(1.0, slmaforce, 1.0, maforce);

  constr_slmaforce = LINALG::CreateVector(*gSlMaActiveForceMap, true);
  LINALG::ExtractMyVector(slmaforce, *constr_slmaforce);
  // consider time integration factor
  constr_slmaforce->Scale(1.0 - combo_.data_.AlphaF());

  str_slmaforce = LINALG::CreateVector(*gSlMaActiveForceMap, true);
  LINALG::ExtractMyVector(str_force, *str_slmaforce);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::PreAsymptoticSwitching::GetGlobalSlMaActiveForceMaps(
    const Epetra_Vector& slforce, const Epetra_Vector& maforce,
    Teuchos::RCP<Epetra_Map>& gSlActiveForceMap, Teuchos::RCP<Epetra_Map>& gMaActiveForceMap) const
{
  // initialize the map pointers
  gSlActiveForceMap = Teuchos::rcp(new Epetra_Map(0, 0, *combo_.data_.CommPtr()));
  gMaActiveForceMap = Teuchos::rcp(new Epetra_Map(0, 0, *combo_.data_.CommPtr()));

  Teuchos::RCP<Epetra_Map> imap = Teuchos::null;

  for (const Teuchos::RCP<CONTACT::CoInterface>& cit : combo_.Get().Interfaces())
  {
    const CONTACT::AUG::Interface& interface = dynamic_cast<AUG::Interface&>(*cit);

    imap = Teuchos::null;
    imap = interface.BuildActiveForceMap(slforce);
    if (not imap.is_null()) gSlActiveForceMap = LINALG::MergeMap(*gSlActiveForceMap, *imap);

    imap = Teuchos::null;
    imap = interface.BuildActiveForceMap(maforce);
    if (not imap.is_null()) gMaActiveForceMap = LINALG::MergeMap(*gMaActiveForceMap, *imap);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
CONTACT::AUG::ComboStrategy::PreAsymptoticSwitching::GetStructuralForceWithoutDbcDofs(
    const CONTACT::ParamsInterface& cparams)
{
  Teuchos::RCP<Epetra_Vector> force_no_dbc_ptr =
      Teuchos::rcp(new Epetra_Vector(*combo_.no_dbc_.slMaMap_));

  const STR::MODELEVALUATOR::Contact& cmodel =
      dynamic_cast<const STR::MODELEVALUATOR::Contact&>(cparams.GetModelEvaluator());

  const std::vector<INPAR::STR::ModelType> without_contact_model(1, cmodel.Type());
  Teuchos::RCP<Epetra_Vector> force_ptr = cmodel.AssembleForceOfModels(&without_contact_model);

  LINALG::Export(*force_ptr, *force_no_dbc_ptr);

  return force_no_dbc_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::PreAsymptoticSwitching::CheckPenetration(std::ostream& os)
{
  // get the overall largest penetration value
  double min_awgap = 0.0;
  combo_.data_.AWGap().MinValue(&min_awgap);
  double max_awgap = 0.0;
  combo_.data_.AWGap().MaxValue(&max_awgap);

  const double max_abs_awgap = std::max(std::abs(min_awgap), std::abs(max_awgap));
  const double penbound = GetPenetrationBound();

#ifdef DEBUG_COMBO_STRATEGY
  IO::cout << __LINE__ << " -- " << __PRETTY_FUNCTION__ << IO::endl;
  IO::cout << "max_awgap = " << min_awgap << " | penbound = " << penbound << IO::endl;
#endif

  const bool pen_check = (max_abs_awgap < penbound);

  // update minimal averaged weighted gap container
  maxabsawgap_.Update(max_abs_awgap);

  os << std::setfill('.') << std::setw(80)
     << "# Checking penetration ( max{ | avg.-w. gaps | } < pen. bound) "
     << (pen_check ? "SUCCEEDED" : "FAILED") << "\n"
     << std::setfill(' ');
  os << "  max abs. gap = " << std::setw(14) << std::setprecision(4) << std::scientific
     << max_abs_awgap << " < " << std::setw(14) << std::setprecision(4) << std::scientific
     << penbound << "\n";
  os << "  ( max. pre-asymptotic = " << maxabsawgap_.pre_asymptotic_ << ", "
     << "max. asymptotic = " << maxabsawgap_.asymptotic_ << " )\n";

  return pen_check;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::ComboStrategy::PreAsymptoticSwitching::GetPenetrationBound() const
{
  double penbound = 1.0e12;
  for (const Teuchos::RCP<CONTACT::CoInterface>& cit : combo_.Interfaces())
  {
    const CONTACT::AUG::Interface& interface = dynamic_cast<CONTACT::AUG::Interface&>(*cit);
    penbound = std::min(penbound, interface.PenetrationBound());
  }

#ifdef DEBUG_COMBO_STRATEGY
  std::cout << __LINE__ << " -- " << __PRETTY_FUNCTION__ << std::endl;
  std::cout << "penbound = " << penbound << std::endl;
#endif

  return penbound;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::PreAsymptoticSwitching::PrintUpdateHead(std::ostream& os) const
{
  os << "--- ComboStrategy::PreAsymptoticSwitching::Update\n";
}
