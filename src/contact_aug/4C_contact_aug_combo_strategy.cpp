/*---------------------------------------------------------------------*/
/*! \file
\brief This strategy allows the combination of an arbitrary number of
       augmented contact solving strategies.

\level 3

*/
/*---------------------------------------------------------------------*/

#include "4C_contact_aug_combo_strategy.hpp"

#include "4C_contact_aug_interface.hpp"
#include "4C_contact_aug_strategy.hpp"
#include "4C_contact_paramsinterface.hpp"
#include "4C_contact_strategy_factory.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_structure_new_solver_factory.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::AbstractStrategy> CONTACT::AUG::ComboStrategy::Create(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data, const Epetra_Map* dof_row_map,
    const Epetra_Map* node_row_map, const Teuchos::ParameterList& params,
    const plain_interface_set& ref_interfaces, const int dim,
    const Teuchos::RCP<const Epetra_Comm>& comm, const int maxdof,
    CONTACT::ParamsInterface* cparams_interface)
{
  const Teuchos::ParameterList& p_combo = params.sublist("AUGMENTED").sublist("COMBO");

  plain_strategy_set strategies(0);
  plain_interface_set strat_interfaces(0);
  plain_lin_solver_set strat_lin_solvers(0);

  unsigned count = 0;
  std::ostringstream strat_count;
  strat_count << "STRATEGY_" << count;
  std::ostringstream lin_solver_count;
  lin_solver_count << "LINEAR_SOLVER_STRATEGY_" << count;
  while (p_combo.isParameter(strat_count.str()))
  {
    const enum INPAR::CONTACT::SolvingStrategy strat_type =
        CORE::UTILS::IntegralValue<enum INPAR::CONTACT::SolvingStrategy>(
            p_combo, strat_count.str());

    create_strategy_interfaces(strat_type, ref_interfaces, strat_interfaces);

    strategies.push_back(STRATEGY::Factory::BuildStrategy(strat_type, params, false, false, maxdof,
        strat_interfaces, dof_row_map, node_row_map, dim, comm, data));

    create_strategy_linear_solvers(
        *strategies.back(), lin_solver_count.str(), params, cparams_interface, strat_lin_solvers);

    /// clear and increase strategy count string
    strat_count.str("");
    lin_solver_count.str("");
    strat_count << "STRATEGY_" << ++count;
    lin_solver_count << "LINEAR_SOLVER_STRATEGY_" << count;
  }

  return Teuchos::rcp(new ComboStrategy(
      data, dof_row_map, node_row_map, params, strategies, strat_lin_solvers, dim, comm, maxdof));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::create_strategy_interfaces(
    const enum INPAR::CONTACT::SolvingStrategy strat_type,
    const plain_interface_set& ref_interfaces, plain_interface_set& strat_interfaces)
{
  strat_interfaces.clear();
  strat_interfaces.reserve(ref_interfaces.size());

  Teuchos::RCP<CONTACT::Interface> newinterface = Teuchos::null;

  for (plain_interface_set::const_iterator cit = ref_interfaces.begin();
       cit != ref_interfaces.end(); ++cit)
  {
    const CONTACT::AUG::Interface& ref_interface = dynamic_cast<CONTACT::AUG::Interface&>(**cit);

    const Teuchos::RCP<CONTACT::AUG::InterfaceDataContainer> idata_ptr =
        ref_interface.shared_interface_data_ptr();
    CONTACT::AUG::InterfaceDataContainer& idata = *idata_ptr;

    /* create a new interface by copying the data pointer from the reference
     * interface */
    newinterface = STRATEGY::Factory::CreateInterface(strat_type, idata.Id(), idata.Comm(),
        idata.Dim(), idata.IMortar(), idata.IsSelfContact(), Teuchos::null, idata_ptr);

    strat_interfaces.push_back(newinterface);
    // reset pointer
    newinterface = Teuchos::null;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::create_strategy_linear_solvers(
    const CONTACT::AbstractStrategy& strategy, const std::string& lin_solver_id_str,
    const Teuchos::ParameterList& params, CONTACT::ParamsInterface* cparams_interface,
    plain_lin_solver_set& strat_lin_solvers)
{
  const Teuchos::ParameterList& p_combo = params.sublist("AUGMENTED").sublist("COMBO");

  int ls_id = p_combo.get<int>(lin_solver_id_str);
  ls_id = (ls_id == -1 ? params.get<int>("LINEAR_SOLVER") : ls_id);
  if (ls_id == -1)
    FOUR_C_THROW(
        "You must specify a reasonable LINEAR_SOLVER ID for the combo "
        "%s either as CONTACT DYNAMIC/LINEAR_SOLVER or as "
        "CONTACT DYNAMIC/AUGMENTED/COMBO/LINEAR_SOLVER_STRATEGY_%c. "
        "However, you provided no LINEAR_SOLVER at all, thus, please fix your "
        "INPUT file. ",
        INPAR::CONTACT::SolvingStrategy2String(strategy.Type()).c_str(), lin_solver_id_str.back());

  if (not cparams_interface)
    FOUR_C_THROW("You have to provide a pointer to the CONTACT::ParamsInterface!");

  DRT::Discretization* str_discret = cparams_interface->Get<DRT::Discretization>();

  strat_lin_solvers.push_back(STR::SOLVER::Factory::build_meshtying_contact_lin_solver(
      *str_discret, strategy.Type(), strategy.SystemType(), ls_id));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::ComboStrategy::ComboStrategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data, const Epetra_Map* dof_row_map,
    const Epetra_Map* node_row_map, const Teuchos::ParameterList& params,
    const plain_strategy_set& strategies, const plain_lin_solver_set& lin_solvers, const int dim,
    const Teuchos::RCP<const Epetra_Comm>& comm, const int maxdof)
    : CONTACT::AbstractStrategy(data, dof_row_map, node_row_map, params, dim, comm, 0.0, maxdof),
      strategies_(strategies),
      lin_solvers_(lin_solvers),
      interface_sets_(0),
      data_(dynamic_cast<CONTACT::AUG::DataContainer&>(*data)),
      no_dbc_(),
      switch_(Switching::Create(*this))
{
  for (plain_strategy_set::const_iterator cit = strategies_.begin(); cit != strategies_.end();
       ++cit)
  {
    const CONTACT::AbstractStrategy& s = **cit;
    const plain_interface_set& sinterfaces = s.ContactInterfaces();
    const unsigned num_interfaces = sinterfaces.size();

    interface_sets_.push_back(plain_interface_set(num_interfaces, Teuchos::null));
    plain_interface_set& interfaces = interface_sets_.back();
    std::copy(sinterfaces.begin(), sinterfaces.end(), interfaces.begin());
  }

  /// extract parameters from the parameter list
  const Teuchos::ParameterList& p_combo = params.sublist("AUGMENTED").sublist("COMBO");
  output_.initScreenOutput(CORE::UTILS::IntegralValue<bool>(p_combo, "PRINT2SCREEN"));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::CONTACT::SolvingStrategy CONTACT::AUG::ComboStrategy::Type() const { return get().Type(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::LINALG::Solver* CONTACT::AUG::ComboStrategy::GetLinearSolver() const
{
  return lin_solvers_.at(switch_->Id()).get();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPreEvaluate(CONTACT::ParamsInterface& cparams)
{
  get().RunPreEvaluate(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPostEvaluate(CONTACT::ParamsInterface& cparams)
{
  const MORTAR::ActionType curr_eval = data_.GetCurrentEvalState();
  IO::cout << MORTAR::ActionType2String(curr_eval) << "\n";
  switch (curr_eval)
  {
    case MORTAR::eval_force:
    {
      run_post_eval_force(cparams);
      break;
    }
    case MORTAR::eval_force_stiff:
    {
      run_post_eval_force_stiff(cparams);
      break;
    }
    case MORTAR::eval_static_constraint_rhs:
    {
      run_post_eval_static_constraint_rhs(cparams);
      break;
    }
    default:
      FOUR_C_THROW("Unexpected MORTAR::ActionType! ( actiontype = %s | %d )",
          MORTAR::ActionType2String(curr_eval).c_str(), curr_eval);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPreSolve(
    const Teuchos::RCP<const Epetra_Vector>& curr_disp, const CONTACT::ParamsInterface& cparams)
{
  get().RunPreSolve(curr_disp, cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::Reset(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& dispnp, const Epetra_Vector& xnew)
{
  get().Reset(cparams, dispnp, xnew);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::IsSaddlePointSystem() const
{
  return get().IsSaddlePointSystem();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::ResetActiveSet() { get().ResetActiveSet(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis)
{
  get().SaveReferenceState(dis);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::ComboStrategy::ConstraintNorm() const { return get().ConstraintNorm(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::ActiveSetConverged() { return get().ActiveSetConverged(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONTACT::AUG::ComboStrategy::ActiveSetSteps() { return get().ActiveSetSteps(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::UpdateActiveSet() { return get().UpdateActiveSet(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::evaluate_rel_mov_predict() { get().evaluate_rel_mov_predict(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::active_set_semi_smooth_converged() const
{
  return get().active_set_semi_smooth_converged();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> CONTACT::AUG::ComboStrategy::get_old_active_row_nodes() const
{
  return get().get_old_active_row_nodes();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> CONTACT::AUG::ComboStrategy::GetOldSlipRowNodes() const
{
  return get().GetOldSlipRowNodes();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> CONTACT::AUG::ComboStrategy::sl_normal_do_f_row_map_ptr(
    const bool& redist) const
{
  return get().sl_normal_do_f_row_map_ptr(redist);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map& CONTACT::AUG::ComboStrategy::SlNormalDoFRowMap(const bool& redist) const
{
  return get().SlNormalDoFRowMap(redist);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> CONTACT::AUG::ComboStrategy::sl_tangential_do_f_row_map_ptr(
    const bool& redist) const
{
  return get().sl_tangential_do_f_row_map_ptr(redist);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map& CONTACT::AUG::ComboStrategy::sl_tangential_do_f_row_map(const bool& redist) const
{
  return get().sl_tangential_do_f_row_map(redist);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::AUG::ComboStrategy::GetRhsBlockPtr(
    const enum CONTACT::VecBlockType& bt) const
{
  return get().GetRhsBlockPtr(bt);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::AUG::ComboStrategy::get_rhs_block_ptr_for_norm_check(
    const enum CONTACT::VecBlockType& bt) const
{
  return get().get_rhs_block_ptr_for_norm_check(bt);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::AUG::ComboStrategy::GetCondensedRhsPtr(
    Epetra_Vector& f, const double& timefac_np) const
{
  return get().GetCondensedRhsPtr(f, timefac_np);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> CONTACT::AUG::ComboStrategy::GetMatrixBlockPtr(
    const enum CONTACT::MatBlockType& bt, const CONTACT::ParamsInterface* cparams) const
{
  return get().GetMatrixBlockPtr(bt, cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix>
CONTACT::AUG::ComboStrategy::get_condensed_matrix_block_ptr(
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& kteff, const double& timefac_np) const
{
  return get().get_condensed_matrix_block_ptr(kteff, timefac_np);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> CONTACT::AUG::ComboStrategy::ConstrRhs() { return get().ConstrRhs(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::Initialize() { FOUR_C_THROW("Unnecessary in this Strategy."); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::EvalConstrRHS() { FOUR_C_THROW("Unnecessary in this Strategy."); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::update_active_set_semi_smooth(const bool firstStepPredictor)
{
  FOUR_C_THROW(
      "Unnecessary in this Strategy. Furthermore, this method is "
      "deprecated in the AUG::Strategy framework and has been replaced by "
      "CONTACT::AUG::update_active_set_semi_smooth( const CONTACT::ParamsInterface& ).");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::DoReadRestart(IO::DiscretizationReader& reader,
    Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr)
{
  get().DoReadRestart(reader, dis, cparams_ptr);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::Update(Teuchos::RCP<const Epetra_Vector> dis)
{
  get().Update(dis);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::was_in_contact_last_iter() const
{
  return get().was_in_contact_last_iter();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::compute_contact_stresses() { get().compute_contact_stresses(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::post_setup(bool redistributed, bool init)
{
  if (redistributed) no_dbc_.Redistribute(data_);

  get().post_setup(redistributed, init);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::post_store_dirichlet_status(
    Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmaps)
{
  no_dbc_.Assemble(*dbcmaps->CondMap(), data_);

  get().post_store_dirichlet_status(dbcmaps);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::GlobalNoDbc::Assemble(
    const Epetra_Map& dbcmap, const CONTACT::AUG::DataContainer& data)
{
  const Epetra_Map& gSlMaDofRowMap = *data.GSlMaDofRowMapPtr();

  Teuchos::RCP<Epetra_Map> gSlMaDbcDofRowMap = CORE::LINALG::IntersectMap(gSlMaDofRowMap, dbcmap);

  slMaMap_ = CORE::LINALG::SplitMap(gSlMaDofRowMap, *gSlMaDbcDofRowMap);

  Reset(*slMaMap_, data);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::GlobalNoDbc::Redistribute(const CONTACT::AUG::DataContainer& data)
{
  slMaMap_ =
      CORE::REBALANCE::RebalanceInAccordanceWithReference(*data.GSlMaDofRowMapPtr(), *slMaMap_);

  Reset(*slMaMap_, data);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::GlobalNoDbc::Reset(
    const Epetra_Map& slMaMap, const CONTACT::AUG::DataContainer& data)
{
  slMap_ = CORE::LINALG::IntersectMap(slMaMap, *data.GSlDofRowMapPtr());
  maMap_ = CORE::LINALG::IntersectMap(slMaMap, *data.GMaDofRowMapPtr());

  slMaForce_ = Teuchos::rcp(new Epetra_Vector(slMaMap, true));
  slForce_ = Teuchos::rcp(new Epetra_Vector(*slMap_, true));
  maForce_ = Teuchos::rcp(new Epetra_Vector(*maMap_, true));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::EvalForce(CONTACT::ParamsInterface& cparams)
{
  get().PreEvalForce(cparams);

  get().AssembleGap();

  get().EvalAugmentedForces();

#ifdef DEBUG_COMBO_STRATEGY
  std::cout << __LINE__ << " -- " << __PRETTY_FUNCTION__ << std::endl;
  std::cout << INPAR::CONTACT::SolvingStrategy2String(Get().Type()) << "\n";
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::run_post_eval_force(CONTACT::ParamsInterface& cparams)
{
  get().AssembleContactRHS();
  get().EvalStrContactRHS();
  get().EvalConstrRHS();

  get().PostEvalForce(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::EvalForceStiff(CONTACT::ParamsInterface& cparams)
{
  EvalForce(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::run_post_eval_force_stiff(CONTACT::ParamsInterface& cparams)
{
  switch_update(cparams);
  run_post_eval_force(cparams);

  get().assemble_contact_stiff();
  get().PostEvalForceStiff(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::eval_static_constraint_rhs(CONTACT::ParamsInterface& cparams)
{
  get().SetCurrentEvalState(cparams);
  get().InitEvalInterface(cparams);

  get().Initialize(cparams.GetActionType());
  get().AssembleGap();

  get().eval_constraint_forces();
  get().ZeroizeLMForces();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::run_post_eval_static_constraint_rhs(
    CONTACT::ParamsInterface& cparams)
{
  //  switch_->Update( cparams );

  get().AssembleContactRHS();
  get().EvalStrContactRHS();
  get().EvalConstrRHS();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::switch_update(CONTACT::ParamsInterface& cparams)
{
  IO::cout(IO::debug) << std::string(40, '*') << "\n";
  IO::cout(IO::debug) << CONTACT_FUNC_NAME << IO::endl;
  IO::cout(IO::debug) << "Correction Type = "
                      << NOX::NLN::CorrectionType2String(cparams.GetCorrectionType()).c_str()
                      << "\n";
  IO::cout(IO::debug) << std::string(40, '*') << "\n";

  /* Do not perform a switch in case of a correction step, since this will lead
   * to an error during the potential backup evaluation. */
  if (cparams.GetCorrectionType() != NOX::NLN::CorrectionType::vague) return;

  switch_->Update(cparams, output_.oscreen());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPostIterate(const CONTACT::ParamsInterface& cparams)
{
  get().RunPostIterate(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPreComputeX(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold, Epetra_Vector& dir_mutable)
{
  get().RunPreComputeX(cparams, xold, dir_mutable);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPostComputeX(const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  get().RunPostComputeX(cparams, xold, dir, xnew);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::run_post_apply_jacobian_inverse(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& rhs, Epetra_Vector& result,
    const Epetra_Vector& xold, const NOX::NLN::Group& grp)
{
  get().run_post_apply_jacobian_inverse(cparams, rhs, result, xold, grp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::remove_condensed_contributions_from_rhs(
    Epetra_Vector& str_rhs) const
{
  get().remove_condensed_contributions_from_rhs(str_rhs);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::CorrectParameters(
    CONTACT::ParamsInterface& cparams, const NOX::NLN::CorrectionType type)
{
  get().CorrectParameters(cparams, type);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::reset_lagrange_multipliers(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xnew)
{
  get().reset_lagrange_multipliers(cparams, xnew);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CONTACT::Interface>>& CONTACT::AUG::ComboStrategy::Interfaces()
{
  return interface_sets_[switch_->Id()];
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::vector<Teuchos::RCP<CONTACT::Interface>>& CONTACT::AUG::ComboStrategy::Interfaces() const
{
  return interface_sets_[switch_->Id()];
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::Strategy& CONTACT::AUG::ComboStrategy::get()
{
  return dynamic_cast<CONTACT::AUG::Strategy&>(*strategies_[switch_->Id()]);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const CONTACT::AUG::Strategy& CONTACT::AUG::ComboStrategy::get() const
{
  return dynamic_cast<const CONTACT::AUG::Strategy&>(*strategies_[switch_->Id()]);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::ComboStrategy::GetPotentialValue(
    const enum NOX::NLN::MeritFunction::MeritFctName mrt_type) const
{
  return get().GetPotentialValue(mrt_type);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::ComboStrategy::get_linearized_potential_value_terms(const Epetra_Vector& dir,
    const enum NOX::NLN::MeritFunction::MeritFctName mrt_type,
    const enum NOX::NLN::MeritFunction::LinOrder linorder,
    const enum NOX::NLN::MeritFunction::LinType lintype) const
{
  return get().get_linearized_potential_value_terms(dir, mrt_type, linorder, lintype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::WriteOutput(IO::DiscretizationWriter& writer) const
{
  return get().WriteOutput(writer);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::evaluate_reference_state() { get().evaluate_reference_state(); }

/*----------------------------------------------------------------------*
 *  *----------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::dyn_redistribute_contact(
    const Teuchos::RCP<const Epetra_Vector>& dis, Teuchos::RCP<const Epetra_Vector> vel,
    const int nlniter)
{
  const bool is_redistributed = get().dyn_redistribute_contact(dis, vel, nlniter);

  // This function must be called manually, since the internal post_setup
  // call will not effect this wrapper class.
  if (is_redistributed) no_dbc_.Redistribute(data_);

  return is_redistributed;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::Output::initScreenOutput(bool print2screen)
{
  if (print2screen)
    oscreen_ = &IO::cout.os(IO::standard);
  else
    oscreen_ = blackhole_.get();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& CONTACT::AUG::ComboStrategy::Output::oscreen() const
{
  if (not oscreen_) FOUR_C_THROW("Call initScreenOutput first");
  return *oscreen_;
}

FOUR_C_NAMESPACE_CLOSE
