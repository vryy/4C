/*---------------------------------------------------------------------*/
/*!
\brief This strategy allows the combination of an arbitrary number of
       augmented contact solving strategies.

\level 3

\maintainer Matthias Mayr
*/
/*---------------------------------------------------------------------*/

#include "contact_aug_combo_strategy.H"
#include "contact_augmented_strategy.H"
#include "contact_augmented_interface.H"

#include "../drt_contact/contact_strategy_factory.H"
#include "../drt_contact/contact_paramsinterface.H"

#include "../drt_structure_new/str_solver_factory.H"

#include "../drt_inpar/inpar_contact.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_utils_parallel.H"

#include "../drt_io/io_pstream.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CoAbstractStrategy> CONTACT::AUG::ComboStrategy::Create(
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
        DRT::INPUT::IntegralValue<enum INPAR::CONTACT::SolvingStrategy>(p_combo, strat_count.str());

    CreateStrategyInterfaces(strat_type, ref_interfaces, strat_interfaces);

    strategies.push_back(STRATEGY::Factory::BuildStrategy(strat_type, params, false, false, maxdof,
        strat_interfaces, dof_row_map, node_row_map, dim, comm, data));

    CreateStrategyLinearSolvers(
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
void CONTACT::AUG::ComboStrategy::CreateStrategyInterfaces(
    const enum INPAR::CONTACT::SolvingStrategy strat_type,
    const plain_interface_set& ref_interfaces, plain_interface_set& strat_interfaces)
{
  strat_interfaces.clear();
  strat_interfaces.reserve(ref_interfaces.size());

  Teuchos::RCP<CONTACT::CoInterface> newinterface = Teuchos::null;

  for (plain_interface_set::const_iterator cit = ref_interfaces.begin();
       cit != ref_interfaces.end(); ++cit)
  {
    const CONTACT::AUG::Interface& ref_interface = dynamic_cast<CONTACT::AUG::Interface&>(**cit);

    const Teuchos::RCP<CONTACT::AUG::IDataContainer> idata_ptr =
        ref_interface.SharedInterfaceDataPtr();
    CONTACT::AUG::IDataContainer& idata = *idata_ptr;

    /* create a new interface by copying the data pointer from the reference
     * interface */
    newinterface = STRATEGY::Factory::CreateInterface(strat_type, idata.Id(), idata.Comm(),
        idata.Dim(), idata.IMortar(), idata.IsSelfContact(), idata.RedundantStorage(),
        Teuchos::null, idata_ptr);

    strat_interfaces.push_back(newinterface);
    // reset pointer
    newinterface = Teuchos::null;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::CreateStrategyLinearSolvers(
    const CONTACT::CoAbstractStrategy& strategy, const std::string& lin_solver_id_str,
    const Teuchos::ParameterList& params, CONTACT::ParamsInterface* cparams_interface,
    plain_lin_solver_set& strat_lin_solvers)
{
  const Teuchos::ParameterList& p_combo = params.sublist("AUGMENTED").sublist("COMBO");

  int ls_id = p_combo.get<int>(lin_solver_id_str);
  ls_id = (ls_id == -1 ? params.get<int>("LINEAR_SOLVER") : ls_id);
  if (ls_id == -1)
    dserror(
        "You must specify a reasonable LINEAR_SOLVER ID for the combo "
        "%s either as CONTACT DYNAMIC/LINEAR_SOLVER or as "
        "CONTACT DYNAMIC/AUGMENTED/COMBO/LINEAR_SOLVER_STRATEGY_%c. "
        "However, you provided no LINEAR_SOLVER at all, thus, please fix your "
        "INPUT file. ",
        INPAR::CONTACT::SolvingStrategy2String(strategy.Type()).c_str(), lin_solver_id_str.back());

  if (not cparams_interface)
    dserror("You have to provide a pointer to the CONTACT::ParamsInterface!");

  DRT::DiscretizationInterface* str_discret =
      cparams_interface->Get<DRT::DiscretizationInterface>();

  strat_lin_solvers.push_back(STR::SOLVER::Factory::BuildMeshtyingContactLinSolver(
      *str_discret, strategy.Type(), strategy.SystemType(), ls_id));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::ComboStrategy::ComboStrategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data, const Epetra_Map* dof_row_map,
    const Epetra_Map* node_row_map, const Teuchos::ParameterList& params,
    const plain_strategy_set& strategies, const plain_lin_solver_set& lin_solvers, const int dim,
    const Teuchos::RCP<const Epetra_Comm>& comm, const int maxdof)
    : CONTACT::CoAbstractStrategy(data, dof_row_map, node_row_map, params, dim, comm, 0.0, maxdof),
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
    const CONTACT::CoAbstractStrategy& s = **cit;
    const plain_interface_set& sinterfaces = s.ContactInterfaces();
    const unsigned num_interfaces = sinterfaces.size();

    interface_sets_.push_back(plain_interface_set(num_interfaces, Teuchos::null));
    plain_interface_set& interfaces = interface_sets_.back();
    std::copy(sinterfaces.begin(), sinterfaces.end(), interfaces.begin());
  }

  /// extract parameters from the parameter list
  const Teuchos::ParameterList& p_combo = params.sublist("AUGMENTED").sublist("COMBO");
  output_.initScreenOutput(DRT::INPUT::IntegralValue<bool>(p_combo, "PRINT2SCREEN"));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::CONTACT::SolvingStrategy CONTACT::AUG::ComboStrategy::Type() const { return Get().Type(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LINALG::Solver* CONTACT::AUG::ComboStrategy::GetLinearSolver() const
{
  return lin_solvers_.at(switch_->Id()).get();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPreEvaluate(CONTACT::ParamsInterface& cparams)
{
  Get().RunPreEvaluate(cparams);
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
      RunPostEvalForce(cparams);
      break;
    }
    case MORTAR::eval_force_stiff:
    {
      RunPostEvalForceStiff(cparams);
      break;
    }
    case MORTAR::eval_static_constraint_rhs:
    {
      RunPostEvalStaticConstraintRHS(cparams);
      break;
    }
    default:
      dserror("Unexpected MORTAR::ActionType! ( actiontype = %s | %d )",
          MORTAR::ActionType2String(curr_eval).c_str(), curr_eval);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPreSolve(
    const Teuchos::RCP<const Epetra_Vector>& curr_disp, const CONTACT::ParamsInterface& cparams)
{
  Get().RunPreSolve(curr_disp, cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::Reset(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& dispnp, const Epetra_Vector& xnew)
{
  Get().Reset(cparams, dispnp, xnew);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::IsSaddlePointSystem() const
{
  return Get().IsSaddlePointSystem();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::ResetActiveSet() { Get().ResetActiveSet(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis)
{
  Get().SaveReferenceState(dis);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::ComboStrategy::ConstraintNorm() const { return Get().ConstraintNorm(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::ActiveSetConverged() { return Get().ActiveSetConverged(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONTACT::AUG::ComboStrategy::ActiveSetSteps() { return Get().ActiveSetSteps(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::UpdateActiveSet() { return Get().UpdateActiveSet(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::EvaluateRelMovPredict() { Get().EvaluateRelMovPredict(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::ActiveSetSemiSmoothConverged() const
{
  return Get().ActiveSetSemiSmoothConverged();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> CONTACT::AUG::ComboStrategy::GetOldActiveRowNodes() const
{
  return Get().GetOldActiveRowNodes();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> CONTACT::AUG::ComboStrategy::GetOldSlipRowNodes() const
{
  return Get().GetOldSlipRowNodes();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> CONTACT::AUG::ComboStrategy::SlNormalDoFRowMapPtr(
    const bool& redist) const
{
  return Get().SlNormalDoFRowMapPtr(redist);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map& CONTACT::AUG::ComboStrategy::SlNormalDoFRowMap(const bool& redist) const
{
  return Get().SlNormalDoFRowMap(redist);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> CONTACT::AUG::ComboStrategy::SlTangentialDoFRowMapPtr(
    const bool& redist) const
{
  return Get().SlTangentialDoFRowMapPtr(redist);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map& CONTACT::AUG::ComboStrategy::SlTangentialDoFRowMap(const bool& redist) const
{
  return Get().SlTangentialDoFRowMap(redist);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::AUG::ComboStrategy::GetRhsBlockPtr(
    const enum DRT::UTILS::VecBlockType& bt) const
{
  return Get().GetRhsBlockPtr(bt);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::AUG::ComboStrategy::GetRhsBlockPtrForNormCheck(
    const enum DRT::UTILS::VecBlockType& bt) const
{
  return Get().GetRhsBlockPtrForNormCheck(bt);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::AUG::ComboStrategy::GetCondensedRhsPtr(
    Epetra_Vector& f, const double& timefac_np) const
{
  return Get().GetCondensedRhsPtr(f, timefac_np);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> CONTACT::AUG::ComboStrategy::GetMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt, const CONTACT::ParamsInterface* cparams) const
{
  return Get().GetMatrixBlockPtr(bt, cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> CONTACT::AUG::ComboStrategy::GetCondensedMatrixBlockPtr(
    Teuchos::RCP<LINALG::SparseMatrix>& kteff, const double& timefac_np) const
{
  return Get().GetCondensedMatrixBlockPtr(kteff, timefac_np);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> CONTACT::AUG::ComboStrategy::ConstrRhs() { return Get().ConstrRhs(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::Initialize() { dserror("Unnecessary in this Strategy."); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::EvalConstrRHS() { dserror("Unnecessary in this Strategy."); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::UpdateActiveSetSemiSmooth(const bool firstStepPredictor)
{
  dserror(
      "Unnecessary in this Strategy. Furthermore, this method is "
      "deprecated in the AUG::Strategy framework and has been replaced by "
      "CONTACT::AUG::UpdateActiveSetSemiSmooth( const CONTACT::ParamsInterface& ).");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::DoReadRestart(IO::DiscretizationReader& reader,
    Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr)
{
  Get().DoReadRestart(reader, dis, cparams_ptr);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::Update(Teuchos::RCP<const Epetra_Vector> dis)
{
  Get().Update(dis);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::WasInContactLastIter() const
{
  return Get().WasInContactLastIter();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::OutputStresses() { Get().OutputStresses(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::PostSetup(bool redistributed, bool init)
{
  if (redistributed) no_dbc_.Redistribute(data_);

  Get().PostSetup(redistributed, init);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::PostStoreDirichletStatus(
    Teuchos::RCP<const LINALG::MapExtractor> dbcmaps)
{
  no_dbc_.Assemble(*dbcmaps->CondMap(), data_);

  Get().PostStoreDirichletStatus(dbcmaps);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::GlobalNoDbc::Assemble(
    const Epetra_Map& dbcmap, const CONTACT::AUG::DataContainer& data)
{
  const Epetra_Map& gSlMaDofRowMap = *data.GSlMaDofRowMapPtr();

  Teuchos::RCP<Epetra_Map> gSlMaDbcDofRowMap = LINALG::IntersectMap(gSlMaDofRowMap, dbcmap);

  slMaMap_ = LINALG::SplitMap(gSlMaDofRowMap, *gSlMaDbcDofRowMap);

  Reset(*slMaMap_, data);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::GlobalNoDbc::Redistribute(const CONTACT::AUG::DataContainer& data)
{
  slMaMap_ =
      DRT::UTILS::RedistributeInAccordanceWithReference(*data.GSlMaDofRowMapPtr(), *slMaMap_);

  Reset(*slMaMap_, data);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::GlobalNoDbc::Reset(
    const Epetra_Map& slMaMap, const CONTACT::AUG::DataContainer& data)
{
  slMap_ = LINALG::IntersectMap(slMaMap, *data.GSlDofRowMapPtr());
  maMap_ = LINALG::IntersectMap(slMaMap, *data.GMaDofRowMapPtr());

  slMaForce_ = Teuchos::rcp(new Epetra_Vector(slMaMap, true));
  slForce_ = Teuchos::rcp(new Epetra_Vector(*slMap_, true));
  maForce_ = Teuchos::rcp(new Epetra_Vector(*maMap_, true));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::EvalForce(CONTACT::ParamsInterface& cparams)
{
  Get().PreEvalForce(cparams);

  Get().AssembleGap();

  Get().EvalAugmentedForces();

#ifdef DEBUG_COMBO_STRATEGY
  std::cout << __LINE__ << " -- " << __PRETTY_FUNCTION__ << std::endl;
  std::cout << INPAR::CONTACT::SolvingStrategy2String(Get().Type()) << "\n";
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPostEvalForce(CONTACT::ParamsInterface& cparams)
{
  Get().AssembleContactRHS();
  Get().EvalStrContactRHS();
  Get().EvalConstrRHS();

  Get().PostEvalForce(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::EvalForceStiff(CONTACT::ParamsInterface& cparams)
{
  EvalForce(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPostEvalForceStiff(CONTACT::ParamsInterface& cparams)
{
  SwitchUpdate(cparams);
  RunPostEvalForce(cparams);

  Get().AssembleContactStiff();
  Get().PostEvalForceStiff(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::EvalStaticConstraintRHS(CONTACT::ParamsInterface& cparams)
{
  Get().SetCurrentEvalState(cparams);
  Get().InitEvalInterface(cparams);

  Get().Initialize(cparams.GetActionType());
  Get().AssembleGap();

  Get().EvalConstraintForces();
  Get().ZeroizeLMForces();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPostEvalStaticConstraintRHS(CONTACT::ParamsInterface& cparams)
{
  //  switch_->Update( cparams );

  Get().AssembleContactRHS();
  Get().EvalStrContactRHS();
  Get().EvalConstrRHS();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::SwitchUpdate(CONTACT::ParamsInterface& cparams)
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
  Get().RunPostIterate(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPreComputeX(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold, Epetra_Vector& dir_mutable)
{
  Get().RunPreComputeX(cparams, xold, dir_mutable);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPostComputeX(const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  Get().RunPostComputeX(cparams, xold, dir, xnew);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RunPostApplyJacobianInverse(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& rhs, Epetra_Vector& result,
    const Epetra_Vector& xold, const NOX::NLN::Group& grp)
{
  Get().RunPostApplyJacobianInverse(cparams, rhs, result, xold, grp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::RemoveCondensedContributionsFromRhs(Epetra_Vector& str_rhs) const
{
  Get().RemoveCondensedContributionsFromRhs(str_rhs);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::CorrectParameters(
    CONTACT::ParamsInterface& cparams, const NOX::NLN::CorrectionType type)
{
  Get().CorrectParameters(cparams, type);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::ResetLagrangeMultipliers(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xnew)
{
  Get().ResetLagrangeMultipliers(cparams, xnew);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CONTACT::CoInterface>>& CONTACT::AUG::ComboStrategy::Interfaces()
{
  return interface_sets_[switch_->Id()];
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::vector<Teuchos::RCP<CONTACT::CoInterface>>& CONTACT::AUG::ComboStrategy::Interfaces()
    const
{
  return interface_sets_[switch_->Id()];
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::Strategy& CONTACT::AUG::ComboStrategy::Get()
{
  return dynamic_cast<CONTACT::AUG::Strategy&>(*strategies_[switch_->Id()]);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const CONTACT::AUG::Strategy& CONTACT::AUG::ComboStrategy::Get() const
{
  return dynamic_cast<const CONTACT::AUG::Strategy&>(*strategies_[switch_->Id()]);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::ComboStrategy::GetPotentialValue(
    const enum NOX::NLN::MeritFunction::MeritFctName mrt_type) const
{
  return Get().GetPotentialValue(mrt_type);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::ComboStrategy::GetLinearizedPotentialValueTerms(const Epetra_Vector& dir,
    const enum NOX::NLN::MeritFunction::MeritFctName mrt_type,
    const enum NOX::NLN::MeritFunction::LinOrder linorder,
    const enum NOX::NLN::MeritFunction::LinType lintype) const
{
  return Get().GetLinearizedPotentialValueTerms(dir, mrt_type, linorder, lintype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::WriteOutput(IO::DiscretizationWriter& writer) const
{
  return Get().WriteOutput(writer);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::ComboStrategy::EvaluateReferenceState(Teuchos::RCP<const Epetra_Vector> vec)
{
  Get().EvaluateReferenceState(vec);
}

/*----------------------------------------------------------------------*
 *  *----------------------------------------------------------------------*/
bool CONTACT::AUG::ComboStrategy::DynRedistributeContact(
    const Teuchos::RCP<const Epetra_Vector>& dis, const int nlniter)
{
  const bool is_redistributed = Get().DynRedistributeContact(dis, nlniter);

  // This function must be called manually, since the internal PostSetup
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
  if (not oscreen_) dserror("Call initScreenOutput first");
  return *oscreen_;
}
