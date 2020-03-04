/*---------------------------------------------------------------------*/
/*! \file
\brief Main abstract class for contact solution strategies

\level 2

\maintainer Matthias Mayr

*/
/*---------------------------------------------------------------------*/
#include "Epetra_SerialComm.h"
#include "contact_abstract_strategy.H"
#include "contact_defines.H"
#include "contact_interface.H"
#include "friction_node.H"
#include "contact_paramsinterface.H"
#include "contact_noxinterface.H"
#include "contact_utils_parallel.H"

#include "../drt_mortar/mortar_defines.H"
#include "../drt_mortar/mortar_utils.H"

#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_colors.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../linalg/linalg_utils_densematrix_communication.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../solver_nonlin_nox/nox_nln_group.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::AbstractStratDataContainer::AbstractStratDataContainer()
    : glmdofrowmap_(Teuchos::null),
      gsnoderowmap_(Teuchos::null),
      gmnoderowmap_(Teuchos::null),
      gsdofrowmap_(Teuchos::null),
      gmdofrowmap_(Teuchos::null),
      gndofrowmap_(Teuchos::null),
      gsmdofrowmap_(Teuchos::null),
      gdisprowmap_(Teuchos::null),
      gactivenodes_(Teuchos::null),
      gactivedofs_(Teuchos::null),
      ginactivenodes_(Teuchos::null),
      ginactivedofs_(Teuchos::null),
      gactiven_(Teuchos::null),
      gactivet_(Teuchos::null),
      gslipnodes_(Teuchos::null),
      gslipdofs_(Teuchos::null),
      gslipt_(Teuchos::null),
      gsdofVertex_(Teuchos::null),
      gsdofEdge_(Teuchos::null),
      gsdofSurf_(Teuchos::null),
      unbalanceEvaluationTime_(0),
      unbalanceNumSlaveElements_(0),
      pglmdofrowmap_(Teuchos::null),
      pgsdofrowmap_(Teuchos::null),
      pgmdofrowmap_(Teuchos::null),
      pgsmdofrowmap_(Teuchos::null),
      pgsdirichtoggle_(Teuchos::null),
      partype_(INPAR::MORTAR::parredist_none),
      initial_elecolmap_(Teuchos::null),
      dmatrix_(Teuchos::null),
      mmatrix_(Teuchos::null),
      g_(Teuchos::null),
      tangrhs_(Teuchos::null),
      inactiverhs_(Teuchos::null),
      strContactRhsPtr_(Teuchos::null),
      constrrhs_(Teuchos::null),
      lindmatrix_(Teuchos::null),
      linmmatrix_(Teuchos::null),
      kteffnew_(Teuchos::null),
      dold_(Teuchos::null),
      mold_(Teuchos::null),
      z_(Teuchos::null),
      zold_(Teuchos::null),
      zincr_(Teuchos::null),
      zuzawa_(Teuchos::null),
      stressnormal_(Teuchos::null),
      stresstangential_(Teuchos::null),
      forcenormal_(Teuchos::null),
      forcetangential_(Teuchos::null),
      stepnp_(-1),
      iter_(-1),
      isincontact_(false),
      wasincontact_(false),
      wasincontactlts_(false),
      isselfcontact_(false),
      friction_(false),
      nonSmoothContact_(false),
      regularized_(false),
      dualquadslavetrafo_(false),
      trafo_(Teuchos::null),
      invtrafo_(Teuchos::null),
      dmatrixmod_(Teuchos::null),
      doldmod_(Teuchos::null),
      inttime_(0.0),
      ivel_(0),
      stype_(INPAR::CONTACT::solution_vague),
      constr_direction_(INPAR::CONTACT::constr_vague)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::CoAbstractStrategy::CoAbstractStrategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr, const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap, const Teuchos::ParameterList& params, const int spatialDim,
    const Teuchos::RCP<const Epetra_Comm>& comm, const double alphaf, const int maxdof)
    : MORTAR::StrategyBase(
          data_ptr, DofRowMap, NodeRowMap, params, spatialDim, comm, alphaf, maxdof),
      glmdofrowmap_(data_ptr->GLmDofRowMapPtr()),
      gsnoderowmap_(data_ptr->GSlNodeRowMapPtr()),
      gmnoderowmap_(data_ptr->GMaNodeRowMapPtr()),
      gsdofrowmap_(data_ptr->GSlDofRowMapPtr()),
      gmdofrowmap_(data_ptr->GMaDofRowMapPtr()),
      gndofrowmap_(data_ptr->GInternalDofRowMapPtr()),
      gsmdofrowmap_(data_ptr->GSlMaDofRowMapPtr()),
      gdisprowmap_(data_ptr->GDispDofRowMapPtr()),
      gactivenodes_(data_ptr->GActiveNodeRowMapPtr()),
      gactivedofs_(data_ptr->GActiveDofRowMapPtr()),
      ginactivenodes_(data_ptr->GInActiveNodeRowMapPtr()),
      ginactivedofs_(data_ptr->GInActiveDofRowMapPtr()),
      gactiven_(data_ptr->GActiveNDofRowMapPtr()),
      gactivet_(data_ptr->GActiveTDofRowMapPtr()),
      gslipnodes_(data_ptr->GSlipNodeRowMapPtr()),
      gslipdofs_(data_ptr->GSlipDofRowMapPtr()),
      gslipt_(data_ptr->GSlipTDofRowMapPtr()),
      gsdofVertex_(data_ptr->GSDofVertexRowMapPtr()),
      gsdofEdge_(data_ptr->GSDofEdgeRowMapPtr()),
      gsdofSurf_(data_ptr->GSDofSurfRowMapPtr()),
      unbalanceEvaluationTime_(data_ptr->UnbalanceTimeFactors()),
      unbalanceNumSlaveElements_(data_ptr->UnbalanceElementFactors()),
      pglmdofrowmap_(data_ptr->PGLmDofRowMapPtr()),
      pgsdofrowmap_(data_ptr->PGSlDofRowMapPtr()),
      pgmdofrowmap_(data_ptr->PGMaDofRowMapPtr()),
      pgsmdofrowmap_(data_ptr->PGSlMaDofRowMapPtr()),
      pgsdirichtoggle_(data_ptr->PGSlDirichToggleDofRowMapPtr()),
      initial_elecolmap_(data_ptr->InitialSlMaEleColMap()),
      dmatrix_(data_ptr->DMatrixPtr()),
      mmatrix_(data_ptr->MMatrixPtr()),
      g_(data_ptr->WGapPtr()),
      tangrhs_(data_ptr->TangRhsPtr()),
      inactiverhs_(data_ptr->InactiveRhsPtr()),
      strcontactrhs_(data_ptr->StrContactRhsPtr()),
      constrrhs_(data_ptr->ConstrRhsPtr()),
      lindmatrix_(data_ptr->DLinMatrixPtr()),
      linmmatrix_(data_ptr->MLinMatrixPtr()),
      kteffnew_(data_ptr->kteffnewMatrixPtr()),
      dold_(data_ptr->OldDMatrixPtr()),
      mold_(data_ptr->OldMMatrixPtr()),
      z_(data_ptr->LmPtr()),
      zold_(data_ptr->OldLmPtr()),
      zincr_(data_ptr->LmIncrPtr()),
      zuzawa_(data_ptr->LmUzawaPtr()),
      stressnormal_(data_ptr->StressNormalPtr()),
      stresstangential_(data_ptr->StressTangentialPtr()),
      forcenormal_(data_ptr->ForceNormalPtr()),
      forcetangential_(data_ptr->ForceTangentialPtr()),
      step_(data_ptr->StepNp()),
      iter_(data_ptr->NlnIter()),
      isincontact_(data_ptr->IsInContact()),
      wasincontact_(data_ptr->WasInContact()),
      wasincontactlts_(data_ptr->WasInContactLastTimeStep()),
      isselfcontact_(data_ptr->IsSelfContact()),
      friction_(data_ptr->IsFriction()),
      nonSmoothContact_(data_ptr->IsNonSmoothContact()),
      regularized_(data_ptr->IsRegularized()),
      dualquadslavetrafo_(data_ptr->IsDualQuadSlaveTrafo()),
      trafo_(data_ptr->TrafoPtr()),
      invtrafo_(data_ptr->InvTrafoPtr()),
      dmatrixmod_(data_ptr->ModifiedDMatrixPtr()),
      doldmod_(data_ptr->OldModifiedDMatrixPtr()),
      inttime_(data_ptr->IntTime()),
      ivel_(data_ptr->MeanInterfaceVels()),
      stype_(data_ptr->SolType()),
      constr_direction_(data_ptr->ConstrDirection()),
      data_ptr_(data_ptr),
      noxinterface_ptr_(Teuchos::null)
{
  // set data container pointer (only PRIVATE direct access!)
  data_ptr_->SolType() =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(params, "STRATEGY");
  data_ptr_->ConstrDirection() = DRT::INPUT::IntegralValue<INPAR::CONTACT::ConstraintDirection>(
      params, "CONSTRAINT_DIRECTIONS");
  data_ptr_->ParType() = DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(
      params.sublist("PARALLEL REDISTRIBUTION"), "PARALLEL_REDIST");

  INPAR::CONTACT::FrictionType ftype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(Params(), "FRICTION");

  // set frictional contact status
  if (ftype != INPAR::CONTACT::friction_none) friction_ = true;

  // set nonsmooth contact status
  if (DRT::INPUT::IntegralValue<int>(Params(), "NONSMOOTH_GEOMETRIES")) nonSmoothContact_ = true;

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::Regularization>(
          Params(), "CONTACT_REGULARIZATION") != INPAR::CONTACT::reg_none)
    regularized_ = true;

  // initialize storage fields for parallel redistribution
  unbalanceEvaluationTime_.clear();
  unbalanceNumSlaveElements_.clear();

  // build the NOX::NLN::CONSTRAINT::Interface::Required object
  noxinterface_ptr_ = Teuchos::rcp(new CONTACT::NoxInterface());
  noxinterface_ptr_->Init(Teuchos::rcp(this, false));
  noxinterface_ptr_->Setup();

  return;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const CONTACT::CoAbstractStrategy& strategy)
{
  strategy.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::CoAbstractStrategy::IsRebalancingNecessary(const bool first_time_step)
{
  // No rebalancing of a serial run, since it makes no sense.
  if (Comm().NumProc() == 1) return false;

  bool perform_rebalancing = false;
  const double max_time_unbalance =
      Params().sublist("PARALLEL REDISTRIBUTION").get<double>("MAX_BALANCE");

  double time_average = 0.0;
  double elements_average = 0.0;
  if (!first_time_step) ComputeAndResetParallelBalanceIndicators(time_average, elements_average);

  switch (WhichParRedist())
  {
    case INPAR::MORTAR::parredist_none:
    {
      break;
    }
    case INPAR::MORTAR::parredist_static:
    {
      // Static redistribution: ONLY at time t=0 or after restart
      if (first_time_step)
      {
        // The user demanded to perform rebalancing, so let's do it.
        perform_rebalancing = true;
      }

      break;
    }
    case INPAR::MORTAR::parredist_dynamic:
    {
      // Dynamic redistribution: whenever system is out of balance

      // This is the first time step (t=0) or restart
      if (first_time_step)
      {
        // Always perform rebalancing in the first time step
        perform_rebalancing = true;
      }

      // This is a regular time step (neither t=0 nor restart)
      else
      {
        /* Decide on redistribution
         *
         * We allow a maximum value of the balance measure in the system as defined in the input
         * parameter MAX_BALANCE, i.e. the maximum local processor workload and the minimum local
         * processor workload for mortar evaluation of all interfaces may not differ by more than
         * (MAX_BALANCE - 1.0)*100%)
         *
         * Moreover, we redistribute if in the majority of iteration steps of the last time step
         * there has been an unbalance in element distribution, i.e. if elements_average >= 0.5)
         */
        if (time_average >= max_time_unbalance || elements_average >= 0.5)
          perform_rebalancing = true;
      }

      break;
    }
  }

  PrintParallelBalanceIndicators(time_average, elements_average, max_time_unbalance);

  return perform_rebalancing;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::ComputeAndResetParallelBalanceIndicators(
    double& time_average, double& elements_average)
{
  dsassert(unbalanceEvaluationTime_.size() > 0, "Vector should have length > 0.");
  dsassert(unbalanceNumSlaveElements_.size() > 0, "Vector should have length > 0.");

  // compute average balance factors of last time step
  for (const auto& time : unbalanceEvaluationTime_) time_average += time;
  time_average /= static_cast<double>(unbalanceEvaluationTime_.size());
  for (const auto& num_elements : unbalanceNumSlaveElements_) elements_average += num_elements;
  elements_average /= static_cast<double>(unbalanceNumSlaveElements_.size());

  // Reset balance factors of last time step
  unbalanceEvaluationTime_.resize(0);
  unbalanceNumSlaveElements_.resize(0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::PrintParallelBalanceIndicators(
    double& time_average, double& elements_average, const double& max_time_unbalance) const
{
  // Screen output only on proc 0
  if (Comm().MyPID() == 0)
  {
    std::cout << "*************** DATA OF PREVIOUS TIME STEP ***************" << std::endl;
    if (time_average > 0)
    {
      std::cout << "Parallel balance (time): " << time_average << " (limit " << max_time_unbalance
                << ")\n"
                << "Parallel balance (eles): " << elements_average << " (limit 0.5)" << std::endl;
    }
    else
      std::cout << "Parallel balance: t=0/restart" << std::endl;
    std::cout << "**********************************************************" << std::endl;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::CoAbstractStrategy::IsUpdateOfGhostingNecessary(
    const INPAR::MORTAR::ExtendGhosting& ghosting_strategy, const bool first_time_step) const
{
  bool enforce_update_of_ghosting = false;
  switch (ghosting_strategy)
  {
    case INPAR::MORTAR::ExtendGhosting::redundant_all:
    case INPAR::MORTAR::ExtendGhosting::redundant_master:
    {
      // this is the first time step (t=0) or restart
      if (first_time_step)
        enforce_update_of_ghosting = true;
      else
        enforce_update_of_ghosting = false;

      break;
    }
    case INPAR::MORTAR::ExtendGhosting::roundrobin:
    case INPAR::MORTAR::ExtendGhosting::binning:
    {
      enforce_update_of_ghosting = true;
      break;
    }
    default:
    {
      dserror("Unknown strategy to extend ghosting if necessary.");
    }
  }

  return enforce_update_of_ghosting;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::CoAbstractStrategy::RedistributeContact(
    Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<const Epetra_Vector> vel)
{
  bool redistributed = false;

  if (CONTACT::UTILS::UseSafeRedistributeAndGhosting(Params()))
    redistributed = RedistributeWithSafeGhosting(*dis, *vel);
  else
  {
    if (Comm().MyPID() == 0)
    {
      std::cout << "+++++++++++++++++++++++++++++++ WARNING +++++++++++++++++++++++++++++++\n"
                << "+++ You're using an outdated contact redistribution implementation, +++\n"
                << "+++ that might deliver an insufficient master-side ghosting.        +++\n"
                << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
                << std::endl;
    }
    redistributed = RedistributeContactOld(dis, vel);
  }

  return redistributed;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::CoAbstractStrategy::RedistributeWithSafeGhosting(
    const Epetra_Vector& displacement, const Epetra_Vector& velocity)
{
  // time measurement
  Comm().Barrier();
  const double t_start = Teuchos::Time::wallTime();

  const INPAR::MORTAR::ExtendGhosting ghosting_strategy =
      Teuchos::getIntegralValue<INPAR::MORTAR::ExtendGhosting>(
          Params().sublist("PARALLEL REDISTRIBUTION"), "GHOSTING_STRATEGY");

  bool first_time_step = false;
  if (unbalanceEvaluationTime_.size() == 0 && unbalanceNumSlaveElements_.size() == 0)
    first_time_step = true;

  const bool perform_rebalancing = IsRebalancingNecessary(first_time_step);
  const bool enforce_ghosting_update =
      IsUpdateOfGhostingNecessary(ghosting_strategy, first_time_step);

  // Prepare for extending the ghosting
  ivel_.resize(Interfaces().size(), 0.0);  // initialize to zero for non-binning strategies
  if (ghosting_strategy == INPAR::MORTAR::ExtendGhosting::binning)
    CalcMeanVelocityForBinning(velocity);

  // Set old and current displacement state (needed for search within redistribution)
  if (perform_rebalancing)
  {
    SetState(MORTAR::state_new_displacement, displacement);
    SetState(MORTAR::state_old_displacement, displacement);
  }

  // Update parallel distribution and ghosting of all interfaces
  for (std::size_t i = 0; i < Interfaces().size(); ++i)
    Interfaces()[i]->UpdateParallelLayoutAndDataStructures(
        perform_rebalancing, enforce_ghosting_update, maxdof_, ivel_[i]);

  // Re-setup strategy to update internal map objects
  if (perform_rebalancing) Setup(true, false);

  // time measurement
  Comm().Barrier();
  double t_end = Teuchos::Time::wallTime() - t_start;
  if (Comm().MyPID() == 0)
    std::cout << "\nTime for parallel redistribution..............." << std::scientific
              << std::setprecision(6) << t_end << " secs\n"
              << std::endl;

  return perform_rebalancing;
}

/*----------------------------------------------------------------------*
 | parallel redistribution                                   popp 09/10 |
 *----------------------------------------------------------------------*/
bool CONTACT::CoAbstractStrategy::RedistributeContactOld(
    Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<const Epetra_Vector> vel)
{
  // decide whether redistribution should be applied or not
  bool first_time_step = false;
  if (unbalanceEvaluationTime_.size() == 0 && unbalanceNumSlaveElements_.size() == 0)
    first_time_step = true;
  const bool doredist = IsRebalancingNecessary(first_time_step);

  // get out of here if simulation is still in balance
  if (!doredist) return false;

  // time measurement
  Comm().Barrier();
  const double t_start = Teuchos::Time::wallTime();

  // Prepare for extending the ghosting
  ivel_.resize(Interfaces().size(), 0.0);  // initialize to zero for non-binning strategies
  if (Teuchos::getIntegralValue<INPAR::MORTAR::ExtendGhosting>(
          Params().sublist("PARALLEL REDISTRIBUTION"), "GHOSTING_STRATEGY") ==
      INPAR::MORTAR::ExtendGhosting::binning)
    CalcMeanVelocityForBinning(*vel);

  /* set old and current displacement state
   * (needed for search within redistribution) */
  SetState(MORTAR::state_new_displacement, *dis);
  SetState(MORTAR::state_old_displacement, *dis);

  // parallel redistribution of all interfaces
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    // redistribute optimally among procs
    Interfaces()[i]->Redistribute();

    // call fill complete again
    Interfaces()[i]->FillComplete(true, maxdof_, ivel_[i]);

    // print new parallel distribution
    if (Comm().MyPID() == 0)
      std::cout << "Interface parallel distribution after rebalancing:" << std::endl;
    Interfaces()[i]->PrintParallelDistribution();

    // re-create binary search tree
    Interfaces()[i]->CreateSearchTree();
  }

  // re-setup strategy with redistributed=TRUE, init=FALSE
  Setup(true, false);

  // time measurement
  Comm().Barrier();
  double t_end = Teuchos::Time::wallTime() - t_start;
  if (Comm().MyPID() == 0)
    std::cout << "\nTime for parallel redistribution..............." << std::scientific
              << std::setprecision(6) << t_end << " secs\n"
              << std::endl;

  return doredist;
}

/*----------------------------------------------------------------------*
 | setup this strategy object                                popp 08/10 |
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::Setup(bool redistributed, bool init)
{
  if (init)
  {
    // set potential global self contact status
    // (this is TRUE if at least one contact interface is a self contact interface)
    bool selfcontact = false;
    for (unsigned i = 0; i < Interfaces().size(); ++i)
      if (Interfaces()[i]->SelfContact()) selfcontact = true;

    if (selfcontact) isselfcontact_ = true;
  }

  // ------------------------------------------------------------------------
  // setup global accessible Epetra_Maps
  // ------------------------------------------------------------------------

  // make sure to remove all existing maps first
  // (do NOT remove map of non-interface dofs after redistribution)
  gsdofrowmap_ = Teuchos::null;
  gmdofrowmap_ = Teuchos::null;
  gsmdofrowmap_ = Teuchos::null;
  glmdofrowmap_ = Teuchos::null;
  gdisprowmap_ = Teuchos::null;
  gsnoderowmap_ = Teuchos::null;
  gmnoderowmap_ = Teuchos::null;
  gactivenodes_ = Teuchos::null;
  gactivedofs_ = Teuchos::null;
  ginactivenodes_ = Teuchos::null;
  ginactivedofs_ = Teuchos::null;
  gactiven_ = Teuchos::null;
  gactivet_ = Teuchos::null;
  if (!redistributed) gndofrowmap_ = Teuchos::null;
  if (init) initial_elecolmap_.clear();
  initial_elecolmap_.resize(0);

  if (friction_)
  {
    gslipnodes_ = Teuchos::null;
    gslipdofs_ = Teuchos::null;
    gslipt_ = Teuchos::null;
  }

  // initialize vertex, edge and surface maps for nonsmooth case
  if (DRT::INPUT::IntegralValue<int>(Params(), "NONSMOOTH_GEOMETRIES"))
  {
    gsdofVertex_ = Teuchos::null;
    gsdofEdge_ = Teuchos::null;
    gsdofSurf_ = Teuchos::null;
  }

  // make numbering of LM dofs consecutive and unique across N interfaces
  int offset_if = 0;

  // merge interface maps to global maps
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    // build Lagrange multiplier dof map
    if (IsSelfContact())
    {
      if (redistributed) dserror("SELF-CONTACT: Parallel redistribution is not supported!");

      CoInterface& inter = *Interfaces()[i];
      Teuchos::RCP<const Epetra_Map> refdofrowmap = Teuchos::null;
      if (inter.SelfContact())
        refdofrowmap = LINALG::MergeMap(inter.SlaveRowDofs(), inter.MasterRowDofs());
      else
        refdofrowmap = inter.SlaveRowDofs();

      Teuchos::RCP<Epetra_Map> selfcontact_lmmap =
          Interfaces()[i]->UpdateLagMultSets(offset_if, redistributed, *refdofrowmap);

      Teuchos::RCP<Epetra_Map>& gsc_refdofmap_ptr = Data().GSelfContactRefDofRowMapPtr();
      Teuchos::RCP<Epetra_Map>& gsc_lmdofmap_ptr = Data().GSelfContactLmDofRowMapPtr();
      gsc_lmdofmap_ptr = LINALG::MergeMap(selfcontact_lmmap, gsc_lmdofmap_ptr);
      gsc_refdofmap_ptr = LINALG::MergeMap(refdofrowmap, gsc_refdofmap_ptr);

      const int loffset_interface = selfcontact_lmmap->NumGlobalElements();
      if (loffset_interface > 0) offset_if += loffset_interface;
    }
    else
    {
      Interfaces()[i]->UpdateLagMultSets(offset_if, redistributed);
      const int loffset_interface = Interfaces()[i]->LagMultDofs()->NumGlobalElements();
      if (loffset_interface > 0) offset_if += loffset_interface;
    }

    // merge interface master, slave maps to global master, slave map
    gsnoderowmap_ = LINALG::MergeMap(SlRowNodesPtr(), Interfaces()[i]->SlaveRowNodes());
    gmnoderowmap_ = LINALG::MergeMap(MaRowNodesPtr(), Interfaces()[i]->MasterRowNodes());
    gsdofrowmap_ = LINALG::MergeMap(SlDoFRowMapPtr(true), Interfaces()[i]->SlaveRowDofs());
    gmdofrowmap_ = LINALG::MergeMap(gmdofrowmap_, Interfaces()[i]->MasterRowDofs());

    // merge active sets and slip sets of all interfaces
    // (these maps are NOT allowed to be overlapping !!!)
    Interfaces()[i]->BuildActiveSet(init);
    gactivenodes_ = LINALG::MergeMap(gactivenodes_, Interfaces()[i]->ActiveNodes(), false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_, Interfaces()[i]->ActiveDofs(), false);

    ginactivenodes_ = LINALG::MergeMap(ginactivenodes_, Interfaces()[i]->InActiveNodes(), false);
    ginactivedofs_ = LINALG::MergeMap(ginactivedofs_, Interfaces()[i]->InActiveDofs(), false);

    gactiven_ = LINALG::MergeMap(gactiven_, Interfaces()[i]->ActiveNDofs(), false);
    gactivet_ = LINALG::MergeMap(gactivet_, Interfaces()[i]->ActiveTDofs(), false);

    // store initial element col map for binning strategy
    initial_elecolmap_.push_back(
        Teuchos::rcp<Epetra_Map>(new Epetra_Map(*Interfaces()[i]->Discret().ElementColMap())));

    // ****************************************************
    // friction
    // ****************************************************
    if (friction_)
    {
      gslipnodes_ = LINALG::MergeMap(gslipnodes_, Interfaces()[i]->SlipNodes(), false);
      gslipdofs_ = LINALG::MergeMap(gslipdofs_, Interfaces()[i]->SlipDofs(), false);
      gslipt_ = LINALG::MergeMap(gslipt_, Interfaces()[i]->SlipTDofs(), false);
    }

    // define maps for nonsmooth case
    if (DRT::INPUT::IntegralValue<int>(Params(), "NONSMOOTH_GEOMETRIES"))
    {
      gsdofVertex_ = LINALG::MergeMap(gsdofVertex_, Interfaces()[i]->SdofVertexRowmap());
      gsdofEdge_ = LINALG::MergeMap(gsdofEdge_, Interfaces()[i]->SdofEdgeRowmap());
      gsdofSurf_ = LINALG::MergeMap(gsdofSurf_, Interfaces()[i]->SdofSurfRowmap());
    }
  }

  // create the global Lagrange multiplier DoF row map
  glmdofrowmap_ = CreateDeterministicLMDofRowMap(*gsdofrowmap_);

  // setup global non-slave-or-master dof map
  // (this is done by splitting from the discretization dof map)
  // (no need to rebuild this map after redistribution)
  if (!redistributed)
  {
    gndofrowmap_ = LINALG::SplitMap(*(ProblemDofs()), SlDoFRowMap(true));
    gndofrowmap_ = LINALG::SplitMap(*gndofrowmap_, *gmdofrowmap_);
  }

  // setup combined global slave and master dof map
  // setup global displacement dof map
  gsmdofrowmap_ = LINALG::MergeMap(SlDoFRowMap(true), *gmdofrowmap_, false);
  gdisprowmap_ = LINALG::MergeMap(*gndofrowmap_, *gsmdofrowmap_, false);

  // TODO: check if necessary!
  // due to boundary modification we have to extend master map to slave dofs
  //  if(DRT::INPUT::IntegralValue<int>(Params(),"NONSMOOTH_GEOMETRIES"))
  //  {
  //    gmdofrowmap_ = LINALG::MergeMap(SlDoFRowMap(true), *gmdofrowmap_, false);
  //  }

  // initialize flags for global contact status
  if (gactivenodes_->NumGlobalElements())
  {
    isincontact_ = true;
    wasincontact_ = true;
    wasincontactlts_ = true;
  }

  // ------------------------------------------------------------------------
  // setup global accessible vectors and matrices
  // ------------------------------------------------------------------------

  // initialize vectors and matrices
  if (!redistributed)
  {
    // setup Lagrange multiplier vectors
    z_ = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
    zincr_ = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
    zold_ = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
    zuzawa_ = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));

    // setup global mortar matrices Dold and Mold
    dold_ = Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true), 1, true, false));
    dold_->Zero();
    dold_->Complete();
    mold_ = Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true), 1, true, false));
    mold_->Zero();
    mold_->Complete(*gmdofrowmap_, SlDoFRowMap(true));
  }

  // In the redistribution case, first check if the vectors and
  // matrices have already been defined, If yes, transform them
  // to the new redistributed maps. If not, initialize them.
  // Moreover, store redistributed quantities into nodes!!!
  else
  {
    // setup Lagrange multiplier vectors
    if (z_ == Teuchos::null)
    {
      z_ = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
    }
    else
    {
      Teuchos::RCP<Epetra_Vector> newz = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
      LINALG::Export(*z_, *newz);
      z_ = newz;
    }

    if (zincr_ == Teuchos::null)
    {
      zincr_ = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
    }
    else
    {
      Teuchos::RCP<Epetra_Vector> newzincr = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
      LINALG::Export(*zincr_, *newzincr);
      zincr_ = newzincr;
    }

    if (zold_ == Teuchos::null)
      zold_ = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
    else
    {
      Teuchos::RCP<Epetra_Vector> newzold = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
      LINALG::Export(*zold_, *newzold);
      zold_ = newzold;
    }

    if (zuzawa_ == Teuchos::null)
      zuzawa_ = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
    else
    {
      Teuchos::RCP<Epetra_Vector> newzuzawa = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
      LINALG::Export(*zuzawa_, *newzuzawa);
      zuzawa_ = newzuzawa;
    }

    // setup global Mortar matrices Dold and Mold
    if (dold_ == Teuchos::null)
    {
      dold_ = Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true), 1, true, false));
      dold_->Zero();
      dold_->Complete();
    }
    else if (dold_->RowMap().NumGlobalElements() > 0)
      dold_ = MORTAR::MatrixRowColTransform(dold_, SlDoFRowMapPtr(true), SlDoFRowMapPtr(true));

    if (mold_ == Teuchos::null)
    {
      mold_ = Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true), 1, true, false));
      mold_->Zero();
      mold_->Complete(*gmdofrowmap_, SlDoFRowMap(true));
    }
    else if (mold_->RowMap().NumGlobalElements() > 0)
      mold_ = MORTAR::MatrixRowColTransform(mold_, SlDoFRowMapPtr(true), gmdofrowmap_);
  }

  // output contact stress vectors
  stressnormal_ = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
  stresstangential_ = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));

  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // These matrices need to be applied to the slave displacements
  // in the cases of dual LM interpolation for tet10/hex20 meshes
  // in 3D or for locally linear Lagrange multipliers for line3 meshes
  // in 2D. Here, the displacement basis functions have been modified
  // in order to assure positivity of the D matrix entries and at
  // the same time biorthogonality. Thus, to scale back the modified
  // discrete displacements \hat{d} to the nodal discrete displacements
  // {d}, we have to apply the transformation matrix T and vice versa
  // with the transformation matrix T^(-1).
  //----------------------------------------------------------------------
  INPAR::MORTAR::ShapeFcn shapefcn =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(), "LM_SHAPEFCN");
  INPAR::MORTAR::LagMultQuad lagmultquad =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(Params(), "LM_QUAD");
  if ((shapefcn == INPAR::MORTAR::shape_dual || shapefcn == INPAR::MORTAR::shape_petrovgalerkin) &&
      (Dim() == 3 || (Dim() == 2 && lagmultquad == INPAR::MORTAR::lagmult_lin)))
    for (int i = 0; i < (int)Interfaces().size(); ++i)
      dualquadslavetrafo_ += Interfaces()[i]->Quadslave();

  //----------------------------------------------------------------------
  // IF SO, COMPUTE TRAFO MATRIX AND ITS INVERSE
  //----------------------------------------------------------------------
  if (Dualquadslavetrafo())
  {
    // for locally linear Lagrange multipliers, consider both slave and master DOFs,
    // and otherwise, only consider slave DOFs
    if (lagmultquad == INPAR::MORTAR::lagmult_lin)
    {
      trafo_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsmdofrowmap_, 10));
      invtrafo_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsmdofrowmap_, 10));
    }
    else
    {
      trafo_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_, 10));
      invtrafo_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_, 10));
    }

    // set of already processed nodes
    // (in order to avoid double-assembly for N interfaces)
    std::set<int> donebefore;

    // for all interfaces
    for (int i = 0; i < (int)Interfaces().size(); ++i)
      Interfaces()[i]->AssembleTrafo(*trafo_, *invtrafo_, donebefore);

    // FillComplete() transformation matrices
    trafo_->Complete();
    invtrafo_->Complete();
  }

  // transform modified old D-matrix in case of friction
  // (ony necessary after parallel redistribution)
  if (redistributed && friction_ && Dualquadslavetrafo())
  {
    if (doldmod_ == Teuchos::null)
    {
      doldmod_ = Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true), 1, true, false));
      doldmod_->Zero();
      doldmod_->Complete();
    }
    else
      doldmod_ =
          MORTAR::MatrixRowColTransform(doldmod_, SlDoFRowMapPtr(true), SlDoFRowMapPtr(true));
  }

  if (init)
  {
    // store interface maps with parallel distribution of underlying
    // problem discretization (i.e. interface maps before parallel
    // redistribution of slave and master sides)
    if (ParRedist())
    {
      for (std::size_t i = 0; i < Interfaces().size(); ++i)
        Interfaces()[i]->StoreUnredistributedMaps();
      if (LMDoFRowMapPtr(true) != Teuchos::null)
        pglmdofrowmap_ = Teuchos::rcp(new Epetra_Map(LMDoFRowMap(true)));
      pgsdofrowmap_ = Teuchos::rcp(new Epetra_Map(SlDoFRowMap(true)));
      pgmdofrowmap_ = Teuchos::rcp(new Epetra_Map(*gmdofrowmap_));
      pgsmdofrowmap_ = Teuchos::rcp(new Epetra_Map(*gsmdofrowmap_));
    }
  }

  PostSetup(redistributed, init);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> CONTACT::CoAbstractStrategy::CreateDeterministicLMDofRowMap(
    const Epetra_Map& gsdofrowmap) const
{
  const unsigned num_my_sdofs = gsdofrowmap.NumMyElements();
  const int* my_sdof_gids = gsdofrowmap.MyGlobalElements();

  std::vector<int> my_lm_gids(num_my_sdofs, -1);

  for (unsigned slid = 0; slid < num_my_sdofs; ++slid)
  {
    const int sgid = my_sdof_gids[slid];

    // find slid of the interface map
    unsigned interface_id = 0;
    int interface_slid = -1;
    for (auto cit = Interfaces().begin(); cit != Interfaces().end(); ++cit, ++interface_id)
    {
      const CoInterface& interface = **cit;
      Teuchos::RCP<const Epetra_Map> sdof_map = interface.SlaveRowDofs();

      interface_slid = sdof_map->LID(sgid);
      if (interface_slid != -1) break;
    }

    if (interface_slid == -1)
      dserror(
          "Couldn't find the global slave dof id #%d in the local interface "
          "maps on proc #%d!",
          sgid, Comm().MyPID());

    // get the corresponding Lagrange Multiplier GID
    const int interface_lmgid = Interfaces()[interface_id]->LagMultDofs()->GID(interface_slid);
    if (interface_lmgid == -1)
      dserror(
          "Couldn't find the corresponding Lagrange multiplier GID! "
          "Note that the UpdateLagMultSets() must be called on each interface "
          "beforehand.");

    my_lm_gids[slid] = interface_lmgid;
  }
  return Teuchos::rcp(
      new Epetra_Map(-1, static_cast<int>(my_lm_gids.size()), &my_lm_gids[0], 0, Comm()));
}


/*----------------------------------------------------------------------*
 | global evaluation method called from time integrator      popp 06/09 |
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
    Teuchos::RCP<LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f, const int timeStep,
    const int nonlinearIteration, bool predictor)
{
  // update step and iteration counters
  step_ = timeStep;
  iter_ = nonlinearIteration;

  // Create timing reports?
  bool doAccurateTimeMeasurements =
      DRT::INPUT::IntegralValue<bool>(Data().SContact(), "TIMING_DETAILS");

  if (doAccurateTimeMeasurements)
  {
    // mortar initialization and evaluation
    Comm().Barrier();
    const double t_start1 = Teuchos::Time::wallTime();
    SetState(MORTAR::state_new_displacement, *dis);
    Comm().Barrier();
    const double t_end1 = Teuchos::Time::wallTime() - t_start1;


    Comm().Barrier();
    const double t_start2 = Teuchos::Time::wallTime();
    //---------------------------------------------------------------
    // For selfcontact the master/slave sets are updated within the -
    // contact search, see SelfBinaryTree.                          -
    // Therefore, we have to initialize the mortar matrices after   -
    // interface evaluations.                                       -
    //---------------------------------------------------------------
    if (IsSelfContact())
    {
      InitEvalInterface();  // evaluate mortar terms (integrate...)
      InitMortar();         // initialize mortar matrices and vectors
      AssembleMortar();     // assemble mortar terms into global matrices
    }
    else
    {
      InitMortar();         // initialize mortar matrices and vectors
      InitEvalInterface();  // evaluate mortar terms (integrate...)
      AssembleMortar();     // assemble mortar terms into global matrices
    }
    Comm().Barrier();
    const double t_end2 = Teuchos::Time::wallTime() - t_start2;

    // evaluate relative movement for friction
    Comm().Barrier();
    const double t_start3 = Teuchos::Time::wallTime();
    if (predictor)
      EvaluateRelMovPredict();
    else
      EvaluateRelMov();

    // update active set
    if (!predictor) UpdateActiveSetSemiSmooth();

    Comm().Barrier();
    const double t_end3 = Teuchos::Time::wallTime() - t_start3;

    // apply contact forces and stiffness
    Comm().Barrier();
    const double t_start4 = Teuchos::Time::wallTime();
    Initialize();          // init lin-matrices
    Evaluate(kt, f, dis);  // assemble lin. matrices, condensation ...
    EvalConstrRHS();       // evaluate the constraint rhs (saddle-point system only)

    Comm().Barrier();
    const double t_end4 = Teuchos::Time::wallTime() - t_start4;

    // only for debugging:
    InterfaceForces();

    if (Comm().MyPID() == 0)
    {
      std::cout << "    -->setstate :\t" << t_end1 << " seconds\n";
      std::cout << "    -->interface eval. :\t" << t_end2 << " seconds\n";
      std::cout << "    -->update active set :\t" << t_end3 << " seconds\n";
      std::cout << "    -->modify global system :\t" << t_end4 << " seconds\n";
    }
  }
  else
  {
    // mortar initialization and evaluation
    SetState(MORTAR::state_new_displacement, *dis);

    //---------------------------------------------------------------
    // For selfcontact the master/slave sets are updated within the -
    // contact search, see SelfBinaryTree.                          -
    // Therefore, we have to initialize the mortar matrices after   -
    // interface evaluations.                                       -
    //---------------------------------------------------------------
    if (IsSelfContact())
    {
      InitEvalInterface();  // evaluate mortar terms (integrate...)
      InitMortar();         // initialize mortar matrices and vectors
      AssembleMortar();     // assemble mortar terms into global matrices
    }
    else
    {
      InitMortar();         // initialize mortar matrices and vectors
      InitEvalInterface();  // evaluate mortar terms (integrate...)
      AssembleMortar();     // assemble mortar terms into global matrices
    }

    // evaluate relative movement for friction
    if (predictor)
      EvaluateRelMovPredict();
    else
      EvaluateRelMov();

    // update active set
    if (!predictor) UpdateActiveSetSemiSmooth();

    // apply contact forces and stiffness
    Initialize();          // init lin-matrices
    Evaluate(kt, f, dis);  // assemble lin. matrices, condensation ...
    EvalConstrRHS();       // evaluate the constraint rhs (saddle-point system only)

    // only for debugging:
    InterfaceForces();
  }

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::SetState(
    const enum MORTAR::StateType& statetype, const Epetra_Vector& vec)
{
  switch (statetype)
  {
    case MORTAR::state_new_displacement:
    case MORTAR::state_old_displacement:
    {
      // set state on interfaces
      for (int i = 0; i < (int)Interfaces().size(); ++i) Interfaces()[i]->SetState(statetype, vec);
      break;
    }
    default:
    {
      dserror(
          "Unsupported state type! (state type = %s)", MORTAR::StateType2String(statetype).c_str());
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | update global master and slave sets (public)               popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::UpdateGlobalSelfContactState()
{
  if (not IsSelfContact()) return;

  // reset global slave / master Epetra Maps
  gsnoderowmap_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
  gsdofrowmap_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
  gmdofrowmap_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
  glmdofrowmap_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));

  // make numbering of LM dofs consecutive and unique across N interfaces
  int offset_if = 0;

  // setup global slave / master Epetra_Maps
  // (this is done by looping over all interfaces and merging)
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    // build Lagrange multiplier dof map
    Interfaces()[i]->UpdateSelfContactLagMultSet(GSelfContactLmMap(), *gsmdofrowmap_);

    // merge interface Lagrange multiplier dof maps to global LM dof map
    glmdofrowmap_ = LINALG::MergeMap(LMDoFRowMapPtr(true), Interfaces()[i]->LagMultDofs());
    offset_if = LMDoFRowMap(true).NumGlobalElements();
    if (offset_if < 0) offset_if = 0;

    // merge interface master, slave maps to global master, slave map
    gsnoderowmap_ = LINALG::MergeMap(SlRowNodesPtr(), Interfaces()[i]->SlaveRowNodes());
    gsdofrowmap_ = LINALG::MergeMap(SlDoFRowMapPtr(true), Interfaces()[i]->SlaveRowDofs());
    gmdofrowmap_ = LINALG::MergeMap(gmdofrowmap_, Interfaces()[i]->MasterRowDofs());
  }

  Teuchos::RCP<Epetra_Vector> tmp_ptr = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));

  {
    const int* oldgids = zincr_->Map().MyGlobalElements();
    for (int i = 0; i < zincr_->Map().NumMyElements(); ++i)
    {
      if (std::abs((*zincr_)[i]) > std::numeric_limits<double>::epsilon())
      {
        const int new_lid = gsdofrowmap_->LID(oldgids[i]);
        if (new_lid == -1)
          dserror(
              "Self contact: The Lagrange multiplier increment vector "
              "could not be transferred consistently.");
        else
          (*tmp_ptr)[new_lid] = (*zincr_)[i];
      }
    }
    zincr_ = Teuchos::rcp(new Epetra_Vector(*tmp_ptr));
  }

  tmp_ptr->PutScalar(0.0);
  {
    const int* oldgids = z_->Map().MyGlobalElements();
    for (int i = 0; i < z_->Map().NumMyElements(); ++i)
    {
      if (std::abs((*z_)[i]) > std::numeric_limits<double>::epsilon())
      {
        const int new_lid = gsdofrowmap_->LID(oldgids[i]);
        if (new_lid == -1)
          dserror(
              "Self contact: The Lagrange multiplier vector "
              "could not be transferred consistently.");
        else
          (*tmp_ptr)[new_lid] = (*z_)[i];
      }
    }
    z_ = tmp_ptr;
  }

  return;
}

/*----------------------------------------------------------------------*
 | Calculate mean. vel. for bin size                         farah 11/13|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::CalcMeanVelocityForBinning(const Epetra_Vector& velocity)
{
  ivel_.clear();
  ivel_.resize(0);

  // for dynamic problems
  if (alphaf_ != 0.0)  // ToDo Improve this check! alphaf_ = 0 does not guarantee Statics.
  {
    // create vector of interface velocities
    for (int i = 0; i < (int)Interfaces().size(); ++i)
    {
      // interface node map
      Teuchos::RCP<Epetra_Vector> interfaceVelocity =
          Teuchos::rcp(new Epetra_Vector(*(Interfaces()[i]->Discret().DofRowMap())));
      LINALG::Export(velocity, *interfaceVelocity);

      double meanVelocity = 0.0;

      int err = interfaceVelocity->MeanValue(&meanVelocity);
      if (err)
        dserror("Calculation of mean velocity for interface %s failed.",
            Interfaces()[i]->Discret().Name().c_str());
      meanVelocity = abs(meanVelocity);

      ivel_.push_back(meanVelocity);
    }
  }
  // static problems
  else
  {
    // TODO: should be fixed for static problems!
    dserror(
        "Binning Strategy is only recommended for dynamic problems! Please use a different "
        "Parallel Strategy!");
  }
  return;
}

/*----------------------------------------------------------------------*
 | initialize + evaluate interface for next Newton step       popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::InitEvalInterface(
    Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr)
{
  // time measurement (on each processor)
  const double t_start = Teuchos::Time::wallTime();

  // get type of parallel strategy
  const Teuchos::ParameterList& mortarParallelRedistParams =
      Params().sublist("PARALLEL REDISTRIBUTION");
  INPAR::MORTAR::ExtendGhosting extendghosting =
      Teuchos::getIntegralValue<INPAR::MORTAR::ExtendGhosting>(
          mortarParallelRedistParams, "GHOSTING_STRATEGY");

  // Evaluation for all interfaces
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    // initialize / reset interfaces
    Interfaces()[i]->Initialize();

    // store required integration time
    inttime_ += Interfaces()[i]->Inttime();

    switch (extendghosting)
    {
      case INPAR::MORTAR::ExtendGhosting::roundrobin:
      {
        // first perform rrloop to detect the required ghosting
        Interfaces()[i]->RoundRobinDetectGhosting();

        // second step --> evaluate
        Interfaces()[i]->Evaluate(0, step_, iter_);
        break;
      }
      case INPAR::MORTAR::ExtendGhosting::binning:
      {
        // required master elements are already ghosted (preparestepcontact) !!!
        // call evaluation
        Interfaces()[i]->Evaluate(0, step_, iter_);
        break;
      }
      case INPAR::MORTAR::ExtendGhosting::redundant_all:
      case INPAR::MORTAR::ExtendGhosting::redundant_master:
      {
        Interfaces()[i]->Evaluate(0, step_, iter_);
        break;
      }
    }
  }  // end interface loop

  // check the parallel distribution
  CheckParallelDistribution(t_start);

  //**********************************************************************
  // OVERVIEW OF PARALLEL MORTAR COUPLING STATUS
  //**********************************************************************
#ifdef CONTACTSTATUS
  // total numbers per processor
  std::vector<int> smpairs(1);
  std::vector<int> smintpairs(1);
  std::vector<int> intcells(1);

  // add numbers of all interfaces
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    smpairs[0] += Interfaces()[i]->SlaveMasterPairs();
    smintpairs[0] += Interfaces()[i]->SlaveMasterIntPairs();
    intcells[0] += Interfaces()[i]->IntegrationCells();
  }

  // vector containing all proc ids
  const int numproc = Comm().NumProc();
  std::vector<int> allproc(numproc);
  for (int i = 0; i < numproc; ++i) allproc[i] = i;

  // global numbers
  std::vector<int> gsmpairs, gsmintpairs, gintcells;
  LINALG::Gather<int>(smpairs, gsmpairs, numproc, &allproc[0], Comm());
  LINALG::Gather<int>(smintpairs, gsmintpairs, numproc, &allproc[0], Comm());
  LINALG::Gather<int>(intcells, gintcells, numproc, &allproc[0], Comm());

  // output to screen
  if (Comm().MyPID() == 0)
  {
    std::cout << "--------------------------------------------------------------------------------"
              << std::endl;
    std::cout << std::setw(10) << "proc ID" << std::setw(16) << "# s/m pairs" << std::setw(16)
              << "# s/m intpairs" << std::setw(16) << "# intcells" << std::endl;
    for (int i = 0; i < numproc; ++i)
    {
      std::cout << std::setw(10) << i << std::setw(16) << gsmpairs[i] << std::setw(16)
                << gsmintpairs[i] << std::setw(16) << gintcells[i] << std::endl;
    }
    std::cout << "--------------------------------------------------------------------------------"
              << std::endl;
  }
#endif  // #ifdef CONTACTSTATUS
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::CheckParallelDistribution(const double& t_start)
{
  const double my_total_time = Teuchos::Time::wallTime() - t_start;
  UpdateParallelDistributionStatus(my_total_time);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::UpdateParallelDistributionStatus(const double& my_total_time)
{
  //**********************************************************************
  // PARALLEL REDISTRIBUTION
  //**********************************************************************
  // don't do this if this is a single processor (serial) job
  if (Comm().NumProc() == 1) return;

  // collect information about participation in coupling evaluation
  // and in parallel distribution of the individual interfaces
  std::vector<int> numloadele((int)Interfaces().size());
  std::vector<int> numcrowele((int)Interfaces().size());
  for (int i = 0; i < (int)Interfaces().size(); ++i)
    Interfaces()[i]->CollectDistributionData(numloadele[i], numcrowele[i]);

  // time measurement (on each processor)
  double t_end_for_minall = my_total_time;
  double t_end_for_maxall = my_total_time;

  // restrict time measurement to procs that own at least some part
  // of the "close" slave interface section(s) on the global level,
  // i.e. restrict to procs that actually have to do some work
  int gnumloadele = 0;
  for (int i = 0; i < (int)numloadele.size(); ++i) gnumloadele += numloadele[i];

  // for non-loaded procs, set time measurement to values 0.0 / 1.0e12,
  // which do not affect the maximum and minimum identification
  if (gnumloadele == 0)
  {
    t_end_for_minall = 1.0e12;
    t_end_for_maxall = 0.0;
  }

  // store time indicator for parallel redistribution
  // (indicator is the maximum local processor time
  // divided by the minimum local processor time)
  double maxall = 0.0;
  double minall = 0.0;
  Comm().MaxAll(&t_end_for_maxall, &maxall, 1);
  Comm().MinAll(&t_end_for_minall, &minall, 1);

  // check for plausibility before storing
  if (maxall == 0.0 && minall == 1.0e12)
    Data().UnbalanceTimeFactors().push_back(1.0);
  else
    Data().UnbalanceTimeFactors().push_back(maxall / minall);

  // obtain info whether there is an unbalance in element distribution
  bool eleunbalance = false;
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    // find out how many close slave elements in total
    int totrowele = 0;
    Comm().SumAll(&numcrowele[i], &totrowele, 1);

    // find out how many procs have work on this interface
    int lhascrowele = 0;
    int ghascrowele = 0;
    if (numcrowele[i] > 0) lhascrowele = 1;
    Comm().SumAll(&lhascrowele, &ghascrowele, 1);

    // minimum number of elements per proc
    int minele = Params().sublist("PARALLEL REDISTRIBUTION").get<int>("MIN_ELEPROC");
    int numproc = Comm().NumProc();

    //--------------------------------------------------------------------
    // check if there is an element unbalance
    //--------------------------------------------------------------------
    // CASE 0: if minimum number of elements per proc is zero, but
    // further procs are still available and more than numproc elements
    if ((minele == 0) && (totrowele > numproc) && (ghascrowele < numproc)) eleunbalance = true;

    // CASE 1: in total too few close slave elements but more than one
    // proc is active (otherwise, i.e. if interface small, we have no choice)
    if ((minele > 0) && (totrowele < ghascrowele * minele) && (ghascrowele > 1))
      eleunbalance = true;

    // CASE 2: in total too many close slave elements, but further procs
    // are still available for redsitribution
    if ((minele > 0) && (totrowele >= (ghascrowele + 1) * minele) && (ghascrowele < numproc))
      eleunbalance = true;
  }

  // obtain global info on element unbalance
  int geleunbalance = 0;
  int leleunbalance = (int)(eleunbalance);
  Comm().SumAll(&leleunbalance, &geleunbalance, 1);
  if (geleunbalance > 0)
    Data().UnbalanceElementFactors().push_back(1);
  else
    Data().UnbalanceElementFactors().push_back(0);

  // debugging output
  // std::cout << "PROC: " << Comm().MyPID() << "\t LOADELE: " << numloadele[0] << "\t ROWELE: " <<
  // numcrowele[0]
  //     << "\t MIN: " << minall << "\t MAX: " << maxall
  //     << "\t tmin: " << t_end_for_minall << "\t tmax: " << t_end_for_maxall
  //     << "\t TUNBALANCE: " << unbalanceEvaluationTime_[(int)unbalanceEvaluationTime_.size()-1]
  //     << "\t EUNBALANCE: " <<
  //     unbalanceNumSlaveElements_[(int)unbalanceNumSlaveElements_.size()-1] << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 | initialize mortar stuff for next Newton step               popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::InitMortar()
{
  // for self contact, slave and master sets may have changed,
  // thus we have to update them before initializing D,M etc.
  UpdateGlobalSelfContactState();

  // initialize Dold and Mold if not done already
  if (dold_ == Teuchos::null)
  {
    dold_ = Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true), 10));
    dold_->Zero();
    dold_->Complete();
  }
  if (mold_ == Teuchos::null)
  {
    mold_ = Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true), 100));
    mold_->Zero();
    mold_->Complete(*gmdofrowmap_, SlDoFRowMap(true));
  }

  // (re)setup global Mortar LINALG::SparseMatrices and Epetra_Vectors
  dmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true), 10));
  mmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true), 100));

  if (constr_direction_ == INPAR::CONTACT::constr_xyz)
    g_ = LINALG::CreateVector(SlDoFRowMap(true), true);
  else if (constr_direction_ == INPAR::CONTACT::constr_ntt)
    g_ = LINALG::CreateVector(SlRowNodes(), true);
  else
    dserror("unknown contact constraint direction");

  // in the case of frictional dual quad 3D, also the modified D matrices are setup
  if (friction_ && Dualquadslavetrafo())
  {
    // initialize Dold and Mold if not done already
    if (doldmod_ == Teuchos::null)
    {
      doldmod_ = Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true), 10));
      doldmod_->Zero();
      doldmod_->Complete();
    }
    // setup of dmatrixmod_
    dmatrixmod_ = Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true), 10));
  }

  return;
}

/*----------------------------------------------------------------------*
 | Assemble mortar stuff for next Newton step                 popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::AssembleMortar()
{
  // for all interfaces
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    // assemble D-, M-matrix and g-vector, store them globally
    Interfaces()[i]->AssembleDM(*dmatrix_, *mmatrix_);
    Interfaces()[i]->AssembleG(*g_);

#ifdef CONTACTFDNORMAL
    // FD check of normal derivatives
    std::cout << " -- CONTACTFDNORMAL- -----------------------------------" << std::endl;
    //    Interfaces()[i]->FDCheckNormalDeriv();
    Interfaces()[i]->FDCheckNormalCPPDeriv();
    std::cout << " -- CONTACTFDNORMAL- -----------------------------------" << std::endl;
#endif  // #ifdef CONTACTFDNORMAL
#ifdef CONTACTFDMORTARD
    // FD check of Mortar matrix D derivatives
    std::cout << " -- CONTACTFDMORTARD -----------------------------------" << std::endl;
    dmatrix_->Complete();
    if (dmatrix_->NormOne()) Interfaces()[i]->FDCheckMortarDDeriv();
    dmatrix_->UnComplete();
    std::cout << " -- CONTACTFDMORTARD -----------------------------------" << std::endl;
#endif  // #ifdef CONTACTFDMORTARD
#ifdef CONTACTFDMORTARM
    // FD check of Mortar matrix M derivatives
    std::cout << " -- CONTACTFDMORTARM -----------------------------------" << std::endl;
    mmatrix_->Complete(*gmdofrowmap_, *gsdofrowmap_);
    if (mmatrix_->NormOne()) Interfaces()[i]->FDCheckMortarMDeriv();
    mmatrix_->UnComplete();
    std::cout << " -- CONTACTFDMORTARM -----------------------------------" << std::endl;
#endif  // #ifdef CONTACTFDMORTARM
  }

  // FillComplete() global Mortar matrices
  dmatrix_->Complete();
  mmatrix_->Complete(*gmdofrowmap_, SlDoFRowMap(true));

  return;
}

/*----------------------------------------------------------------------*
 | evaluate reference state                               gitterle 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::EvaluateReferenceState(Teuchos::RCP<const Epetra_Vector> vec)
{
  // flag for initialization of contact with nodal gaps
  bool initcontactbygap = DRT::INPUT::IntegralValue<int>(Params(), "INITCONTACTBYGAP");

  // only do something for frictional case
  // or for initialization of initial contact set with nodal gap
  if (!friction_ and !initcontactbygap) return;

  // set state and do mortar calculation
  SetState(MORTAR::state_new_displacement, *vec);
  InitMortar();
  InitEvalInterface();
  AssembleMortar();

  // (1) GAP INITIALIZATION CASE
  // initialize init contact with nodal gap
  if (initcontactbygap)
  {
    // merge interface maps to global maps
    for (int i = 0; i < (int)Interfaces().size(); ++i)
    {
      // merge active sets and slip sets of all interfaces
      // (these maps are NOT allowed to be overlapping !!!)
      Interfaces()[i]->BuildActiveSet(true);
      gactivenodes_ = LINALG::MergeMap(gactivenodes_, Interfaces()[i]->ActiveNodes(), false);
      gactivedofs_ = LINALG::MergeMap(gactivedofs_, Interfaces()[i]->ActiveDofs(), false);
      gactiven_ = LINALG::MergeMap(gactiven_, Interfaces()[i]->ActiveNDofs(), false);
      gactivet_ = LINALG::MergeMap(gactivet_, Interfaces()[i]->ActiveTDofs(), false);

      if (friction_)
      {
        gslipnodes_ = LINALG::MergeMap(gslipnodes_, Interfaces()[i]->SlipNodes(), false);
        gslipdofs_ = LINALG::MergeMap(gslipdofs_, Interfaces()[i]->SlipDofs(), false);
        gslipt_ = LINALG::MergeMap(gslipt_, Interfaces()[i]->SlipTDofs(), false);
      }
    }

    // initialize flags for global contact status
    if (gactivenodes_->NumGlobalElements())
    {
      isincontact_ = true;
      wasincontact_ = true;
      wasincontactlts_ = true;
    }

    // error if no nodes are initialized to active
    if (gactivenodes_->NumGlobalElements() == 0)
      dserror("ERROR: No active nodes: Choose bigger value for INITCONTACTGAPVALUE!");
  }

  // (2) FRICTIONAL CONTACT CASE
  // do some friction stuff
  if (friction_)
  {
    // store contact state to contact nodes (active or inactive)
    StoreNodalQuantities(MORTAR::StrategyBase::activeold);

    // store D and M to old ones
    StoreDM("old");

    // store nodal entries from D and M to old ones
    StoreToOld(MORTAR::StrategyBase::dm);

    // store nodal normals
    StoreToOld(MORTAR::StrategyBase::n_old);

    // transform dold_ in the case of dual quadratic 3d
    if (Dualquadslavetrafo())
    {
      Teuchos::RCP<LINALG::SparseMatrix> tempold =
          LINALG::MLMultiply(*dold_, false, *invtrafo_, false, false, false, true);
      doldmod_ = tempold;
    }

    // evaluate relative movement
    // needed because it is not called in the predictor of the
    // lagrange multiplier strategy
    EvaluateRelMov();
  }

  // reset unbalance factors for redistribution
  // (since the interface has been evaluated once above)
  unbalanceEvaluationTime_.resize(0);
  unbalanceNumSlaveElements_.resize(0);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate relative movement of contact bodies           gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::EvaluateRelMov()
{
  // only for fricional contact
  if (!friction_) return;

  // transformation of slave displacement dofs
  // Dmod       ---->   D * T^(-1)
  if (Dualquadslavetrafo())
  {
    Teuchos::RCP<LINALG::SparseMatrix> temp =
        LINALG::MLMultiply(*dmatrix_, false, *invtrafo_, false, false, false, true);
    dmatrixmod_ = temp;
  }

  // vector of slave coordinates xs
  Teuchos::RCP<Epetra_Vector> xsmod = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));

  for (int i = 0; i < (int)Interfaces().size(); ++i) Interfaces()[i]->AssembleSlaveCoord(xsmod);

  // in case of 3D dual quadratic case, slave coordinates xs are modified
  if (Dualquadslavetrafo()) invtrafo_->Apply(*xsmod, *xsmod);

  // ATTENTION: for EvaluateRelMov() we need the vector xsmod in
  // fully overlapping layout. Thus, export here. First, allreduce
  // slave dof row map to obtain fully overlapping slave dof map.
  Teuchos::RCP<Epetra_Map> fullsdofs = LINALG::AllreduceEMap(SlDoFRowMap(true));
  Teuchos::RCP<Epetra_Vector> xsmodfull = Teuchos::rcp(new Epetra_Vector(*fullsdofs));
  LINALG::Export(*xsmod, *xsmodfull);
  xsmod = xsmodfull;

  // evaluation of obj. invariant slip increment
  // do the evaluation on the interface
  // loop over all slave row nodes on the current interface
  if (DRT::INPUT::IntegralValue<int>(Params(), "GP_SLIP_INCR") == false)
    for (int i = 0; i < (int)Interfaces().size(); ++i)
      Interfaces()[i]->EvaluateRelMov(xsmod, dmatrixmod_, doldmod_);

  return;
}

/*----------------------------------------------------------------------*
 | call appropriate evaluate for contact evaluation           popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::Evaluate(Teuchos::RCP<LINALG::SparseOperator>& kteff,
    Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis)
{
  // treat frictional and frictionless cases differently
  if (friction_)
    EvaluateFriction(kteff, feff);
  else
    EvaluateContact(kteff, feff);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate matrix of normals (for VelocityUpdate)            popp 10/11|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> CONTACT::CoAbstractStrategy::EvaluateNormals(
    Teuchos::RCP<Epetra_Vector> dis)
{
  // set displacement state and evaluate nodal normals
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    Interfaces()[i]->SetState(MORTAR::state_new_displacement, *dis);
    Interfaces()[i]->EvaluateNodalNormals();
  }

  // create empty global matrix
  // (rectangular: rows=snodes, cols=sdofs)
  Teuchos::RCP<LINALG::SparseMatrix> normals =
      Teuchos::rcp(new LINALG::SparseMatrix(SlRowNodes(), 3));

  // assemble nodal normals
  for (int i = 0; i < (int)Interfaces().size(); ++i) Interfaces()[i]->AssembleNormals(*normals);

  // complete global matrix
  // (rectangular: rows=snodes, cols=sdofs)
  normals->Complete(SlDoFRowMap(true), SlRowNodes());

  return normals;
}

/*----------------------------------------------------------------------*
 |  Store Lagrange mulitpliers and disp. jumps into CNode     popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::StoreNodalQuantities(MORTAR::StrategyBase::QuantityType type)
{
  // loop over all interfaces
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    // get global quantity to be stored in nodes
    Teuchos::RCP<Epetra_Vector> vectorglobal = Teuchos::null;

    // start type switch
    switch (type)
    {
      case MORTAR::StrategyBase::lmold:
      {
        vectorglobal = LagrMultOld();
        break;
      }
      case MORTAR::StrategyBase::lmcurrent:
      case MORTAR::StrategyBase::lmupdate:
      {
        vectorglobal = LagrMult();
        break;
      }
      case MORTAR::StrategyBase::lmuzawa:
      {
        vectorglobal = LagrMultUzawa();
        break;
      }
      case MORTAR::StrategyBase::activeold:
      case MORTAR::StrategyBase::slipold:
      {
        break;
      }
      default:
        dserror("ERROR: StoreNodalQuantities: Unknown state std::string variable!");
        break;
    }  // switch

    // slave dof and node map of the interface
    // columnmap for current or updated LM
    // rowmap for remaining cases
    Teuchos::RCP<Epetra_Map> sdofmap, snodemap;
    if (type == MORTAR::StrategyBase::lmupdate or type == MORTAR::StrategyBase::lmcurrent)
    {
      sdofmap = Interfaces()[i]->SlaveColDofs();
      snodemap = Interfaces()[i]->SlaveColNodes();
    }
    else
    {
      sdofmap = Interfaces()[i]->SlaveRowDofs();
      snodemap = Interfaces()[i]->SlaveRowNodes();
    }

    // export global quantity to current interface slave dof map (column or row)
    Teuchos::RCP<Epetra_Vector> vectorinterface = Teuchos::null;
    vectorinterface = Teuchos::rcp(new Epetra_Vector(*sdofmap));
    if (vectorglobal != Teuchos::null)  // necessary for case "activeold" and wear
      LINALG::Export(*vectorglobal, *vectorinterface);

    // loop over all slave nodes (column or row) on the current interface
    for (int j = 0; j < snodemap->NumMyElements(); ++j)
    {
      int gid = snodemap->GID(j);
      DRT::Node* node = Interfaces()[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      // be aware of problem dimension
      const int dim = Dim();
      const int numdof = cnode->NumDof();
      if (dim != numdof) dserror("ERROR: Inconsisteny Dim <-> NumDof");

      // find indices for DOFs of current node in Epetra_Vector
      // and extract this node's quantity from vectorinterface
      std::vector<int> locindex(dim);

      for (int dof = 0; dof < dim; ++dof)
      {
        locindex[dof] = (vectorinterface->Map()).LID(cnode->Dofs()[dof]);
        if (locindex[dof] < 0) dserror("ERROR: StoreNodalQuantites: Did not find dof in map");

        switch (type)
        {
          case MORTAR::StrategyBase::lmcurrent:
          {
            cnode->MoData().lm()[dof] = (*vectorinterface)[locindex[dof]];
            break;
          }
          case MORTAR::StrategyBase::lmold:
          {
            cnode->MoData().lmold()[dof] = (*vectorinterface)[locindex[dof]];
            break;
          }
          case MORTAR::StrategyBase::lmuzawa:
          {
            cnode->MoData().lmuzawa()[dof] = (*vectorinterface)[locindex[dof]];
            break;
          }
          case MORTAR::StrategyBase::lmupdate:
          {
#ifndef CONTACTPSEUDO2D
            // throw a dserror if node is Active and DBC
            if (cnode->IsDbc() && cnode->Active())
              dserror("ERROR: Slave node %i is active AND carries D.B.C.s!", cnode->Id());
#endif  // #ifndef CONTACTPSEUDO2D

            // store updated LM into node
            cnode->MoData().lm()[dof] = (*vectorinterface)[locindex[dof]];
            break;
          }
          case MORTAR::StrategyBase::activeold:
          {
            cnode->CoData().ActiveOld() = cnode->Active();
            break;
          }
          case MORTAR::StrategyBase::slipold:
          {
            if (!friction_) dserror("ERROR: Slip just for friction problems!");

            FriNode* fnode = dynamic_cast<FriNode*>(cnode);
            fnode->FriData().SlipOld() = fnode->FriData().Slip();
            break;
          }
          default:
            dserror("ERROR: StoreNodalQuantities: Unknown state std::string variable!");
            break;
        }  // switch
      }
    }  // end slave loop
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Output vector of normal/tang. contact stresses        gitterle 08/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::ComputeContactStresses()
{
  // reset contact stress class variables
  stressnormal_ = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
  stresstangential_ = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));

  // loop over all interfaces
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    // loop over all slave row nodes on the current interface
    for (int j = 0; j < Interfaces()[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = Interfaces()[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = Interfaces()[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      // be aware of problem dimension
      int dim = Dim();
      int numdof = cnode->NumDof();
      if (dim != numdof) dserror("ERROR: Inconsisteny Dim <-> NumDof");

      double nn[3];
      double nt1[3];
      double nt2[3];
      double lmn = 0.0;
      double lmt1 = 0.0;
      double lmt2 = 0.0;

      for (int j = 0; j < 3; ++j)
      {
        nn[j] = cnode->MoData().n()[j];
        nt1[j] = cnode->CoData().txi()[j];
        nt2[j] = cnode->CoData().teta()[j];
        lmn += nn[j] * cnode->MoData().lm()[j];
        lmt1 += nt1[j] * cnode->MoData().lm()[j];
        lmt2 += nt2[j] * cnode->MoData().lm()[j];
      }

      // find indices for DOFs of current node in Epetra_Vector
      // and put node values (normal and tangential stress components) at these DOFs

      std::vector<int> locindex(dim);

      // normal stress components
      for (int dof = 0; dof < dim; ++dof)
      {
        locindex[dof] = (stressnormal_->Map()).LID(cnode->Dofs()[dof]);
        (*stressnormal_)[locindex[dof]] = -lmn * nn[dof];
      }

      // tangential stress components
      for (int dof = 0; dof < dim; ++dof)
      {
        locindex[dof] = (stresstangential_->Map()).LID(cnode->Dofs()[dof]);
        (*stresstangential_)[locindex[dof]] = -lmt1 * nt1[dof] - lmt2 * nt2[dof];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Store dirichlet B.C. status into CNode                    popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::StoreDirichletStatus(
    Teuchos::RCP<const LINALG::MapExtractor> dbcmaps)
{
  // loop over all interfaces
  for (unsigned i = 0; i < Interfaces().size(); ++i)
  {
    // currently this only works safely for 1 interface
    // if (i>0) dserror("ERROR: StoreDirichletStatus: Double active node check needed for n
    // interfaces!");

    // loop over all slave row nodes on the current interface
    for (int j = 0; j < Interfaces()[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = Interfaces()[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = Interfaces()[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      // check if this node's dofs are in dbcmap
      for (int k = 0; k < cnode->NumDof(); ++k)
      {
        int currdof = cnode->Dofs()[k];
        int lid = (dbcmaps->CondMap())->LID(currdof);

        // store dbc status if found
        if (lid >= 0 && cnode->DbcDofs()[k] == false) cnode->SetDbc() = true;

        // check compatibility of contact symmetry condition and displacement dirichlet conditions
        if (lid < 0 && cnode->DbcDofs()[k] == true)
        {
          std::cout << "node " << cnode->Id() << " at: " << cnode->X()[0] << " " << cnode->X()[1]
                    << " " << cnode->X()[2] << std::endl;
          std::cout << "dbcdofs: " << cnode->DbcDofs()[0] << cnode->DbcDofs()[1]
                    << cnode->DbcDofs()[2] << std::endl;
          dserror("Inconsistency in structure Dirichlet conditions and Mortar symmetry conditions");
        }
      }
    }
  }
  // create old style dirichtoggle vector (supposed to go away)
  pgsdirichtoggle_ = LINALG::CreateVector(SlDoFRowMap(true), true);
  Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*(dbcmaps->CondMap())));
  temp->PutScalar(1.0);
  LINALG::Export(*temp, *pgsdirichtoggle_);

  PostStoreDirichletStatus(dbcmaps);

  return;
}

/*----------------------------------------------------------------------*
 |  Store D and M last coverged step <-> current step         popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::StoreDM(const std::string& state)
{
  // store Dold and Mold matrix in D and M
  if (state == "current")
  {
    dmatrix_ = dold_;
    mmatrix_ = mold_;
  }

  // store D and M matrix in Dold and Mold
  else if (state == "old")
  {
    dold_ = dmatrix_;
    mold_ = mmatrix_;
    if (friction_ && Dualquadslavetrafo()) doldmod_ = dmatrixmod_;
  }

  // unknown conversion
  else
  {
    dserror("ERROR: StoreDM: Unknown conversion requested!");
  }

  return;
}

/*----------------------------------------------------------------------*
 | Store nodal quant. to old ones (last conv. time step)  gitterle 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::StoreToOld(MORTAR::StrategyBase::QuantityType type)
{
  // loop over all interfaces
  for (int i = 0; i < (int)Interfaces().size(); ++i) Interfaces()[i]->StoreToOld(type);

  return;
}

/*----------------------------------------------------------------------*
 |  Update and output contact at end of time step             popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::Update(Teuchos::RCP<const Epetra_Vector> dis)
{
  // store Lagrange multipliers, D and M
  // (we need this for interpolation of the next generalized mid-point)
  // in the case of self contact, the size of z may have changed
  if (IsSelfContact()) zold_ = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));

  zold_->Scale(1.0, *z_);
  StoreNodalQuantities(MORTAR::StrategyBase::lmold);
  StoreDM("old");

  // store contact state to contact nodes (active or inactive)
  StoreNodalQuantities(MORTAR::StrategyBase::activeold);

  // old displacements in nodes
  // (this is NOT only needed for friction but also for calculating
  // the auxiliary positions in binarytree contact search)
  SetState(MORTAR::state_old_displacement, *dis);

  // reset active set status for next time step
  ResetActiveSet();

  // update flag for global contact status of last time step
  if (gactivenodes_->NumGlobalElements())
  {
    wasincontact_ = true;
    wasincontactlts_ = true;
  }
  else
  {
    wasincontact_ = false;
    wasincontactlts_ = false;
  }

  //----------------------------------------friction: store history values
  // in the case of frictional contact we have to store several
  // information and quantities at the end of a time step (converged
  // state) which is needed in the next time step as history
  // information / quantities.
  if (friction_)
  {
    // store contact state to friction nodes (slip or stick)
    StoreNodalQuantities(MORTAR::StrategyBase::slipold);

    // store nodal entries of D and M to old ones
    StoreToOld(MORTAR::StrategyBase::dm);

    // store nodal entries form penalty contact tractions to old ones
    StoreToOld(MORTAR::StrategyBase::pentrac);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  write restart information for contact                     popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::DoWriteRestart(
    std::map<std::string, Teuchos::RCP<Epetra_Vector>>& restart_vectors, bool forcedrestart) const
{
  // initalize
  Teuchos::RCP<Epetra_Vector> activetoggle = Teuchos::rcp(new Epetra_Vector(SlRowNodes()));
  Teuchos::RCP<Epetra_Vector> sliptoggle = Teuchos::null;

  // write toggle
  restart_vectors["activetoggle"] = activetoggle;
  if (friction_)
  {
    sliptoggle = Teuchos::rcp(new Epetra_Vector(SlRowNodes()));
    restart_vectors["sliptoggle"] = sliptoggle;
  }

  // loop over all interfaces
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < Interfaces()[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = Interfaces()[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = Interfaces()[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);
      int dof = (activetoggle->Map()).LID(gid);

      if (forcedrestart)
      {
        // set value active / inactive in toggle vector
        if (cnode->CoData().ActiveOld()) (*activetoggle)[dof] = 1;
      }
      else
      {
        // set value active / inactive in toggle vector
        if (cnode->Active()) (*activetoggle)[dof] = 1;
      }

      // set value slip / stick in the toggle vector
      if (friction_)
      {
        CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(cnode);
        if (forcedrestart)
        {
          if (frinode->FriData().SlipOld()) (*sliptoggle)[dof] = 1;
        }
        else
        {
          if (frinode->FriData().Slip()) (*sliptoggle)[dof] = 1;
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  read restart information for contact                      popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::DoReadRestart(IO::DiscretizationReader& reader,
    Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr)
{
  // check whether this is a restart with contact of a previously
  // non-contact simulation run (if yes, we have to be careful not
  // to try to read certain, in this case non-existing, vectors
  // such as the activetoggle or sliptoggle vectors, but rather
  // initialize the restart active and slip sets as being empty)
  bool restartwithcontact = DRT::INPUT::IntegralValue<int>(Params(), "RESTART_WITH_CONTACT");

  // set restart displacement state
  SetState(MORTAR::state_new_displacement, *dis);
  SetState(MORTAR::state_old_displacement, *dis);

  // evaluate interface and restart mortar quantities
  // in the case of SELF CONTACT, also re-setup master/slave maps
  InitMortar();
  InitEvalInterface(cparams_ptr);
  AssembleMortar();

  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // D         ---->   D * T^(-1)
  //----------------------------------------------------------------------
  if (Dualquadslavetrafo())
  {
    // modify dmatrix_
    Teuchos::RCP<LINALG::SparseMatrix> temp =
        LINALG::MLMultiply(*dmatrix_, false, *invtrafo_, false, false, false, true);
    dmatrix_ = temp;
  }

  // read restart information on active set and slip set (leave sets empty
  // if this is a restart with contact of a non-contact simulation run)
  Teuchos::RCP<Epetra_Vector> activetoggle = Teuchos::rcp(new Epetra_Vector(SlRowNodes(), true));
  if (!restartwithcontact) reader.ReadVector(activetoggle, "activetoggle");

  // friction
  Teuchos::RCP<Epetra_Vector> sliptoggle;
  Teuchos::RCP<Epetra_Vector> weightedwear;

  if (friction_)
  {
    sliptoggle = Teuchos::rcp(new Epetra_Vector(SlRowNodes()));
    if (!restartwithcontact) reader.ReadVector(sliptoggle, "sliptoggle");
  }

  // store restart information on active set and slip set
  // into nodes, therefore first loop over all interfaces
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < (Interfaces()[i]->SlaveRowNodes())->NumMyElements(); ++j)
    {
      int gid = (Interfaces()[i]->SlaveRowNodes())->GID(j);
      int dof = (activetoggle->Map()).LID(gid);

      if ((*activetoggle)[dof] == 1)
      {
        DRT::Node* node = Interfaces()[i]->Discret().gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %", gid);
        CoNode* cnode = dynamic_cast<CoNode*>(node);

        // set value active / inactive in cnode
        cnode->Active() = true;

        if (friction_)
        {
          // set value stick / slip in cnode
          if ((*sliptoggle)[dof] == 1)
            dynamic_cast<CONTACT::FriNode*>(cnode)->FriData().Slip() = true;
        }
      }
    }
  }

  // read restart information on Lagrange multipliers
  z_ = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
  zold_ = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
  if (!restartwithcontact)
    if (!(DRT::Problem::Instance()->StructuralDynamicParams().get<std::string>("INT_STRATEGY") ==
                "Standard" &&
            IsPenalty()))
    {
      reader.ReadVector(LagrMult(), "lagrmultold");
      reader.ReadVector(LagrMultOld(), "lagrmultold");
    }

  // Lagrange multiplier increment is always zero (no restart value to be read)
  zincr_ = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
  // store restart information on Lagrange multipliers into nodes
  StoreNodalQuantities(MORTAR::StrategyBase::lmcurrent);
  StoreNodalQuantities(MORTAR::StrategyBase::lmold);

  // only for Uzawa augmented strategy
  // TODO: this should be moved to contact_penalty_strategy
  if (stype_ == INPAR::CONTACT::solution_uzawa)
  {
    zuzawa_ = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
    if (!restartwithcontact) reader.ReadVector(LagrMultUzawa(), "lagrmultold");
    StoreNodalQuantities(MORTAR::StrategyBase::lmuzawa);
  }

  // store restart Mortar quantities
  StoreDM("old");

  if (friction_)
  {
    StoreNodalQuantities(MORTAR::StrategyBase::activeold);
    StoreToOld(MORTAR::StrategyBase::dm);
  }

  // (re)setup active global Epetra_Maps
  gactivenodes_ = Teuchos::null;
  gactivedofs_ = Teuchos::null;
  gactiven_ = Teuchos::null;
  gactivet_ = Teuchos::null;
  gslipnodes_ = Teuchos::null;
  gslipdofs_ = Teuchos::null;
  gslipt_ = Teuchos::null;

  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    Interfaces()[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_, Interfaces()[i]->ActiveNodes(), false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_, Interfaces()[i]->ActiveDofs(), false);
    gactiven_ = LINALG::MergeMap(gactiven_, Interfaces()[i]->ActiveNDofs(), false);
    gactivet_ = LINALG::MergeMap(gactivet_, Interfaces()[i]->ActiveTDofs(), false);
    if (friction_)
    {
      gslipnodes_ = LINALG::MergeMap(gslipnodes_, Interfaces()[i]->SlipNodes(), false);
      gslipdofs_ = LINALG::MergeMap(gslipdofs_, Interfaces()[i]->SlipDofs(), false);
      gslipt_ = LINALG::MergeMap(gslipt_, Interfaces()[i]->SlipTDofs(), false);
    }
  }

  // update flags for global contact status
  if (gactivenodes_->NumGlobalElements())
  {
    isincontact_ = true;
    wasincontact_ = true;
    wasincontactlts_ = true;
  }

  // evaluate relative movement (jump)
  // needed because it is not called in the predictor of the
  // lagrange multiplier strategy
  EvaluateRelMov();

  // reset unbalance factors for redistribution
  // (during restart the interface has been evaluated once)
  unbalanceEvaluationTime_.resize(0);
  unbalanceNumSlaveElements_.resize(0);

  return;
}

/*----------------------------------------------------------------------*
 |  Compute interface forces (for debugging only)             popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::InterfaceForces(bool output)
{
  // check chosen output option
  INPAR::CONTACT::EmOutputType emtype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::EmOutputType>(Params(), "EMOUTPUT");

  // get out of here if no output wanted
  if (emtype == INPAR::CONTACT::output_none) return;

  // compute discrete slave and master interface forces
  Teuchos::RCP<Epetra_Vector> fcslavetemp = Teuchos::rcp(new Epetra_Vector(dmatrix_->RowMap()));
  Teuchos::RCP<Epetra_Vector> fcmastertemp = Teuchos::rcp(new Epetra_Vector(mmatrix_->DomainMap()));

  // for self contact, slave and master sets may have changed,
  // thus we have to export z to new D and M dimensions
  if (IsSelfContact())
  {
    Teuchos::RCP<Epetra_Vector> zexp = Teuchos::rcp(new Epetra_Vector(dmatrix_->RowMap()));
    if (dmatrix_->RowMap().NumGlobalElements()) LINALG::Export(*z_, *zexp);
    dmatrix_->Multiply(true, *zexp, *fcslavetemp);
    mmatrix_->Multiply(true, *zexp, *fcmastertemp);
  }
  // if there is no self contact everything is ok
  else
  {
    dmatrix_->Multiply(true, *z_, *fcslavetemp);
    mmatrix_->Multiply(true, *z_, *fcmastertemp);
  }

  // export the interface forces to full dof layout
  Teuchos::RCP<Epetra_Vector> fcslave = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  Teuchos::RCP<Epetra_Vector> fcmaster = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fcslavetemp, *fcslave);
  LINALG::Export(*fcmastertemp, *fcmaster);

  // contact forces and moments
  std::vector<double> gfcs(3);
  std::vector<double> ggfcs(3);
  std::vector<double> gfcm(3);
  std::vector<double> ggfcm(3);
  std::vector<double> gmcs(3);
  std::vector<double> ggmcs(3);
  std::vector<double> gmcm(3);
  std::vector<double> ggmcm(3);

  std::vector<double> gmcsnew(3);
  std::vector<double> ggmcsnew(3);
  std::vector<double> gmcmnew(3);
  std::vector<double> ggmcmnew(3);

  // weighted gap vector
  Teuchos::RCP<Epetra_Vector> gapslave = Teuchos::rcp(new Epetra_Vector(dmatrix_->RowMap()));
  Teuchos::RCP<Epetra_Vector> gapmaster = Teuchos::rcp(new Epetra_Vector(mmatrix_->DomainMap()));

  // loop over all interfaces
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < Interfaces()[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = Interfaces()[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = Interfaces()[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      std::vector<double> nodeforce(3);
      std::vector<double> position(3);

      // forces and positions
      for (int d = 0; d < Dim(); ++d)
      {
        int dofid = (fcslavetemp->Map()).LID(cnode->Dofs()[d]);
        if (dofid < 0) dserror("ERROR: ContactForces: Did not find slave dof in map");
        nodeforce[d] = (*fcslavetemp)[dofid];
        gfcs[d] += nodeforce[d];
        position[d] = cnode->xspatial()[d];
      }

      // moments
      std::vector<double> nodemoment(3);
      nodemoment[0] = position[1] * nodeforce[2] - position[2] * nodeforce[1];
      nodemoment[1] = position[2] * nodeforce[0] - position[0] * nodeforce[2];
      nodemoment[2] = position[0] * nodeforce[1] - position[1] * nodeforce[0];
      for (int d = 0; d < 3; ++d) gmcs[d] += nodemoment[d];

      // weighted gap
      Epetra_SerialDenseVector posnode(Dim());
      std::vector<int> lm(Dim());
      std::vector<int> lmowner(Dim());
      for (int d = 0; d < Dim(); ++d)
      {
        posnode[d] = cnode->xspatial()[d];
        lm[d] = cnode->Dofs()[d];
        lmowner[d] = cnode->Owner();
      }
      LINALG::Assemble(*gapslave, posnode, lm, lmowner);
    }

    // loop over all master nodes on the current interface
    for (int j = 0; j < Interfaces()[i]->MasterRowNodes()->NumMyElements(); ++j)
    {
      int gid = Interfaces()[i]->MasterRowNodes()->GID(j);
      DRT::Node* node = Interfaces()[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      std::vector<double> nodeforce(3);
      std::vector<double> position(3);

      // forces and positions
      for (int d = 0; d < Dim(); ++d)
      {
        int dofid = (fcmastertemp->Map()).LID(cnode->Dofs()[d]);
        if (dofid < 0) dserror("ERROR: ContactForces: Did not find master dof in map");
        nodeforce[d] = -(*fcmastertemp)[dofid];
        gfcm[d] += nodeforce[d];
        position[d] = cnode->xspatial()[d];
      }

      // moments
      std::vector<double> nodemoment(3);
      nodemoment[0] = position[1] * nodeforce[2] - position[2] * nodeforce[1];
      nodemoment[1] = position[2] * nodeforce[0] - position[0] * nodeforce[2];
      nodemoment[2] = position[0] * nodeforce[1] - position[1] * nodeforce[0];
      for (int d = 0; d < 3; ++d) gmcm[d] += nodemoment[d];

      // weighted gap
      Epetra_SerialDenseVector posnode(Dim());
      std::vector<int> lm(Dim());
      std::vector<int> lmowner(Dim());
      for (int d = 0; d < Dim(); ++d)
      {
        posnode[d] = cnode->xspatial()[d];
        lm[d] = cnode->Dofs()[d];
        lmowner[d] = cnode->Owner();
      }
      LINALG::Assemble(*gapmaster, posnode, lm, lmowner);
    }
  }

  // weighted gap
  Teuchos::RCP<Epetra_Vector> gapslavefinal = Teuchos::rcp(new Epetra_Vector(dmatrix_->RowMap()));
  Teuchos::RCP<Epetra_Vector> gapmasterfinal = Teuchos::rcp(new Epetra_Vector(mmatrix_->RowMap()));
  dmatrix_->Multiply(false, *gapslave, *gapslavefinal);
  mmatrix_->Multiply(false, *gapmaster, *gapmasterfinal);
  Teuchos::RCP<Epetra_Vector> gapfinal = Teuchos::rcp(new Epetra_Vector(dmatrix_->RowMap()));
  gapfinal->Update(1.0, *gapslavefinal, 0.0);
  gapfinal->Update(-1.0, *gapmasterfinal, 1.0);

  // again, for alternative moment: lambda x gap
  // loop over all interfaces
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < Interfaces()[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = Interfaces()[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = Interfaces()[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      std::vector<double> lm(3);
      std::vector<double> nodegaps(3);
      std::vector<double> nodegapm(3);

      // LMs and gaps
      for (int d = 0; d < Dim(); ++d)
      {
        int dofid = (fcslavetemp->Map()).LID(cnode->Dofs()[d]);
        if (dofid < 0) dserror("ERROR: ContactForces: Did not find slave dof in map");
        nodegaps[d] = (*gapslavefinal)[dofid];
        nodegapm[d] = (*gapmasterfinal)[dofid];
        lm[d] = cnode->MoData().lm()[d];
      }

      // moments
      std::vector<double> nodemoments(3);
      std::vector<double> nodemomentm(3);
      nodemoments[0] = nodegaps[1] * lm[2] - nodegaps[2] * lm[1];
      nodemoments[1] = nodegaps[2] * lm[0] - nodegaps[0] * lm[2];
      nodemoments[2] = nodegaps[0] * lm[1] - nodegaps[1] * lm[0];
      nodemomentm[0] = nodegapm[1] * lm[2] - nodegapm[2] * lm[1];
      nodemomentm[1] = nodegapm[2] * lm[0] - nodegapm[0] * lm[2];
      nodemomentm[2] = nodegapm[0] * lm[1] - nodegapm[1] * lm[0];
      for (int d = 0; d < 3; ++d)
      {
        gmcsnew[d] += nodemoments[d];
        gmcmnew[d] -= nodemomentm[d];
      }
    }
  }

  // summing up over all processors
  for (int i = 0; i < 3; ++i)
  {
    Comm().SumAll(&gfcs[i], &ggfcs[i], 1);
    Comm().SumAll(&gfcm[i], &ggfcm[i], 1);
    Comm().SumAll(&gmcs[i], &ggmcs[i], 1);
    Comm().SumAll(&gmcm[i], &ggmcm[i], 1);
    Comm().SumAll(&gmcsnew[i], &ggmcsnew[i], 1);
    Comm().SumAll(&gmcmnew[i], &ggmcmnew[i], 1);
  }

  // print interface results to file
  if (emtype == INPAR::CONTACT::output_file || emtype == INPAR::CONTACT::output_both)
  {
    // do this at end of time step only (output==true)!
    // processor 0 does all the work
    if (output && Comm().MyPID() == 0)
    {
      FILE* MyFile = NULL;
      std::ostringstream filename;
      const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
      filename << filebase << ".interface";
      MyFile = fopen(filename.str().c_str(), "at+");

      if (MyFile)
      {
        for (int i = 0; i < 3; i++) fprintf(MyFile, "%g\t", ggfcs[i]);
        for (int i = 0; i < 3; i++) fprintf(MyFile, "%g\t", ggfcm[i]);
        for (int i = 0; i < 3; i++) fprintf(MyFile, "%g\t", ggmcs[i]);
        for (int i = 0; i < 3; i++) fprintf(MyFile, "%g\t", ggmcm[i]);
        // for (int i=0; i<3; i++) fprintf(MyFile, "%g\t", gsfgh[i]);
        // for (int i=0; i<3; i++) fprintf(MyFile, "%g\t", gsmgh[i]);
        fprintf(MyFile, "\n");
        fclose(MyFile);
      }
      else
        dserror("ERROR: File for writing meshtying forces could not be opened.");
    }
  }

  // print interface results to screen
  if (emtype == INPAR::CONTACT::output_screen || emtype == INPAR::CONTACT::output_both)
  {
    // do this during Newton steps only (output==false)!
    // processor 0 does all the work
    if (!output && Comm().MyPID() == 0)
    {
      double snorm = sqrt(ggfcs[0] * ggfcs[0] + ggfcs[1] * ggfcs[1] + ggfcs[2] * ggfcs[2]);
      double mnorm = sqrt(ggfcm[0] * ggfcm[0] + ggfcm[1] * ggfcm[1] + ggfcm[2] * ggfcm[2]);
      printf("Slave Contact Force:   % e  % e  % e \tNorm: % e\n", ggfcs[0], ggfcs[1], ggfcs[2],
          snorm);
      printf("Master Contact Force:  % e  % e  % e \tNorm: % e\n", ggfcm[0], ggfcm[1], ggfcm[2],
          mnorm);
      printf("Slave Contact Moment:  % e  % e  % e\n", ggmcs[0], ggmcs[1], ggmcs[2]);
      // printf("Slave Contact Moment:  % e  % e  % e\n",ggmcsnew[0],ggmcsnew[1],ggmcsnew[2]);
      printf("Master Contact Moment: % e  % e  % e\n", ggmcm[0], ggmcm[1], ggmcm[2]);
      // printf("Master Contact Moment: % e  % e  % e\n",ggmcmnew[0],ggmcmnew[1],ggmcmnew[2]);
      fflush(stdout);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  print interfaces (public)                                mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::Print(std::ostream& os) const
{
  if (Comm().MyPID() == 0)
  {
    os << "--------------------------------- CONTACT::CoAbstractStrategy\n"
       << "Contact interfaces: " << (int)Interfaces().size() << std::endl
       << "-------------------------------------------------------------\n";
  }
  Comm().Barrier();
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    std::cout << *(Interfaces()[i]);
  }
  Comm().Barrier();

  return;
}

/*----------------------------------------------------------------------*
 | print active set information                               popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::PrintActiveSet() const
{
  // output message
  Comm().Barrier();
  if (Comm().MyPID() == 0)
  {
    printf("\nActive contact set--------------------------------------------------------------\n");
    fflush(stdout);
  }

  //**********************************************************************
  // detailed active set output
  //**********************************************************************

#ifdef CONTACTASOUTPUT
  // create storage for local and global data
  std::vector<int> lnid, gnid;
  std::vector<double> llmn, glmn;
  std::vector<double> lgap, ggap;

  std::vector<double> Xposl, Xposg;
  std::vector<double> Yposl, Yposg;
  std::vector<double> Zposl, Zposg;

  std::vector<double> xposl, xposg;
  std::vector<double> yposl, yposg;
  std::vector<double> zposl, zposg;

  // introduce integer variable status
  // (0=inactive, 1=active, 2=slip, 3=stick)
  // (this is necessary as all data will be written by proc 0, but
  // the knowledge of the above status ONLY exists on the owner
  // processor of the respective node. Thus this information must
  // also be communicated to proc 0 in addition to the actual data!)
  std::vector<int> lsta, gsta;

  // some more storage for local and global friction data
  std::vector<double> llmt, glmt;
  std::vector<double> ljtx, gjtx;
  std::vector<double> ljte, gjte;
  std::vector<double> lwear, gwear;

  // loop over all interfaces
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    // if (i>0) dserror("ERROR: PrintActiveSet: Double active node check needed for n interfaces!");

    // loop over all slave row nodes on the current interface
    for (int j = 0; j < Interfaces()[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      // gid of current node
      int gid = Interfaces()[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = Interfaces()[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);

      //--------------------------------------------------------------------
      // FRICTIONLESS CASE
      //--------------------------------------------------------------------
      if (!friction_)
      {
        // cast to CoNode
        CoNode* cnode = dynamic_cast<CoNode*>(node);

        // compute weighted gap
        double wgap = (*g_)[g_->Map().LID(gid)];

        double Xpos = cnode->X()[0];
        double Ypos = cnode->X()[1];
        double Zpos = cnode->X()[2];

        double xpos = cnode->xspatial()[0];
        double ypos = cnode->xspatial()[1];
        double zpos = cnode->xspatial()[2];

        // compute normal part of Lagrange multiplier
        double nz = 0.0;
        for (int k = 0; k < 3; ++k) nz += cnode->MoData().n()[k] * cnode->MoData().lm()[k];

        // store node id
        lnid.push_back(gid);

        // store relevant data
        llmn.push_back(nz);
        lgap.push_back(wgap);
        Xposl.push_back(Xpos);
        Yposl.push_back(Ypos);
        Zposl.push_back(Zpos);
        xposl.push_back(xpos);
        yposl.push_back(ypos);
        zposl.push_back(zpos);

        // store status (0=inactive, 1=active, 2=slip, 3=stick)
        if (cnode->Active())
          lsta.push_back(1);
        else
          lsta.push_back(0);
      }

      //--------------------------------------------------------------------
      // FRICTIONAL CASE
      //--------------------------------------------------------------------
      else
      {
        // cast to CoNode and FriNode
        CoNode* cnode = dynamic_cast<CoNode*>(node);
        FriNode* frinode = dynamic_cast<FriNode*>(cnode);

        // compute weighted gap
        double wgap = (*g_)[g_->Map().LID(gid)];

        // compute normal part of Lagrange multiplier
        double nz = 0.0;
        for (int k = 0; k < 3; ++k) nz += frinode->MoData().n()[k] * frinode->MoData().lm()[k];

        // compute tangential parts of Lagrange multiplier and jumps and wear
        double txiz = 0.0;
        double tetaz = 0.0;
        double jumptxi = 0.0;
        double jumpteta = 0.0;
        double wear = 0.0;

        for (int k = 0; k < Dim(); ++k)
        {
          txiz += frinode->CoData().txi()[k] * frinode->MoData().lm()[k];
          tetaz += frinode->CoData().teta()[k] * frinode->MoData().lm()[k];
          jumptxi += frinode->CoData().txi()[k] * frinode->FriData().jump()[k];
          jumpteta += frinode->CoData().teta()[k] * frinode->FriData().jump()[k];
        }

        // total tangential component
        double tz = sqrt(txiz * txiz + tetaz * tetaz);

        // check for dimensions
        if (Dim() == 2 && abs(jumpteta) > 0.0001) dserror("Error: Jumpteta should be zero for 2D");

        // store node id
        lnid.push_back(gid);

        // store relevant data
        llmn.push_back(nz);
        lgap.push_back(wgap);
        llmt.push_back(tz);
        ljtx.push_back(jumptxi);
        ljte.push_back(jumpteta);
        lwear.push_back(wear);

        // store status (0=inactive, 1=active, 2=slip, 3=stick)
        if (cnode->Active())
        {
          if (frinode->FriData().Slip())
            lsta.push_back(2);
          else
            lsta.push_back(3);
        }
        else
        {
          lsta.push_back(0);
        }
      }
    }
  }

  // we want to gather data from on all procs
  std::vector<int> allproc(Comm().NumProc());
  for (int i = 0; i < Comm().NumProc(); ++i) allproc[i] = i;

  // communicate all data to proc 0
  LINALG::Gather<int>(lnid, gnid, (int)allproc.size(), &allproc[0], Comm());
  LINALG::Gather<double>(llmn, glmn, (int)allproc.size(), &allproc[0], Comm());
  LINALG::Gather<double>(lgap, ggap, (int)allproc.size(), &allproc[0], Comm());
  LINALG::Gather<int>(lsta, gsta, (int)allproc.size(), &allproc[0], Comm());

  LINALG::Gather<double>(Xposl, Xposg, (int)allproc.size(), &allproc[0], Comm());
  LINALG::Gather<double>(Yposl, Yposg, (int)allproc.size(), &allproc[0], Comm());
  LINALG::Gather<double>(Zposl, Zposg, (int)allproc.size(), &allproc[0], Comm());

  LINALG::Gather<double>(xposl, xposg, (int)allproc.size(), &allproc[0], Comm());
  LINALG::Gather<double>(yposl, yposg, (int)allproc.size(), &allproc[0], Comm());
  LINALG::Gather<double>(zposl, zposg, (int)allproc.size(), &allproc[0], Comm());

  // communicate some more data to proc 0 for friction
  if (friction_)
  {
    LINALG::Gather<double>(llmt, glmt, (int)allproc.size(), &allproc[0], Comm());
    LINALG::Gather<double>(ljtx, gjtx, (int)allproc.size(), &allproc[0], Comm());
    LINALG::Gather<double>(ljte, gjte, (int)allproc.size(), &allproc[0], Comm());
    LINALG::Gather<double>(lwear, gwear, (int)allproc.size(), &allproc[0], Comm());
  }

  // output is solely done by proc 0
  if (Comm().MyPID() == 0)
  {
    //--------------------------------------------------------------------
    // FRICTIONLESS CASE
    //--------------------------------------------------------------------
    if (!friction_)
    {
      // loop over all nodes
      for (int k = 0; k < (int)gnid.size(); ++k)
      {
        // print nodes of inactive set *************************************
        if (gsta[k] == 0)
        {
          printf("INACTIVE: %d \t wgap: % e \t lm: % e \t Xref: % e \t Yref: % e \t Zref: % e \n",
              gnid[k], ggap[k], glmn[k], Xposg[k], Yposg[k], Zposg[k]);
          fflush(stdout);
        }

        // print nodes of active set ***************************************
        else if (gsta[k] == 1)
        {
          printf("ACTIVE:   %d \t wgap: % e \t lm: % e \t Xref: % e \t Yref: % e \t Zref: % e \n",
              gnid[k], ggap[k], glmn[k], Xposg[k], Yposg[k], Zposg[k]);
          fflush(stdout);
        }

        // invalid status **************************************************
        else
          dserror("ERROR: Invalid node status %i for frictionless case", gsta[k]);
      }
    }

    //--------------------------------------------------------------------
    // FRICTIONAL CASE
    //--------------------------------------------------------------------
    else
    {
#ifdef CONTACTEXPORT
      // export variables
      double sum_jumpx = 0.0;
      double sum_jumpe = 0.0;
      double sum_jumpall = 0.0;
      double iter_slip = 0.0;
#endif

      // loop over all nodes
      for (int k = 0; k < (int)gnid.size(); ++k)
      {
        // print nodes of slip set **************************************
        if (gsta[k] == 2)
        {
          printf("SLIP:  %d \t lm_n: % e \t lm_t: % e \t jump1: % e \t jump2: % e \t wear: % e \n",
              gnid[k], glmn[k], glmt[k], gjtx[k], gjte[k], gwear[k]);
          fflush(stdout);
#ifdef CONTACTEXPORT
          // preparation for output
          sum_jumpx += gjtx[k];
          sum_jumpe += gjte[k];
          sum_jumpall += sqrt(gjtx[k] * gjtx[k] + gjte[k] * gjte[k]);
          iter_slip = iter_slip + 1.0;
#endif
        }

        // print nodes of stick set *************************************
        else if (gsta[k] == 3)
        {
          printf("STICK: %d \t lm_n: % e \t lm_t: % e \t jump1: % e \t jump2: % e \t wear: % e \n",
              gnid[k], glmn[k], glmt[k], gjtx[k], gjte[k], gwear[k]);
          fflush(stdout);
        }

        // print nodes of inactive set *************************************
        else if (gsta[k] == 0)
        {
          // do nothing
        }

        // invalid status **************************************************
        else
          dserror("ERROR: Invalid node status %i for frictional case", gsta[k]);
      }

#ifdef CONTACTEXPORT
      // export averaged slip increments to xxx.jump
      double sum_jumpx_final = 0.0;
      double sum_jumpe_final = 0.0;
      double sum_jumpall_final = 0.0;

      if (iter_slip > 0.0)
      {
        sum_jumpx_final = sum_jumpx / iter_slip;
        sum_jumpe_final = sum_jumpe / iter_slip;
        sum_jumpall_final = sum_jumpall / iter_slip;
      }

      FILE* MyFile = NULL;
      std::ostringstream filename;
      const std::string filebase =
          DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix();
      filename << filebase << ".jump";
      MyFile = fopen(filename.str().c_str(), "at+");
      if (MyFile)
      {
        // fprintf(MyFile,valuename.c_str());
        fprintf(MyFile, "%g\t", sum_jumpx_final);
        fprintf(MyFile, "%g\t", sum_jumpe_final);
        fprintf(MyFile, "%g\n", sum_jumpall_final);
        fclose(MyFile);
      }
      else
        dserror("ERROR: File for Output could not be opened.");
#endif
    }
  }

#else
  //**********************************************************************
  // reduced active set output
  //**********************************************************************

  // counters
  int activenodes = 0;
  int gactivenodes = 0;
  int inactivenodes = 0;
  int ginactivenodes = 0;
  int slipnodes = 0;
  int gslipnodes = 0;

  // counters for non-smooth contact
  int surfacenodes = 0;
  int gsurfacenodes = 0;
  int edgenodes = 0;
  int gedgenodes = 0;
  int cornernodes = 0;
  int gcornernodes = 0;

  // nonsmooth contact active?
  bool nonsmooth = DRT::INPUT::IntegralValue<int>(Params(), "NONSMOOTH_GEOMETRIES");

  // loop over all interfaces
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < Interfaces()[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = Interfaces()[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = Interfaces()[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);

      // increase active counters
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      if (cnode->Active())
        activenodes += 1;
      else
        inactivenodes += 1;

      // increase friction counters
      if (friction_)
      {
        FriNode* frinode = dynamic_cast<FriNode*>(cnode);
        if (cnode->Active() && frinode->FriData().Slip()) slipnodes += 1;
      }

      // get nonsmooth contact states
      if (nonsmooth)
      {
        if (cnode->Active() && cnode->IsOnEdge() && !cnode->IsOnCorner()) edgenodes += 1;
        if (cnode->Active() && cnode->IsOnCorner()) cornernodes += 1;
        if (cnode->Active() && !cnode->IsOnCornerEdge()) surfacenodes += 1;
      }
    }
  }

  // sum among all processors
  Comm().SumAll(&activenodes, &gactivenodes, 1);
  Comm().SumAll(&inactivenodes, &ginactivenodes, 1);
  Comm().SumAll(&slipnodes, &gslipnodes, 1);
  Comm().SumAll(&edgenodes, &gedgenodes, 1);
  Comm().SumAll(&cornernodes, &gcornernodes, 1);
  Comm().SumAll(&surfacenodes, &gsurfacenodes, 1);

  // print active set information
  if (Comm().MyPID() == 0)
  {
    if (nonsmooth)
    {
      std::cout << BLUE2_LIGHT << "Total ACTIVE SURFACE nodes:\t" << gsurfacenodes << END_COLOR
                << std::endl;
      std::cout << BLUE2_LIGHT << "Total    ACTIVE EDGE nodes:\t" << gedgenodes << END_COLOR
                << std::endl;
      std::cout << BLUE2_LIGHT << "Total  ACTIVE CORNER nodes:\t" << gcornernodes << END_COLOR
                << std::endl;
      std::cout << RED_LIGHT << "Total       INACTIVE nodes:\t" << ginactivenodes << END_COLOR
                << std::endl;
    }
    else if (friction_)
    {
      std::cout << BLUE2_LIGHT << "Total     SLIP nodes:\t" << gslipnodes << END_COLOR << std::endl;
      std::cout << BLUE2_LIGHT << "Total    STICK nodes:\t" << gactivenodes - gslipnodes
                << END_COLOR << std::endl;
      std::cout << RED_LIGHT << "Total INACTIVE nodes:\t" << ginactivenodes << END_COLOR
                << std::endl;
    }
    else
    {
      std::cout << BLUE2_LIGHT << "Total   ACTIVE nodes:\t" << gactivenodes << END_COLOR
                << std::endl;
      std::cout << RED_LIGHT << "Total INACTIVE nodes:\t" << ginactivenodes << END_COLOR
                << std::endl;
    }
  }
#endif  // #ifdef CONTACTASOUTPUT
  // output line
  Comm().Barrier();
  if (Comm().MyPID() == 0)
  {
    printf("--------------------------------------------------------------------------------\n\n");
    fflush(stdout);
  }

  return;
}

/*----------------------------------------------------------------------*
 | Visualization of contact segments with gmsh                popp 08/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::VisualizeGmsh(const int step, const int iter)
{
  // visualization with gmsh
  for (int i = 0; i < (int)Interfaces().size(); ++i) Interfaces()[i]->VisualizeGmsh(step, iter);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::CollectMapsForPreconditioner(
    Teuchos::RCP<Epetra_Map>& MasterDofMap, Teuchos::RCP<Epetra_Map>& SlaveDofMap,
    Teuchos::RCP<Epetra_Map>& InnerDofMap, Teuchos::RCP<Epetra_Map>& ActiveDofMap) const
{
  InnerDofMap = gndofrowmap_;   // global internal dof row map
  ActiveDofMap = gactivedofs_;  // global active slave dof row map

  // check if parallel redistribution is used
  // if parallel redistribution is activated, then use (original) maps before redistribution
  // otherwise we use just the standard master/slave maps
  if (pgsdofrowmap_ != Teuchos::null)
    SlaveDofMap = pgsdofrowmap_;
  else
    SlaveDofMap = gsdofrowmap_;
  if (pgmdofrowmap_ != Teuchos::null)
    MasterDofMap = pgmdofrowmap_;
  else
    MasterDofMap = gmdofrowmap_;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::Reset(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& dispnp, const Epetra_Vector& xnew)
{
  SetState(MORTAR::state_new_displacement, dispnp);
  ResetLagrangeMultipliers(cparams, xnew);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::Evaluate(CONTACT::ParamsInterface& cparams,
    const std::vector<Teuchos::RCP<const Epetra_Vector>>* eval_vec,
    const std::vector<Teuchos::RCP<Epetra_Vector>>* eval_vec_mutable)
{
  PreEvaluate(cparams);

  const enum MORTAR::ActionType& act = cparams.GetActionType();
  switch (act)
  {
    // -------------------------------------------------------------------
    // evaluate only the weighted gap
    // -------------------------------------------------------------------
    case MORTAR::eval_weighted_gap:
    {
      EvalWeightedGap(cparams);
      break;
    }
    // -------------------------------------------------------------------
    // evaluate only the contact forces / contact right hand side
    // -------------------------------------------------------------------
    case MORTAR::eval_force:
    {
      EvalForce(cparams);
      break;
    }
    // -------------------------------------------------------------------
    // evaluate the contact matrix blocks and the rhs contributions
    // -------------------------------------------------------------------
    case MORTAR::eval_force_stiff:
    {
      EvalForceStiff(cparams);
      break;
    }
    // -------------------------------------------------------------------
    // run before an evaluate call in the STR::ModelEvaluator class
    // -------------------------------------------------------------------
    case MORTAR::eval_run_pre_evaluate:
    {
      RunPreEvaluate(cparams);
      break;
    }
    // -------------------------------------------------------------------
    // run after an evaluate call in the STR::ModelEvaluator class
    // -------------------------------------------------------------------
    case MORTAR::eval_run_post_evaluate:
    {
      RunPostEvaluate(cparams);
      break;
    }
    // -------------------------------------------------------------------
    // reset internal stored solution quantities (e.g.
    // displacement state, Lagrange multi.)
    // -------------------------------------------------------------------
    case MORTAR::eval_reset:
    {
      if (not eval_vec) dserror("Missing evaluation vectors!");
      if (eval_vec->size() != 2)
        dserror(
            "The \"MORTAR::eval_reset\" action expects \n"
            "exactly 2 evaluation vector pointers! But you \n"
            "passed %i vector pointers!",
            eval_vec->size());
      const Epetra_Vector& dispnp = *((*eval_vec)[0]);
      const Epetra_Vector& xnew = *((*eval_vec)[1]);
      Reset(cparams, dispnp, xnew);

      break;
    }
    // -------------------------------------------------------------------
    // recover internal stored solution quantities (e.g. Lagrange multi.)
    // -------------------------------------------------------------------
    case MORTAR::eval_run_post_compute_x:
    {
      if (not eval_vec) dserror("Missing evaluation vectors!");

      if (eval_vec->size() != 3)
        dserror(
            "The \"MORTAR::eval_recover\" action expects \n"
            "exactly 3 evaluation vector pointers! But you \n"
            "passed %i vector pointers!",
            eval_vec->size());

      const Teuchos::RCP<const Epetra_Vector>& xold_ptr = (*eval_vec)[0];
      if (xold_ptr.is_null() or !xold_ptr.is_valid_ptr()) dserror("xold_ptr is undefined!");

      const Teuchos::RCP<const Epetra_Vector>& dir_ptr = (*eval_vec)[1];
      if (dir_ptr.is_null() or !dir_ptr.is_valid_ptr()) dserror("dir_ptr is undefined!");

      const Teuchos::RCP<const Epetra_Vector>& xnew_ptr = (*eval_vec)[2];
      if (xnew_ptr.is_null() or !xnew_ptr.is_valid_ptr()) dserror("xnew_ptr is undefined!");

      RunPostComputeX(cparams, *xold_ptr, *dir_ptr, *xnew_ptr);

      break;
    }
    case MORTAR::eval_run_pre_compute_x:
    {
      if (not eval_vec) dserror("Missing constant evaluation vectors!");
      if (not eval_vec_mutable) dserror("Missing mutable evaluation vectors!");

      if (eval_vec->size() != 1)
        dserror(
            "The \"MORTAR::eval_augment_direction\" action expects \n"
            "exactly 1 constant evaluation vector pointer! But you \n"
            "passed %i vector pointers!",
            eval_vec->size());

      if (eval_vec_mutable->size() != 1)
        dserror(
            "The \"MORTAR::eval_augment_direction\" action expects \n"
            "exactly 1 mutable evaluation vector pointer! But you \n"
            "passed %i vector pointers!",
            eval_vec->size());

      const Teuchos::RCP<const Epetra_Vector>& xold_ptr = eval_vec->front();
      if (xold_ptr.is_null()) dserror("Missing xold vector!");

      const Teuchos::RCP<Epetra_Vector>& dir_mutable_ptr = eval_vec_mutable->front();
      if (dir_mutable_ptr.is_null()) dserror("Missing dir_mutable vector!");

      RunPreComputeX(cparams, *xold_ptr, *dir_mutable_ptr);

      break;
    }
    case MORTAR::eval_run_post_iterate:
    {
      RunPostIterate(cparams);

      break;
    }
    case MORTAR::eval_run_post_apply_jacobian_inverse:
    {
      const Epetra_Vector* rhs = cparams.Get<const Epetra_Vector>(0);
      Epetra_Vector* result = cparams.Get<Epetra_Vector>(1);
      const Epetra_Vector* xold = cparams.Get<const Epetra_Vector>(2);
      const NOX::NLN::Group* grp = cparams.Get<const NOX::NLN::Group>(3);

      RunPostApplyJacobianInverse(cparams, *rhs, *result, *xold, *grp);

      break;
    }
    case MORTAR::eval_correct_parameters:
    {
      const NOX::NLN::CorrectionType* type = cparams.Get<const NOX::NLN::CorrectionType>(0);

      CorrectParameters(cparams, *type);

      break;
    }
    case MORTAR::eval_wgap_gradient_error:
    {
      EvalWeightedGapGradientError(cparams);

      break;
    }
    case MORTAR::eval_static_constraint_rhs:
    {
      EvalStaticConstraintRHS(cparams);

      break;
    }
    case MORTAR::remove_condensed_contributions_from_str_rhs:
    {
      if (not eval_vec_mutable) dserror("The mutable evaluation vector is expected!");
      if (eval_vec_mutable->size() < 1)
        dserror(
            "The eval_vec_mutable is supposed to have at least a length"
            " of 1!");

      Epetra_Vector& str_rhs = *eval_vec_mutable->front();
      RemoveCondensedContributionsFromRhs(str_rhs);

      break;
    }
    case MORTAR::eval_run_pre_solve:
    {
      if (not eval_vec) dserror("The read-only evaluation vector is expected!");
      if (eval_vec->size() < 1)
        dserror(
            "The eval_vec is supposed to have at least a length"
            " of 1!");

      const Teuchos::RCP<const Epetra_Vector>& curr_disp = eval_vec->front();
      RunPreSolve(curr_disp, cparams);

      break;
    }
    // -------------------------------------------------------------------
    // no suitable action could be found
    // -------------------------------------------------------------------
    default:
    {
      dserror("Unsupported action type: %i | %s", act, ActionType2String(act).c_str());
      break;
    }
  }

  PostEvaluate(cparams);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::EvalWeightedGap(CONTACT::ParamsInterface& cparams)
{
  dserror(
      "Not yet implemented! See the XCONTACT::Strategy for an "
      "example.");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::EvalForce(CONTACT::ParamsInterface& cparams)
{
  dserror(
      "Not yet implemented! See the CONTACT::AUG::Strategy for an "
      "example.");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::EvalForceStiff(CONTACT::ParamsInterface& cparams)
{
  dserror(
      "Not yet implemented! See the CONTACT::AUG::Strategy for an "
      "example.");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::EvalStaticConstraintRHS(CONTACT::ParamsInterface& cparams)
{
  dserror(
      "Not yet implemented! See the CONTACT::AUG::Strategy for an "
      "example.");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::RemoveCondensedContributionsFromRhs(Epetra_Vector& str_rhs) const
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::RunPreEvaluate(CONTACT::ParamsInterface& cparams)
{
  /* Not yet implemented! See the CONTACT::AUG framework for an
   * example. */
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::RunPostEvaluate(CONTACT::ParamsInterface& cparams)
{
  /* Not yet implemented! See the CONTACT::AUG::ComboStrategy for an
   * example. */
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::RunPostComputeX(const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  dserror(
      "Not yet implemented! See the CONTACT::AUG::Strategy for an "
      "example.");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::RunPreComputeX(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold, Epetra_Vector& dir_mutable)
{
  // do nothing
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::RunPostIterate(const CONTACT::ParamsInterface& cparams)
{
  // do nothing
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::RunPreSolve(
    const Teuchos::RCP<const Epetra_Vector>& curr_disp, const CONTACT::ParamsInterface& cparams)
{
  // do nothing
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::RunPostApplyJacobianInverse(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& rhs, Epetra_Vector& result,
    const Epetra_Vector& xold, const NOX::NLN::Group& grp)
{
  // do nothing
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::EvalWeightedGapGradientError(CONTACT::ParamsInterface& cparams)
{
  dserror(
      "Not yet implemented! See the CONTACT::AUG::Strategy for an "
      "example.");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::ResetLagrangeMultipliers(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xnew)
{
  dserror(
      "Not yet implemented! See the CONTACT::AUG::Strategy for an "
      "example.");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::CorrectParameters(
    CONTACT::ParamsInterface& cparams, const NOX::NLN::CorrectionType type)
{
  /* do nothing */
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::CoAbstractStrategy::IsSaddlePointSystem() const
{
  if ((stype_ == INPAR::CONTACT::solution_lagmult) and
      SystemType() == INPAR::CONTACT::system_saddlepoint)
  {
    if (IsInContact() or WasInContact() or WasInContactLastTimeStep()) return true;
  }
  return false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::CoAbstractStrategy::IsCondensedSystem() const
{
  if (stype_ == INPAR::CONTACT::solution_lagmult and
      SystemType() != INPAR::CONTACT::system_saddlepoint)
  {
    if (IsInContact() or WasInContact() or WasInContactLastTimeStep()) return true;
  }
  return false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::FillMapsForPreconditioner(
    std::vector<Teuchos::RCP<Epetra_Map>>& maps) const
{
  /* FixMe This function replaces the deprecated CollectMapsForPreconditioner(),
   * the old version can be deleted, as soon as the contact uses the new
   * structure framework. */

  if (maps.size() != 4) dserror("The vector size has to be 4!");
  /* check if parallel redistribution is used
   * if parallel redistribution is activated, then use (original) maps
   * before redistribution otherwise we use just the standard master/slave
   * maps */

  // (0) masterDofMap
  if (pgmdofrowmap_ != Teuchos::null)
    maps[0] = pgmdofrowmap_;
  else
    maps[0] = gmdofrowmap_;

  // (1) slaveDofMap
  if (pgsdofrowmap_ != Teuchos::null)
    maps[1] = pgsdofrowmap_;
  else
    maps[1] = gsdofrowmap_;

  // (2) innerDofMap
  maps[2] = gndofrowmap_;

  // (3) activeDofMap
  maps[3] = gactivedofs_;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::CoAbstractStrategy::computePreconditioner(
    const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams)
{
  dserror("Not implemented!");
  return false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::CoAbstractStrategy::GetLagrMultNp(
    const bool& redist) const
{
  if ((redist) or not ParRedist()) return z_;

  Teuchos::RCP<Epetra_Vector> z_unredist = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(false)));
  LINALG::Export(*z_, *z_unredist);
  return z_unredist;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::CoAbstractStrategy::GetLagrMultN(
    const bool& redist) const
{
  if ((redist) or not ParRedist()) return zold_;

  Teuchos::RCP<Epetra_Vector> zold_unredist = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(false)));
  LINALG::Export(*zold_, *zold_unredist);
  return zold_unredist;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::CoAbstractStrategy::GetPotentialValue(
    const enum NOX::NLN::MeritFunction::MeritFctName mrt_type) const
{
  dserror("The currently active strategy \"%s\" does not support this method!",
      INPAR::CONTACT::SolvingStrategy2String(Type()).c_str());
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::CoAbstractStrategy::GetLinearizedPotentialValueTerms(const Epetra_Vector& dir,
    const enum NOX::NLN::MeritFunction::MeritFctName mrt_type,
    const enum NOX::NLN::MeritFunction::LinOrder linorder,
    const enum NOX::NLN::MeritFunction::LinType lintype) const
{
  dserror("The currently active strategy \"%s\" does not support this method!",
      INPAR::CONTACT::SolvingStrategy2String(Type()).c_str());
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::PostprocessQuantitiesPerInterface(
    Teuchos::RCP<Teuchos::ParameterList> outputParams)
{
  using Teuchos::RCP;

  // Evaluate slave and master forces
  {
    RCP<Epetra_Vector> fcslave = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true), true));
    RCP<Epetra_Vector> fcmaster = Teuchos::rcp(new Epetra_Vector(MaDoFRowMap(true), true));

    // Mortar matrices might not be initialized, e.g. in the initial state. If so, keep zero vector.
    if (!DMatrix().is_null()) DMatrix()->Multiply(true, *zold_, *fcslave);
    if (!MMatrix().is_null()) MMatrix()->Multiply(true, *zold_, *fcmaster);

    // Append data to parameter list
    outputParams->set<RCP<const Epetra_Vector>>("interface traction", zold_);
    outputParams->set<RCP<const Epetra_Vector>>("slave forces", fcslave);
    outputParams->set<RCP<const Epetra_Vector>>("master forces", fcmaster);
  }

  // Postprocess contact stresses
  {
    ComputeContactStresses();

    // Append data to parameter list
    outputParams->set<RCP<const Epetra_Vector>>("norcontactstress", stressnormal_);
    outputParams->set<RCP<const Epetra_Vector>>("tancontactstress", stresstangential_);
  }

  for (auto& interface : Interfaces()) interface->PostprocessQuantities(*outputParams);
}
