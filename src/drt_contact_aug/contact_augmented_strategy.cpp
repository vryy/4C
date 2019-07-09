/*----------------------------------------------------------------------*/
/*!
\brief Augmented contact solving strategy with standard Lagrangian
       multipliers.

\level 3

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/
#include "contact_augmented_strategy.H"
#include "contact_augmented_interface.H"
#include "contact_aug_potential.H"
#include "contact_aug_steepest_ascent_strategy.H"
#include "contact_aug_active_set.H"
#include "contact_aug_timemonitor.H"
#include "contact_aug_parallel_distribution_controller.H"

#include "../drt_mortar/mortar_matrix_transform.H"

#include "../drt_contact/contact_node.H"
#include "../drt_contact/contact_lagrange_strategy.H"
#include "../drt_contact/contact_defines.H"
#include "../drt_contact/contact_paramsinterface.H"
#include "../drt_contact/contact_utils.H"
#include "../drt_mortar/mortar_utils.H"

#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/epetra_utils.H"

#include "../drt_io/io_pstream.H"
#include "../drt_io/io.H"

#include "../drt_inpar/inpar_contact.H"

#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_utils.H"

#include <unordered_map>
#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::DataContainer::DataContainer()
    : wasincontactlastiter_(false),
      isactivesetconverged_(false),
      printlinearconservation_(false),
      printangularconservation_(false),
      is_semi_smooth_newton_(false),
      matrix_maps_valid_(false),
      vector_maps_valid_(false),
      cn_(-1.0),
      eval_state_(MORTAR::eval_none),
      ghosting_strategy_(INPAR::MORTAR::ghosting_redundant),
      var_type_(INPAR::CONTACT::var_unknown),
      fd_check_type_(INPAR::CONTACT::FDCheck::off),
      potentialPtr_(Teuchos::null),
      mat_row_col_transformer_(Teuchos::null),
      BMatrixPtr_(Teuchos::null),
      dGLmLinMatrixPtr_(Teuchos::null),
      dGGLinMatrixPtr_(Teuchos::null),
      dLmNWGapLinMatrixPtr_(Teuchos::null),
      dLmTLmTMatrixPtr_(Teuchos::null),
      dLmTLmTLinMatrixPtr_(Teuchos::null),
      inactiveLinMatrixPtr_(Teuchos::null),
      inactiveDiagMatrixPtr_(Teuchos::null),
      aPtr_(Teuchos::null),
      kappaPtr_(Teuchos::null),
      lmNPtr_(Teuchos::null),
      aWGapPtr_(Teuchos::null),
      wGapAllPtr_(Teuchos::null),
      dLmTLmTRhsPtr_(Teuchos::null),
      slForceLmPtr_(Teuchos::null),
      slForceGPtr_(Teuchos::null),
      maForceLmPtr_(Teuchos::null),
      maForceGPtr_(Teuchos::null),
      cnPtr_(Teuchos::null),
      gsndofrowmapPtr_(Teuchos::null),
      gstdofrowmapPtr_(Teuchos::null),
      gOldActiveSlaveNodesPtr_(Teuchos::null)
{
  // empty constructor

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::DataContainer::InitMatrixRowColTransformer()
{
  MORTAR::MatrixRowColTransformer::plain_block_map_pairs redistributed_row(4);
  redistributed_row[DRT::UTILS::block_displ_displ] = &GSlMaDofRowMapPtr();
  redistributed_row[DRT::UTILS::block_displ_lm] = &GDispDofRowMapPtr();
  redistributed_row[DRT::UTILS::block_lm_displ] = &GLmDofRowMapPtr();
  redistributed_row[DRT::UTILS::block_lm_lm] = &GLmDofRowMapPtr();

  MORTAR::MatrixRowColTransformer::plain_block_map_pairs redistributed_col(4);
  redistributed_col[DRT::UTILS::block_displ_displ] = &GSlMaDofRowMapPtr();
  redistributed_col[DRT::UTILS::block_displ_lm] = &GLmDofRowMapPtr();
  redistributed_col[DRT::UTILS::block_lm_displ] = &GDispDofRowMapPtr();
  redistributed_col[DRT::UTILS::block_lm_lm] = &GLmDofRowMapPtr();

  MORTAR::MatrixRowColTransformer::plain_block_map_pairs unredistributed_row(4);
  unredistributed_row[DRT::UTILS::block_displ_displ] = &PGSlMaDofRowMapPtr();
  unredistributed_row[DRT::UTILS::block_displ_lm] = &ProbDofsPtr();
  unredistributed_row[DRT::UTILS::block_lm_displ] = &PGLmDofRowMapPtr();
  unredistributed_row[DRT::UTILS::block_lm_lm] = &PGLmDofRowMapPtr();


  MORTAR::MatrixRowColTransformer::plain_block_map_pairs unredistributed_col(4);
  unredistributed_col[DRT::UTILS::block_displ_displ] = &PGSlMaDofRowMapPtr();
  unredistributed_col[DRT::UTILS::block_displ_lm] = &PGLmDofRowMapPtr();
  unredistributed_col[DRT::UTILS::block_lm_displ] = &ProbDofsPtr();
  unredistributed_col[DRT::UTILS::block_lm_lm] = &PGLmDofRowMapPtr();

  mat_row_col_transformer_->Init(
      redistributed_row, redistributed_col, unredistributed_row, unredistributed_col);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::DataContainer::InitSubDataContainer(
    const INPAR::CONTACT::SolvingStrategy strat_type)
{
  switch (strat_type)
  {
    case INPAR::CONTACT::solution_steepest_ascent:
    case INPAR::CONTACT::solution_steepest_ascent_sp:
      // do not initialize it twice (combo strategy)
      if (sa_data_ptr_.is_null())
        sa_data_ptr_ = Teuchos::RCP<CONTACT::AUG::STEEPESTASCENT::DataContainer>(
            new CONTACT::AUG::STEEPESTASCENT::DataContainer());
      break;
    default:
      dserror(
          "There is no known sub-data container for the given "
          "strategy type! ( strat_type = %s )",
          INPAR::CONTACT::SolvingStrategy2String(strat_type).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::Strategy::Strategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap, const Teuchos::ParameterList& params,
    const plain_interface_set& interfaces, int dim, const Teuchos::RCP<const Epetra_Comm>& comm,
    int maxdof)
    : ::CONTACT::CoAbstractStrategy(
          data_ptr, DofRowMap, NodeRowMap, params, dim, comm, 0.0, maxdof),
      augDataPtr_(Teuchos::rcp_dynamic_cast<CONTACT::AUG::DataContainer>(data_ptr, true)),
      augData_(*augDataPtr_)
{
  // store values of the parameter list
  const Teuchos::ParameterList& p_aug = params.sublist("AUGMENTED");

  Data().PrintLinearMomConservation() =
      DRT::INPUT::IntegralValue<bool>(p_aug, "PRINT_LINEAR_CONSERVATION");

  Data().PrintAngularMomConservation() =
      DRT::INPUT::IntegralValue<bool>(p_aug, "PRINT_ANGULAR_CONSERVATION");

  Data().AddInactivForceContributions() =
      DRT::INPUT::IntegralValue<bool>(p_aug, "ADD_INACTIVE_FORCE_CONTRIBUTIONS");

  Data().SetIsSemiSmoothNewton(DRT::INPUT::IntegralValue<bool>(params, "SEMI_SMOOTH_NEWTON"));

  Data().SetConstantCn(Params().get<double>("SEMI_SMOOTH_CN"));

  Data().SetGhostingStrategy(DRT::INPUT::IntegralValue<INPAR::MORTAR::GhostingStrategy>(
      Params().sublist("PARALLEL REDISTRIBUTION"), "GHOSTING_STRATEGY"));

  Data().SetVariationalApproachType(DRT::INPUT::IntegralValue<INPAR::CONTACT::VariationalApproach>(
      p_aug, "VARIATIONAL_APPROACH"));

  if (Data().AddInactivForceContributions() and
      Data().VariationalApproachType() != INPAR::CONTACT::var_complete)
    dserror(
        "The \"ADD_INACTIVE_FORCE_CONTRIBUTIONS\" option is only supported "
        "by the complete variational approach.");

  Data().SetFDCheckType(Teuchos::getIntegralValue<INPAR::CONTACT::FDCheck>(p_aug, "FD_CHECK"));

  Data().SetPotential(Teuchos::rcp(new Potential(*this, *augDataPtr_)));

  Data().SetMatrixRowColTransformer(Teuchos::rcp(new MORTAR::MatrixRowColTransformer(4)));

  const int par_redist_interval = p_aug.get<int>("PARALLEL_REDIST_INTERVAL");
  Data().SetPDController(
      Teuchos::rcp(new ParallelDistributionController(*this, Data(), par_redist_interval)));

  interface_.reserve(interfaces.size());
  for (plain_interface_set::const_iterator cit = interfaces.begin(); cit != interfaces.end(); ++cit)
  {
    const Teuchos::RCP<CONTACT::CoInterface>& interface = *cit;
    // cast to augmented interfaces just as sanity check
    interface_.push_back(Teuchos::rcp_dynamic_cast<CONTACT::AUG::Interface>(interface, true));
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::Strategy::IsSaddlePointSystem() const
{
  return (IsInContact() or WasInContact() or WasInContactLastTimeStep());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::plain_interface_set& CONTACT::AUG::Strategy::Interfaces() { return interface_; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const CONTACT::AUG::plain_interface_set& CONTACT::AUG::Strategy::Interfaces() const
{
  return interface_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::PostSetup(bool redistributed, bool init)
{
  if (init or redistributed)
  {
    // reassemble the global slave normal/tangential dof row maps
    AssembleGlobalSlNTDofRowMaps();

    AssembleGlobalEleMaps();
  }

  if (init)
  {
    // initialize cn
    InitializeCn(Data().ConstantCn());

    if (ParRedist()) Data().InitMatrixRowColTransformer();
  }

  // just used for the redistributed case
  if (redistributed)
  {
    // redistribute the cn-vector
    RedistributeCn();

    // redistribute the global augmented old active slave nodes map
    if ((not Data().GOldActiveSlaveNodesPtr().is_null()) and
        (Data().GOldActiveSlaveNodes().NumGlobalElements() > 0))
    {
      Data().GOldActiveSlaveNodesPtr() = DRT::UTILS::RedistributeInAccordanceWithReference(
          SlRowNodes(), Data().GOldActiveSlaveNodes());
    }
  }

  // in both cases the maps change and we have to re-build all matrices
  Data().SetMatrixMapsValid(false);
  Data().SetVectorMapsValid(false);

  // setup the potential class with the current maps
  Data().Potential().Setup();

  // setup the row column transformer object
  if (ParRedist()) Data().MatrixRowColTransformer().Setup();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::Update(Teuchos::RCP<const Epetra_Vector> dis)
{
  CONTACT::CoAbstractStrategy::Update(dis);
  InitializeCn(Data().ConstantCn());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AssembleGlobalSlNTDofRowMaps()
{
  Data().GSlNormalDofRowMapPtr() = Teuchos::null;
  Data().GSlTangentialDofRowMapPtr() = Teuchos::null;

  for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end(); ++cit)
  {
    const CONTACT::AUG::Interface& interface = dynamic_cast<CONTACT::AUG::Interface&>(**cit);

    Data().GSlNormalDofRowMapPtr() =
        LINALG::MergeMap(Data().GSlNormalDofRowMapPtr(), interface.SlaveRowNDofs());
    Data().GSlTangentialDofRowMapPtr() =
        LINALG::MergeMap(Data().GSlTangentialDofRowMapPtr(), interface.SlaveRowTDofs());
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AssembleGlobalEleMaps()
{
  Data().GSeleRowMapPtr() = Teuchos::null;
  Data().GSeleColMapPtr() = Teuchos::null;

  Data().GMeleRowMapPtr() = Teuchos::null;
  Data().GMeleColMapPtr() = Teuchos::null;

  for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end(); ++cit)
  {
    Data().GSeleRowMapPtr() = LINALG::MergeMap(Data().GSeleRowMapPtr(), (*cit)->SlaveRowElements());
    Data().GSeleColMapPtr() = LINALG::MergeMap(Data().GSeleColMapPtr(), (*cit)->SlaveColElements());

    Data().GMeleRowMapPtr() =
        LINALG::MergeMap(Data().GMeleRowMapPtr(), (*cit)->MasterRowElements());
    Data().GMeleColMapPtr() =
        LINALG::MergeMap(Data().GMeleColMapPtr(), (*cit)->MasterColElements());
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::InitializeCn(const double cn_init)
{
  if (cn_init <= 0.0) dserror("The initial CN value must be greater than zero!");

  if (Data().CnPtr().is_null() or Data().Cn().GlobalLength() == 0)
    Data().CnPtr() = LINALG::CreateVector(SlRowNodes(), true);

  // set all nodal cn-values to the input value
  Data().Cn().PutScalar(cn_init);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::RedistributeCn()
{
  // redistribute the cn-vector
  Teuchos::RCP<Epetra_Vector> newcn = Teuchos::rcp(new Epetra_Vector(SlRowNodes()));
  LINALG::Export(Data().Cn(), *newcn);
  Data().CnPtr() = newcn;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::DoReadRestart(IO::DiscretizationReader& reader,
    Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr)
{
  CONTACT::CoAbstractStrategy::DoReadRestart(reader, dis, cparams_ptr);
  PostSetup(false, false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::InitMortar()
{
  // for self contact, slave and master sets may have changed,
  // thus we have to update them before initializing Dn, Mn etc.
  UpdateGlobalSelfContactState();

  // (re)setup global Mortar LINALG::SparseMatrices and Epetra_Vectors
  Data().DMatrixPtr() = Teuchos::null;
  Data().MMatrixPtr() = Teuchos::null;
  Data().BMatrixPtr() = Teuchos::rcp(new LINALG::SparseMatrix(
      Data().GSlNormalDofRowMap(), 100, false, false, LINALG::SparseMatrix::FE_MATRIX));

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AssembleMortar()
{
  // for all interfaces
  for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end(); ++cit)
  {
    const CONTACT::AUG::Interface& interface = dynamic_cast<CONTACT::AUG::Interface&>(**cit);

    interface.AssembleBMatrix(Data().BMatrix());
  }

  Data().BMatrix().Complete(SlMaDoFRowMap(true), Data().GSlNormalDofRowMap());

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::SplitMortar()
{
  TEUCHOS_FUNC_TIME_MONITOR(CONTACT_FUNC_NAME);

  Data().DMatrixPtr() =
      ExtractMatrix(Data().BMatrix(), Data().GActiveNDofRowMap(), *Data().GSlDofRowMapPtr());

  Data().MMatrixPtr() =
      ExtractMatrix(Data().BMatrix(), Data().GActiveNDofRowMap(), *Data().GMaDofRowMapPtr());

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::ZeroizeStiffnessState()
{
  // *** zeroize existing global matrices ***
  Data().DGLmLinMatrix().Reset();
  Data().DGGLinMatrix().Reset();

  Data().DLmNWGapLinMatrix().Reset();
  Data().DLmTLmTMatrix().Reset();
  Data().DLmTLmTLinMatrix().Reset();

  Data().InactiveLinMatrix().Reset();
  Data().InactiveDDMatrix().Reset();
  Data().InactiveDiagMatrix().PutScalar(0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::CreateStiffnessState(const Epetra_Map& gAugInactiveSlaveDofs)
{
  // *** (re)setup global matrices ***
  Data().DGLmLinMatrixPtr() = Teuchos::rcp(new LINALG::SparseMatrix(
      SlMaDoFRowMap(true), 100, false, false, LINALG::SparseMatrix::FE_MATRIX));
  Data().DGGLinMatrixPtr() = Teuchos::rcp(new LINALG::SparseMatrix(
      SlMaDoFRowMap(true), 100, false, false, LINALG::SparseMatrix::FE_MATRIX));

  Data().DLmNWGapLinMatrixPtr() =
      Teuchos::rcp(new LINALG::SparseMatrix(Data().GActiveNDofRowMap(), 100, false, false));
  Data().DLmTLmTMatrixPtr() =
      Teuchos::rcp(new LINALG::SparseMatrix(Data().GActiveTDofRowMap(), 100, false, false));
  Data().DLmTLmTLinMatrixPtr() =
      Teuchos::rcp(new LINALG::SparseMatrix(Data().GActiveTDofRowMap(), 100, false, false));

  Data().InactiveDiagMatrixPtr() = LINALG::CreateVector(gAugInactiveSlaveDofs, true);
  Data().InactiveLinMatrixPtr() =
      Teuchos::rcp(new LINALG::SparseMatrix(gAugInactiveSlaveDofs, 100, false, false));
  Data().InactiveDDMatrixPtr() = Teuchos::rcp(new LINALG::SparseMatrix(
      SlDoFRowMap(true), 100, false, false, LINALG::SparseMatrix::FE_MATRIX));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::ZeroizeRhsState()
{
  // *** zeroize existing global matrices ***
  Data().LmN().PutScalar(0.0);
  Data().AWGap().PutScalar(0.0);

  Data().DLmTLmTRhs().PutScalar(0.0);
  Data().InactiveRhs().PutScalar(0.0);

  Data().AVec().PutScalar(0.0);
  Data().KappaVec().PutScalar(0.0);
  Data().WGap().PutScalar(0.0);
  Data().WGapAllSlNodes().PutScalar(0.0);

  Data().SlForceLm().PutScalar(0.0);
  Data().SlForceLmInactive().PutScalar(0.0);
  Data().SlForceG().PutScalar(0.0);
  Data().MaForceLm().PutScalar(0.0);
  Data().MaForceG().PutScalar(0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::CreateRhsState(const Epetra_Map& gAugInactiveSlaveDofs)
{
  // *** (re)setup global augmented Epetra_Vectors ***
  Data().LmNPtr() = Teuchos::rcp(new Epetra_Vector(Data().GActiveNDofRowMap(), true));
  Data().AWGapPtr() = Teuchos::rcp(new Epetra_Vector(Data().GActiveNDofRowMap(), true));
  Data().DLmTLmTRhsPtr() = Teuchos::rcp(new Epetra_Vector(Data().GActiveTDofRowMap(), true));
  Data().InactiveRhsPtr() = Teuchos::rcp(new Epetra_Vector(gAugInactiveSlaveDofs, true));

  Data().AVecPtr() = Teuchos::rcp(new Epetra_Vector(SlRowNodes(), true));
  Data().KappaVecPtr() = Teuchos::rcp(new Epetra_Vector(Data().GActiveNodeRowMap(), true));
  Data().WGapPtr() = Teuchos::rcp(new Epetra_Vector(Data().GActiveNDofRowMap(), true));
  Data().WGapAllSlNodesPtr() = Teuchos::rcp(new Epetra_Vector(Data().GSlNormalDofRowMap(), true));

  Data().SlForceLmPtr() = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true), true));
  Data().SlForceLmInactivePtr() = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true), true));
  Data().SlForceGPtr() = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)), true);
  Data().MaForceLmPtr() = Teuchos::rcp(new Epetra_Vector(MaDoFRowMap(true), true));
  Data().MaForceGPtr() = Teuchos::rcp(new Epetra_Vector(MaDoFRowMap(true)), true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::Initialize(enum MORTAR::ActionType actiontype)
{
  // get inactive slave dofs
  Teuchos::RCP<Epetra_Map> gAugInactiveSlaveDofs = Teuchos::null;
  LINALG::SplitMap(SlDoFRowMap(true), Data().GActiveDofRowMap());

  switch (actiontype)
  {
    case MORTAR::eval_force_stiff:
    {
      if (Data().MatrixMapsValid())
      {
        ZeroizeStiffnessState();
      }
      else
      {
        gAugInactiveSlaveDofs = LINALG::SplitMap(SlDoFRowMap(true), Data().GActiveDofRowMap());
        CreateStiffnessState(*gAugInactiveSlaveDofs);
        Data().SetMatrixMapsValid(true);
      }
    }
    // no break
    case MORTAR::eval_force:
    case MORTAR::eval_static_constraint_rhs:
    {
      if (Data().VectorMapsValid())
      {
        ZeroizeRhsState();
      }
      else
      {
        if (gAugInactiveSlaveDofs.is_null())
          gAugInactiveSlaveDofs = LINALG::SplitMap(SlDoFRowMap(true), Data().GActiveDofRowMap());
        CreateRhsState(*gAugInactiveSlaveDofs);
        Data().SetVectorMapsValid(true);
      }

      break;
    }
    default:
    {
      dserror(
          "Unsupported action type detected: %s", MORTAR::ActionType2String(actiontype).c_str());
      exit(EXIT_FAILURE);
    }
  }

  if (Data().IsFriction())
    dserror(
        "AugmentedLagrangeStrategy::Initialize: "
        "Frictional case is not yet considered!");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::EvalForceStiff(CONTACT::ParamsInterface& cparams)
{
  // call the evaluate force routine
  EvalForce(cparams);

  // --- Assemble stiffness matrix ---------------------------------------
  AssembleContactStiff();

  PostEvalForceStiff(cparams);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::PostEvalForceStiff(CONTACT::ParamsInterface& cparams)
{
  if (not IsInContact() and not WasInContact() and not WasInContactLastTimeStep()) return;

  // --- DEBUGGING -------------------------------------------------------
  switch (Data().FDCheckType())
  {
    // Finite Difference check at Gauss-point level
    case INPAR::CONTACT::FDCheck::gauss_point:
    {
      AugFDCheckGP(cparams);
      break;
    }
    // Finite Difference check at global level
    case INPAR::CONTACT::FDCheck::global:
    {
      AugFDCheckGlobal(cparams);
      break;
    }
    case INPAR::CONTACT::FDCheck::off:
    {
      /* do nothing */
      break;
    }
    default:
    {
      dserror("Unknown FDCheck type. (enum=%d)", Data().FDCheckType());
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::InitEvalInterface(Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr)
{
  TEUCHOS_FUNC_TIME_MONITOR(CONTACT_FUNC_NAME);

  // set variational approach
  cparams_ptr->SetVariationalApproachType(Data().VariationalApproachType());
  Data().PDController().setup(*cparams_ptr);

  // get type of parallel strategy
  INPAR::MORTAR::GhostingStrategy strat = Data().GhostingStrategy();

  // Evaluation for all interfaces
  for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end(); ++cit)
  {
    CONTACT::AUG::Interface& interface = dynamic_cast<CONTACT::AUG::Interface&>(**cit);

    // initialize / reset interfaces
    interface.Initialize();

    // store required integration time
    Data().IntTime() += interface.Inttime();
    switch (strat)
    {
      /*----------------------------------------------------------*
       |  Fully redundant ghosting of master side                 |
       *----------------------------------------------------------*/
      case INPAR::MORTAR::ghosting_redundant:
      {
        EvalInterface(interface, 0, cparams_ptr);

        break;
      }
      default:
      {
        dserror(
            "ERROR: Augmented Lagrange strategy supports only "
            "\"ghosting_redundant\" as \"GHOSTING_STRATEGY\".");
        break;
      }
    }
  }  // end interface loop

  Data().PDController().check(*cparams_ptr);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::SpreadGlobalSeleEvalTimesToInterfaces()
{
  // Evaluation for all interfaces
  for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end(); ++cit)
  {
    CONTACT::AUG::Interface& interface = dynamic_cast<CONTACT::AUG::Interface&>(**cit);

    interface.StoreSeleEvalTimes(*Data().GSeleEvalTimesPtr());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::CheckParallelDistribution(const GlobalTimeMonitor& global_timer)
{
  const double my_total_time = global_timer.getMyTotalTime();
  UpdateParallelDistributionStatus(my_total_time);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::AUG::Strategy::DynRedistributeContact(
    const Teuchos::RCP<const Epetra_Vector>& dis, const int nlniter)
{
  return Data().PDController().redistribute(dis, nlniter);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::EvalInterface(CONTACT::AUG::Interface& interface, const int rriter,
    const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  const enum MORTAR::ActionType atype = cparams_ptr->GetActionType();
  switch (atype)
  {
    case MORTAR::eval_force:
    case MORTAR::eval_force_stiff:
    {
      // evaluate averaged weighted gap
      interface.Evaluate(rriter, cparams_ptr);

      // evaluate remaining entities and linearization
      interface.RedEvaluate(cparams_ptr);

      break;
    }
    case MORTAR::eval_static_constraint_rhs:
    {
      interface.EvalActiveContributions(rriter, cparams_ptr);
      interface.RedEvaluate(cparams_ptr);

      break;
    }
    case MORTAR::eval_wgap_gradient_error:
    {
      interface.EvalActiveContributions(rriter, cparams_ptr);

      break;
    }
    default:
    {
      dserror("What shall be integrated? (enum=%d | \"%s\")", atype,
          MORTAR::ActionType2String(atype).c_str());

      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::RunPostComputeX(const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  /* Since the augmented Lagrangian strategy does not support any kind
   * of condensation, we use this routine just to store the Lagrange
   * multiplier increment. */
  Epetra_Vector zincr_exp(LMDoFRowMap(true));
  LINALG::Export(dir, zincr_exp);
  int err = zincr_exp.ReplaceMap(SlDoFRowMap(true));
  if (err) dserror("ReplaceMap failed with error code %d.", err);

  // get the current step length
  const double stepLength = cparams.GetStepLength();
  // ---------------------------------------------------------------------
  /* store the SCALED Lagrange multiplier increment in the contact
   * strategy */
  // ---------------------------------------------------------------------
  CATCH_EPETRA_ERROR(Data().LmIncrPtr()->Update(stepLength, zincr_exp, 0.0));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::Reset(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& dispnp, const Epetra_Vector& xnew)
{
  Data().SetCurrentEvalState(MORTAR::eval_none);
  CONTACT::CoAbstractStrategy::Reset(cparams, dispnp, xnew);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::ResetLagrangeMultipliers(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xnew)
{
  /* Since the augmented Lagrangian strategy does not support any kind
   * of condensation, we do not have to check if it is a saddle
   * point system. */
  Epetra_Vector& z = *Data().LmPtr();

  z.PutScalar(0.0);

  if (z.ReplaceMap(LMDoFRowMap(true))) dserror("ReplaceMap failed!");

  LINALG::Export(xnew, z);

  if (z.ReplaceMap(SlDoFRowMap(true))) dserror("ReplaceMap failed!");

  // ---------------------------------------------------------------------
  // store the new Lagrange multiplier in the nodes
  // ---------------------------------------------------------------------
  StoreNodalQuantities(MORTAR::StrategyBase::lmupdate);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::PreEvalForce(CONTACT::ParamsInterface& cparams)
{
  // set current evaluation action type
  SetCurrentEvalState(cparams);

  /*--------------------------------------------------------------*
   | For self-contact the master/slave sets are updated within the|
   | contact search, see SelfBinaryTree.                          |
   | Therefore, we have to initialize the mortar matrices after   |
   | interface evaluations.                                       |
   *--------------------------------------------------------------*/
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr = Teuchos::rcp(&cparams, false);
  if (IsSelfContact())
  {
    // evaluate mortar terms (integrate...)
    InitEvalInterface(cparams_ptr);
    // initialize mortar matrices and vectors
    InitMortar();
    // assemble mortar terms into global matrices
    AssembleMortar();
  }
  else
  {
    // initialize mortar matrices and vectors
    InitMortar();
    // evaluate mortar terms (integrate...)
    InitEvalInterface(cparams_ptr);
    // assemble mortar terms into global matrices
    AssembleMortar();
  }

  if (cparams.IsPredictor())
  {
    // evaluate relative movement for friction
    EvaluateRelMovPredict();
  }
  else
    EvaluateRelMov();

  // update active set
  UpdateActiveSetSemiSmooth(cparams);

  /* Split the Dn/Mn matrices to get only the active rows
   * (only necessary for the augmented Lagrangian formulation) */
  SplitMortar();

  // initialize all rhs vectors and linearization matrices
  Initialize(cparams.GetActionType());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::EvalForce(CONTACT::ParamsInterface& cparams)
{
  // --- Prepare the evaluation and integrate all quantities ------------------
  PreEvalForce(cparams);

  // --- Assemble the gap vectors ---------------------------------------------
  AssembleGap();

  // --- compute the augmented forces -----------------------------------------
  EvalAugmentedForces();

  // --- Assemble the right hand side terms -----------------------------------
  AssembleContactRHS();

  /* Evaluate structural and constraint rhs. This is also necessary, if the
   * rhs did not change during the predictor step, but a redistribution was
   * executed! */
  EvalStrContactRHS();  // update structural contact rhs
  EvalConstrRHS();      // update the constrRHS

  PostEvalForce(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::PostEvalForce(CONTACT::ParamsInterface& cparams)
{
  // Check linear and angular momentum conservation
  if (cparams.GetActionType() == MORTAR::eval_force or  // only one per Newton
      cparams.GetNlnIter() == 0)                        // predictor
    CheckConservationLaws(cparams);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AssembleGap()
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return;
  /*--------------------------------------------------------------------*
   | Assembly                                                           |
   *--------------------------------------------------------------------*
   | --> weighted gap                                                   |
   | --> averaged weighted gap                                          |
   *--------------------------------------------------------------------*/
  for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end(); ++cit)
  {
    CONTACT::AUG::Interface& interface = dynamic_cast<CONTACT::AUG::Interface&>(**cit);

    interface.AssembleActiveGapVectors(Data().AWGap(), Data().WGap());
    interface.AssembleGapVectorOfAllSlNodes(Data().WGapAllSlNodes());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::AUG::Strategy::AssembleContactRHS()
{
  TEUCHOS_FUNC_TIME_MONITOR(CONTACT_FUNC_NAME);

  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return false;
  /*--------------------------------------------------------------------*
   |             ASSEMBLE THE CONTACT RIGHT HAND SIDE                   |
   *--------------------------------------------------------------------*/
  /*--------------------------------------------------------------------*
   | Assembly                                                           |
   *--------------------------------------------------------------------*
   | --> normal Lagrange multiplier                                     |
   | --> tangential constraint right hand side for the frictionless case|
   | --> normal and tangential inactive rhs                             |
   *--------------------------------------------------------------------*/
  for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end(); ++cit)
  {
    CONTACT::AUG::Interface& interface = dynamic_cast<CONTACT::AUG::Interface&>(**cit);

    // --- augmented Lagrange formulation --------------------------------
    // --- FORCE BALANCE -------------------------------------------------
    interface.AssembleLmNVector(Data().LmN());

    // --- CONSTRAINTS ---------------------------------------------------
    // active - normal direction
    // --> wGapRhs_

    // tributary area of inactive and active nodes
    interface.AssembleAugAVector(Data().AVec(), Data().KappaVec());

    // active - tangential direction
    interface.AssembleDLmTLmTRhs(Data().DLmTLmTRhs());

    // inactive - all directions
    interface.AssembleAugInactiveRhs(Data().InactiveRhs(), Data().Cn(), InactiveScaleFactor());
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AssembleContactStiff()
{
  TEUCHOS_FUNC_TIME_MONITOR(CONTACT_FUNC_NAME);

  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return;

  /*--------------------------------------------------------------------*
   |             ASSEMBLE THE TANGENTIAL STIFFNESS MATRIX               |
   *--------------------------------------------------------------------*/
  // --- augmented Lagrange formulation ----------------------------------
  for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end(); ++cit)
  {
    const CONTACT::AUG::Interface& interface = dynamic_cast<CONTACT::AUG::Interface&>(**cit);

    // --- Force Balance ------------------------------------------------
    // linearization w.r.t. displ.
    interface.AssembleDGLmLinMatrix(Data().DGLmLinMatrix());
    interface.AssembleDGGLinMatrix(Data().DGGLinMatrix(), Data().Cn());

    // --- Constraints --------------------------------------------------
    // linearization w.r.t. LM
    interface.AssembleDLmTLmTMatrix(Data().DLmTLmTMatrix());
    interface.AssembleAugInactiveDiagMatrix(
        Data().InactiveDiagMatrix(), Data().Cn(), InactiveScaleFactor());

    // linearization w.r.t. displ.
    // active - normal direction
    interface.AssembleDLmNWGapLinMatrix(Data().DLmNWGapLinMatrix());

    // active - tangential direction
    interface.AssembleDLmTLmTLinMatrix(Data().DLmTLmTLinMatrix());

    /*--------------------------------------------------------------------*
     | The linearization of the nodal area w.r.t. the displ. for inactive |
     | nodes can help to prevent oscillations of the active set, because  |
     | the reduction of the inactive lm-value is decelerated.             |
     |                                                                    |
     | In general, the following relation does NOT hold anymore:          |
     |                    Delta(z_{n,i}^{k+1}) = - z_{n,i}^{k}            |
     |                          z_{n,i}^{k+1}  =   0                      |
     *--------------------------------------------------------------------*/
    interface.AssembleAugInactiveLinMatrix(
        Data().InactiveLinMatrix(), Data().Cn(), InactiveScaleFactor());

    interface.AssembleInactiveDDMatrix(
        Data().InactiveDDMatrix(), Data().Cn(), InactiveScaleFactor());
  }

  // --- START - FillComplete matrices ----------------------------------
  // domainmap: columnmap | rangemap: rowmap
  // get inactive slave dofs
  Teuchos::RCP<Epetra_Map> gAugInactiveSlaveDofs =
      LINALG::SplitMap(SlDoFRowMap(true), Data().GActiveDofRowMap());

  // --- Force Balance --------------------------------------------------
  // linearization w.r.t. displ.
  Data().DGLmLinMatrix().Complete(SlMaDoFRowMap(true), SlMaDoFRowMap(true));
  Data().DGGLinMatrix().Complete(SlMaDoFRowMap(true), SlMaDoFRowMap(true));

  // --- Constraints ----------------------------------------------------
  // linearization w.r.t. LM
  Data().DLmTLmTMatrix().Complete(Data().GActiveTDofRowMap(), Data().GActiveTDofRowMap());

  // linearization w.r.t. displ.
  Data().DLmNWGapLinMatrix().Complete(SlMaDoFRowMap(true), Data().GActiveNDofRowMap());
  Data().DLmTLmTLinMatrix().Complete(SlDoFRowMap(true), Data().GActiveTDofRowMap());
  Data().InactiveLinMatrix().Complete(SlDoFRowMap(true), *gAugInactiveSlaveDofs);
  Data().InactiveDDMatrix().Complete(SlDoFRowMap(true), SlDoFRowMap(true));

  // --- END - FillComplete matrices ------------------------------------

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::UpdateActiveSetSemiSmooth(const CONTACT::ParamsInterface& cparams)
{
  ActiveSet active_set(*this);
  active_set.Update(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::EvalStrContactRHS()
{
  if (!IsInContact() and !WasInContact() and !WasInContactLastTimeStep())
  {
    Data().StrContactRhsPtr() = Teuchos::null;
    return;
  }
  Data().StrContactRhsPtr() = Teuchos::rcp(new Epetra_Vector(*ProblemDofs(), true));


  // For self contact, slave and master sets may have changed,
  if (IsSelfContact())
    dserror("ERROR: Augmented Lagrange Formulation: Self contact is not yet considered!");

  // --- add contact force terms ----------------------------------------------
  // *** Slave side ***
  Epetra_Vector augfs(Data().SlForceLm());
  CATCH_EPETRA_ERROR(augfs.Update(-1.0, Data().SlForceG(), 1.0));
  if (Data().AddInactivForceContributions())
    CATCH_EPETRA_ERROR(augfs.Update(1.0, Data().SlForceLmInactive(), 1.0));

  Epetra_Vector augfs_exp(*ProblemDofs());
  LINALG::Export(augfs, augfs_exp);
  Data().StrContactRhs().Scale(-1.0, augfs_exp);

  // Master side
  Epetra_Vector augfm(Data().MaForceLm());
  CATCH_EPETRA_ERROR(augfm.Update(-1.0, Data().MaForceG(), 1.0));

  Epetra_Vector augfm_exp(*ProblemDofs());
  LINALG::Export(augfm, augfm_exp);
  CATCH_EPETRA_ERROR(Data().StrContactRhs().Update(-1.0, augfm_exp, 1.0));

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::EvalConstrRHS()
{
  if (!IsInContact() and !WasInContact() and !WasInContactLastTimeStep())
  {
    // (re)setup the vector
    Data().ConstrRhsPtr() = Teuchos::null;
    return;
  }

  // initialize constraint r.h.s. (still with wrong map)
  Teuchos::RCP<Epetra_Vector> augConstrRhs =
      Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true), true));

  // We solve for the incremental Lagrange multiplier dz_. Hence,
  // we can keep the contact force terms on the right-hand side!
  AddContributionsToConstrRHS(*augConstrRhs);

  // replace row map
  augConstrRhs->ReplaceMap(LMDoFRowMap(true));

  // export and set constraint rhs vector
  if (ParRedist())
  {
    Data().ConstrRhsPtr() = Teuchos::rcp(new Epetra_Vector(LMDoFRowMap(false)));
    LINALG::Export(*augConstrRhs, *Data().ConstrRhsPtr());
  }
  else
    Data().ConstrRhsPtr() = augConstrRhs;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::EvalStaticConstraintRHS(CONTACT::ParamsInterface& cparams)
{
  SetCurrentEvalState(cparams);
  InitEvalInterface(cparams);

  // --- Assemble the gap vectors ---------------------------------------------
  Initialize(cparams.GetActionType());
  AssembleGap();

  // --- Evaluate only the forces coming from the constraints -----------------
  EvalConstraintForces();

  // --- Zeroize all contributions from the LM values -------------------------
  ZeroizeLMForces();

  // --- Assemble the right hand side terms -----------------------------------
  AssembleContactRHS();

  // --- update the rhs w.r.t. all constraint contributions -------------------
  EvalStrContactRHS();
  EvalConstrRHS();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AddContributionsToConstrRHS(Epetra_Vector& augConstrRhs) const
{
  // Add active constraints in normal direction:
  LINALG::AssembleMyVector(0.0, augConstrRhs, 1.0, *Data().WGapPtr());

  if (IO::cout.RequestedOutputLevel() >= IO::debug)
  {
    double wgap_nrm2 = 0.0;
    Data().WGapPtr()->Norm2(&wgap_nrm2);
    IO::cout << __FUNCTION__ << " [wgap-norm2] = " << wgap_nrm2 << IO::endl;
  }

  // Add inactive constraints
  LINALG::AssembleMyVector(1.0, augConstrRhs, 1.0, *Data().InactiveRhsPtr());

  // Add tangential frictionless constraints
  LINALG::AssembleMyVector(1.0, augConstrRhs, 1.0, *Data().DLmTLmTRhsPtr());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::CheckConservationLaws(CONTACT::ParamsInterface& cparams)
{
  // get forces due to Lagrange multipliers of the slave side
  const Teuchos::RCP<Epetra_Vector>& augfs_lm = Data().SlForceLmPtr();

  // get forces due to the regularization term of the slave side
  const Teuchos::RCP<Epetra_Vector>& augfs_g = Data().SlForceGPtr();

  // get forces due to Lagrange multipliers of the master side
  const Teuchos::RCP<Epetra_Vector>& augfm_lm = Data().MaForceLmPtr();

  // get forces due to the regularization term of the master side
  const Teuchos::RCP<Epetra_Vector>& augfm_g = Data().MaForceGPtr();

  // complete force on the slave side
  Epetra_Vector augfs(Data().SlForceLm());
  CATCH_EPETRA_ERROR(augfs.Update(-1.0, Data().SlForceG(), 1.0));
  if (Data().AddInactivForceContributions())
    CATCH_EPETRA_ERROR(augfs.Update(1.0, Data().SlForceLmInactive(), 1.0));

  // compelet force on the master side
  Epetra_Vector augfm(Data().MaForceLm());
  CATCH_EPETRA_ERROR(augfm.Update(-1.0, Data().MaForceG(), 1.0));

  /*-------------------------------*
   | LINEAR MOMENTUM CONSERVATION  |
   *-------------------------------*/
  if (Data().PrintLinearMomConservation())
  {
    const int probdim = Dim();

    double lssum[3] = {0.0, 0.0, 0.0};  // local slave sum
    double gssum[3] = {0.0, 0.0, 0.0};  // global slave sum
    double lmsum[3] = {0.0, 0.0, 0.0};  // local master sum
    double gmsum[3] = {0.0, 0.0, 0.0};  // global master sum
    double gcsum[3] = {0.0, 0.0, 0.0};  // global complete sum
    // slave
    for (int i = 0; i < augfs_lm->MyLength(); ++i) lssum[i % probdim] += (*augfs_lm)[i];

    Comm().SumAll(lssum, gssum, probdim);
    // master
    for (int i = 0; i < augfm_lm->MyLength(); ++i) lmsum[i % probdim] += (*augfm_lm)[i];

    Comm().SumAll(lmsum, gmsum, probdim);
    // complete balance check
    for (int i = 0; i < probdim; ++i) gcsum[i] = gssum[i] + gmsum[i];

    if (Comm().MyPID() == 0)
    {
      std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<\n";
      std::cout << ">>      Linear Momentum Conservation      <<\n";
      std::cout << ">>  comp.-wise in x-, y- and z-direction  <<\n";
      std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<\n";
      std::cout << ">>      Standard terms (lm)               <<" << std::endl;
      std::cout << "SLAVE:   " << std::setw(14) << gssum[0] << "," << std::setw(14) << gssum[1]
                << "," << std::setw(14) << gssum[2] << std::endl;

      std::cout << "MASTER:  " << std::setw(14) << gmsum[0] << "," << std::setw(14) << gmsum[1]
                << "," << std::setw(14) << gmsum[2] << std::endl;

      std::cout << "Balance: " << std::setw(14) << gcsum[0] << "," << std::setw(14) << gcsum[1]
                << "," << std::setw(14) << gcsum[2] << std::endl;
      std::cout << "--------------------------------------------" << std::endl;
      std::cout << ">>      Regularization terms (awgap)      <<" << std::endl;
    }

    for (int i = 0; i < probdim; ++i)
    {
      if (std::abs(gcsum[i]) >
          1.0e-11 * std::max(1.0, std::max(std::abs(gssum[i]), std::abs(gmsum[i]))))
        dserror("Conservation of linear momentum is not fulfilled!");
    }

    // slave
    std::fill(lssum, lssum + probdim, 0.0);
    for (int i = 0; i < augfs_g->MyLength(); ++i) lssum[i % probdim] += (*augfs_g)[i];

    Comm().SumAll(lssum, gssum, probdim);

    // master
    std::fill(lmsum, lmsum + probdim, 0.0);
    for (int i = 0; i < augfm_g->MyLength(); ++i) lmsum[i % probdim] += (*augfm_g)[i];

    Comm().SumAll(lmsum, gmsum, probdim);

    // complete balance check
    for (int i = 0; i < probdim; ++i) gcsum[i] = gssum[i] + gmsum[i];

    if (Comm().MyPID() == 0)
    {
      std::cout << "SLAVE:   " << std::setw(14) << gssum[0] << "," << std::setw(14) << gssum[1]
                << "," << std::setw(14) << gssum[2] << std::endl;

      std::cout << "MASTER:  " << std::setw(14) << gmsum[0] << "," << std::setw(14) << gmsum[1]
                << "," << std::setw(14) << gmsum[2] << std::endl;

      std::cout << "Balance: " << std::setw(14) << gcsum[0] << "," << std::setw(14) << gcsum[1]
                << "," << std::setw(14) << gcsum[2] << std::endl;
      std::cout << "--------------------------------------------" << std::endl;
      std::cout << ">>      Complete                          <<" << std::endl;
    }

    for (int i = 0; i < probdim; ++i)
    {
      if (std::abs(gcsum[i]) >
          1.0e-11 * std::max(1.0, std::max(std::abs(gssum[i]), std::abs(gmsum[i]))))
        dserror("Conservation of linear momentum is not fulfilled!");
    }

    // slave
    std::fill(lssum, lssum + probdim, 0.0);
    for (int i = 0; i < augfs.MyLength(); ++i) lssum[i % probdim] += augfs[i];

    Comm().SumAll(lssum, gssum, probdim);

    // master
    std::fill(lmsum, lmsum + probdim, 0.0);
    for (int i = 0; i < augfm.MyLength(); ++i) lmsum[i % probdim] += augfm[i];

    Comm().SumAll(lmsum, gmsum, probdim);

    // complete balance check
    for (int i = 0; i < probdim; ++i) gcsum[i] = gssum[i] + gmsum[i];

    if (Comm().MyPID() == 0)
    {
      std::cout << "SLAVE:   " << std::setw(14) << gssum[0] << "," << std::setw(14) << gssum[1]
                << "," << std::setw(14) << gssum[2] << std::endl;

      std::cout << "MASTER:  " << std::setw(14) << gmsum[0] << "," << std::setw(14) << gmsum[1]
                << "," << std::setw(14) << gmsum[2] << std::endl;

      std::cout << "Balance: " << std::setw(14) << gcsum[0] << "," << std::setw(14) << gcsum[1]
                << "," << std::setw(14) << gcsum[2] << std::endl;
    }

    for (int i = 0; i < probdim; ++i)
    {
      if (std::abs(gcsum[i]) >
          1.0e-11 * std::max(1.0, std::max(std::abs(gssum[i]), std::abs(gmsum[i]))))
        dserror("Conservation of linear momentum is not fulfilled!");
    }
  }
  /*-------------------------------*
   | ANGULAR MOMENTUM CONSERVATION |
   *-------------------------------*/
  if (Data().PrintAngularMomConservation())
  {
    if (Comm().MyPID() == 0)
    {
      std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<\n";
      std::cout << ">>      Angular Momentum Conservation     <<\n";
      std::cout << ">>  comp.-wise in x-, y- and z-direction  <<\n";
      std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
    }

    unsigned icount = 0;
    const int nln_iter = cparams.GetNlnIter();

    LINALG::SerialDenseMatrix conservation_data(18, 1, false);
    for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end();
         ++cit, ++icount)
    {
      const CoInterface& interface = **cit;

      if (Comm().MyPID() == 0)
      {
        std::cout << ">>----- Interface " << std::setw(2) << icount;
        std::cout << " ---------------------<<" << std::endl;
        std::cout << ">>      Standard terms (lm)               <<" << std::endl;
      }
      interface.EvalResultantMoment(*augfs_lm, *augfm_lm, &conservation_data);
      CONTACT::UTILS::WriteConservationDataToFile(Comm().MyPID(), icount, nln_iter,
          conservation_data, cparams.GetOutputFilePath(), "lm_terms");

      if (Comm().MyPID() == 0)
      {
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << ">>      Regularization terms (awgap)      <<" << std::endl;
      }
      interface.EvalResultantMoment(*augfs_g, *augfm_g, &conservation_data);
      CONTACT::UTILS::WriteConservationDataToFile(Comm().MyPID(), icount, nln_iter,
          conservation_data, cparams.GetOutputFilePath(), "regularization_terms");

      if (Comm().MyPID() == 0)
      {
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << ">>      Complete                          <<" << std::endl;
      }
      interface.EvalResultantMoment(augfs, augfm, &conservation_data);
      CONTACT::UTILS::WriteConservationDataToFile(Comm().MyPID(), icount, nln_iter,
          conservation_data, cparams.GetOutputFilePath(), "complete");
    }
    if (Comm().MyPID() == 0)
      std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
  }
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::EvalAugmentedForces()
{
  if (!IsInContact() and !WasInContact() and !WasInContactLastTimeStep()) return;

  // augmented force terms
  enum ForceTerm
  {
    lm_term = 0,
    g_term = 1
  };

  const Epetra_Vector& zn = Data().Potential().GetZnActive();
  Epetra_Vector awgapn(Data().GActiveNDofRowMap(), true);
  LINALG::ExtractMyVector(Data().AWGap(), awgapn);

  // scale the averaged weighted gap elementwise by cn
  MultiplyElementwise(Data().Cn(), Data().GActiveNodeRowMap(), awgapn, false);

  double* values[2] = {NULL, NULL};

  values[lm_term] = zn.Values();
  values[g_term] = awgapn.Values();

  const Epetra_MultiVector zn_awgapn(View, Data().GActiveNDofRowMap(), values, 2);

  /*----------------- SLAVE SIDE ---------------------------------------------*/
  double* f_values[2] = {NULL, NULL};
  f_values[lm_term] = Data().SlForceLm().Values();
  f_values[g_term] = Data().SlForceG().Values();

  // interface forces on the slave side
  Epetra_MultiVector slForces(View, SlDoFRowMap(true), f_values, 2);

  Data().DMatrix().Multiply(true, zn_awgapn, slForces);

  for (auto& inter : interface_)
  {
    const Interface& aug_inter = dynamic_cast<const Interface&>(*inter);
    aug_inter.Add_Var_A_GG(*slForces(g_term), Data().Cn());
    aug_inter.Assemble_SlForceLmInactive(
        Data().SlForceLmInactive(), Data().Cn(), InactiveScaleFactor());
  }

  /*----------------- MASTER SIDE --------------------------------------------*/
  f_values[lm_term] = Data().MaForceLm().Values();
  f_values[g_term] = Data().MaForceG().Values();

  // interface forces on the slave side
  Epetra_MultiVector maForces(View, MaDoFRowMap(true), f_values, 2);

  Data().MMatrix().Multiply(true, zn_awgapn, maForces);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::EvalConstraintForces()
{
  if (!IsInContact() and !WasInContact() and !WasInContactLastTimeStep()) return;

  Epetra_Vector awgapn(Data().GActiveNDofRowMap(), true);
  LINALG::ExtractMyVector(Data().AWGap(), awgapn);

  // scale the averaged weighted gap elementwise by cn
  MultiplyElementwise(Data().Cn(), Data().GActiveNodeRowMap(), awgapn, false);

  Data().DMatrix().Multiply(true, awgapn, Data().SlForceG());

  Data().MMatrix().Multiply(true, awgapn, Data().MaForceG());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::ZeroizeLMForces()
{
  // set the lagrange multiplier forces to zero
  Data().SlForceLm().PutScalar(0.0);
  Data().MaForceLm().PutScalar(0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AugForces(Epetra_Vector& gAugFs_lm, Epetra_Vector& gAugFs_g,
    Epetra_Vector& gAugFm_lm, Epetra_Vector& gAugFm_g) const
{
  if (!IsInContact()) return;

  /****************** SLAVE SIDE ****************************************/
  // *** standard Lagrange multiplier fraction ***
  // Export
  LINALG::Export(*Data().SlForceLmPtr(), gAugFs_lm);

  // *** regularization fraction ***
  // Export
  LINALG::Export(*Data().SlForceGPtr(), gAugFs_g);
  gAugFs_g.Scale(-1.0);

  /****************** MASTER SIDE ***************************************/
  // *** standard lagrange multiplier fraction ***
  // Export
  LINALG::Export(*Data().MaForceLmPtr(), gAugFm_lm);

  // *** regularization fraction ***
  // Export
  LINALG::Export(*Data().MaForceGPtr(), gAugFm_g);
  gAugFm_g.Scale(-1.0);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::ComputeContactStresses()
{
  // reset contact stress class variables
  Data().StressNormalPtr() = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
  Data().StressTangentialPtr() = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));

  // loop over all interfaces
  for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end(); ++cit)
  {
    const CoInterface& interface = **cit;

    // loop over all slave row nodes on the current interface
    for (int j = 0; j < interface.SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface.SlaveRowNodes()->GID(j);
      DRT::Node* node = interface.Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      // be aware of problem dimension
      int dim = Dim();
      int numdof = cnode->NumDof();
      if (dim != numdof) dserror("ERROR: Inconsisteny Dim <-> NumDof");

      // get nodal normal and tangential directions
      double* nn = cnode->MoData().n();
      double* nt1 = cnode->CoData().txi();
      double* nt2 = cnode->CoData().teta();
      double lmn = cnode->MoData().lm()[0];
      double lmt1 = cnode->MoData().lm()[1];
      double lmt2 = cnode->MoData().lm()[2];

      // find indices for DOFs of current node in Epetra_Vector
      // and put node values (normal and tangential stress components) at these DOFs

      std::vector<int> locindex(dim);

      // normal stress components
      for (int dof = 0; dof < dim; ++dof)
      {
        locindex[dof] = (Data().StressNormalPtr()->Map()).LID(cnode->Dofs()[dof]);
        (*Data().StressNormalPtr())[locindex[dof]] = -lmn * nn[dof];
      }

      // tangential stress components
      for (int dof = 0; dof < dim; ++dof)
      {
        locindex[dof] = (Data().StressTangentialPtr()->Map()).LID(cnode->Dofs()[dof]);
        (*Data().StressTangentialPtr())[locindex[dof]] = -lmt1 * nt1[dof] - lmt2 * nt2[dof];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::WriteOutput(IO::DiscretizationWriter& writer) const
{
  Teuchos::RCP<Epetra_Vector> augfs_lm = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  Teuchos::RCP<Epetra_Vector> augfs_g = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  Teuchos::RCP<Epetra_Vector> augfm_lm = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  Teuchos::RCP<Epetra_Vector> augfm_g = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));

  // evaluate augmented contact forces
  AugForces(*augfs_lm, *augfs_g, *augfm_lm, *augfm_g);

  // contact forces on slave and master side
  writer.WriteVector("norslaveforcelm", augfs_lm);
  writer.WriteVector("norslaveforceg", augfs_g);
  writer.WriteVector("normasterforcelm", augfm_lm);
  writer.WriteVector("normasterforceg", augfm_g);

  Epetra_Vector str_row_node_owners(*ProblemNodes(), false);
  str_row_node_owners.PutScalar(-1.0);

  for (auto& cinterface : Interfaces())
  {
    const CONTACT::AUG::Interface& interface =
        dynamic_cast<const CONTACT::AUG::Interface&>(*cinterface);

    Teuchos::RCP<Epetra_Vector> irow_node_owners =
        interface.CollectRowNodeOwners(writer.GetDiscret());

    LINALG::Export(*irow_node_owners, str_row_node_owners);
  }

  writer.WriteVector("contactowner", Teuchos::rcpFromRef(str_row_node_owners), IO::nodevector);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::AUG::Strategy::GetRhsBlockPtr(
    const enum DRT::UTILS::VecBlockType& bt) const
{
  // if there are no active contact contributions
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return Teuchos::null;

  // get the desired vector and return it (read-only)
  Teuchos::RCP<const Epetra_Vector> vec_ptr = Teuchos::null;
  switch (bt)
  {
    case DRT::UTILS::block_displ:
    {
      vec_ptr = Data().StrContactRhsPtr();

      if (IO::cout.RequestedOutputLevel() >= IO::debug)
      {
        double vec_nrm2 = 0.0;
        vec_ptr->Norm2(&vec_nrm2);
        IO::cout << __FUNCTION__ << " [DRT::UTILS::block_displ] = " << vec_nrm2 << IO::endl;
      }

      break;
    }
    case DRT::UTILS::block_constraint:
    {
      vec_ptr = Data().ConstrRhsPtr();
      if (IO::cout.RequestedOutputLevel() >= IO::debug)
      {
        double vec_nrm2 = 0.0;
        vec_ptr->Norm2(&vec_nrm2);
        IO::cout << __FUNCTION__ << " [DRT::UTILS::block_constraint] = " << vec_nrm2 << IO::endl;
      }

      break;
    }
    default:
    {
      dserror("Unknown STR::VecBlockType!");
      break;
    }
  }

  return vec_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> CONTACT::AUG::Strategy::GetMatrixBlockPtr(
    const enum DRT::UTILS::MatBlockType& bt, const CONTACT::ParamsInterface* cparams) const
{
  TEUCHOS_FUNC_TIME_MONITOR(CONTACT_FUNC_NAME);

  // if there are no active contact contributions
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return Teuchos::null;

  Teuchos::RCP<LINALG::SparseMatrix> mat_ptr = Teuchos::null;
  switch (bt)
  {
    case DRT::UTILS::block_displ_displ:
    {
      Teuchos::RCP<LINALG::SparseMatrix>& kdd_ptr = mat_ptr;
      kdd_ptr = Teuchos::rcp(new LINALG::SparseMatrix(SlMaDoFRowMap(true), 100, false, true));

      // build matrix kdd
      AddContributionsToMatrixBlockDisplDispl(*kdd_ptr);
      kdd_ptr->Complete(SlMaDoFRowMap(true), SlMaDoFRowMap(true));

      break;
    }
    case DRT::UTILS::block_displ_lm:
    {
      Teuchos::RCP<LINALG::SparseMatrix>& kdz_ptr = mat_ptr;
      kdz_ptr =
          Teuchos::rcp(new LINALG::SparseMatrix(*Data().GDispDofRowMapPtr(), 100, false, true));

      // build constraint matrix kdz
      AddContributionsToMatrixBlockDisplLm(*kdz_ptr);
      kdz_ptr->Complete(SlDoFRowMap(true), *Data().GDispDofRowMapPtr());

      // transform constraint matrix kzd to lmdofmap (MatrixColTransform)
      static Teuchos::RCP<Epetra_Map> newcolmap = Teuchos::null;
      MORTAR::ReplaceColumnAndDomainMap(*kdz_ptr, *LMDoFRowMapPtr(true), &newcolmap);

      break;
    }
    case DRT::UTILS::block_lm_displ:
    {
      Teuchos::RCP<LINALG::SparseMatrix>& kzd_ptr = mat_ptr;
      kzd_ptr = Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true), 100, false, true));

      // build constraint matrix kzd
      AddContributionsToMatrixBlockLmDispl(*kzd_ptr);
      int err = kzd_ptr->ReplaceRowMap(*LMDoFRowMapPtr(true));
      if (err) dserror("ReplaceMap failed on kzd_ptr! (err = %d)", err);

      kzd_ptr->Complete(*Data().GDispDofRowMapPtr(), *LMDoFRowMapPtr(true));

      break;
    }
    case DRT::UTILS::block_lm_lm:
    {
      Teuchos::RCP<LINALG::SparseMatrix>& kzz_ptr = mat_ptr;
      kzz_ptr = Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true), 100, false, true));

      // build constraint matrix kzz
      AddContributionsToMatrixBlockLmLm(*kzz_ptr);

      // replace row map
      int err = kzz_ptr->ReplaceRowMap(*LMDoFRowMapPtr(true));
      if (err) dserror("ReplaceMap failed on kzz_ptr! (err = %d)", err);

      kzz_ptr->Complete(SlDoFRowMap(true), *LMDoFRowMapPtr(true));

      // transform constraint matrix kzz to lmdofmap (columns, 2nd step)
      static Teuchos::RCP<Epetra_Map> newcolmap = Teuchos::null;
      MORTAR::ReplaceColumnAndDomainMap(*kzz_ptr, *LMDoFRowMapPtr(true), &newcolmap);

      break;
    }
    default:
    {
      dserror("Unknown STR::MatBlockType!");
      break;
    }
  }

  // transform parallel row/column distribution
  // (only necessary in the parallel redistribution case)
  if (ParRedist())
  {
    MORTAR::MatrixRowColTransformer& transformer = Data().MatrixRowColTransformer();
    mat_ptr = transformer.RedistributedToUnredistributed(bt, *mat_ptr);
  }

  return mat_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AddContributionsToMatrixBlockDisplDispl(
    LINALG::SparseMatrix& kdd, const CONTACT::ParamsInterface* cparams) const
{
  kdd.Add(*Data().DGLmLinMatrixPtr(), false, -1.0, 1.0);
  kdd.Add(*Data().DGGLinMatrixPtr(), false, 1.0, 1.0);

  if (Data().AddInactivForceContributions())
    kdd.Add(*Data().InactiveDDMatrixPtr(), false, -1.0, 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AddContributionsToMatrixBlockDisplLm(LINALG::SparseMatrix& kdz) const
{
  kdz.Add(*Data().DMatrixPtr(), true, -1.0, 1.0);
  kdz.Add(*Data().MMatrixPtr(), true, -1.0, 1.0);

  if (Data().AddInactivForceContributions())
    kdz.Add(*Data().InactiveLinMatrixPtr(), true, -1.0, 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AddContributionsToMatrixBlockLmDispl(LINALG::SparseMatrix& kzd) const
{
  kzd.Add(*Data().DLmNWGapLinMatrixPtr(), false, 1.0, 1.0);
  kzd.Add(*Data().DLmTLmTLinMatrixPtr(), false, 1.0, 1.0);
  kzd.Add(*Data().InactiveLinMatrixPtr(), false, 1.0, 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AddContributionsToMatrixBlockLmLm(LINALG::SparseMatrix& kzz) const
{
  if (LINALG::InsertMyRowDiagonalIntoUnfilledMatrix(kzz, *Data().InactiveDiagMatrixPtr()))
  {
    Epetra_Vector kzz_diag = Epetra_Vector(kzz.RangeMap(), true);
    // extract the diagonal and avoid to replace already set values
    kzz.ExtractDiagonalCopy(kzz_diag);
    LINALG::AssembleMyVector(1.0, kzz_diag, 1.0, *Data().InactiveDiagMatrixPtr());

    // if the matrix is filled, we try to replace the diagonal
    if (kzz.ReplaceDiagonalValues(kzz_diag)) dserror("ReplaceDiagonalValues failed!");
  }

  kzz.Add(*Data().DLmTLmTMatrixPtr(), false, 1.0, 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::Strategy::ConstraintNorm() const
{
  double nrm2 = 0.0;
  Data().ConstrRhsPtr()->Norm2(&nrm2);
  return nrm2;
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
const Epetra_Vector& CONTACT::AUG::Strategy::GetWeightedGap(const enum MapType type) const
{
  switch (type)
  {
    case MapType::active_slave_nodes:
    {
      if (Data().WGapPtr().is_null()) dserror("The weighted gap vector is not initialized!");

      return *Data().WGapPtr();
    }
    case MapType::all_slave_nodes:
    {
      if (Data().WGapAllSlNodesPtr().is_null())
        dserror("The weighted gap of all slave nodes is not initialized!");

      return *Data().WGapAllSlNodesPtr();
    }
    default:
    {
      dserror("Unknown weighted gap type! (type=%d)", type);
      exit(EXIT_FAILURE);
    }
  }
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::SparseMatrix> CONTACT::AUG::Strategy::GetWeightedGapGradient(
    const enum WGapGradientType grad_type, const enum MapType map_type) const
{
  switch (grad_type)
  {
    case WGapGradientType::force_balance:
    {
      switch (map_type)
      {
        case MapType::all_slave_nodes:
        {
          if (Data().BMatrixPtr().is_null())
            dserror("The modified weighted gap gradient is not initialized!");

          return Data().BMatrixPtr();
        }
        case MapType::active_slave_nodes:
        {
          Teuchos::RCP<LINALG::SparseMatrix> wgap_grad =
              Teuchos::rcp(new LINALG::SparseMatrix(*Data().GActiveNDofRowMapPtr(), 100));

          if (Data().DMatrixPtr().is_null() or Data().MMatrixPtr().is_null())
            dserror("D-Matrix or/and M-Matrix are NULL!");

          wgap_grad->Add(*Data().DMatrixPtr(), false, 1.0, 0.0);
          wgap_grad->Add(*Data().MMatrixPtr(), false, 1.0, 0.0);

          wgap_grad->Complete(*Data().GSlMaDofRowMapPtr(), *Data().GActiveNDofRowMapPtr());

          return wgap_grad;
        }
        default:
          dserror("Unsupported MapType!");
          exit(EXIT_FAILURE);
      }
      break;
    }
    case WGapGradientType::constraint_enforcement:
    {
      switch (map_type)
      {
        case MapType::all_slave_nodes:
        {
          // If the complete variational approach is active, there is no need
          // to re-assemble any quantities. Due to the fact that the system is
          // symmetric, the B-Matrix can be returned.
          if (Data().VariationalApproachType() == INPAR::CONTACT::var_complete)
          {
            return GetWeightedGapGradient(
                WGapGradientType::force_balance, MapType::all_slave_nodes);
          }

          Teuchos::RCP<LINALG::SparseMatrix> wgap_grad = Teuchos::rcp(new LINALG::SparseMatrix(
              *Data().GSlNormalDofRowMapPtr(), 100, false, false, LINALG::SparseMatrix::FE_MATRIX));

          for (plain_interface_set::const_iterator cit = interface_.begin();
               cit != interface_.end(); ++cit)
          {
            const CONTACT::AUG::Interface& interface =
                dynamic_cast<CONTACT::AUG::Interface&>(**cit);

            interface.AssembleDLmNWGapLinMatrix(*wgap_grad, map_type);
          }

          wgap_grad->Complete(*Data().GSlMaDofRowMapPtr(), *Data().GSlNormalDofRowMapPtr());

          return wgap_grad;
        }
        case MapType::active_slave_nodes:
        {
          if (Data().DLmNWGapLinMatrixPtr().is_null())
            dserror("The active weighted gap gradient is not initialized!");

          return Data().DLmNWGapLinMatrixPtr();
        }
        default:
          dserror("Unsupported MapType!");
          exit(EXIT_FAILURE);
      }
    }
    default:
      dserror("Unsupported WGapGradientType!");
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::EvalWeightedGapGradientError(CONTACT::ParamsInterface& cparams)
{
  const int num_my_activenodes = Data().GActiveNodeRowMapPtr()->NumMyElements();
  const int num_gactivenodes = Data().GActiveNodeRowMapPtr()->NumGlobalElements();
  if (num_gactivenodes == 0) return;

  std::unordered_map<int, Deriv1stMap> grad_error_ma;
  grad_error_ma.reserve(num_my_activenodes);
  std::unordered_map<int, Deriv1stMap> grad_error_jac;
  grad_error_jac.reserve(num_my_activenodes);

  cparams.SetUnorderedMap(&grad_error_ma, 0);
  cparams.SetUnorderedMap(&grad_error_jac, 1);

  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr = Teuchos::rcpFromRef(cparams);
  InitEvalInterface(cparams_ptr);

  // --- Nodal gradient error due to the convective parametric master coordinate
  {
    std::vector<std::pair<int, double>>& grad_error_ma_per_node =
        Data().NodalGradientError_MaProj();

    L2ErrorNormPerNode(grad_error_ma, grad_error_ma_per_node);
    IO::cout(IO::standard) << "Nodal gradient error: projected parametric "
                              "master coordinate\n";
    for (auto& ge_per_node : grad_error_ma_per_node)
      IO::cout(IO::standard) << "GID #" << ge_per_node.first << ", e = " << ge_per_node.second
                             << IO::endl;
  }

  // --- Nodal gradient error due to the slave jacobian determinant
  {
    std::vector<std::pair<int, double>>& grad_error_jac_per_node =
        Data().NodalGradientError_Jacobian();

    L2ErrorNormPerNode(grad_error_jac, grad_error_jac_per_node);
    IO::cout(IO::standard) << "Nodal gradient error: jacobian\n";
    for (auto& ge_per_node : grad_error_jac_per_node)
      IO::cout(IO::standard) << "GID #" << ge_per_node.first << ", e = " << ge_per_node.second
                             << IO::endl;
  }

  // --- total gradient error (master proj + slave jacobian)
  {
    std::unordered_map<int, Deriv1stMap>* error_map_vec[2] = {&grad_error_ma, &grad_error_jac};
    double my_total_error = MyTotalSquareError(error_map_vec, 2);
    double& total_error = Data().TotalGradientError();
    Comm().SumAll(&my_total_error, &total_error, 1);

    total_error /= num_gactivenodes;
    total_error = std::sqrt(total_error);

    IO::cout(IO::standard) << "total_error = " << total_error << IO::endl;
  }

  cparams.ClearAll(GEN::AnyDataContainer::DataType::unordered_map);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::MyTotalSquareError(
    const std::unordered_map<int, Deriv1stMap>* const* error_map_vec, const unsigned num_vecs)
{
  double total_error = 0.0;
  for (unsigned i = 0; i < num_vecs; ++i)
    for (auto& error_map_i : (*error_map_vec[i]))
      for (auto& error_map_ij : error_map_i.second)
        total_error += error_map_ij.second * error_map_ij.second;

  return total_error;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::L2ErrorNormPerNode(const std::unordered_map<int, Deriv1stMap>& error_map,
    std::vector<std::pair<int, double>>& error_norm_per_node)
{
  error_norm_per_node.clear();
  error_norm_per_node.resize(error_map.size(), std::pair<int, double>());

  unsigned count = 0;
  for (auto& error_map_i : error_map)
  {
    std::pair<int, double>& nodal_error_norm = error_norm_per_node[count++];
    nodal_error_norm.first = error_map_i.first;
    for (auto& error_map_ij : error_map_i.second)
      nodal_error_norm.second += error_map_ij.second * error_map_ij.second;

    // compute L2-norm for each nodal contribution
    nodal_error_norm.second = std::sqrt(nodal_error_norm.second);
    //    std::cout << "nodal_error_norm of #" << nodal_error_norm.first << " = " <<
    //        nodal_error_norm.second << std::endl;
  }

  if (count != error_map.size()) dserror("Size <--> count mismatch!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::Strategy::GetPotentialValue(
    const enum NOX::NLN::MeritFunction::MeritFctName mrt_type) const
{
  Potential& pot = Data().MutablePotential();
  pot.Compute();

  switch (mrt_type)
  {
    case NOX::NLN::MeritFunction::mrtfct_lagrangian:
    {
      return Data().Potential().Get(POTENTIAL::Type::lagrangian, POTENTIAL::SetType::all);
    }
    case NOX::NLN::MeritFunction::mrtfct_lagrangian_active:
    {
      return Data().Potential().Get(POTENTIAL::Type::lagrangian, POTENTIAL::SetType::active);
    }
    case NOX::NLN::MeritFunction::mrtfct_infeasibility_two_norm:
    {
      return Data().Potential().Get(
          POTENTIAL::Type::infeasibility_measure, POTENTIAL::SetType::all);
    }
    case NOX::NLN::MeritFunction::mrtfct_infeasibility_two_norm_active:
    {
      return Data().Potential().Get(
          POTENTIAL::Type::infeasibility_measure, POTENTIAL::SetType::active);
    }
    default:
    {
      dserror(
          "The specified merit function type is not yet supported. "
          "type = %s | %d",
          NOX::NLN::MeritFunction::MeritFuncName2String(mrt_type).c_str(), mrt_type);
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double CONTACT::AUG::Strategy::GetLinearizedPotentialValueTerms(const Epetra_Vector& dir,
    const enum NOX::NLN::MeritFunction::MeritFctName mrt_type,
    const enum NOX::NLN::MeritFunction::LinOrder linorder,
    const enum NOX::NLN::MeritFunction::LinType lintype) const
{
  switch (mrt_type)
  {
    case NOX::NLN::MeritFunction::mrtfct_lagrangian:
    {
      return GetLinearizedPotentialModelTerms(
          dir, POTENTIAL::Type::lagrangian, POTENTIAL::SetType::all, linorder, lintype);
    }
    case NOX::NLN::MeritFunction::mrtfct_lagrangian_active:
    {
      return GetLinearizedPotentialModelTerms(
          dir, POTENTIAL::Type::lagrangian, POTENTIAL::SetType::active, linorder, lintype);
    }
    case NOX::NLN::MeritFunction::mrtfct_infeasibility_two_norm:
    {
      return GetLinearizedPotentialModelTerms(
          dir, POTENTIAL::Type::infeasibility_measure, POTENTIAL::SetType::all, linorder, lintype);
    }
    case NOX::NLN::MeritFunction::mrtfct_infeasibility_two_norm_active:
    {
      return GetLinearizedPotentialModelTerms(dir, POTENTIAL::Type::infeasibility_measure,
          POTENTIAL::SetType::active, linorder, lintype);
    }
    default:
    {
      dserror(
          "The specified merit function type is not yet supported. "
          "type = %s | %d",
          NOX::NLN::MeritFunction::MeritFuncName2String(mrt_type).c_str(), mrt_type);
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::Strategy::GetLinearizedPotentialModelTerms(const Epetra_Vector& dir,
    const enum POTENTIAL::Type pottype, const enum POTENTIAL::SetType potset,
    const enum NOX::NLN::MeritFunction::LinOrder linorder,
    const enum NOX::NLN::MeritFunction::LinType lintype) const
{
  double linval = 0.0;

  Potential& pot = Data().MutablePotential();
  pot.ComputeLin(dir);

  switch (linorder)
  {
    case NOX::NLN::MeritFunction::linorder_first:
    {
      GetLinearizedPotentialModelTerms_1stOrder(pottype, potset, lintype, linval);

      break;
    }
    case NOX::NLN::MeritFunction::linorder_second:
    {
      GetLinearizedPotentialModelTerms_2ndOrder(pottype, potset, lintype, linval);

      break;
    }
    case NOX::NLN::MeritFunction::linorder_all:
    {
      GetLinearizedPotentialModelTerms_1stOrder(pottype, potset, lintype, linval);
      GetLinearizedPotentialModelTerms_2ndOrder(pottype, potset, lintype, linval);

      break;
    }
  }

  return linval;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::GetLinearizedPotentialModelTerms_1stOrder(
    const enum POTENTIAL::Type pottype, const enum POTENTIAL::SetType potset,
    const enum NOX::NLN::MeritFunction::LinType lintype, double& linval) const
{
  const CONTACT::AUG::Potential& pot = Data().Potential();

  switch (lintype)
  {
    case NOX::NLN::MeritFunction::lin_wrt_primary_dofs:
    {
      linval += pot.GetLin(pottype, potset, POTENTIAL::LinTerm::wrt_d);

      break;
    }
    case NOX::NLN::MeritFunction::lin_wrt_lagrange_multiplier_dofs:
    {
      linval += pot.GetLin(pottype, potset, POTENTIAL::LinTerm::wrt_z);

      break;
    }
    case NOX::NLN::MeritFunction::lin_wrt_mixed_dofs:
    {
      break;
    }
    case NOX::NLN::MeritFunction::lin_wrt_all_dofs:
    {
      GetLinearizedPotentialModelTerms_1stOrder(
          pottype, potset, NOX::NLN::MeritFunction::lin_wrt_primary_dofs, linval);
      GetLinearizedPotentialModelTerms_1stOrder(
          pottype, potset, NOX::NLN::MeritFunction::lin_wrt_lagrange_multiplier_dofs, linval);

      break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::GetLinearizedPotentialModelTerms_2ndOrder(
    const enum POTENTIAL::Type pottype, const enum POTENTIAL::SetType potset,
    const enum NOX::NLN::MeritFunction::LinType lintype, double& linval) const
{
  const CONTACT::AUG::Potential& pot = Data().Potential();

  switch (lintype)
  {
    case NOX::NLN::MeritFunction::lin_wrt_primary_dofs:
    {
      break;
    }
    case NOX::NLN::MeritFunction::lin_wrt_lagrange_multiplier_dofs:
    {
      linval += pot.GetLin(pottype, potset, POTENTIAL::LinTerm::wrt_z_and_z);

      break;
    }
    case NOX::NLN::MeritFunction::lin_wrt_mixed_dofs:
    {
      linval += pot.GetLin(pottype, potset, POTENTIAL::LinTerm::wrt_d_and_z);

      break;
    }
    case NOX::NLN::MeritFunction::lin_wrt_all_dofs:
    {
      GetLinearizedPotentialModelTerms_2ndOrder(
          pottype, potset, NOX::NLN::MeritFunction::lin_wrt_primary_dofs, linval);
      GetLinearizedPotentialModelTerms_2ndOrder(
          pottype, potset, NOX::NLN::MeritFunction::lin_wrt_lagrange_multiplier_dofs, linval);
      GetLinearizedPotentialModelTerms_2ndOrder(
          pottype, potset, NOX::NLN::MeritFunction::lin_wrt_mixed_dofs, linval);

      break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::SplitStateVector(const Epetra_Vector& full_state,
    Teuchos::RCP<Epetra_Vector>& displ_state_slma_ptr,
    Teuchos::RCP<Epetra_Vector>& z_state_active_ptr,
    Teuchos::RCP<Epetra_Vector>& z_state_inactive_ptr) const
{
  if (displ_state_slma_ptr.is_null())
    displ_state_slma_ptr = Teuchos::rcp(new Epetra_Vector(*Data().GSlMaDofRowMapPtr()));

  if (z_state_active_ptr.is_null())
    z_state_active_ptr = Teuchos::rcp(new Epetra_Vector(*Data().GActiveNDofRowMapPtr()));

  if (z_state_inactive_ptr.is_null())
  {
    Teuchos::RCP<Epetra_Map> ginactivendofs =
        LINALG::SplitMap(*Data().GSlNormalDofRowMapPtr(), *Data().GActiveNDofRowMapPtr());
    z_state_inactive_ptr = Teuchos::rcp(new Epetra_Vector(*ginactivendofs));
  }

  SplitStateVector(full_state, *displ_state_slma_ptr, *z_state_active_ptr, *z_state_inactive_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::SplitStateVector(const Epetra_Vector& full_state,
    Epetra_Vector& displ_state_slma, Epetra_Vector& z_state_active,
    Epetra_Vector& z_state_inactive) const
{
  // extract slave/master part of the displacement increment
  Epetra_Vector displ_exp(*Data().GDispDofRowMapPtr());
  LINALG::Export(full_state, displ_exp);
  LINALG::ExtractMyVector(displ_exp, displ_state_slma);

  // extract active/inactive part of the solution vector
  Epetra_Vector z_exp(LMDoFRowMap(true));
  LINALG::Export(full_state, z_exp);
  CATCH_EPETRA_ERROR(z_exp.ReplaceMap(SlDoFRowMap(true)));

  LINALG::ExtractMyVector(z_exp, z_state_active);
  LINALG::ExtractMyVector(z_exp, z_state_inactive);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::Strategy::CharacteristicInterfaceElementLength(
    const enum CONTACT::AUG::SideType stype) const
{
  //  find the maximal characteristic interface element length
  double my_max_ih = -1.0;

  for (const Teuchos::RCP<CONTACT::CoInterface>& iptr : interface_)
  {
    const Interface& interface = dynamic_cast<const Interface&>(*iptr);
    my_max_ih = std::max(interface.MyCharacteristicElementLength(stype), my_max_ih);
  }

  double max_ih = -1;
  Comm().MaxAll(&my_max_ih, &max_ih, 1);

  return max_ih;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::SetCurrentEvalState(const CONTACT::ParamsInterface& cparams)
{
  Data().SetCurrentEvalState(cparams.GetActionType());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::EvaluateReferenceState(Teuchos::RCP<const Epetra_Vector> vec)
{
  // do nothing for the augmented strategy
}
