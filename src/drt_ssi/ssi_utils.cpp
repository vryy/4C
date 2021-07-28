/*----------------------------------------------------------------------*/
/*! \file
 \brief Utility methods for SSI

 \level 1


 *------------------------------------------------------------------------------------------------*/

#include "ssi_utils.H"

#include "ssi_monolithic.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_inpar/inpar_s2i.H"

#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_utils_gid_vector.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../drt_lib/drt_matchingoctree.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                        AN,JH 09/2014 */
/* Function for checking that the different time steps are a
 multiplicative of each other                                           */

int SSI::UTILS::CheckTimeStepping(double dt1, double dt2)
{
  const double workdt1 = std::min(dt1, dt2);
  const double workdt2 = std::max(dt1, dt2);
  int i = 0;

  while (true)
  {
    i++;
    const double t1 = i * workdt1;

    if (std::abs(t1 - workdt2) < 10E-10)
      break;
    else if (t1 > workdt2)
      dserror("Chosen time steps %f and %f are not a multiplicative of each other", dt1, dt2);
  }
  return i;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                        AN,JH 10/2014 */
// Modification of time parameter list for problem with different time step size

void SSI::UTILS::ChangeTimeParameter(const Epetra_Comm& comm, Teuchos::ParameterList& ssiparams,
    Teuchos::ParameterList& scatradyn, Teuchos::ParameterList& sdyn)
{
  bool difftimestep = DRT::INPUT::IntegralValue<int>(ssiparams, "DIFFTIMESTEPSIZE");

  if (difftimestep)  // Create subproblems with different time steps
  {
    // Check correct choice of time stepping for single fields
    double scatrastep = scatradyn.get<double>("TIMESTEP");
    double solidstep = sdyn.get<double>("TIMESTEP");

    SSI::UTILS::CheckTimeStepping(scatrastep, solidstep);

    // modify global time step size
    ssiparams.set<double>("TIMESTEP", std::min(scatrastep, solidstep));
  }
  else
  {
    // -------------------------------------------------------------------
    // overrule certain parameters for coupled problems
    // -------------------------------------------------------------------
    // the default time step size
    scatradyn.set<double>("TIMESTEP", ssiparams.get<double>("TIMESTEP"));
    sdyn.set<double>("TIMESTEP", ssiparams.get<double>("TIMESTEP"));
    // maximum simulation time
    scatradyn.set<double>("MAXTIME", ssiparams.get<double>("MAXTIME"));
    sdyn.set<double>("MAXTIME", ssiparams.get<double>("MAXTIME"));
    // maximum number of timesteps
    scatradyn.set<int>("NUMSTEP", ssiparams.get<int>("NUMSTEP"));
    sdyn.set<int>("NUMSTEP", ssiparams.get<int>("NUMSTEP"));
  }

  // Check correct input of restart. Code relies that both time value RESTARTEVRYTIME and
  // RESULTSEVRYTIME are given if restart from time is applied
  double restarttime = ssiparams.get<double>("RESTARTEVRYTIME");
  double updatetime = ssiparams.get<double>("RESULTSEVRYTIME");
  if ((updatetime > 0.0) or (restarttime > 0.0))
  {
    if (updatetime <= 0.0 and restarttime <= 0.0)
    {
      dserror(
          "If time controlled output and restart is desired, both parameters RESTARTEVRYTIME and "
          "RESULTSEVRYTIME has to be set");
    }
  }

  // set restart params
  int scatrarestart;
  int structurerestart;

  if (restarttime > 0.0)
  {
    scatrarestart = SSI::UTILS::CheckTimeStepping(scatradyn.get<double>("TIMESTEP"), restarttime);
    structurerestart = SSI::UTILS::CheckTimeStepping(sdyn.get<double>("TIMESTEP"), restarttime);
  }
  else
  {
    int restart = ssiparams.get<int>("RESTARTEVRY");
    scatrarestart = restart;
    structurerestart = restart;
  }

  // set output params
  int scatraupres;
  int structureupres;

  if (updatetime > 0.0)
  {
    scatraupres = SSI::UTILS::CheckTimeStepping(scatradyn.get<double>("TIMESTEP"), updatetime);
    structureupres = SSI::UTILS::CheckTimeStepping(sdyn.get<double>("TIMESTEP"), updatetime);
  }
  else
  {
    int update = ssiparams.get<int>("RESULTSEVRY");
    scatraupres = update;
    structureupres = update;
  }

  // restart
  scatradyn.set<int>("RESTARTEVRY", scatrarestart);
  sdyn.set<int>("RESTARTEVRY", structurerestart);
  // solution output
  scatradyn.set<int>("RESULTSEVRY", scatraupres);
  sdyn.set<int>("RESULTSEVRY", structureupres);

  if (comm.MyPID() == 0)
  {
    std::cout << "====================== Overview of chosen time stepping: "
                 "==============================\n"
              << "\t Timestep scatra:           " << scatradyn.get<double>("TIMESTEP") << "\n"
              << "\t Timestep structure:        " << sdyn.get<double>("TIMESTEP") << "\n"
              << "\t Result step scatra:        " << scatradyn.get<int>("RESULTSEVRY") << "\n"
              << "\t Result step structure:     " << sdyn.get<int>("RESULTSEVRY") << "\n"
              << "\t Restart step scatra:       " << scatradyn.get<int>("RESTARTEVRY") << "\n"
              << "\t Restart step structure:    " << sdyn.get<int>("RESTARTEVRY") << "\n"
              << "================================================================================="
                 "=======\n \n";
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<std::pair<DRT::Condition* const, DRT::Condition* const>>
SSI::UTILS::BuildSlaveMasterPairing(const std::vector<DRT::Condition*>& conditions)
{
  // sort conditions by condition ID and slave/master using a map
  std::map<const int, DRT::Condition* const> slave_conditions;
  std::map<const int, DRT::Condition* const> master_conditions;
  for (auto* condition : conditions)
  {
    const auto pair =
        std::pair<const int, DRT::Condition* const>(condition->GetInt("ConditionID"), condition);
    switch (condition->GetInt("interface side"))
    {
      case INPAR::S2I::side_slave:
        slave_conditions.insert(pair);
        break;
      case INPAR::S2I::side_master:
        master_conditions.insert(pair);
        break;
      default:
        dserror("Coupling side must either be slave or master");
    }
  }

  // safety check
  if (master_conditions.size() != slave_conditions.size())
    dserror("must provide both, master and slave side conditions");

  // create slave-master pair for matching condition IDs
  std::vector<std::pair<DRT::Condition* const, DRT::Condition* const>> condition_pairs;
  for (const auto& slave_condition : slave_conditions)
  {
    for (const auto master_condition : master_conditions)
    {
      if (slave_condition.first == master_condition.first)
      {
        condition_pairs.emplace_back(
            std::make_pair(slave_condition.second, master_condition.second));
        break;
      }
    }
  }

  return condition_pairs;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<ADAPTER::Coupling> SSI::UTILS::SetupInterfaceCouplingAdapterStructure(
    Teuchos::RCP<DRT::Discretization> structdis, bool meshtying_3_domain_intersection,
    const std::string& conditionname_coupling,
    const std::string& conditionname_3_domain_intersection)
{
  auto interfacecouplingadapter = Teuchos::rcp(new ADAPTER::Coupling());

  // extract ssi coupling conditions from discretization and assign to slave/master side maps
  std::vector<DRT::Condition*> conditions(0, nullptr);
  structdis->GetCondition(conditionname_coupling, conditions);
  auto condition_pairs = BuildSlaveMasterPairing(conditions);

  // Setup coupling adapter from ssi conditions
  if (meshtying_3_domain_intersection)
  {
    // strategy:
    // - remove all nodes in conditionname_3_domain_intersection from conditions except for one
    // arbitrary condition
    // - setup coupling adapter with modified nodes

    // get 3 domain intersection conditions
    std::vector<DRT::Condition*> conditions_3_domain_intersection(0, nullptr);
    structdis->GetCondition(conditionname_3_domain_intersection, conditions_3_domain_intersection);
    auto condition_3_domain_pairs = BuildSlaveMasterPairing(conditions_3_domain_intersection);

    // safety check
    if (condition_3_domain_pairs.size() != 2)
      dserror("Currently, exactly 2 coupling pairs with 3-domain-meshtying supported");
    // both master conditions must be the same
    DRT::UTILS::HaveSameNodes(
        condition_3_domain_pairs[0].second, condition_3_domain_pairs[1].second, true);

    // select an arbitrary 3 domain meshtying condition (0) and find other conditions
    // (arbitrary_coupling_id) with same nodes
    const auto& condition_3_domain_pair = condition_3_domain_pairs[0];
    int arbitrary_coupling_id = -1;
    for (const auto& condition_pair : condition_pairs)
    {
      // do all nodes match on this proc?
      int my_match_slave_nodes = 1;
      int my_match_master_nodes = 1;

      for (int gid_line_slave : *condition_3_domain_pair.first->Nodes())
        if (!condition_pair.first->ContainsNode(gid_line_slave)) my_match_slave_nodes = 0;

      for (int gid_line_master : *condition_3_domain_pair.second->Nodes())
        if (!condition_pair.second->ContainsNode(gid_line_master)) my_match_master_nodes = 0;

      // do all nodes match on all procs? -> other condition with matching nodes is found
      int match_slave_nodes;
      int match_master_nodes;
      structdis->Comm().SumAll(&my_match_slave_nodes, &match_slave_nodes, 1);
      structdis->Comm().SumAll(&my_match_master_nodes, &match_master_nodes, 1);
      const int numproc = structdis->Comm().NumProc();
      if (match_slave_nodes == numproc and match_master_nodes == numproc)
      {
        arbitrary_coupling_id = condition_pair.first->GetInt("ConditionID");
        break;
      }
    }

    // safety check
    if (arbitrary_coupling_id == -1)
      dserror("cannot find surface condition with nodes from line condition");

    // vectors for global ids of slave and master interface nodes for each condition
    std::vector<std::vector<int>> islavenodegidvec_cond;
    std::vector<std::vector<int>> imasternodegidvec_cond;

    // loop over slave conditions and build vector of nodes for slave and master condition with same
    // coupling ID
    for (const auto& condition_pair : condition_pairs)
    {
      std::vector<int> islavenodegidvec;
      std::vector<int> imasternodegidvec;

      // Build GID vector of nodes on this proc, sort, and remove duplicates
      DRT::UTILS::AddOwnedNodeGIDVector(
          structdis, *condition_pair.first->Nodes(), islavenodegidvec);

      DRT::UTILS::AddOwnedNodeGIDVector(
          structdis, *condition_pair.second->Nodes(), imasternodegidvec);

      // remove all nodes from line conditions except on conditions with arbitrary_coupling_id
      if (condition_pair.first->GetInt("ConditionID") != arbitrary_coupling_id)
      {
        for (const auto& condition_3_domain_intersection : conditions_3_domain_intersection)
        {
          DRT::UTILS::RemoveNodeGIDsFromVector(
              structdis, *condition_3_domain_intersection->Nodes(), islavenodegidvec);
          DRT::UTILS::RemoveNodeGIDsFromVector(
              structdis, *condition_3_domain_intersection->Nodes(), imasternodegidvec);
        }
      }

      imasternodegidvec_cond.push_back(imasternodegidvec);
      islavenodegidvec_cond.push_back(islavenodegidvec);
    }

    interfacecouplingadapter->SetupCoupling(*structdis, *structdis, imasternodegidvec_cond,
        islavenodegidvec_cond, DRT::Problem::Instance()->NDim(), true, 1.0e-8);
  }
  else
  {
    // vectors for global ids of slave and master interface nodes for all conditions
    std::vector<int> inodegidvec_master;
    std::vector<int> inodegidvec_slave;

    // Build GID vector of nodes on this proc, sort, and remove duplicates
    for (const auto& condition_pair : condition_pairs)
    {
      DRT::UTILS::AddOwnedNodeGIDVector(
          structdis, *condition_pair.first->Nodes(), inodegidvec_slave);
      DRT::UTILS::AddOwnedNodeGIDVector(
          structdis, *condition_pair.second->Nodes(), inodegidvec_master);
    }

    interfacecouplingadapter->SetupCoupling(*structdis, *structdis, inodegidvec_master,
        inodegidvec_slave, DRT::Problem::Instance()->NDim(), true, 1.0e-8);
  }

  return interfacecouplingadapter;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<ADAPTER::Coupling>
SSI::UTILS::SetupInterfaceCouplingAdapterStructure3DomainIntersection(
    Teuchos::RCP<DRT::Discretization> structdis,
    const std::string& conditionname_3_domain_intersection)
{
  // strategy:
  // - setup coupling adapter with all nodes from conditions_3_domain_intersection except for nodes
  // from arbitrary condition (see SetupInterfaceCouplingAdapterStructure)

  std::vector<int> ilinenodegidvec_master;
  std::vector<int> ilinenodegidvec_slave;

  // Build slave-master pairs
  std::vector<DRT::Condition*> conditions_3_domain_intersection(0, nullptr);
  structdis->GetCondition(conditionname_3_domain_intersection, conditions_3_domain_intersection);
  auto condition_3_domain_pairs = BuildSlaveMasterPairing(conditions_3_domain_intersection);

  // safety check
  if (condition_3_domain_pairs.size() != 2)
    dserror("Currently, exactly 2 coupling pairs with 3-domain-meshtying supported");
  // both master conditions must be the same
  DRT::UTILS::HaveSameNodes(
      condition_3_domain_pairs[0].second, condition_3_domain_pairs[1].second, true);

  // Build GID vector of nodes for master and slave side except for arbitrary condition (0)
  const auto& condition_3_domain_pair = condition_3_domain_pairs[1];
  DRT::UTILS::AddOwnedNodeGIDVector(
      structdis, *condition_3_domain_pair.first->Nodes(), ilinenodegidvec_slave);
  DRT::UTILS::AddOwnedNodeGIDVector(
      structdis, *condition_3_domain_pair.second->Nodes(), ilinenodegidvec_master);

  // setup coupling adapter based on master and slave side GID node vectors
  auto coupling_line_structure = Teuchos::rcp(new ADAPTER::Coupling());
  coupling_line_structure->SetupCoupling(*structdis, *structdis, ilinenodegidvec_master,
      ilinenodegidvec_slave, DRT::Problem::Instance()->NDim(), true, 1.0e-8);

  return coupling_line_structure;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::ParameterList SSI::UTILS::CloneScaTraManifoldParams(
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& sublist_manifold_params)
{
  auto* scatra_manifold_params = new Teuchos::ParameterList(scatraparams);

  switch (DRT::INPUT::IntegralValue<INPAR::SCATRA::InitialField>(
      sublist_manifold_params, "INITIALFIELD"))
  {
    case INPAR::SCATRA::initfield_zero_field:
    {
      scatra_manifold_params->set<std::string>("INITIALFIELD", "zero_field");
      scatra_manifold_params->set<int>("INITFUNCNO", -1);
      break;
    }
    case INPAR::SCATRA::initfield_field_by_function:
    {
      scatra_manifold_params->set<std::string>("INITIALFIELD", "field_by_function");
      scatra_manifold_params->set<int>(
          "INITFUNCNO", sublist_manifold_params.get<int>("INITFUNCNO"));
      break;
    }
    case INPAR::SCATRA::initfield_field_by_condition:
    {
      scatra_manifold_params->set<std::string>("INITIALFIELD", "field_by_condition");
      scatra_manifold_params->set<int>("INITFUNCNO", -1);
      break;
    }
    default:
      dserror("Initial field type on manifold not supported.");
      break;
  }

  if (DRT::INPUT::IntegralValue<INPAR::SCATRA::OutputScalarType>(scatraparams, "OUTPUTSCALARS") !=
      INPAR::SCATRA::outputscalars_none)
    scatra_manifold_params->set<bool>("output_file_name_discretization", true);

  scatra_manifold_params->set<std::string>("ADAPTIVE_TIMESTEPPING", "No");

  return *scatra_manifold_params;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::ParameterList SSI::UTILS::ModifyScaTraParams(const Teuchos::ParameterList& scatraparams)
{
  auto* scatraparams_mutable = new Teuchos::ParameterList(scatraparams);

  if (DRT::INPUT::IntegralValue<INPAR::SCATRA::OutputScalarType>(scatraparams, "OUTPUTSCALARS") !=
      INPAR::SCATRA::outputscalars_none)
    scatraparams_mutable->set<bool>("output_file_name_discretization", true);

  return *scatraparams_mutable;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::UTILS::SSIMatrices::SSIMatrices(Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps,
    const LINALG::MatrixType ssi_matrixtype, const LINALG::MatrixType scatra_matrixtype,
    const bool is_scatra_manifold)
    : is_scatra_manifold_(is_scatra_manifold),
      scatra_matrixtype_(scatra_matrixtype),
      scatra_dofrowmap_(ssi_maps->ScaTraDofRowMap()),
      scatramanifold_dofrowmap_(Teuchos::null),
      structure_dofrowmap_(ssi_maps->StructureDofRowMap()),
      system_matrix_(Teuchos::null),
      scatra_matrix_(Teuchos::null),
      scatramanifold_structure_matrix_(Teuchos::null),
      scatra_structure_matrix_(Teuchos::null),
      structure_scatra_matrix_(Teuchos::null),
      structure_matrix_(Teuchos::null),
      manifold_matrix_(Teuchos::null),
      scatra_scatramanifold_matrix_(Teuchos::null),
      scatramanifold_scatra_matrix_(Teuchos::null)
{
  // fill maps related to scalar transport manifold if relevant
  if (is_scatra_manifold_) scatramanifold_dofrowmap_ = ssi_maps->ScaTraManifoldDofRowMap();

  InitializeSystemMatrix(ssi_maps, ssi_matrixtype);

  InitializeMainDiagMatrices(ssi_maps);

  InitializeOffDiagMatrices(ssi_maps);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::InitializeMainDiagMatrices(
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps)
{
  structure_matrix_ = SetupSparseMatrix(structure_dofrowmap_);

  switch (scatra_matrixtype_)
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      scatra_matrix_ = SetupBlockMatrix(ssi_maps->BlockMapScaTra(), ssi_maps->BlockMapScaTra());
      if (is_scatra_manifold_)
        manifold_matrix_ = SetupBlockMatrix(
            ssi_maps->BlockMapScaTraManifold(), ssi_maps->BlockMapScaTraManifold());

      break;
    }

    case LINALG::MatrixType::sparse:
    {
      scatra_matrix_ = SetupSparseMatrix(scatra_dofrowmap_);

      if (is_scatra_manifold_) manifold_matrix_ = SetupSparseMatrix(scatramanifold_dofrowmap_);

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::InitializeOffDiagMatrices(
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps)
{
  switch (scatra_matrixtype_)
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      scatra_structure_matrix_ =
          SetupBlockMatrix(ssi_maps->BlockMapScaTra(), ssi_maps->BlockMapStructure());

      structure_scatra_matrix_ =
          SetupBlockMatrix(ssi_maps->BlockMapStructure(), ssi_maps->BlockMapScaTra());

      if (is_scatra_manifold_)
      {
        scatramanifold_structure_matrix_ =
            SetupBlockMatrix(ssi_maps->BlockMapScaTraManifold(), ssi_maps->BlockMapStructure());
        scatramanifold_scatra_matrix_ =
            SetupBlockMatrix(ssi_maps->BlockMapScaTraManifold(), ssi_maps->BlockMapScaTra());
        scatra_scatramanifold_matrix_ =
            SetupBlockMatrix(ssi_maps->BlockMapScaTra(), ssi_maps->BlockMapScaTraManifold());
      }

      break;
    }

    case LINALG::MatrixType::sparse:
    {
      scatra_structure_matrix_ = SetupSparseMatrix(scatra_dofrowmap_);
      structure_scatra_matrix_ = SetupSparseMatrix(structure_dofrowmap_);

      if (is_scatra_manifold_)
      {
        scatramanifold_structure_matrix_ = SetupSparseMatrix(scatramanifold_dofrowmap_);
        scatramanifold_scatra_matrix_ = SetupSparseMatrix(scatramanifold_dofrowmap_);
        scatra_scatramanifold_matrix_ = SetupSparseMatrix(scatra_dofrowmap_);
      }

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::CompleteScaTraManifoldScaTraMatrix()
{
  switch (scatra_matrixtype_)
  {
    case LINALG::MatrixType::sparse:
      ScaTraManifoldScaTraMatrix()()->Complete(*scatra_dofrowmap_, *scatramanifold_dofrowmap_);
      break;
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
      ScaTraManifoldScaTraMatrix()()->Complete();
      break;
    default:
      dserror("Not supported LINALG::MatrixType!");
      break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::CompleteScaTraManifoldStructureMatrix()
{
  switch (scatra_matrixtype_)
  {
    case LINALG::MatrixType::sparse:
      ScaTraManifoldStructureMatrix()->Complete(*structure_dofrowmap_, *scatramanifold_dofrowmap_);
      break;
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
      ScaTraManifoldStructureMatrix()->Complete();
      break;
    default:
      dserror("Not supported LINALG::MatrixType!");
      break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::CompleteScaTraScaTraManifoldMatrix()
{
  switch (scatra_matrixtype_)
  {
    case LINALG::MatrixType::sparse:
      ScaTraScaTraManifoldMatrix()->Complete(*scatramanifold_dofrowmap_, *scatra_dofrowmap_);
      break;
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
      ScaTraScaTraManifoldMatrix()->Complete();
      break;
    default:
      dserror("Not supported LINALG::MatrixType!");
      break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::CompleteScaTraStructureMatrix()
{
  switch (scatra_matrixtype_)
  {
    case LINALG::MatrixType::sparse:
      ScaTraStructureMatrix()->Complete(*structure_dofrowmap_, *scatra_dofrowmap_);
      break;
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
      ScaTraStructureMatrix()->Complete();
      break;
    default:
      dserror("Not supported LINALG::MatrixType!");
      break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::CompleteStructureScaTraMatrix()
{
  switch (scatra_matrixtype_)
  {
    case LINALG::MatrixType::sparse:
      StructureScaTraMatrix()->Complete(*scatra_dofrowmap_, *structure_dofrowmap_);
      break;
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
      StructureScaTraMatrix()->Complete();
      break;
    default:
      dserror("Not supported LINALG::MatrixType!");
      break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::ClearMatrices()
{
  system_matrix_->Zero();
  scatra_matrix_->Zero();
  scatra_structure_matrix_->Zero();
  structure_scatra_matrix_->Zero();
  structure_matrix_->Zero();

  if (is_scatra_manifold_)
  {
    scatramanifold_structure_matrix_->Zero();
    manifold_matrix_->Zero();
    scatra_scatramanifold_matrix_->Zero();
    scatramanifold_scatra_matrix_->Zero();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::UTILS::SSIVectors::SSIVectors(
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, const bool is_scatra_manifold)
    : increment_(LINALG::CreateVector(*(ssi_maps->MapsSubProblems()->FullMap()), true)),
      is_scatra_manifold_(is_scatra_manifold),
      manifold_residual_(is_scatra_manifold
                             ? LINALG::CreateVector(*(ssi_maps->ScaTraManifoldDofRowMap()), true)
                             : Teuchos::null),
      residual_(LINALG::CreateVector(*(ssi_maps->MapsSubProblems()->FullMap()), true)),
      scatra_residual_(LINALG::CreateVector(*(ssi_maps->ScaTraDofRowMap()), true)),
      structure_residual_(LINALG::CreateVector(*(ssi_maps->StructureDofRowMap()), true))
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SSIVectors::ClearIncrement() { increment_->PutScalar(0.0); }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SSIVectors::ClearResiduals()
{
  residual_->PutScalar(0.0);
  scatra_residual_->PutScalar(0.0);
  structure_residual_->PutScalar(0.0);
  if (is_scatra_manifold_) manifold_residual_->PutScalar(0.0);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::InitializeSystemMatrix(
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, const LINALG::MatrixType ssi_matrixtype)
{
  switch (ssi_matrixtype)
  {
    case LINALG::MatrixType::block_field:
    {
      system_matrix_ =
          SetupBlockMatrix(ssi_maps->BlockMapSystemMatrix(), ssi_maps->BlockMapSystemMatrix());
      break;
    }

    case LINALG::MatrixType::sparse:
    {
      system_matrix_ = SetupSparseMatrix(ssi_maps->MapSystemMatrix());
      break;
    }

    default:
    {
      dserror("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> SSI::UTILS::SSIMatrices::SetupBlockMatrix(
    Teuchos::RCP<const LINALG::MultiMapExtractor> row_map,
    Teuchos::RCP<const LINALG::MultiMapExtractor> col_map)
{
  const int expected_entries_per_row = 81;
  const bool explicitdirichlet = false;
  const bool savegraph = true;

  return Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
      *col_map, *row_map, expected_entries_per_row, explicitdirichlet, savegraph));
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> SSI::UTILS::SSIMatrices::SetupSparseMatrix(
    const Teuchos::RCP<const Epetra_Map> row_map)
{
  const int expected_entries_per_row = 27;
  const bool explicitdirichlet = false;
  const bool savegraph = true;

  return Teuchos::rcp(
      new LINALG::SparseMatrix(*row_map, expected_entries_per_row, explicitdirichlet, savegraph));
}

/* ----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::UTILS::SSIMaps::SSIMaps(const SSI::SSIMono& ssi_mono_algorithm)
    : block_maps_sub_problems_(),
      block_map_system_matrix_(Teuchos::null),
      map_system_matrix_(Teuchos::null),
      maps_sub_problems_(Teuchos::null),
      scatra_matrixtype_(ssi_mono_algorithm.ScaTraField()->MatrixType()),
      scatra_manifold_matrixtype_(ssi_mono_algorithm.IsScaTraManifold()
                                      ? ssi_mono_algorithm.ScaTraManifold()->MatrixType()
                                      : LINALG::MatrixType::undefined),
      ssi_matrixtype_(ssi_mono_algorithm.MatrixType())
{
  std::vector<Teuchos::RCP<const Epetra_Map>> partial_maps(
      ssi_mono_algorithm.IsScaTraManifold() ? 3 : 2, Teuchos::null);
  Teuchos::RCP<const Epetra_Map> merged_map;

  partial_maps[GetProblemPosition(Subproblem::scalar_transport)] =
      Teuchos::rcp(new Epetra_Map(*ssi_mono_algorithm.ScaTraField()->DofRowMap()));
  partial_maps[GetProblemPosition(Subproblem::structure)] =
      Teuchos::rcp(new Epetra_Map(*ssi_mono_algorithm.StructureField()->DofRowMap()));
  if (ssi_mono_algorithm.IsScaTraManifold())
  {
    partial_maps[GetProblemPosition(Subproblem::manifold)] =
        Teuchos::rcp(new Epetra_Map(*ssi_mono_algorithm.ScaTraManifold()->DofRowMap()));
    auto temp_map = LINALG::MergeMap(partial_maps[0], partial_maps[1], false);
    merged_map = LINALG::MergeMap(temp_map, partial_maps[2], false);
  }
  else
    merged_map = LINALG::MergeMap(partial_maps[0], partial_maps[1], false);

  maps_sub_problems_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*merged_map, partial_maps));
  // check global map extractor
  maps_sub_problems_->CheckForValidMapExtractor();

  switch (ssi_matrixtype_)
  {
    case LINALG::MatrixType::block_field:
    {
      auto block_map_structure = Teuchos::rcp(new LINALG::MultiMapExtractor(
          *ssi_mono_algorithm.StructureField()->Discretization()->DofRowMap(),
          std::vector<Teuchos::RCP<const Epetra_Map>>(
              1, ssi_mono_algorithm.StructureField()->DofRowMap())));

      block_map_structure->CheckForValidMapExtractor();

      block_maps_sub_problems_.insert(std::make_pair(Subproblem::structure, block_map_structure));

      switch (scatra_matrixtype_)
      {
        case LINALG::MatrixType::sparse:
        {
          auto block_map_scatra = Teuchos::rcp(new LINALG::MultiMapExtractor(
              *ssi_mono_algorithm.ScaTraField()->Discretization()->DofRowMap(),
              std::vector<Teuchos::RCP<const Epetra_Map>>(
                  1, ssi_mono_algorithm.ScaTraField()->DofRowMap())));

          block_map_scatra->CheckForValidMapExtractor();

          block_maps_sub_problems_.insert(
              std::make_pair(Subproblem::scalar_transport, block_map_scatra));

          if (ssi_mono_algorithm.IsScaTraManifold())
          {
            auto block_map_scatra_manifold = Teuchos::rcp(new LINALG::MultiMapExtractor(
                *ssi_mono_algorithm.ScaTraManifold()->Discretization()->DofRowMap(),
                std::vector<Teuchos::RCP<const Epetra_Map>>(
                    1, ssi_mono_algorithm.ScaTraManifold()->DofRowMap())));

            block_map_scatra_manifold->CheckForValidMapExtractor();

            block_maps_sub_problems_.insert(
                std::make_pair(Subproblem::manifold, block_map_scatra_manifold));
          }
          break;
        }
        case LINALG::MatrixType::block_condition:
        case LINALG::MatrixType::block_condition_dof:
        {
          block_maps_sub_problems_.insert(std::make_pair(Subproblem::scalar_transport,
              Teuchos::rcpFromRef(ssi_mono_algorithm.ScaTraField()->BlockMaps())));

          if (ssi_mono_algorithm.IsScaTraManifold())
          {
            block_maps_sub_problems_.insert(std::make_pair(Subproblem::manifold,
                Teuchos::rcpFromRef(ssi_mono_algorithm.ScaTraManifold()->BlockMaps())));
          }
          else
          {
            block_maps_sub_problems_.insert(std::make_pair(Subproblem::manifold, Teuchos::null));
          }

          break;
        }
        default:
        {
          dserror("Invalid matrix type associated with scalar transport field!");
          break;
        }
      }

      CreateAndCheckBlockMapsSubProblems(ssi_mono_algorithm);

      break;
    }

    case LINALG::MatrixType::sparse:
    {
      block_maps_sub_problems_.insert(std::make_pair(Subproblem::structure, Teuchos::null));
      block_maps_sub_problems_.insert(std::make_pair(Subproblem::scalar_transport, Teuchos::null));
      block_maps_sub_problems_.insert(std::make_pair(Subproblem::manifold, Teuchos::null));
      break;
    }

    default:
    {
      dserror("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }

  map_system_matrix_ = maps_sub_problems_->FullMap();
}

/* ----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::UTILS::SSIMeshTyingMaps::SSIMeshTyingMaps(
    Teuchos::RCP<ADAPTER::Coupling> interface_coupling_adapter_structure,
    Teuchos::RCP<ADAPTER::Coupling> interface_coupling_adapter_structure_3domain_intersection,
    const bool meshtying_3_domain_intersection, Teuchos::RCP<DRT::Discretization> structure_dis)
    : maps_coup_struct_(Teuchos::null), maps_coup_struct_3_domain_intersection_(Teuchos::null)
{
  // setup map with interior and master side structural dofs
  Teuchos::RCP<Epetra_Map> map_structure_interior_master(Teuchos::null);
  if (meshtying_3_domain_intersection)
  {
    auto combined_slave_dof_maps =
        LINALG::MultiMapExtractor::MergeMaps({interface_coupling_adapter_structure->SlaveDofMap(),
            interface_coupling_adapter_structure_3domain_intersection->SlaveDofMap()});

    // set up map for interior and master-side structural degrees of freedom
    map_structure_interior_master =
        LINALG::SplitMap(*structure_dis->DofRowMap(), *combined_slave_dof_maps);
  }
  else
  {
    // set up map for interior and master-side structural degrees of freedom
    map_structure_interior_master = LINALG::SplitMap(
        *structure_dis->DofRowMap(), *interface_coupling_adapter_structure->SlaveDofMap());
  }

  // setup map extractor with interior, master, and slave side structural dofs
  {
    // set up structural map extractor holding interior and interface maps of degrees of freedom
    std::vector<Teuchos::RCP<const Epetra_Map>> maps_surf(0, Teuchos::null);
    maps_surf.emplace_back(LINALG::SplitMap(
        *map_structure_interior_master, *interface_coupling_adapter_structure->MasterDofMap()));
    maps_surf.emplace_back(interface_coupling_adapter_structure->SlaveDofMap());
    maps_surf.emplace_back(interface_coupling_adapter_structure->MasterDofMap());
    maps_coup_struct_ =
        Teuchos::rcp(new LINALG::MultiMapExtractor(*structure_dis->DofRowMap(), maps_surf));
    maps_coup_struct_->CheckForValidMapExtractor();

    if (meshtying_3_domain_intersection)
    {
      std::vector<Teuchos::RCP<const Epetra_Map>> maps_line(0, Teuchos::null);
      maps_line.emplace_back(LINALG::SplitMap(
          *interface_coupling_adapter_structure_3domain_intersection->MasterDofMap(),
          *interface_coupling_adapter_structure_3domain_intersection->MasterDofMap()));
      maps_line.emplace_back(
          interface_coupling_adapter_structure_3domain_intersection->SlaveDofMap());
      maps_line.emplace_back(
          interface_coupling_adapter_structure_3domain_intersection->MasterDofMap());
      maps_coup_struct_3_domain_intersection_ =
          Teuchos::rcp(new LINALG::MultiMapExtractor(*structure_dis->DofRowMap(), maps_line));
      maps_coup_struct_3_domain_intersection_->CheckForValidMapExtractor();
    }
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> SSI::UTILS::SSIMeshTyingMaps::MapStructureInterior() const
{
  return maps_coup_struct_->Map(0);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> SSI::UTILS::SSIMeshTyingMaps::MapStructureMaster() const
{
  return maps_coup_struct_->Map(2);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> SSI::UTILS::SSIMeshTyingMaps::MapStructureSlave() const
{
  return maps_coup_struct_->Map(1);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> SSI::UTILS::SSIMeshTyingMaps::MapStructureMaster3DomainIntersection()
    const
{
  return maps_coup_struct_3_domain_intersection_->Map(2);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> SSI::UTILS::SSIMeshTyingMaps::MapStructureSlave3DomainIntersection()
    const
{
  return maps_coup_struct_3_domain_intersection_->Map(1);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<int>> SSI::UTILS::SSIMaps::GetBlockPositions(Subproblem subproblem) const
{
  dsassert(ssi_matrixtype_ != LINALG::MatrixType::sparse, "Sparse matrices have just one block");

  Teuchos::RCP<std::vector<int>> block_position = Teuchos::rcp(new std::vector<int>(0));

  switch (subproblem)
  {
    case Subproblem::structure:
    {
      if (scatra_matrixtype_ == LINALG::MatrixType::sparse)
        block_position->emplace_back(1);
      else
        block_position->emplace_back(BlockMapScaTra()->NumMaps());
      break;
    }
    case Subproblem::scalar_transport:
    {
      if (scatra_matrixtype_ == LINALG::MatrixType::sparse)
        block_position->emplace_back(0);
      else
      {
        for (int i = 0; i < BlockMapScaTra()->NumMaps(); ++i) block_position->emplace_back(i);
      }
      break;
    }
    case Subproblem::manifold:
    {
      if (scatra_manifold_matrixtype_ == LINALG::MatrixType::sparse)
        block_position->emplace_back(2);
      else
      {
        auto scatra_manifold_num_block_maps = BlockMapScaTraManifold()->NumMaps();

        for (int i = 0; i < scatra_manifold_num_block_maps; ++i)
          block_position->emplace_back(BlockMapScaTra()->NumMaps() + 1 + i);
      }
      break;
    }
    default:
    {
      dserror("Unknown type of subproblem");
      break;
    }
  }

  return block_position;
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
int SSI::UTILS::SSIMaps::GetProblemPosition(Subproblem subproblem)
{
  int position = -1;

  switch (subproblem)
  {
    case Subproblem::structure:
    {
      position = 1;
      break;
    }
    case Subproblem::scalar_transport:
    {
      position = 0;
      break;
    }
    case Subproblem::manifold:
    {
      position = 2;
      break;
    }
    default:
    {
      dserror("Unknown type of sub problem");
      break;
    }
  }

  return position;
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::SSIMaps::CreateAndCheckBlockMapsSubProblems(const SSI::SSIMono& ssi_mono_algorithm)
{
  const int num_blocks_systemmatrix =
      BlockMapScaTra()->NumMaps() + BlockMapStructure()->NumMaps() +
      (ssi_mono_algorithm.IsScaTraManifold() ? BlockMapScaTraManifold()->NumMaps() : 0);

  std::vector<Teuchos::RCP<const Epetra_Map>> partial_maps_system_matrix(
      num_blocks_systemmatrix, Teuchos::null);

  for (int i = 0; i < BlockMapScaTra()->NumMaps(); ++i)

  {
    auto block_positions_scatra = GetBlockPositions(Subproblem::scalar_transport);
    partial_maps_system_matrix[block_positions_scatra->at(i)] = BlockMapScaTra()->Map(i);
  }

  partial_maps_system_matrix.at(GetBlockPositions(Subproblem::structure)->at(0)) =
      BlockMapStructure()->FullMap();

  if (ssi_mono_algorithm.IsScaTraManifold())
  {
    for (int i = 0; i < BlockMapScaTraManifold()->NumMaps(); ++i)
    {
      auto block_positions_manifold = GetBlockPositions(Subproblem::manifold);
      partial_maps_system_matrix[block_positions_manifold->at(i)] =
          BlockMapScaTraManifold()->Map(i);
    }
  }

  block_map_system_matrix_ = Teuchos::rcp(
      new LINALG::MultiMapExtractor(*maps_sub_problems_->FullMap(), partial_maps_system_matrix));

  block_map_system_matrix_->CheckForValidMapExtractor();
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::MultiMapExtractor> SSI::UTILS::SSIMaps::BlockMapScaTra() const
{
  return block_maps_sub_problems_.at(Subproblem::scalar_transport);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::MultiMapExtractor> SSI::UTILS::SSIMaps::BlockMapScaTraManifold() const
{
  return block_maps_sub_problems_.at(Subproblem::manifold);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::MultiMapExtractor> SSI::UTILS::SSIMaps::BlockMapStructure() const
{
  return block_maps_sub_problems_.at(Subproblem::structure);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> SSI::UTILS::SSIMaps::ScaTraDofRowMap() const
{
  return MapsSubProblems()->Map(SSIMaps::GetProblemPosition(Subproblem::scalar_transport));
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> SSI::UTILS::SSIMaps::ScaTraManifoldDofRowMap() const
{
  return MapsSubProblems()->Map(SSIMaps::GetProblemPosition(Subproblem::manifold));
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> SSI::UTILS::SSIMaps::StructureDofRowMap() const
{
  return MapsSubProblems()->Map(SSIMaps::GetProblemPosition(Subproblem::structure));
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
SSI::UTILS::SSISlaveSideConverter::SSISlaveSideConverter(
    Teuchos::RCP<ADAPTER::Coupling> icoup_structure,
    Teuchos::RCP<ADAPTER::Coupling> icoup_structure_3_domain_intersection,
    bool meshtying_3_domain_intersection)
    : icoup_structure_slave_converter_(
          Teuchos::rcp(new ADAPTER::CouplingSlaveConverter(*icoup_structure))),
      icoup_structure_slave_converter_3_domain_intersection_(
          meshtying_3_domain_intersection ? Teuchos::rcp(new ADAPTER::CouplingSlaveConverter(
                                                *icoup_structure_3_domain_intersection))
                                          : Teuchos::null)
{
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::CheckConsistencyOfSSIInterfaceContactCondition(
    const std::vector<DRT::Condition*>& conditionsToBeTested,
    Teuchos::RCP<DRT::Discretization>& structdis)
{
  // get conditions to check against
  std::vector<DRT::Condition*> s2ikinetics_conditions;
  structdis->GetCondition("S2IKinetics", s2ikinetics_conditions);
  std::vector<DRT::Condition*> contactconditions;
  structdis->GetCondition("Contact", contactconditions);

  // loop over all ssi conditions and check them
  for (const auto* conditionToBeTested : conditionsToBeTested)
  {
    std::vector<DRT::Condition*> InterfaceS2IConditions;
    std::vector<DRT::Condition*> InterfaceContactConditions;

    const int S2IKineticsID = conditionToBeTested->GetInt("S2IKineticsID");
    const int contactconditionID = conditionToBeTested->GetInt("ContactConditionID");

    if (S2IKineticsID != contactconditionID)
    {
      dserror(
          "For the 'SSIInterfaceContact' condition we have to demand, that the 'S2ICouplingID' and "
          "the 'ContactConditionID' have the same value as the contact strategy factory relies on "
          "this to set the scatra interface parameters correctly.");
    }

    // loop over all scatra-scatra interface conditions and add them to the vector, if IDs match
    for (auto* s2ikinetics_cond : s2ikinetics_conditions)
    {
      if (s2ikinetics_cond->GetInt("ConditionID") != S2IKineticsID) continue;

      InterfaceS2IConditions.push_back(s2ikinetics_cond);
    }

    // loop over all contact conditions and add them to the vector, if IDs match
    for (auto* contactcondition : contactconditions)
    {
      if (contactcondition->GetInt("Interface ID") != contactconditionID) continue;

      InterfaceContactConditions.push_back(contactcondition);
    }

    if (InterfaceContactConditions.empty())
      dserror("Did not find 'Contact' condition as defined in 'SSIInterfaceContact' condition!");

    if (InterfaceS2IConditions.empty())
      dserror(
          "Did not find 'S2ICoupling' condition as defined in 'SSIInterfaceContact' condition!");

    // now get the nodes
    std::vector<int> InterfaceS2INodes;
    std::vector<int> InterfaceContactNodes;
    DRT::UTILS::FindConditionedNodes(*structdis, InterfaceS2IConditions, InterfaceS2INodes);
    DRT::UTILS::FindConditionedNodes(*structdis, InterfaceContactConditions, InterfaceContactNodes);

    // and compare whether same nodes are defined
    for (const auto InterfaceS2INode : InterfaceS2INodes)
    {
      bool foundit(false);

      for (const auto InterfaceContactNode : InterfaceContactNodes)
      {
        if (InterfaceS2INode == InterfaceContactNode) foundit = true;
      }

      if (!foundit)
      {
        dserror(
            "The following conditions are defined on a non-consistent set of nodes!\n"
            "The Conditions defined from 'SSIInterfaceContact' condition with the condition-ID: "
            "%i:\n"
            "'S2ICoupling' conditions with ID: %i;\n"
            "'Contact' conditions with ID: %i;\n"
            "The last two conditions are NOT defined on the same node-sets which is not "
            "reasonable. Check your Input-File!",
            conditionToBeTested->GetInt("ConditionID"), S2IKineticsID, contactconditionID);
      }
    }
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
SSI::UTILS::SSIStructureMeshTying::SSIStructureMeshTying(const std::string& conditionname_coupling,
    const std::string& conditionname_3_domain_intersection,
    const bool meshtying_3_domain_intersection, Teuchos::RCP<DRT::Discretization> struct_dis)
    : icoup_structure_(Teuchos::null),
      icoup_structure_3_domain_intersection_(Teuchos::null),
      meshtying_3_domain_intersection_(meshtying_3_domain_intersection),
      slave_side_converter_(Teuchos::null),
      ssi_meshtyingmaps_(Teuchos::null),
      meshtying_handler_()
{
  icoup_structure_ = SSI::UTILS::SetupInterfaceCouplingAdapterStructure(struct_dis,
      meshtying_3_domain_intersection, conditionname_coupling, conditionname_3_domain_intersection);

  if (meshtying_3_domain_intersection)
  {
    icoup_structure_3_domain_intersection_ =
        SSI::UTILS::SetupInterfaceCouplingAdapterStructure3DomainIntersection(
            struct_dis, conditionname_3_domain_intersection);
  }

  slave_side_converter_ = Teuchos::rcp(new SSI::UTILS::SSISlaveSideConverter(
      icoup_structure_, icoup_structure_3_domain_intersection_, meshtying_3_domain_intersection));

  ssi_meshtyingmaps_ = Teuchos::rcp(new SSI::UTILS::SSIMeshTyingMaps(icoup_structure_,
      icoup_structure_3_domain_intersection_, meshtying_3_domain_intersection, struct_dis));
}
