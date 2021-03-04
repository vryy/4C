/*----------------------------------------------------------------------*/
/*! \file
 \brief Utility methods for SSI

 \level 1


 *------------------------------------------------------------------------------------------------*/

#include "ssi_utils.H"
#include "ssi_monolithic.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/adapter_coupling.H"

#include "../drt_inpar/inpar_s2i.H"

#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_utils_gid_vector.H"
#include "../drt_lib/drt_utils_vector.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

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
void SSI::UTILS::CheckConsistencyWithS2IMeshtyingCondition(
    const std::vector<DRT::Condition*>& conditionsToBeTested,
    Teuchos::RCP<DRT::Discretization>& structdis)
{
  std::vector<DRT::Condition*> s2iconditions;
  structdis->GetCondition("S2ICoupling", s2iconditions);

  // loop over all conditions to be tested and check for a consistent initialization of the s2i
  // conditions
  for (const auto& conditionToBeTested : conditionsToBeTested)
  {
    if (conditionToBeTested->GType() != DRT::Condition::Surface) continue;
    bool matchingconditions(false);
    bool isslave(true);
    const int s2icouplingid = conditionToBeTested->GetInt("S2ICouplingID");
    const auto* side = conditionToBeTested->Get<std::string>("Side");
    // check interface side
    if (*side == "Slave")
      isslave = true;
    else if (*side == "Master")
      isslave = false;
    else
    {
      dserror(
          "Interface side of tested condition not recognized, has to be either 'Slave' or "
          "'Master'");
    }

    // loop over all s2i conditions to find the one that is matching the current ssi condition
    for (const auto& s2icondition : s2iconditions)
    {
      const int s2iconditionid = s2icondition->GetInt("ConditionID");
      // only do further checks if Ids match
      if (s2icouplingid != s2iconditionid) continue;

      // check the interface side
      switch (s2icondition->GetInt("interface side"))
      {
        case INPAR::S2I::side_slave:
        {
          if (isslave)
            matchingconditions = DRT::UTILS::HaveSameNodes(conditionToBeTested, s2icondition);

          break;
        }
        case INPAR::S2I::side_master:
        {
          if (!isslave)
            matchingconditions = DRT::UTILS::HaveSameNodes(conditionToBeTested, s2icondition);

          break;
        }
        default:
        {
          dserror("interface side of 'S2iCondition' has to be either 'Slave' or 'Master'");
          break;
        }
      }
    }

    if (!matchingconditions)
    {
      dserror(
          "Did not find 'S2ICoupling' condition with ID: %i and interface side: %s as defined in "
          "the condition to be tested",
          s2icouplingid, side->c_str());
    }
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
    if (*condition->Get<std::string>("Side") == "Slave")
      slave_conditions.insert(pair);
    else if (*condition->Get<std::string>("Side") == "Master")
      master_conditions.insert(pair);
    else
      dserror("Coupling side must either be slave or master");
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
        condition_pairs.emplace_back(std::pair<DRT::Condition* const, DRT::Condition* const>(
            slave_condition.second, master_condition.second));
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
    if (!DRT::UTILS::HaveSameNodes(
            condition_3_domain_pairs[0].second, condition_3_domain_pairs[1].second))
      dserror("All master conditions with 3-domain-meshtying must be the same");

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
  if (!DRT::UTILS::HaveSameNodes(
          condition_3_domain_pairs[0].second, condition_3_domain_pairs[1].second))
    dserror("All master conditions with 3-domain-meshtying must be the same");

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
    const Teuchos::ParameterList& sublist_manifold_params, const Epetra_Comm& comm)
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

  scatra_manifold_params->set<std::string>("OUTPUTSCALARS", "none");
  scatra_manifold_params->set<std::string>("ADAPTIVE_TIMESTEPPING", "No");

  return *scatra_manifold_params;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::MultiMapExtractor> SSI::UTILS::CreateManifoldMultiMapExtractor(
    Teuchos::RCP<DRT::Discretization> dis)
{
  std::vector<DRT::Condition*> conditions(0, nullptr);
  dis->GetCondition("SSISurfaceManifold", conditions);
  if (conditions.empty()) dserror("Condition SSISurfaceManifold not found");

  // Build GID vector of nodes on this proc, sort, and remove duplicates
  std::vector<int> condition_node_vec;
  for (auto* condition : conditions)
  {
    if (condition->GetInt("coupling id") != 1)
      dserror("'coupling id' in 'SSISurfaceManifold' must be set to 1");

    DRT::UTILS::AddOwnedNodeGIDVector(dis, *condition->Nodes(), condition_node_vec);
  }
  DRT::UTILS::SortAndRemoveDuplicateVectorElements(condition_node_vec);

  // Build GID vector of dofs
  std::vector<int> condition_dof_vec;
  for (int condition_node : condition_node_vec)
  {
    const auto* curr_node = dis->gNode(condition_node);
    for (int j = 0; j < dis->NumDof(0, curr_node); ++j)
      condition_dof_vec.emplace_back(dis->Dof(0, curr_node, j));
  }
  // maps of conditioned dofs and other dofs
  const auto condition_dof_map = Teuchos::rcp(
      new const Epetra_Map(-1, condition_dof_vec.size(), &condition_dof_vec[0], 0, dis->Comm()));
  const auto non_condition_dof_map = LINALG::SplitMap(*dis->DofRowMap(), *condition_dof_map);

  return Teuchos::rcp(
      new LINALG::MapExtractor(*dis->DofRowMap(), non_condition_dof_map, condition_dof_map));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::UTILS::SSIMatrices::SSIMatrices(
    const SSI::SSIMono& ssi_mono_algorithm, Teuchos::RCP<const Epetra_Map> interface_map_scatra)
{
  InitializeSystemMatrix(ssi_mono_algorithm);

  InitializeOffDiagMatrices(ssi_mono_algorithm, interface_map_scatra);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::MultiMapExtractor> SSI::UTILS::SSIMatrices::GetScaTraInterfaceBlockMap(
    const SSI::SSIMono& ssi_mono_algorithm,
    Teuchos::RCP<const Epetra_Map> interface_map_scatra) const
{
  Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_scatra_interface(Teuchos::null);

  switch (ssi_mono_algorithm.ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::sparse:
    {
      // do nothing
      break;
    }
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      if (ssi_mono_algorithm.SSIInterfaceMeshtying())
      {
        const int maps_systemmatrix_scatra = ssi_mono_algorithm.MapsScatra()->NumMaps();

        // build block map for scatra interface by merging slave and master side for each block
        std::vector<Teuchos::RCP<const Epetra_Map>> partial_blockmapscatrainterface(
            maps_systemmatrix_scatra, Teuchos::null);
        for (int iblockmap = 0; iblockmap < maps_systemmatrix_scatra; ++iblockmap)
        {
          partial_blockmapscatrainterface.at(iblockmap) = LINALG::MultiMapExtractor::MergeMaps(
              {ssi_mono_algorithm.MeshtyingStrategyS2I()->BlockMapsSlave().Map(iblockmap),
                  ssi_mono_algorithm.MeshtyingStrategyS2I()->BlockMapsMaster().Map(iblockmap)});
        }
        block_map_scatra_interface = Teuchos::rcp(
            new LINALG::MultiMapExtractor(*interface_map_scatra, partial_blockmapscatrainterface));
        block_map_scatra_interface->CheckForValidMapExtractor();
      }
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  return block_map_scatra_interface;
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::InitializeOffDiagMatrices(
    const SSI::SSIMono& ssi_mono_algorithm, Teuchos::RCP<const Epetra_Map> interface_map_scatra)
{
  switch (ssi_mono_algorithm.ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      scatrastructuredomain_ =
          SetupBlockMatrix(Teuchos::rcpFromRef(ssi_mono_algorithm.ScaTraField()->BlockMaps()),
              ssi_mono_algorithm.MapStructure());

      structurescatradomain_ = SetupBlockMatrix(ssi_mono_algorithm.MapStructure(),
          Teuchos::rcpFromRef(ssi_mono_algorithm.ScaTraField()->BlockMaps()));

      if (ssi_mono_algorithm.IsScaTraManifold())
      {
        // structure dofs on manifold discretization
        const auto map_structure_manifold = Teuchos::rcp(new LINALG::MultiMapExtractor(
            *ssi_mono_algorithm.MapStructureOnScaTraManifold()->Map(0),
            std::vector<Teuchos::RCP<const Epetra_Map>>(
                1, ssi_mono_algorithm.MapStructureOnScaTraManifold()->Map(0))));

        scatramanifoldstructuredomain_ =
            SetupBlockMatrix(Teuchos::rcpFromRef(ssi_mono_algorithm.ScaTraManifold()->BlockMaps()),
                map_structure_manifold);
      }

      if (ssi_mono_algorithm.SSIInterfaceMeshtying())
      {
        auto block_map_scatra_interface =
            GetScaTraInterfaceBlockMap(ssi_mono_algorithm, interface_map_scatra);

        scatrastructureinterface_ =
            SetupBlockMatrix(block_map_scatra_interface, ssi_mono_algorithm.MapStructure());
      }

      break;
    }

    case LINALG::MatrixType::sparse:
    {
      scatrastructuredomain_ = SetupSparseMatrix(ssi_mono_algorithm.ScaTraField()->DofRowMap());
      structurescatradomain_ = SetupSparseMatrix(ssi_mono_algorithm.StructureField()->DofRowMap());

      if (ssi_mono_algorithm.IsScaTraManifold())
      {
        scatramanifoldstructuredomain_ =
            SetupSparseMatrix(ssi_mono_algorithm.ScaTraManifold()->DofRowMap());
      }

      if (ssi_mono_algorithm.SSIInterfaceMeshtying())
        scatrastructureinterface_ = SetupSparseMatrix(interface_map_scatra);

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
SSI::UTILS::SSIVectors::SSIVectors(const SSI::SSIMono& ssi_mono_algorithm)
    : increment_(Teuchos::null), residual_(Teuchos::null)
{
  increment_ = LINALG::CreateVector(*(ssi_mono_algorithm.DofRowMap()), true);
  residual_ = LINALG::CreateVector(*(ssi_mono_algorithm.DofRowMap()), true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SSIVectors::ClearIncrement() { increment_->PutScalar(0.0); }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SSIVectors::ClearResiduals() { residual_->PutScalar(0.0); }

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::InitializeSystemMatrix(const SSI::SSIMono& ssi_mono_algorithm)
{
  switch (ssi_mono_algorithm.MatrixType())
  {
    case LINALG::MatrixType::block_field:
    {
      systemmatrix_ = SetupBlockMatrix(
          ssi_mono_algorithm.MapsSystemMatrix(), ssi_mono_algorithm.MapsSystemMatrix());
      break;
    }

    case LINALG::MatrixType::sparse:
    {
      systemmatrix_ = SetupSparseMatrix(ssi_mono_algorithm.DofRowMap());
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
  std::vector<DRT::Condition*> s2iconditions;
  structdis->GetCondition("S2ICoupling", s2iconditions);
  std::vector<DRT::Condition*> contactconditions;
  structdis->GetCondition("Contact", contactconditions);

  // loop over all ssi conditions and check them
  for (const auto* conditionToBeTested : conditionsToBeTested)
  {
    std::vector<DRT::Condition*> InterfaceS2IConditions;
    std::vector<DRT::Condition*> InterfaceContactConditions;

    const int s2iconditionID = conditionToBeTested->GetInt("S2ICouplingID");
    const int contactconditionID = conditionToBeTested->GetInt("ContactConditionID");

    if (s2iconditionID != contactconditionID)
    {
      dserror(
          "For the 'SSIInterfaceContact' condition we have to demand, that the 'S2ICouplingID' and "
          "the 'ContactConditionID' have the same value as the contact strategy factory relies on "
          "this to set the scatra interface parameters correctly.");
    }

    // loop over all scatra-scatra interface conditions and add them to the vector, if IDs match
    for (auto* s2icondition : s2iconditions)
    {
      if (s2icondition->GetInt("ConditionID") != s2iconditionID) continue;

      InterfaceS2IConditions.push_back(s2icondition);
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
            conditionToBeTested->GetInt("ConditionID"), s2iconditionID, contactconditionID);
      }
    }
  }
}
