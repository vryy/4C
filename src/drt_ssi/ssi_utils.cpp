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

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../linalg/linalg_utils_sparse_algebra_create.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                        AN,JH 09/2014 */
/* Function for checking that the different time steps are a
 multiplicative of each other                                           */

int SSI::UTILS::CheckTimeStepping(double dt1, double dt2)
{
  double workdt1 = std::min(dt1, dt2);
  double workdt2 = std::max(dt1, dt2);
  double t1 = 0.0;
  int i = 0;

  while (true)
  {
    i++;
    t1 = i * workdt1;

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
Teuchos::RCP<ADAPTER::Coupling> SSI::UTILS::SetupInterfaceCouplingAdapterStructure(
    Teuchos::RCP<DRT::Discretization> structdis, bool meshtying_3_domain_intersection,
    const std::string& conditionname_coupling,
    const std::string& conditionname_3_domain_intersection)
{
  auto interfacecouplingadapter = Teuchos::rcp(new ADAPTER::Coupling());

  // extract ssi coupling conditions from discretization and assign to slave/master side maps
  std::vector<DRT::Condition*> conditions(0, nullptr);
  structdis->GetCondition(conditionname_coupling, conditions);
  std::map<const int, DRT::Condition* const> slaveconditions;
  std::map<const int, DRT::Condition* const> masterconditions;
  for (auto* condition : conditions)
  {
    const int condid = condition->GetInt("ConditionID");
    if (condid < 0)
      dserror("Invalid condition ID %i for SSI Interface Meshtying Condition!", condid);

    if (*condition->Get<std::string>("Side") == "Slave")
    {
      if (slaveconditions.find(condid) == slaveconditions.end())
        slaveconditions.insert(std::pair<const int, DRT::Condition* const>(condid, condition));
      else
      {
        dserror(
            "Cannot have multiple slave-side ssi coupling conditions with the same ID %i!", condid);
      }
    }
    else if (*condition->Get<std::string>("Side") == "Master")
    {
      if (masterconditions.find(condid) == masterconditions.end())
        masterconditions.insert(std::pair<const int, DRT::Condition* const>(condid, condition));
      else
      {
        dserror("Cannot have multiple master-side ssi coupling conditions with the same ID %i!",
            condid);
      }
    }
    else
      dserror("Invalid ssi coupling condition!");
  }

  // Setup coupling adapter from ssi conditions
  if (meshtying_3_domain_intersection)
  {
    // vectors for global ids of slave and master interface nodes for each condition
    std::vector<std::vector<int>> islavenodegidvec_cond;
    std::vector<std::vector<int>> imasternodegidvec_cond;

    // loop over slave conditions and build vector of nodes for slave and master condition with same
    // ID
    for (auto& slavecondition : slaveconditions)
    {
      std::vector<int> islavenodegidvec;
      std::vector<int> imasternodegidvec;

      // Build GID vector of nodes on this proc, sort, and remove duplicates
      BuildNodeGidVector(structdis, slavecondition.second->Nodes(), islavenodegidvec);
      SortAndEraseNodeGidVector(islavenodegidvec);

      auto mastercondition = masterconditions.find(slavecondition.first);
      if (mastercondition != masterconditions.end())
      {
        BuildNodeGidVector(structdis, mastercondition->second->Nodes(), imasternodegidvec);
      }
      SortAndEraseNodeGidVector(imasternodegidvec);

      // remove nodes from slave side line condition to avoid non-unique slave-master relation
      std::vector<DRT::Condition*> conditions_3_domain_intersection(0, nullptr);
      structdis->GetCondition(
          conditionname_3_domain_intersection, conditions_3_domain_intersection);

      for (const auto& condition_3_domain_intersection : conditions_3_domain_intersection)
      {
        if (slavecondition.first != condition_3_domain_intersection->GetInt("SSISurfMasterID") and
            condition_3_domain_intersection->GType() == DRT::Condition::Line)
        {
          RemoveNodesGidVector(
              structdis, condition_3_domain_intersection->Nodes(), islavenodegidvec);
          RemoveNodesGidVector(
              structdis, condition_3_domain_intersection->Nodes(), imasternodegidvec);
        }
      }

      imasternodegidvec_cond.push_back(imasternodegidvec);
      islavenodegidvec_cond.push_back(islavenodegidvec);
    }

    interfacecouplingadapter->SetupCoupling(*structdis, *structdis, imasternodegidvec_cond,
        islavenodegidvec_cond, DRT::Problem::Instance()->NDim(), false, 1.e-8);
  }
  else
  {
    // vectors for global ids of slave and master interface nodes for all conditions
    std::vector<int> inodegidvec_master;
    std::vector<int> inodegidvec_slave;

    // Build GID vector of nodes on this proc, sort, and remove duplicates
    for (const auto& slavecondition : slaveconditions)
      BuildNodeGidVector(structdis, slavecondition.second->Nodes(), inodegidvec_slave);
    SortAndEraseNodeGidVector(inodegidvec_slave);

    for (const auto& mastercondition : masterconditions)
      BuildNodeGidVector(structdis, mastercondition.second->Nodes(), inodegidvec_master);
    SortAndEraseNodeGidVector(inodegidvec_master);


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
  std::vector<int> ilinenodegidvec_master;
  std::vector<int> ilinenodegidvec_slave;

  std::vector<DRT::Condition*> conditionsssi(0, nullptr);
  structdis->GetCondition(conditionname_3_domain_intersection, conditionsssi);

  // Build GID vector of nodes on this proc, sort, and remove duplicates
  for (const auto& condition : conditionsssi)
  {
    if (condition->GType() != DRT::Condition::Line) dserror("only line allowed");
    const std::vector<int>* const inodegids = condition->Nodes();

    if (*condition->Get<std::string>("Side") == "Slave")
      BuildNodeGidVector(structdis, inodegids, ilinenodegidvec_slave);
    else
      BuildNodeGidVector(structdis, inodegids, ilinenodegidvec_master);
  }

  SortAndEraseNodeGidVector(ilinenodegidvec_master);
  SortAndEraseNodeGidVector(ilinenodegidvec_slave);

  auto coupling_line_structure = Teuchos::rcp(new ADAPTER::Coupling());
  coupling_line_structure->SetupCoupling(*structdis, *structdis, ilinenodegidvec_master,
      ilinenodegidvec_slave, DRT::Problem::Instance()->NDim(), true, 1.0e-8);

  return coupling_line_structure;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SortAndEraseNodeGidVector(std::vector<int>& nodegidvec)
{
  std::sort(nodegidvec.begin(), nodegidvec.end());
  nodegidvec.erase(unique(nodegidvec.begin(), nodegidvec.end()), nodegidvec.end());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::BuildNodeGidVector(Teuchos::RCP<DRT::Discretization> dis,
    const std::vector<int>* nodegids, std::vector<int>& nodegidvec)
{
  for (const int nodegid : *nodegids)
  {
    // insert global ID of current node into associated vector only if node is owned by current
    // processor need to make sure that node is stored on current processor, otherwise cannot
    // resolve "->Owner()"
    if (dis->HaveGlobalNode(nodegid) and dis->gNode(nodegid)->Owner() == dis->Comm().MyPID())
      nodegidvec.push_back(nodegid);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::RemoveNodesGidVector(Teuchos::RCP<DRT::Discretization> dis,
    const std::vector<int>* nodegidstoremove, std::vector<int>& nodegidvec)
{
  for (int h = static_cast<int>(nodegidvec.size()) - 1; h >= 0; --h)
  {
    for (int nodegid : *nodegidstoremove)
    {
      if (nodegidvec.at(h) == nodegid)
      {
        if (!dis->HaveGlobalNode(nodegid) or dis->gNode(nodegid)->Owner() != dis->Comm().MyPID())
          continue;

        nodegidvec.erase(nodegidvec.begin() + h);
        break;
      }
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::UTILS::SSIMatrices::SSIMatrices(
    const SSI::SSIMono& ssi_mono_algorithm, Teuchos::RCP<Epetra_Map> interface_map_scatra)
{
  // first: build interface maps
  Teuchos::RCP<LINALG::MultiMapExtractor> blockmapscatrainterface(Teuchos::null);
  switch (ssi_mono_algorithm.ScaTraField()->MatrixType())
  {
    // one single main-diagonal matrix block associated with scalar transport field
    case LINALG::MatrixType::sparse:
    {
      if (ssi_mono_algorithm.SSIInterfaceMeshtying())
      {
        blockmapscatrainterface = Teuchos::rcp(new LINALG::MultiMapExtractor(*interface_map_scatra,
            std::vector<Teuchos::RCP<const Epetra_Map>>(1, interface_map_scatra)));
        blockmapscatrainterface->CheckForValidMapExtractor();
      }
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
        blockmapscatrainterface = Teuchos::rcp(
            new LINALG::MultiMapExtractor(*interface_map_scatra, partial_blockmapscatrainterface));
        blockmapscatrainterface->CheckForValidMapExtractor();
      }
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // second: initialze system matrix
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

  // third: initialize scatra-structure block and structure-scatra block of global system matrix
  switch (ssi_mono_algorithm.ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      // initialize scatra-structure block
      scatrastructuredomain_ =
          SetupBlockMatrix(Teuchos::rcpFromRef(ssi_mono_algorithm.ScaTraField()->BlockMaps()),
              ssi_mono_algorithm.MapStructure());
      structurescatradomain_ = SetupBlockMatrix(ssi_mono_algorithm.MapStructure(),
          Teuchos::rcpFromRef(ssi_mono_algorithm.ScaTraField()->BlockMaps()));

      if (ssi_mono_algorithm.SSIInterfaceMeshtying())
      {
        scatrastructureinterface_ =
            SetupBlockMatrix(blockmapscatrainterface, ssi_mono_algorithm.MapStructure());
      }

      break;
    }

    case LINALG::MatrixType::sparse:
    {
      scatrastructuredomain_ = SetupSparseMatrix(ssi_mono_algorithm.ScaTraField()->DofRowMap());
      structurescatradomain_ = SetupSparseMatrix(ssi_mono_algorithm.StructureField()->DofRowMap());
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