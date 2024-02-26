/*----------------------------------------------------------------------*/
/*! \file
\level 1


\brief BACI implementation of main class to control all solid meshtying
*/
/*----------------------------------------------------------------------*/

#include "baci_contact_meshtying_manager.hpp"

#include "baci_contact_meshtying_defines.hpp"
#include "baci_contact_meshtying_lagrange_strategy.hpp"
#include "baci_contact_meshtying_penalty_strategy.hpp"
#include "baci_contact_meshtying_poro_lagrange_strategy.hpp"
#include "baci_global_data.hpp"
#include "baci_inpar_contact.hpp"
#include "baci_inpar_mortar.hpp"
#include "baci_io.hpp"
#include "baci_linalg_utils_sparse_algebra_manipulation.hpp"
#include "baci_mortar_defines.hpp"
#include "baci_mortar_element.hpp"
#include "baci_mortar_interface.hpp"
#include "baci_mortar_node.hpp"
#include "baci_mortar_utils.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 03/08|
 *----------------------------------------------------------------------*/
CONTACT::MtManager::MtManager(DRT::Discretization& discret, double alphaf) : MORTAR::ManagerBase()
{
  // overwrite base class communicator
  comm_ = Teuchos::rcp(discret.Comm().Clone());

  // create some local variables (later to be stored in strategy)
  const int spatialDim = GLOBAL::Problem::Instance()->NDim();
  if (spatialDim != 2 && spatialDim != 3) dserror("Meshtying problem must be 2D or 3D.");

  std::vector<Teuchos::RCP<MORTAR::Interface>> interfaces;
  Teuchos::ParameterList mtparams;

  // read and check meshtying input parameters
  if (Comm().MyPID() == 0)
    std::cout << "Checking meshtying input parameters..........." << std::endl;

  ReadAndCheckInput(mtparams, discret);
  if (Comm().MyPID() == 0) std::cout << "done!" << std::endl;

  // check for FillComplete of discretization
  if (!discret.Filled()) dserror("Discretization of underlying problem is not fillcomplete.");

  // let's check for meshtying boundary conditions in discret
  // and detect groups of matching conditions
  // for each group, create a contact interface and store it
  if (Comm().MyPID() == 0)
    std::cout << "Building meshtying interface(s)..............." << std::endl;

  std::vector<DRT::Condition*> contactconditions(0);
  discret.GetCondition("Mortar", contactconditions);

  // there must be more than one meshtying condition
  if (contactconditions.size() < 2) dserror("Not enough contact conditions in discretization");

  // find all pairs of matching meshtying conditions
  // there is a maximum of (conditions / 2) groups
  std::vector<int> foundgroups(0);
  int numgroupsfound = 0;

  // get nurbs information
  const bool nurbs = mtparams.get<bool>("NURBS");

  // maximum dof number in discretization
  // later we want to create NEW Lagrange multiplier degrees of
  // freedom, which of course must not overlap with displacement dofs
  const int maxdof = discret.DofRowMap()->MaxAllGID();

  for (unsigned i = 0; i < contactconditions.size(); ++i)
  {
    // initialize vector for current group of conditions and temp condition
    std::vector<DRT::Condition*> currentgroup(0);
    DRT::Condition* tempcond = nullptr;

    // try to build meshtying group around this condition
    currentgroup.push_back(contactconditions[i]);
    const auto groupid1 = *currentgroup[0]->Get<int>("Interface ID");
    bool foundit = false;

    for (unsigned j = 0; j < contactconditions.size(); ++j)
    {
      if (j == i) continue;  // do not detect contactconditions[i] again
      tempcond = contactconditions[j];
      const auto groupid2 = *tempcond->Get<int>("Interface ID");
      if (groupid1 != groupid2) continue;  // not in the group
      foundit = true;                      // found a group entry
      currentgroup.push_back(tempcond);    // store it in currentgroup
    }

    // now we should have found a group of conds
    if (!foundit) dserror("Cannot find matching contact condition for id %d", groupid1);

    // see whether we found this group before
    bool foundbefore = false;
    for (int j = 0; j < numgroupsfound; ++j)
    {
      if (groupid1 == foundgroups[j])
      {
        foundbefore = true;
        break;
      }
    }

    // if we have processed this group before, do nothing
    if (foundbefore) continue;

    // we have not found this group before, process it
    foundgroups.push_back(groupid1);
    ++numgroupsfound;

    // find out which sides are Master and Slave
    bool hasslave = false;
    bool hasmaster = false;
    std::vector<const std::string*> sides(currentgroup.size());
    std::vector<bool> isslave(currentgroup.size());

    for (unsigned j = 0; j < sides.size(); ++j)
    {
      sides[j] = currentgroup[j]->Get<std::string>("Side");
      if (*sides[j] == "Slave")
      {
        hasslave = true;
        isslave[j] = true;
      }
      else if (*sides[j] == "Master")
      {
        hasmaster = true;
        isslave[j] = false;
      }
      else
      {
        dserror("MtManager: Unknown mortar side qualifier!");
      }
    }

    if (!hasslave) dserror("Slave side missing in contact condition group!");
    if (!hasmaster) dserror("Master side missing in contact condition group!");

    // find out which sides are initialized as Active
    std::vector<const std::string*> active(currentgroup.size());
    std::vector<bool> isactive(currentgroup.size());

    for (unsigned j = 0; j < sides.size(); ++j)
    {
      active[j] = currentgroup[j]->Get<std::string>("Initialization");
      if (*sides[j] == "Slave")
      {
        // slave sides must be initialized as "Active"
        if (*active[j] == "Active")
          isactive[j] = true;
        else if (*active[j] == "Inactive")
          dserror(" Slave side must be active for meshtying!");
        else
          dserror("Unknown initialization qualifier for slave side of mortar meshtying interface!");
      }
      else if (*sides[j] == "Master")
      {
        // master sides must NOT be initialized as "Active" as this makes no sense
        if (*active[j] == "Active")
          dserror("Master side cannot be active!");
        else if (*active[j] == "Inactive")
          isactive[j] = false;
        else
          dserror(
              "Unknown initialization qualifier for master side of mortar meshtying interface!");
      }
      else
      {
        dserror("MtManager: Unknown contact side qualifier!");
      }
    }

    // create an empty meshtying interface and store it in this Manager
    interfaces.push_back(MORTAR::Interface::Create(groupid1, Comm(), spatialDim, mtparams));

    // get it again
    Teuchos::RCP<MORTAR::Interface> interface = interfaces.back();

    // note that the nodal IDs are unique because they come from
    // one global problem discretization containing all nodes of the
    // contact interface
    // We rely on this fact, therefore it is not possible to
    // do meshtying between two distinct discretizations here

    //-------------------------------------------------- process nodes
    for (unsigned j = 0; j < currentgroup.size(); ++j)
    {
      // get all nodes and add them
      const std::vector<int>* nodeids = currentgroup[j]->GetNodes();
      if (!nodeids) dserror("Condition does not have Node Ids");
      for (unsigned k = 0; k < nodeids->size(); ++k)
      {
        int gid = (*nodeids)[k];
        // do only nodes that I have in my discretization
        if (!discret.NodeColMap()->MyGID(gid)) continue;
        DRT::Node* node = discret.gNode(gid);
        if (!node) dserror("Cannot find node with gid %", gid);

        // create Node object
        Teuchos::RCP<MORTAR::Node> mtnode = Teuchos::rcp(new MORTAR::Node(
            node->Id(), node->X(), node->Owner(), discret.Dof(0, node), isslave[j]));
        //-------------------
        // get nurbs weight!
        if (nurbs) MORTAR::UTILS::PrepareNURBSNode(node, mtnode);

        // get edge and corner information:
        std::vector<DRT::Condition*> contactcornercond(0);
        discret.GetCondition("mrtrcorner", contactcornercond);
        for (unsigned j = 0; j < contactcornercond.size(); j++)
        {
          if (contactcornercond.at(j)->ContainsNode(node->Id()))
          {
            mtnode->SetOnCorner() = true;
          }
        }
        std::vector<DRT::Condition*> contactedgecond(0);
        discret.GetCondition("mrtredge", contactedgecond);
        for (unsigned j = 0; j < contactedgecond.size(); j++)
        {
          if (contactedgecond.at(j)->ContainsNode(node->Id()))
          {
            mtnode->SetOnEdge() = true;
          }
        }

        // Check, if this node (and, in case, which dofs) are in the contact symmetry condition
        std::vector<DRT::Condition*> contactSymconditions(0);
        discret.GetCondition("mrtrsym", contactSymconditions);

        for (unsigned j = 0; j < contactSymconditions.size(); j++)
          if (contactSymconditions.at(j)->ContainsNode(node->Id()))
          {
            const std::vector<int>* onoff =
                contactSymconditions.at(j)->Get<std::vector<int>>("onoff");
            for (unsigned k = 0; k < onoff->size(); k++)
              if (onoff->at(k) == 1) mtnode->DbcDofs()[k] = true;
          }

        // note that we do not have to worry about double entries
        // as the AddNode function can deal with this case!
        interface->AddMortarNode(mtnode);
      }
    }

    //----------------------------------------------- process elements
    int ggsize = 0;
    for (unsigned j = 0; j < currentgroup.size(); ++j)
    {
      // get elements from condition j of current group
      std::map<int, Teuchos::RCP<DRT::Element>>& currele = currentgroup[j]->Geometry();

      // elements in a boundary condition have a unique id
      // but ids are not unique among 2 distinct conditions
      // due to the way elements in conditions are build.
      // We therefore have to give the second, third,... set of elements
      // different ids. ids do not have to be continous, we just add a large
      // enough number ggsize to all elements of cond2, cond3,... so they are
      // different from those in cond1!!!
      // note that elements in ele1/ele2 already are in column (overlapping) map
      int lsize = static_cast<int>(currele.size());
      int gsize = 0;
      Comm().SumAll(&lsize, &gsize, 1);


      for (const auto& element : currele)
      {
        Teuchos::RCP<DRT::Element> ele = element.second;
        Teuchos::RCP<MORTAR::Element> mtele = Teuchos::rcp(new MORTAR::Element(ele->Id() + ggsize,
            ele->Owner(), ele->Shape(), ele->NumNode(), ele->NodeIds(), isslave[j], nurbs));
        //------------------------------------------------------------------
        // get knotvector, normal factor and zero-size information for nurbs
        if (nurbs) MORTAR::UTILS::PrepareNURBSElement(discret, ele, mtele, spatialDim);

        interface->AddMortarElement(mtele);
      }

      ggsize += gsize;  // update global element counter
    }

    /* -------------------- finalize the meshtying interface construction
     *
     * If this is the final parallel distribution, we need to assign degrees of freedom during
     * during FillComplete(). If parallel redistribution is enabled, there will be another call to
     * FillComplete(), so we skip this expensive operation here and do it later. DOFs have to be
     * assigned only once!
     */
    {
      const INPAR::MORTAR::ParallelRedist parallelRedist =
          Teuchos::getIntegralValue<INPAR::MORTAR::ParallelRedist>(
              mtparams.sublist("PARALLEL REDISTRIBUTION"), "PARALLEL_REDIST");
      bool isFinalDistribution = false;
      if (parallelRedist == INPAR::MORTAR::ParallelRedist::redist_none or comm_->NumProc() == 1)
        isFinalDistribution = true;

      interface->FillComplete(isFinalDistribution, maxdof);
    }
  }
  if (Comm().MyPID() == 0) std::cout << "done!" << std::endl;

  //**********************************************************************
  // create the solver strategy object
  // and pass all necessary data to it
  if (Comm().MyPID() == 0)
    std::cout << "Building meshtying strategy object............" << std::endl;

  const GLOBAL::ProblemType problemtype = GLOBAL::Problem::Instance()->GetProblemType();

  INPAR::CONTACT::SolvingStrategy stype =
      INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(mtparams, "STRATEGY");
  if (stype == INPAR::CONTACT::solution_lagmult)
  {
    // finally we should use another criteria to decide which strategy
    if (problemtype != GLOBAL::ProblemType::poroelast && problemtype != GLOBAL::ProblemType::fpsi &&
        problemtype != GLOBAL::ProblemType::fpsi_xfem && problemtype != GLOBAL::ProblemType::fps3i)
    {
      strategy_ = Teuchos::rcp(new MtLagrangeStrategy(discret.DofRowMap(), discret.NodeRowMap(),
          mtparams, interfaces, spatialDim, comm_, alphaf, maxdof));
    }
    else
    {
      strategy_ = Teuchos::rcp(new PoroMtLagrangeStrategy(discret.DofRowMap(), discret.NodeRowMap(),
          mtparams, interfaces, spatialDim, comm_, alphaf, maxdof));
    }
  }
  else if (stype == INPAR::CONTACT::solution_penalty or stype == INPAR::CONTACT::solution_uzawa)
    strategy_ = Teuchos::rcp(new MtPenaltyStrategy(discret.DofRowMap(), discret.NodeRowMap(),
        mtparams, interfaces, spatialDim, comm_, alphaf, maxdof));
  else
    dserror("Unrecognized strategy");

  if (Comm().MyPID() == 0) std::cout << "done!" << std::endl;
  //**********************************************************************

  //**********************************************************************
  // parallel redistribution of all interfaces
  GetStrategy().RedistributeMeshtying();
  //**********************************************************************

  // create binary search tree
  for (auto& interface : interfaces) interface->CreateSearchTree();

  return;
}


/*----------------------------------------------------------------------*
 |  read and check input parameters (public)                  popp 04/08|
 *----------------------------------------------------------------------*/
bool CONTACT::MtManager::ReadAndCheckInput(
    Teuchos::ParameterList& mtparams, const DRT::Discretization& discret)
{
  // read parameter lists from GLOBAL::Problem
  const Teuchos::ParameterList& mortar = GLOBAL::Problem::Instance()->MortarCouplingParams();
  const Teuchos::ParameterList& meshtying = GLOBAL::Problem::Instance()->ContactDynamicParams();
  const Teuchos::ParameterList& wearlist = GLOBAL::Problem::Instance()->WearParams();

  // read Problem Type and Problem Dimension from GLOBAL::Problem
  const GLOBAL::ProblemType problemtype = GLOBAL::Problem::Instance()->GetProblemType();
  const int spatialDim = GLOBAL::Problem::Instance()->NDim();
  CORE::FE::ShapeFunctionType distype = GLOBAL::Problem::Instance()->SpatialApproximationType();

  // get mortar information
  std::vector<DRT::Condition*> mtcond(0);
  std::vector<DRT::Condition*> ccond(0);

  discret.GetCondition("Mortar", mtcond);
  discret.GetCondition("Contact", ccond);

  bool onlymeshtying = false;
  bool meshtyingandcontact = false;

  // check for case
  if (mtcond.size() != 0 and ccond.size() != 0) meshtyingandcontact = true;

  if (mtcond.size() != 0 and ccond.size() == 0) onlymeshtying = true;

  // *********************************************************************
  // invalid parallel strategies
  // *********************************************************************
  const Teuchos::ParameterList& mortarParallelRedistParams =
      mortar.sublist("PARALLEL REDISTRIBUTION");

  if (Teuchos::getIntegralValue<INPAR::MORTAR::ExtendGhosting>(mortarParallelRedistParams,
          "GHOSTING_STRATEGY") == INPAR::MORTAR::ExtendGhosting::roundrobin)
    dserror(
        "Extending the ghosting via a Round-Robin loop is not implemented for mortar meshtying.");

  // *********************************************************************
  // invalid parameter combinations
  // *********************************************************************
  if (INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") ==
          INPAR::CONTACT::solution_penalty &&
      meshtying.get<double>("PENALTYPARAM") <= 0.0)
    dserror("Penalty parameter eps <= 0, must be greater than 0");

  if (INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") ==
          INPAR::CONTACT::solution_uzawa &&
      meshtying.get<double>("PENALTYPARAM") <= 0.0)
    dserror("Penalty parameter eps <= 0, must be greater than 0");

  if (INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") ==
          INPAR::CONTACT::solution_uzawa &&
      meshtying.get<int>("UZAWAMAXSTEPS") < 2)
    dserror("Maximum number of Uzawa / Augmentation steps must be at least 2");

  if (INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") ==
          INPAR::CONTACT::solution_uzawa &&
      meshtying.get<double>("UZAWACONSTRTOL") <= 0.0)
    dserror("Constraint tolerance for Uzawa / Augmentation scheme must be greater than 0");

  if (onlymeshtying && INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(meshtying, "FRICTION") !=
                           INPAR::CONTACT::friction_none)
    dserror("Friction law supplied for mortar meshtying");

  if (INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") ==
          INPAR::CONTACT::solution_lagmult &&
      INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
          INPAR::MORTAR::shape_standard &&
      (INPUT::IntegralValue<INPAR::CONTACT::SystemType>(meshtying, "SYSTEM") ==
              INPAR::CONTACT::system_condensed ||
          INPUT::IntegralValue<INPAR::CONTACT::SystemType>(meshtying, "SYSTEM") ==
              INPAR::CONTACT::system_condensed_lagmult))
    dserror("Condensation of linear system only possible for dual Lagrange multipliers");

  if (Teuchos::getIntegralValue<INPAR::MORTAR::ParallelRedist>(mortarParallelRedistParams,
          "PARALLEL_REDIST") == INPAR::MORTAR::ParallelRedist::redist_dynamic and
      onlymeshtying)
    dserror("Dynamic parallel redistribution not possible for meshtying");

  if (Teuchos::getIntegralValue<INPAR::MORTAR::ParallelRedist>(mortarParallelRedistParams,
          "PARALLEL_REDIST") != INPAR::MORTAR::ParallelRedist::redist_none &&
      mortarParallelRedistParams.get<int>("MIN_ELEPROC") < 0)
    dserror(
        "ERROR: Minimum number of elements per processor for parallel redistribution must be >= 0");

  if (INPUT::IntegralValue<INPAR::MORTAR::ConsistentDualType>(mortar, "LM_DUAL_CONSISTENT") !=
          INPAR::MORTAR::consistent_none &&
      INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") !=
          INPAR::CONTACT::solution_lagmult &&
      INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") !=
          INPAR::MORTAR::shape_standard)
    dserror(
        "ERROR: Consistent dual shape functions in boundary elements only for Lagrange multiplier "
        "strategy.");

  if (INPUT::IntegralValue<INPAR::MORTAR::ConsistentDualType>(mortar, "LM_DUAL_CONSISTENT") !=
          INPAR::MORTAR::consistent_none &&
      INPUT::IntegralValue<INPAR::MORTAR::IntType>(mortar, "INTTYPE") ==
          INPAR::MORTAR::inttype_elements &&
      (INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
              INPAR::MORTAR::shape_dual ||
          INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
              INPAR::MORTAR::shape_petrovgalerkin))

    // *********************************************************************
    // not (yet) implemented combinations
    // *********************************************************************
    if (INPUT::IntegralValue<int>(mortar, "CROSSPOINTS") == true && spatialDim == 3)
      dserror("Crosspoints / edge node modification not yet implemented for 3D");

  if (INPUT::IntegralValue<int>(mortar, "CROSSPOINTS") == true &&
      INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar, "LM_QUAD") ==
          INPAR::MORTAR::lagmult_lin)
    dserror("Crosspoints and linear LM interpolation for quadratic FE not yet compatible");

  if (INPUT::IntegralValue<int>(mortar, "CROSSPOINTS") == true &&
      Teuchos::getIntegralValue<INPAR::MORTAR::ParallelRedist>(mortarParallelRedistParams,
          "PARALLEL_REDIST") != INPAR::MORTAR::ParallelRedist::redist_none)
    dserror("Crosspoints and parallel redistribution not yet compatible");

  if (INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
          INPAR::MORTAR::shape_petrovgalerkin and
      onlymeshtying)
    dserror("Petrov-Galerkin approach makes no sense for meshtying");

  // *********************************************************************
  // 3D quadratic mortar (choice of interpolation and testing fcts.)
  // *********************************************************************
  if (INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar, "LM_QUAD") ==
          INPAR::MORTAR::lagmult_pwlin &&
      INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
          INPAR::MORTAR::shape_dual)
    dserror(
        "ERROR: No pwlin approach (for LM) implemented for quadratic meshtying with DUAL shape "
        "fct.");

  // *********************************************************************
  // element-based vs. segment-based mortar integration
  // *********************************************************************
  INPAR::MORTAR::IntType inttype = INPUT::IntegralValue<INPAR::MORTAR::IntType>(mortar, "INTTYPE");

  if (inttype == INPAR::MORTAR::inttype_elements && mortar.get<int>("NUMGP_PER_DIM") <= 0)
    dserror("Invalid Gauss point number NUMGP_PER_DIM for element-based integration.");

  if (inttype == INPAR::MORTAR::inttype_elements_BS && mortar.get<int>("NUMGP_PER_DIM") <= 0)
    dserror(
        "ERROR: Invalid Gauss point number NUMGP_PER_DIM for element-based integration with "
        "boundary segmentation."
        "\nPlease note that the value you have to provide only applies to the element-based "
        "integration"
        "\ndomain, while pre-defined default values will be used in the segment-based boundary "
        "domain.");

  if ((inttype == INPAR::MORTAR::inttype_elements ||
          inttype == INPAR::MORTAR::inttype_elements_BS) &&
      mortar.get<int>("NUMGP_PER_DIM") <= 1)
    dserror("Invalid Gauss point number NUMGP_PER_DIM for element-based integration.");

  // *********************************************************************
  // warnings
  // *********************************************************************
  if (mortar.get<double>("SEARCH_PARAM") == 0.0 && Comm().MyPID() == 0)
    std::cout << ("Warning: Meshtying search called without inflation of bounding volumes\n")
              << std::endl;

  // get parameter lists
  mtparams.setParameters(mortar);
  mtparams.setParameters(meshtying);
  mtparams.setParameters(wearlist);

  // *********************************************************************
  // predefined params for meshtying and contact
  // *********************************************************************
  if (meshtyingandcontact)
  {
    // set options for mortar coupling
    mtparams.set<std::string>("SEARCH_ALGORITHM", "Binarytree");
    mtparams.set<double>("SEARCH_PARAM", 0.3);
    mtparams.set<std::string>("SEARCH_USE_AUX_POS", "no");
    mtparams.set<std::string>("LM_SHAPEFCN", "dual");
    mtparams.set<std::string>("SYSTEM", "condensed");
    mtparams.set<bool>("NURBS", false);
    mtparams.set<int>("NUMGP_PER_DIM", -1);
    mtparams.set<std::string>("STRATEGY", "LagrangianMultipliers");
    mtparams.set<std::string>("INTTYPE", "segments");
    mtparams.sublist("PARALLEL REDISTRIBUTION").set<std::string>("REDUNDANT_STORAGE", "Master");
    mtparams.sublist("PARALLEL REDISTRIBUTION").set<std::string>("PARALLEL_REDIST", "static");
  }
  // *********************************************************************
  // smooth interfaces
  // *********************************************************************
  // NURBS PROBLEM?
  switch (distype)
  {
    case CORE::FE::ShapeFunctionType::nurbs:
    {
      mtparams.set<bool>("NURBS", true);
      break;
    }
    default:
    {
      mtparams.set<bool>("NURBS", false);
      break;
    }
  }

  // *********************************************************************
  // poroelastic meshtying
  // *********************************************************************
  if ((problemtype == GLOBAL::ProblemType::poroelast || problemtype == GLOBAL::ProblemType::fpsi ||
          problemtype == GLOBAL::ProblemType::fpsi_xfem) &&
      (INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") !=
              INPAR::MORTAR::shape_dual &&
          INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") !=
              INPAR::MORTAR::shape_petrovgalerkin))
    dserror("POROCONTACT: Only dual and petrovgalerkin shape functions implemented yet!");

  if ((problemtype == GLOBAL::ProblemType::poroelast || problemtype == GLOBAL::ProblemType::fpsi ||
          problemtype == GLOBAL::ProblemType::fpsi_xfem) &&
      Teuchos::getIntegralValue<INPAR::MORTAR::ParallelRedist>(mortarParallelRedistParams,
          "PARALLEL_REDIST") != INPAR::MORTAR::ParallelRedist::redist_none)
    dserror(
        "POROCONTACT: Parallel Redistribution not implemented yet!");  // Since we use Pointers to
                                                                       // Parent Elements, which are
                                                                       // not copied to other procs!

  if ((problemtype == GLOBAL::ProblemType::poroelast || problemtype == GLOBAL::ProblemType::fpsi ||
          problemtype == GLOBAL::ProblemType::fpsi_xfem) &&
      INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") !=
          INPAR::CONTACT::solution_lagmult)
    dserror("POROCONTACT: Use Lagrangean Strategy for poro meshtying!");

  if ((problemtype == GLOBAL::ProblemType::poroelast || problemtype == GLOBAL::ProblemType::fpsi ||
          problemtype == GLOBAL::ProblemType::fpsi_xfem) &&
      INPUT::IntegralValue<INPAR::CONTACT::SystemType>(meshtying, "SYSTEM") !=
          INPAR::CONTACT::system_condensed_lagmult)
    dserror("POROCONTACT: Just lagrange multiplier should be condensed for poro meshtying!");

  if ((problemtype == GLOBAL::ProblemType::poroelast || problemtype == GLOBAL::ProblemType::fpsi ||
          problemtype == GLOBAL::ProblemType::fpsi_xfem) &&
      (spatialDim != 3) && (spatialDim != 2))
  {
    const Teuchos::ParameterList& porodyn = GLOBAL::Problem::Instance()->PoroelastDynamicParams();
    if (INPUT::IntegralValue<int>(porodyn, "CONTACTNOPEN"))
      dserror("POROCONTACT: PoroMeshtying with no penetration just tested for 3d (and 2d)!");
  }

  mtparams.setName("CONTACT DYNAMIC / MORTAR COUPLING");

  // no parallel redistribution in the serial case
  if (Comm().NumProc() == 1)
    mtparams.sublist("PARALLEL REDISTRIBUTION").set<std::string>("PARALLEL_REDIST", "None");

  return true;
}

/*----------------------------------------------------------------------*
 |  write restart information for meshtying (public)          popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtManager::WriteRestart(IO::DiscretizationWriter& output, bool forcedrestart)
{
  output.WriteVector("mt_lagrmultold", GetStrategy().LagrMultOld());

  return;
}

/*----------------------------------------------------------------------*
 |  read restart information for meshtying (public)           popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtManager::ReadRestart(IO::DiscretizationReader& reader,
    Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<Epetra_Vector> zero)
{
  // this is meshtying, thus we need zeros for restart
  // let strategy object do all the work
  GetStrategy().DoReadRestart(reader, zero);

  return;
}

/*----------------------------------------------------------------------*
 |  write interface tractions for postprocessing (public)     popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtManager::PostprocessQuantities(IO::DiscretizationWriter& output)
{
  // evaluate interface tractions
  Teuchos::RCP<Epetra_Map> problem = GetStrategy().ProblemDofs();
  Teuchos::RCP<Epetra_Vector> traction =
      Teuchos::rcp(new Epetra_Vector(*(GetStrategy().LagrMultOld())));
  Teuchos::RCP<Epetra_Vector> tractionexp = Teuchos::rcp(new Epetra_Vector(*problem));
  CORE::LINALG::Export(*traction, *tractionexp);

  // evaluate slave and master forces
  Teuchos::RCP<Epetra_Vector> fcslave =
      Teuchos::rcp(new Epetra_Vector(GetStrategy().DMatrix()->RowMap()));
  Teuchos::RCP<Epetra_Vector> fcmaster =
      Teuchos::rcp(new Epetra_Vector(GetStrategy().MMatrix()->DomainMap()));
  Teuchos::RCP<Epetra_Vector> fcslaveexp = Teuchos::rcp(new Epetra_Vector(*problem));
  Teuchos::RCP<Epetra_Vector> fcmasterexp = Teuchos::rcp(new Epetra_Vector(*problem));
  GetStrategy().DMatrix()->Multiply(true, *traction, *fcslave);
  GetStrategy().MMatrix()->Multiply(true, *traction, *fcmaster);
  CORE::LINALG::Export(*fcslave, *fcslaveexp);
  CORE::LINALG::Export(*fcmaster, *fcmasterexp);

  // write to output
  output.WriteVector("interfacetraction", tractionexp);
  output.WriteVector("slaveforces", fcslaveexp);
  output.WriteVector("masterforces", fcmasterexp);

  return;
}

/*-----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::MtManager::PostprocessQuantitiesPerInterface(
    Teuchos::RCP<Teuchos::ParameterList> outputParams)
{
  GetStrategy().PostprocessQuantitiesPerInterface(outputParams);
}

BACI_NAMESPACE_CLOSE
