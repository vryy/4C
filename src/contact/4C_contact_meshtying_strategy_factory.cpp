/*---------------------------------------------------------------------*/
/*! \file
\brief Factory to create the desired meshtying strategy.


\level 3

*/
/*---------------------------------------------------------------------*/

#include "4C_contact_meshtying_strategy_factory.hpp"

#include "4C_contact_abstract_strategy.hpp"
#include "4C_contact_meshtying_abstract_strategy.hpp"
#include "4C_contact_meshtying_lagrange_strategy.hpp"
#include "4C_contact_meshtying_penalty_strategy.hpp"
#include "4C_contact_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io.hpp"
#include "4C_io_pstream.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_mortar_element.hpp"
#include "4C_mortar_node.hpp"
#include "4C_mortar_utils.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::FactoryMT::Setup()
{
  check_init();
  MORTAR::STRATEGY::Factory::Setup();

  set_is_setup();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::FactoryMT::read_and_check_input(Teuchos::ParameterList& params) const
{
  // read parameter lists from GLOBAL::Problem
  const Teuchos::ParameterList& mortar = GLOBAL::Problem::Instance()->mortar_coupling_params();
  const Teuchos::ParameterList& meshtying = GLOBAL::Problem::Instance()->contact_dynamic_params();
  const Teuchos::ParameterList& wearlist = GLOBAL::Problem::Instance()->WearParams();

  // read Problem Type and Problem Dimension from GLOBAL::Problem
  const GLOBAL::ProblemType problemtype = GLOBAL::Problem::Instance()->GetProblemType();
  int dim = GLOBAL::Problem::Instance()->NDim();
  CORE::FE::ShapeFunctionType distype = GLOBAL::Problem::Instance()->spatial_approximation_type();

  // get mortar information
  std::vector<CORE::Conditions::Condition*> mtcond(0);
  std::vector<CORE::Conditions::Condition*> ccond(0);

  discret().GetCondition("Mortar", mtcond);
  discret().GetCondition("Contact", ccond);

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
    FOUR_C_THROW(
        "Extending the ghosting via a Round-Robin loop is not implemented for mortar meshtying.");

  // *********************************************************************
  // invalid parameter combinations
  // *********************************************************************
  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") ==
          INPAR::CONTACT::solution_penalty &&
      meshtying.get<double>("PENALTYPARAM") <= 0.0)
    FOUR_C_THROW("Penalty parameter eps = 0, must be greater than 0");

  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") ==
          INPAR::CONTACT::solution_uzawa &&
      meshtying.get<double>("PENALTYPARAM") <= 0.0)
    FOUR_C_THROW("Penalty parameter eps = 0, must be greater than 0");

  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") ==
          INPAR::CONTACT::solution_uzawa &&
      meshtying.get<int>("UZAWAMAXSTEPS") < 2)
    FOUR_C_THROW("Maximum number of Uzawa / Augmentation steps must be at least 2");

  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") ==
          INPAR::CONTACT::solution_uzawa &&
      meshtying.get<double>("UZAWACONSTRTOL") <= 0.0)
    FOUR_C_THROW("Constraint tolerance for Uzawa / Augmentation scheme must be greater than 0");

  if (onlymeshtying && CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(
                           meshtying, "FRICTION") != INPAR::CONTACT::friction_none)
    FOUR_C_THROW("Friction law supplied for mortar meshtying");

  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") ==
          INPAR::CONTACT::solution_lagmult &&
      CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
          INPAR::MORTAR::shape_standard &&
      (CORE::UTILS::IntegralValue<INPAR::CONTACT::SystemType>(meshtying, "SYSTEM") ==
              INPAR::CONTACT::system_condensed ||
          CORE::UTILS::IntegralValue<INPAR::CONTACT::SystemType>(meshtying, "SYSTEM") ==
              INPAR::CONTACT::system_condensed_lagmult))
    FOUR_C_THROW("Condensation of linear system only possible for dual Lagrange multipliers");

  if (Teuchos::getIntegralValue<INPAR::MORTAR::ParallelRedist>(mortarParallelRedistParams,
          "PARALLEL_REDIST") == INPAR::MORTAR::ParallelRedist::redist_dynamic and
      onlymeshtying)
    FOUR_C_THROW("Dynamic parallel redistribution not possible for meshtying");

  if (Teuchos::getIntegralValue<INPAR::MORTAR::ParallelRedist>(mortarParallelRedistParams,
          "PARALLEL_REDIST") != INPAR::MORTAR::ParallelRedist::redist_none &&
      mortarParallelRedistParams.get<int>("MIN_ELEPROC") < 0)
    FOUR_C_THROW(
        "ERROR: Minimum number of elements per processor for parallel redistribution must be >= 0");

  if (CORE::UTILS::IntegralValue<INPAR::MORTAR::ConsistentDualType>(mortar, "LM_DUAL_CONSISTENT") !=
          INPAR::MORTAR::consistent_none &&
      CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") !=
          INPAR::CONTACT::solution_lagmult &&
      CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") !=
          INPAR::MORTAR::shape_standard)
    FOUR_C_THROW(
        "ERROR: Consistent dual shape functions in boundary elements only for Lagrange multiplier "
        "strategy.");

  if (CORE::UTILS::IntegralValue<INPAR::MORTAR::ConsistentDualType>(mortar, "LM_DUAL_CONSISTENT") !=
          INPAR::MORTAR::consistent_none &&
      CORE::UTILS::IntegralValue<INPAR::MORTAR::IntType>(mortar, "INTTYPE") ==
          INPAR::MORTAR::inttype_elements &&
      (CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
              INPAR::MORTAR::shape_dual ||
          CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
              INPAR::MORTAR::shape_petrovgalerkin))

    // *********************************************************************
    // not (yet) implemented combinations
    // *********************************************************************
    if (CORE::UTILS::IntegralValue<int>(mortar, "CROSSPOINTS") == true && dim == 3)
      FOUR_C_THROW("Crosspoints / edge node modification not yet implemented for 3D");

  if (CORE::UTILS::IntegralValue<int>(mortar, "CROSSPOINTS") == true &&
      CORE::UTILS::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar, "LM_QUAD") ==
          INPAR::MORTAR::lagmult_lin)
    FOUR_C_THROW("Crosspoints and linear LM interpolation for quadratic FE not yet compatible");

  if (CORE::UTILS::IntegralValue<int>(mortar, "CROSSPOINTS") == true &&
      Teuchos::getIntegralValue<INPAR::MORTAR::ParallelRedist>(mortarParallelRedistParams,
          "PARALLEL_REDIST") != INPAR::MORTAR::ParallelRedist::redist_none)
    FOUR_C_THROW("Crosspoints and parallel redistribution not yet compatible");

  if (CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
          INPAR::MORTAR::shape_petrovgalerkin and
      onlymeshtying)
    FOUR_C_THROW("Petrov-Galerkin approach makes no sense for meshtying");

  // *********************************************************************
  // 3D quadratic mortar (choice of interpolation and testing fcts.)
  // *********************************************************************
  if (CORE::UTILS::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar, "LM_QUAD") ==
          INPAR::MORTAR::lagmult_pwlin &&
      CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
          INPAR::MORTAR::shape_dual)
    FOUR_C_THROW(
        "ERROR: No pwlin approach (for LM) implemented for quadratic meshtying with DUAL shape "
        "fct.");

  // *********************************************************************
  // element-based vs. segment-based mortar integration
  // *********************************************************************
  INPAR::MORTAR::IntType inttype =
      CORE::UTILS::IntegralValue<INPAR::MORTAR::IntType>(mortar, "INTTYPE");

  if (inttype == INPAR::MORTAR::inttype_elements && mortar.get<int>("NUMGP_PER_DIM") <= 0)
    FOUR_C_THROW("Invalid Gauss point number NUMGP_PER_DIM for element-based integration.");

  if (inttype == INPAR::MORTAR::inttype_elements_BS && mortar.get<int>("NUMGP_PER_DIM") <= 0)
    FOUR_C_THROW(
        "ERROR: Invalid Gauss point number NUMGP_PER_DIM for element-based integration with "
        "boundary segmentation."
        "\nPlease note that the value you have to provide only applies to the element-based "
        "integration"
        "\ndomain, while pre-defined default values will be used in the segment-based boundary "
        "domain.");

  if ((inttype == INPAR::MORTAR::inttype_elements ||
          inttype == INPAR::MORTAR::inttype_elements_BS) &&
      mortar.get<int>("NUMGP_PER_DIM") <= 1)
    FOUR_C_THROW("Invalid Gauss point number NUMGP_PER_DIM for element-based integration.");

  // *********************************************************************
  // warnings
  // *********************************************************************
  if (mortar.get<double>("SEARCH_PARAM") == 0.0 && comm().MyPID() == 0)
    std::cout << ("Warning: Meshtying search called without inflation of bounding volumes\n")
              << std::endl;

  // get parameter lists
  params.setParameters(mortar);
  params.setParameters(meshtying);
  params.setParameters(wearlist);

  // *********************************************************************
  // predefined params for meshtying and contact
  // *********************************************************************
  if (meshtyingandcontact)
  {
    // set options for mortar coupling
    params.set<std::string>("SEARCH_ALGORITHM", "Binarytree");
    params.set<double>("SEARCH_PARAM", 0.3);
    params.set<std::string>("SEARCH_USE_AUX_POS", "no");
    params.set<std::string>("LM_SHAPEFCN", "dual");
    params.set<std::string>("SYSTEM", "condensed");
    params.set<bool>("NURBS", false);
    params.set<int>("NUMGP_PER_DIM", -1);
    params.set<std::string>("STRATEGY", "LagrangianMultipliers");
    params.set<std::string>("INTTYPE", "segments");
    params.sublist("PARALLEL REDISTRIBUTION").set<std::string>("REDUNDANT_STORAGE", "Master");
    params.sublist("PARALLEL REDISTRIBUTION").set<std::string>("PARALLEL_REDIST", "static");
  }
  // *********************************************************************
  // smooth interfaces
  // *********************************************************************
  // NURBS PROBLEM?
  switch (distype)
  {
    case CORE::FE::ShapeFunctionType::nurbs:
    {
      params.set<bool>("NURBS", true);
      break;
    }
    default:
    {
      params.set<bool>("NURBS", false);
      break;
    }
  }

  // *********************************************************************
  // poroelastic meshtying
  // *********************************************************************
  if ((problemtype == GLOBAL::ProblemType::poroelast || problemtype == GLOBAL::ProblemType::fpsi ||
          problemtype == GLOBAL::ProblemType::fpsi_xfem) &&
      (CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") !=
              INPAR::MORTAR::shape_dual &&
          CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") !=
              INPAR::MORTAR::shape_petrovgalerkin))
    FOUR_C_THROW("POROCONTACT: Only dual and petrovgalerkin shape functions implemented yet!");

  if ((problemtype == GLOBAL::ProblemType::poroelast || problemtype == GLOBAL::ProblemType::fpsi ||
          problemtype == GLOBAL::ProblemType::fpsi_xfem) &&
      Teuchos::getIntegralValue<INPAR::MORTAR::ParallelRedist>(mortarParallelRedistParams,
          "PARALLEL_REDIST") != INPAR::MORTAR::ParallelRedist::redist_none)
    FOUR_C_THROW(
        "POROCONTACT: Parallel Redistribution not implemented yet!");  // Since we use Pointers to
                                                                       // Parent Elements, which are
                                                                       // not copied to other procs!

  if ((problemtype == GLOBAL::ProblemType::poroelast || problemtype == GLOBAL::ProblemType::fpsi ||
          problemtype == GLOBAL::ProblemType::fpsi_xfem) &&
      CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") !=
          INPAR::CONTACT::solution_lagmult)
    FOUR_C_THROW("POROCONTACT: Use Lagrangean Strategy for poro meshtying!");

  if ((problemtype == GLOBAL::ProblemType::poroelast || problemtype == GLOBAL::ProblemType::fpsi ||
          problemtype == GLOBAL::ProblemType::fpsi_xfem) &&
      CORE::UTILS::IntegralValue<INPAR::CONTACT::SystemType>(meshtying, "SYSTEM") !=
          INPAR::CONTACT::system_condensed_lagmult)
    FOUR_C_THROW("POROCONTACT: Just lagrange multiplier should be condensed for poro meshtying!");

  if ((problemtype == GLOBAL::ProblemType::poroelast || problemtype == GLOBAL::ProblemType::fpsi ||
          problemtype == GLOBAL::ProblemType::fpsi_xfem) &&
      (dim != 3) && (dim != 2))
  {
    const Teuchos::ParameterList& porodyn = GLOBAL::Problem::Instance()->poroelast_dynamic_params();
    if (CORE::UTILS::IntegralValue<int>(porodyn, "CONTACTNOPEN"))
      FOUR_C_THROW("POROCONTACT: PoroMeshtying with no penetration just tested for 3d (and 2d)!");
  }

  params.setName("CONTACT DYNAMIC / MORTAR COUPLING");

  // no parallel redistribution in the serial case
  if (comm().NumProc() == 1)
    params.sublist("PARALLEL REDISTRIBUTION").set<std::string>("PARALLEL_REDIST", "None");

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::FactoryMT::BuildInterfaces(const Teuchos::ParameterList& mtparams,
    std::vector<Teuchos::RCP<MORTAR::Interface>>& interfaces, bool& poroslave,
    bool& poromaster) const
{
  int dim = GLOBAL::Problem::Instance()->NDim();

  // start building interfaces
  if (comm().MyPID() == 0)
  {
    std::cout << "Building contact interface(s)...............";
    fflush(stdout);
  }

  std::vector<CORE::Conditions::Condition*> contactconditions(0);
  discret().GetCondition("Mortar", contactconditions);

  // there must be more than one meshtying condition
  if ((int)contactconditions.size() < 2)
    FOUR_C_THROW("Not enough contact conditions in discretization");

  // find all pairs of matching meshtying conditions
  // there is a maximum of (conditions / 2) groups
  std::vector<int> foundgroups(0);
  int numgroupsfound = 0;

  // get nurbs information
  const bool nurbs = mtparams.get<bool>("NURBS");

  // maximum dof number in discretization
  // later we want to create NEW Lagrange multiplier degrees of
  // freedom, which of course must not overlap with displacement dofs
  int maxdof = discret().dof_row_map()->MaxAllGID();

  for (int i = 0; i < (int)contactconditions.size(); ++i)
  {
    // initialize vector for current group of conditions and temp condition
    std::vector<CORE::Conditions::Condition*> currentgroup(0);
    CORE::Conditions::Condition* tempcond = nullptr;

    // try to build meshtying group around this condition
    currentgroup.push_back(contactconditions[i]);
    const int groupid1 = currentgroup[0]->parameters().Get<int>("Interface ID");
    bool foundit = false;

    for (int j = 0; j < (int)contactconditions.size(); ++j)
    {
      if (j == i) continue;  // do not detect contactconditions[i] again
      tempcond = contactconditions[j];
      const int groupid2 = tempcond->parameters().Get<int>("Interface ID");

      if (groupid1 != groupid2) continue;  // not in the group
      foundit = true;                      // found a group entry
      currentgroup.push_back(tempcond);    // store it in currentgroup
    }

    // now we should have found a group of conds
    if (!foundit) FOUR_C_THROW("Cannot find matching contact condition for id %d", groupid1);

    // see whether we found this group before
    bool foundbefore = false;
    for (int j = 0; j < numgroupsfound; ++j)
      if (groupid1 == foundgroups[j])
      {
        foundbefore = true;
        break;
      }

    // if we have processed this group before, do nothing
    if (foundbefore) continue;

    // we have not found this group before, process it
    foundgroups.push_back(groupid1);
    ++numgroupsfound;

    // find out which sides are Master and Slave
    bool hasslave = false;
    bool hasmaster = false;
    std::vector<const std::string*> sides((int)currentgroup.size());
    std::vector<bool> isslave((int)currentgroup.size());

    for (int j = 0; j < (int)sides.size(); ++j)
    {
      sides[j] = &currentgroup[j]->parameters().Get<std::string>("Side");
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
        FOUR_C_THROW("MtManager: Unknown contact side qualifier!");
      }
    }

    if (!hasslave) FOUR_C_THROW("Slave side missing in contact condition group!");
    if (!hasmaster) FOUR_C_THROW("Master side missing in contact condition group!");

    // find out which sides are initialized as Active
    std::vector<const std::string*> active((int)currentgroup.size());
    std::vector<bool> isactive((int)currentgroup.size());

    for (int j = 0; j < (int)sides.size(); ++j)
    {
      active[j] = &currentgroup[j]->parameters().Get<std::string>("Initialization");
      if (*sides[j] == "Slave")
      {
        // slave sides must be initialized as "Active"
        if (*active[j] == "Active")
          isactive[j] = true;
        else if (*active[j] == "Inactive")
          FOUR_C_THROW("Slave side must be active for meshtying!");
        else
          FOUR_C_THROW("Unknown contact init qualifier!");
      }
      else if (*sides[j] == "Master")
      {
        // master sides must NOT be initialized as "Active" as this makes no sense
        if (*active[j] == "Active")
          FOUR_C_THROW("Master side cannot be active!");
        else if (*active[j] == "Inactive")
          isactive[j] = false;
        else
          FOUR_C_THROW("Unknown contact init qualifier!");
      }
      else
      {
        FOUR_C_THROW("MtManager: Unknown contact side qualifier!");
      }
    }

    // create an empty meshtying interface and store it in this Manager
    // (for structural meshtying we currently choose redundant master storage)
    interfaces.push_back(MORTAR::Interface::Create(groupid1, comm(), dim, mtparams));

    // get it again
    Teuchos::RCP<MORTAR::Interface> interface = interfaces[(int)interfaces.size() - 1];

    // note that the nodal ids are unique because they come from
    // one global problem discretization conatining all nodes of the
    // contact interface
    // We rely on this fact, therefore it is not possible to
    // do meshtying between two distinct discretizations here

    //-------------------------------------------------- process nodes
    for (int j = 0; j < (int)currentgroup.size(); ++j)
    {
      // get all nodes and add them
      const std::vector<int>* nodeids = currentgroup[j]->GetNodes();
      if (!nodeids) FOUR_C_THROW("Condition does not have Node Ids");
      for (int k = 0; k < (int)(*nodeids).size(); ++k)
      {
        int gid = (*nodeids)[k];
        // do only nodes that I have in my discretization
        if (!discret_ptr_->NodeColMap()->MyGID(gid)) continue;
        CORE::Nodes::Node* node = discret().gNode(gid);
        if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);

        // create Node object
        Teuchos::RCP<MORTAR::Node> mtnode = Teuchos::rcp(new MORTAR::Node(
            node->Id(), node->X(), node->Owner(), discret().Dof(0, node), isslave[j]));
        //-------------------
        // get nurbs weight!
        if (nurbs)
        {
          MORTAR::UTILS::prepare_nurbs_node(node, mtnode);
        }

        // get edge and corner information:
        std::vector<CORE::Conditions::Condition*> contactcornercond(0);
        discret().GetCondition("mrtrcorner", contactcornercond);
        for (unsigned j = 0; j < contactcornercond.size(); j++)
        {
          if (contactcornercond.at(j)->ContainsNode(node->Id()))
          {
            mtnode->SetOnCorner() = true;
          }
        }
        std::vector<CORE::Conditions::Condition*> contactedgecond(0);
        discret().GetCondition("mrtredge", contactedgecond);
        for (unsigned j = 0; j < contactedgecond.size(); j++)
        {
          if (contactedgecond.at(j)->ContainsNode(node->Id()))
          {
            mtnode->SetOnEdge() = true;
          }
        }

        // Check, if this node (and, in case, which dofs) are in the contact symmetry condition
        std::vector<CORE::Conditions::Condition*> contactSymconditions(0);
        discret().GetCondition("mrtrsym", contactSymconditions);

        for (unsigned j = 0; j < contactSymconditions.size(); j++)
          if (contactSymconditions.at(j)->ContainsNode(node->Id()))
          {
            const auto& onoff =
                contactSymconditions.at(j)->parameters().Get<std::vector<int>>("onoff");
            for (unsigned k = 0; k < onoff.size(); k++)
              if (onoff.at(k) == 1) mtnode->DbcDofs()[k] = true;
          }

        // note that we do not have to worry about double entries
        // as the AddNode function can deal with this case!
        interface->AddMortarNode(mtnode);
      }
    }

    //----------------------------------------------- process elements
    int ggsize = 0;
    for (int j = 0; j < (int)currentgroup.size(); ++j)
    {
      // get elements from condition j of current group
      std::map<int, Teuchos::RCP<CORE::Elements::Element>>& currele = currentgroup[j]->Geometry();

      // elements in a boundary condition have a unique id
      // but ids are not unique among 2 distinct conditions
      // due to the way elements in conditions are build.
      // We therefore have to give the second, third,... set of elements
      // different ids. ids do not have to be continous, we just add a large
      // enough number ggsize to all elements of cond2, cond3,... so they are
      // different from those in cond1!!!
      // note that elements in ele1/ele2 already are in column (overlapping) map
      int lsize = (int)currele.size();
      int gsize = 0;
      comm().SumAll(&lsize, &gsize, 1);


      std::map<int, Teuchos::RCP<CORE::Elements::Element>>::iterator fool;
      for (fool = currele.begin(); fool != currele.end(); ++fool)
      {
        Teuchos::RCP<CORE::Elements::Element> ele = fool->second;
        Teuchos::RCP<MORTAR::Element> mtele = Teuchos::rcp(new MORTAR::Element(ele->Id() + ggsize,
            ele->Owner(), ele->Shape(), ele->num_node(), ele->NodeIds(), isslave[j], nurbs));
        //------------------------------------------------------------------
        // get knotvector, normal factor and zero-size information for nurbs
        if (nurbs)
        {
          MORTAR::UTILS::prepare_nurbs_element(*discret_ptr_, ele, mtele, dim);
        }

        interface->AddMortarElement(mtele);
      }  // for (fool=ele1.start(); fool != ele1.end(); ++fool)

      ggsize += gsize;  // update global element counter
    }

    //-------------------- finalize the meshtying interface construction
    interface->fill_complete(true, maxdof);

  }  // for (int i=0; i<(int)contactconditions.size(); ++i)
  if (comm().MyPID() == 0) std::cout << "done!" << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::MtAbstractStrategy> MORTAR::STRATEGY::FactoryMT::BuildStrategy(
    const Teuchos::ParameterList& params, const bool& poroslave, const bool& poromaster,
    const int& dof_offset, std::vector<Teuchos::RCP<MORTAR::Interface>>& interfaces) const
{
  const INPAR::CONTACT::SolvingStrategy stype =
      CORE::UTILS::IntegralValue<enum INPAR::CONTACT::SolvingStrategy>(params, "STRATEGY");
  Teuchos::RCP<CONTACT::AbstractStratDataContainer> data_ptr = Teuchos::null;

  return BuildStrategy(stype, params, poroslave, poromaster, dof_offset, interfaces,
      discret().dof_row_map(), discret().NodeRowMap(), dim(), comm_ptr(), data_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::MtAbstractStrategy> MORTAR::STRATEGY::FactoryMT::BuildStrategy(
    const INPAR::CONTACT::SolvingStrategy stype, const Teuchos::ParameterList& params,
    const bool& poroslave, const bool& poromaster, const int& dof_offset,
    std::vector<Teuchos::RCP<MORTAR::Interface>>& interfaces, const Epetra_Map* dof_row_map,
    const Epetra_Map* node_row_map, const int dim, const Teuchos::RCP<const Epetra_Comm>& comm_ptr,
    Teuchos::RCP<MORTAR::StratDataContainer> data_ptr)
{
  Teuchos::RCP<CONTACT::MtAbstractStrategy> strategy_ptr = Teuchos::null;

  //**********************************************************************
  // create the solver strategy object
  // and pass all necessary data to it
  if (comm_ptr->MyPID() == 0)
  {
    std::cout << "Building meshtying strategy object............";
    fflush(stdout);
  }

  // Set dummy parameter. The correct parameter will be read directly from time integrator. We still
  // need to pass an argument as long as we want to support the same strategy contructor as the old
  // time integration.
  const double dummy = -1.0;

  if (stype == INPAR::CONTACT::solution_lagmult)
  {
    strategy_ptr = Teuchos::rcp(new CONTACT::MtLagrangeStrategy(
        dof_row_map, node_row_map, params, interfaces, dim, comm_ptr, dummy, dof_offset));
  }
  else if (stype == INPAR::CONTACT::solution_penalty or stype == INPAR::CONTACT::solution_uzawa)
    strategy_ptr = Teuchos::rcp(new CONTACT::MtPenaltyStrategy(
        dof_row_map, node_row_map, params, interfaces, dim, comm_ptr, dummy, dof_offset));
  else
    FOUR_C_THROW("Unrecognized strategy");

  if (comm_ptr->MyPID() == 0) std::cout << "done!" << std::endl;
  return strategy_ptr;
}

FOUR_C_NAMESPACE_CLOSE
