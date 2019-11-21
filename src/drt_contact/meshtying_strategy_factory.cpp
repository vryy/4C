/*---------------------------------------------------------------------*/
/*! \file
\brief Factory to create the desired meshtying strategy.

\maintainer Matthias Mayr

\level 3

*/
/*---------------------------------------------------------------------*/

#include "contact_abstract_strategy.H"
#include "contact_utils.H"
#include "meshtying_abstract_strategy.H"
#include "meshtying_lagrange_strategy.H"
#include "meshtying_penalty_strategy.H"
#include "meshtying_strategy_factory.H"

#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/inpar_contact.H"

#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_node.H"
#include "../drt_mortar/mortar_utils.H"

#include "../drt_structure_new/str_timint_basedataglobalstate.H"
#include "../drt_structure_new/str_utils.H"

#include "../drt_structure_xstructure/xstr_multi_discretization_wrapper.H"

#include "../linalg/linalg_utils.H"

#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::FactoryMT::Setup()
{
  CheckInit();
  MORTAR::STRATEGY::Factory::Setup();

  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::FactoryMT::ReadAndCheckInput(Teuchos::ParameterList& params) const
{
  // read parameter lists from DRT::Problem
  const Teuchos::ParameterList& mortar = DRT::Problem::Instance()->MortarCouplingParams();
  const Teuchos::ParameterList& meshtying = DRT::Problem::Instance()->ContactDynamicParams();
  const Teuchos::ParameterList& wearlist = DRT::Problem::Instance()->WearParams();

  // read Problem Type and Problem Dimension from DRT::Problem
  const ProblemType problemtype = DRT::Problem::Instance()->GetProblemType();
  int dim = DRT::Problem::Instance()->NDim();
  ShapeFunctionType distype = DRT::Problem::Instance()->SpatialApproximationType();

  // get mortar information
  std::vector<DRT::Condition*> mtcond(0);
  std::vector<DRT::Condition*> ccond(0);

  Discret().GetCondition("Mortar", mtcond);
  Discret().GetCondition("Contact", ccond);

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
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::RedundantStorage>(
          mortarParallelRedistParams, "REDUNDANT_STORAGE") == INPAR::MORTAR::redundant_master and
      DRT::INPUT::IntegralValue<INPAR::MORTAR::GhostingStrategy>(
          mortarParallelRedistParams, "GHOSTING_STRATEGY") != INPAR::MORTAR::ghosting_redundant)
    dserror(
        "ERROR: Redundant storage only reasonable in combination with parallel strategy: "
        "ghosting_redundant !");

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::RedundantStorage>(
          mortarParallelRedistParams, "REDUNDANT_STORAGE") == INPAR::MORTAR::redundant_all and
      DRT::INPUT::IntegralValue<INPAR::MORTAR::GhostingStrategy>(
          mortarParallelRedistParams, "GHOSTING_STRATEGY") != INPAR::MORTAR::ghosting_redundant)
    dserror(
        "ERROR: Redundant storage only reasonable in combination with parallel strategy: "
        "redundant_ghosting !");

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::GhostingStrategy>(
          mortarParallelRedistParams, "GHOSTING_STRATEGY") == INPAR::MORTAR::roundrobinghost)
    dserror("ERROR: Round-Robin strategies not for mortar meshtying!");

  // *********************************************************************
  // invalid parameter combinations
  // *********************************************************************
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") ==
          INPAR::CONTACT::solution_penalty &&
      meshtying.get<double>("PENALTYPARAM") <= 0.0)
    dserror("ERROR: Penalty parameter eps = 0, must be greater than 0");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") ==
          INPAR::CONTACT::solution_uzawa &&
      meshtying.get<double>("PENALTYPARAM") <= 0.0)
    dserror("ERROR: Penalty parameter eps = 0, must be greater than 0");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") ==
          INPAR::CONTACT::solution_uzawa &&
      meshtying.get<int>("UZAWAMAXSTEPS") < 2)
    dserror("ERROR: Maximum number of Uzawa / Augmentation steps must be at least 2");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") ==
          INPAR::CONTACT::solution_uzawa &&
      meshtying.get<double>("UZAWACONSTRTOL") <= 0.0)
    dserror("ERROR: Constraint tolerance for Uzawa / Augmentation scheme must be greater than 0");

  if (onlymeshtying && DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(
                           meshtying, "FRICTION") != INPAR::CONTACT::friction_none)
    dserror("ERROR: Friction law supplied for mortar meshtying");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") ==
          INPAR::CONTACT::solution_lagmult &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
          INPAR::MORTAR::shape_standard &&
      (DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(meshtying, "SYSTEM") ==
              INPAR::CONTACT::system_condensed ||
          DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(meshtying, "SYSTEM") ==
              INPAR::CONTACT::system_condensed_lagmult))
    dserror("ERROR: Condensation of linear system only possible for dual Lagrange multipliers");

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(
          mortarParallelRedistParams, "PARALLEL_REDIST") == INPAR::MORTAR::parredist_dynamic and
      onlymeshtying)
    dserror("ERROR: Dynamic parallel redistribution not possible for meshtying");

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(
          mortarParallelRedistParams, "PARALLEL_REDIST") != INPAR::MORTAR::parredist_none &&
      mortarParallelRedistParams.get<int>("MIN_ELEPROC") < 0)
    dserror(
        "ERROR: Minimum number of elements per processor for parallel redistribution must be >= 0");

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ConsistentDualType>(mortar, "LM_DUAL_CONSISTENT") !=
          INPAR::MORTAR::consistent_none &&
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") !=
          INPAR::CONTACT::solution_lagmult &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") !=
          INPAR::MORTAR::shape_standard)
    dserror(
        "ERROR: Consistent dual shape functions in boundary elements only for Lagrange multiplier "
        "strategy.");

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ConsistentDualType>(mortar, "LM_DUAL_CONSISTENT") !=
          INPAR::MORTAR::consistent_none &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(mortar, "INTTYPE") ==
          INPAR::MORTAR::inttype_elements &&
      (DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
              INPAR::MORTAR::shape_dual ||
          DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
              INPAR::MORTAR::shape_petrovgalerkin))

    // *********************************************************************
    // not (yet) implemented combinations
    // *********************************************************************
    if (DRT::INPUT::IntegralValue<int>(mortar, "CROSSPOINTS") == true && dim == 3)
      dserror("ERROR: Crosspoints / edge node modification not yet implemented for 3D");

  if (DRT::INPUT::IntegralValue<int>(mortar, "CROSSPOINTS") == true &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar, "LM_QUAD") ==
          INPAR::MORTAR::lagmult_lin)
    dserror("ERROR: Crosspoints and linear LM interpolation for quadratic FE not yet compatible");

  if (DRT::INPUT::IntegralValue<int>(mortar, "CROSSPOINTS") == true &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(
          mortarParallelRedistParams, "PARALLEL_REDIST") != INPAR::MORTAR::parredist_none)
    dserror("ERROR: Crosspoints and parallel redistribution not yet compatible");

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
          INPAR::MORTAR::shape_petrovgalerkin and
      onlymeshtying)
    dserror("ERROR: Petrov-Galerkin approach makes no sense for meshtying");

  // *********************************************************************
  // 3D quadratic mortar (choice of interpolation and testing fcts.)
  // *********************************************************************
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar, "LM_QUAD") ==
          INPAR::MORTAR::lagmult_pwlin &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
          INPAR::MORTAR::shape_dual)
    dserror(
        "ERROR: No pwlin approach (for LM) implemented for quadratic meshtying with DUAL shape "
        "fct.");

  // *********************************************************************
  // element-based vs. segment-based mortar integration
  // *********************************************************************
  INPAR::MORTAR::IntType inttype =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(mortar, "INTTYPE");

  if (inttype == INPAR::MORTAR::inttype_elements && mortar.get<int>("NUMGP_PER_DIM") <= 0)
    dserror("ERROR: Invalid Gauss point number NUMGP_PER_DIM for element-based integration.");

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
    dserror("ERROR: Invalid Gauss point number NUMGP_PER_DIM for element-based integration.");

  // *********************************************************************
  // warnings
  // *********************************************************************
  if (mortar.get<double>("SEARCH_PARAM") == 0.0 && Comm().MyPID() == 0)
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
    case ShapeFunctionType::shapefunction_nurbs:
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
  if ((problemtype == prb_poroelast || problemtype == prb_fpsi || problemtype == prb_fpsi_xfem) &&
      (DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") !=
              INPAR::MORTAR::shape_dual &&
          DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") !=
              INPAR::MORTAR::shape_petrovgalerkin))
    dserror("POROCONTACT: Only dual and petrovgalerkin shape functions implemented yet!");

  if ((problemtype == prb_poroelast || problemtype == prb_fpsi || problemtype == prb_fpsi_xfem) &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(
          mortarParallelRedistParams, "PARALLEL_REDIST") != INPAR::MORTAR::parredist_none)
    dserror(
        "POROCONTACT: Parallel Redistribution not implemented yet!");  // Since we use Pointers to
                                                                       // Parent Elements, which are
                                                                       // not copied to other procs!

  if ((problemtype == prb_poroelast || problemtype == prb_fpsi || problemtype == prb_fpsi_xfem) &&
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying, "STRATEGY") !=
          INPAR::CONTACT::solution_lagmult)
    dserror("POROCONTACT: Use Lagrangean Strategy for poro meshtying!");

  if ((problemtype == prb_poroelast || problemtype == prb_fpsi || problemtype == prb_fpsi_xfem) &&
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(meshtying, "SYSTEM") !=
          INPAR::CONTACT::system_condensed_lagmult)
    dserror("POROCONTACT: Just lagrange multiplier should be condensed for poro meshtying!");

  if ((problemtype == prb_poroelast || problemtype == prb_fpsi || problemtype == prb_fpsi_xfem) &&
      (dim != 3) && (dim != 2))
  {
    const Teuchos::ParameterList& porodyn = DRT::Problem::Instance()->PoroelastDynamicParams();
    if (DRT::INPUT::IntegralValue<int>(porodyn, "CONTACTNOPEN"))
      dserror("POROCONTACT: PoroMeshtying with no penetration just tested for 3d (and 2d)!");
  }

  params.setName("CONTACT DYNAMIC / MORTAR COUPLING");

  // no parallel redistribution in the serial case
  if (Comm().NumProc() == 1)
    params.sublist("PARALLEL REDISTRIBUTION").set<std::string>("PARALLEL_REDIST", "None");

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::FactoryMT::BuildInterfaces(const Teuchos::ParameterList& mtparams,
    std::vector<Teuchos::RCP<MORTAR::MortarInterface>>& interfaces, bool& poroslave,
    bool& poromaster) const
{
  int dim = DRT::Problem::Instance()->NDim();

  // start building interfaces
  if (Comm().MyPID() == 0)
  {
    std::cout << "Building contact interface(s)...............";
    fflush(stdout);
  }

  std::vector<DRT::Condition*> contactconditions(0);
  Discret().GetCondition("Mortar", contactconditions);

  // there must be more than one meshtying condition
  if ((int)contactconditions.size() < 2)
    dserror("ERROR: Not enough contact conditions in discretization");

  // find all pairs of matching meshtying conditions
  // there is a maximum of (conditions / 2) groups
  std::vector<int> foundgroups(0);
  int numgroupsfound = 0;

  // get nurbs information
  const bool nurbs = mtparams.get<bool>("NURBS");

  // maximum dof number in discretization
  // later we want to create NEW Lagrange multiplier degrees of
  // freedom, which of course must not overlap with displacement dofs
  int maxdof = Discret().DofRowMap()->MaxAllGID();

  for (int i = 0; i < (int)contactconditions.size(); ++i)
  {
    // initialize vector for current group of conditions and temp condition
    std::vector<DRT::Condition*> currentgroup(0);
    DRT::Condition* tempcond = NULL;

    // try to build meshtying group around this condition
    currentgroup.push_back(contactconditions[i]);
    const std::vector<int>* group1v = currentgroup[0]->Get<std::vector<int>>("Interface ID");
    if (!group1v) dserror("ERROR: Contact Conditions does not have value 'Interface ID'");
    int groupid1 = (*group1v)[0];
    bool foundit = false;

    for (int j = 0; j < (int)contactconditions.size(); ++j)
    {
      if (j == i) continue;  // do not detect contactconditions[i] again
      tempcond = contactconditions[j];
      const std::vector<int>* group2v = tempcond->Get<std::vector<int>>("Interface ID");
      if (!group2v) dserror("ERROR: Contact Conditions does not have value 'Interface ID'");
      int groupid2 = (*group2v)[0];
      if (groupid1 != groupid2) continue;  // not in the group
      foundit = true;                      // found a group entry
      currentgroup.push_back(tempcond);    // store it in currentgroup
    }

    // now we should have found a group of conds
    if (!foundit) dserror("ERROR: Cannot find matching contact condition for id %d", groupid1);

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
        dserror("ERROR: MtManager: Unknown contact side qualifier!");
      }
    }

    if (!hasslave) dserror("ERROR: Slave side missing in contact condition group!");
    if (!hasmaster) dserror("ERROR: Master side missing in contact condition group!");

    // find out which sides are initialized as Active
    std::vector<const std::string*> active((int)currentgroup.size());
    std::vector<bool> isactive((int)currentgroup.size());

    for (int j = 0; j < (int)sides.size(); ++j)
    {
      active[j] = currentgroup[j]->Get<std::string>("Initialization");
      if (*sides[j] == "Slave")
      {
        // slave sides must be initialized as "Active"
        if (*active[j] == "Active")
          isactive[j] = true;
        else if (*active[j] == "Inactive")
          dserror("ERROR: Slave side must be active for meshtying!");
        else
          dserror("ERROR: Unknown contact init qualifier!");
      }
      else if (*sides[j] == "Master")
      {
        // master sides must NOT be initialized as "Active" as this makes no sense
        if (*active[j] == "Active")
          dserror("ERROR: Master side cannot be active!");
        else if (*active[j] == "Inactive")
          isactive[j] = false;
        else
          dserror("ERROR: Unknown contact init qualifier!");
      }
      else
      {
        dserror("ERROR: MtManager: Unknown contact side qualifier!");
      }
    }

    // create an empty meshtying interface and store it in this Manager
    // (for structural meshtying we currently choose redundant master storage)
    INPAR::MORTAR::RedundantStorage redundant =
        DRT::INPUT::IntegralValue<INPAR::MORTAR::RedundantStorage>(
            mtparams.sublist("PARALLEL REDISTRIBUTION"), "REDUNDANT_STORAGE");
    //    if (redundant != INPAR::MORTAR::redundant_master)
    //      dserror("ERROR: MtManager: Meshtying requires redundant master storage");
    interfaces.push_back(
        MORTAR::MortarInterface::Create(groupid1, Comm(), dim, mtparams, redundant));

    // get it again
    Teuchos::RCP<MORTAR::MortarInterface> interface = interfaces[(int)interfaces.size() - 1];

    // note that the nodal ids are unique because they come from
    // one global problem discretization conatining all nodes of the
    // contact interface
    // We rely on this fact, therefore it is not possible to
    // do meshtying between two distinct discretizations here

    //-------------------------------------------------- process nodes
    for (int j = 0; j < (int)currentgroup.size(); ++j)
    {
      // get all nodes and add them
      const std::vector<int>* nodeids = currentgroup[j]->Nodes();
      if (!nodeids) dserror("ERROR: Condition does not have Node Ids");
      for (int k = 0; k < (int)(*nodeids).size(); ++k)
      {
        int gid = (*nodeids)[k];
        // do only nodes that I have in my discretization
        if (!Teuchos::rcp_dynamic_cast<DRT::Discretization>(discret_ptr_)->NodeColMap()->MyGID(gid))
          continue;
        DRT::Node* node = Discret().gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %", gid);

        // create MortarNode object
        Teuchos::RCP<MORTAR::MortarNode> mtnode =
            Teuchos::rcp(new MORTAR::MortarNode(node->Id(), node->X(), node->Owner(),
                Discret().NumDof(0, node), Discret().Dof(0, node), isslave[j]));
        //-------------------
        // get nurbs weight!
        if (nurbs)
        {
          MORTAR::UTILS::PrepareNURBSNode(node, mtnode);
        }

        // get edge and corner information:
        std::vector<DRT::Condition*> contactcornercond(0);
        Discret().GetCondition("mrtrcorner", contactcornercond);
        for (unsigned j = 0; j < contactcornercond.size(); j++)
        {
          if (contactcornercond.at(j)->ContainsNode(node->Id()))
          {
            mtnode->SetOnCorner() = true;
          }
        }
        std::vector<DRT::Condition*> contactedgecond(0);
        Discret().GetCondition("mrtredge", contactedgecond);
        for (unsigned j = 0; j < contactedgecond.size(); j++)
        {
          if (contactedgecond.at(j)->ContainsNode(node->Id()))
          {
            mtnode->SetOnEdge() = true;
          }
        }

        // Check, if this node (and, in case, which dofs) are in the contact symmetry condition
        std::vector<DRT::Condition*> contactSymconditions(0);
        Discret().GetCondition("mrtrsym", contactSymconditions);

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
    for (int j = 0; j < (int)currentgroup.size(); ++j)
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
      int lsize = (int)currele.size();
      int gsize = 0;
      Comm().SumAll(&lsize, &gsize, 1);


      std::map<int, Teuchos::RCP<DRT::Element>>::iterator fool;
      for (fool = currele.begin(); fool != currele.end(); ++fool)
      {
        Teuchos::RCP<DRT::Element> ele = fool->second;
        Teuchos::RCP<MORTAR::MortarElement> mtele =
            Teuchos::rcp(new MORTAR::MortarElement(ele->Id() + ggsize, ele->Owner(), ele->Shape(),
                ele->NumNode(), ele->NodeIds(), isslave[j], nurbs));
        //------------------------------------------------------------------
        // get knotvector, normal factor and zero-size information for nurbs
        if (nurbs)
        {
          MORTAR::UTILS::PrepareNURBSElement(
              *Teuchos::rcp_dynamic_cast<DRT::Discretization>(discret_ptr_), ele, mtele, dim);
        }

        interface->AddMortarElement(mtele);
      }  // for (fool=ele1.start(); fool != ele1.end(); ++fool)

      ggsize += gsize;  // update global element counter
    }

    //-------------------- finalize the meshtying interface construction
    interface->FillComplete(true, maxdof);

  }  // for (int i=0; i<(int)contactconditions.size(); ++i)
  if (Comm().MyPID() == 0) std::cout << "done!" << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::MtAbstractStrategy> MORTAR::STRATEGY::FactoryMT::BuildStrategy(
    const Teuchos::ParameterList& params, const bool& poroslave, const bool& poromaster,
    const int& dof_offset, std::vector<Teuchos::RCP<MORTAR::MortarInterface>>& interfaces) const
{
  const INPAR::CONTACT::SolvingStrategy stype =
      DRT::INPUT::IntegralValue<enum INPAR::CONTACT::SolvingStrategy>(params, "STRATEGY");
  Teuchos::RCP<CONTACT::AbstractStratDataContainer> data_ptr = Teuchos::null;

  return BuildStrategy(stype, params, poroslave, poromaster, dof_offset, interfaces,
      Discret().DofRowMap(), Discret().NodeRowMap(), Dim(), CommPtr(), data_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::MtAbstractStrategy> MORTAR::STRATEGY::FactoryMT::BuildStrategy(
    const INPAR::CONTACT::SolvingStrategy stype, const Teuchos::ParameterList& params,
    const bool& poroslave, const bool& poromaster, const int& dof_offset,
    std::vector<Teuchos::RCP<MORTAR::MortarInterface>>& interfaces, const Epetra_Map* dof_row_map,
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

  const double alphaf = STR::TIMINT::GetTimIntFactor();

  if (stype == INPAR::CONTACT::solution_lagmult)
  {
    strategy_ptr = Teuchos::rcp(new CONTACT::MtLagrangeStrategy(
        dof_row_map, node_row_map, params, interfaces, dim, comm_ptr, alphaf, dof_offset));
  }
  else if (stype == INPAR::CONTACT::solution_penalty or stype == INPAR::CONTACT::solution_uzawa)
    strategy_ptr = Teuchos::rcp(new CONTACT::MtPenaltyStrategy(
        dof_row_map, node_row_map, params, interfaces, dim, comm_ptr, alphaf, dof_offset));
  else
    dserror("ERROR: Unrecognized strategy");

  if (comm_ptr->MyPID() == 0) std::cout << "done!" << std::endl;
  return strategy_ptr;
}
