/*---------------------------------------------------------------------*/
/*!
\file meshtying_strategy_factory.cpp

\brief Factory to create the desired meshtying strategy.

\maintainer Alexander Seitz

\level 3

*/
/*---------------------------------------------------------------------*/

#include "meshtying_strategy_factory.H"

#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_node.H"

#include "../drt_structure_new/str_timint_basedataglobalstate.H"
#include "../drt_structure_xstructure/xstr_multi_discretization_wrapper.H"

#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/inpar_contact.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../linalg/linalg_utils.H"

#include "../drt_io/io_pstream.H"
#include "../drt_io/io.H"

#include <Teuchos_ParameterList.hpp>

#include "contact_utils.H"
#include "../drt_mortar/mortar_utils.H"

#include "contact_abstract_strategy.H"
#include "meshtying_abstract_strategy.H"
#include "meshtying_lagrange_strategy.H"
#include "meshtying_penalty_strategy.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::FactoryMT::Setup()
{
  CheckInit();
  MORTAR::STRATEGY::Factory::Setup();
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::FactoryMT::CheckDimension() const
{
  if (Dim() != 2 && Dim() != 3) dserror("ERROR: Contact problem must be 2D or 3D");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::FactoryMT::ReadAndCheckInput(Teuchos::ParameterList& mtparams) const
{
  // read parameter lists from DRT::Problem
  const Teuchos::ParameterList& mortar = DRT::Problem::Instance()->MortarCouplingParams();
  const Teuchos::ParameterList& meshtying = DRT::Problem::Instance()->ContactDynamicParams();
  const Teuchos::ParameterList& wearlist = DRT::Problem::Instance()->WearParams();

  // read Problem Type and Problem Dimension from DRT::Problem
  const PROBLEM_TYP problemtype = DRT::Problem::Instance()->ProblemType();
  int dim = DRT::Problem::Instance()->NDim();
  std::string distype = DRT::Problem::Instance()->SpatialApproximation();

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
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::RedundantStorage>(mortar, "REDUNDANT_STORAGE") ==
          INPAR::MORTAR::redundant_master and
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ParallelStrategy>(mortar, "PARALLEL_STRATEGY") !=
          INPAR::MORTAR::ghosting_redundant)
    dserror(
        "ERROR: Redundant storage only reasonable in combination with parallel strategy: "
        "ghosting_redundant !");

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::RedundantStorage>(mortar, "REDUNDANT_STORAGE") ==
          INPAR::MORTAR::redundant_all and
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ParallelStrategy>(mortar, "PARALLEL_STRATEGY") !=
          INPAR::MORTAR::ghosting_redundant)
    dserror(
        "ERROR: Redundant storage only reasonable in combination with parallel strategy: "
        "redundant_ghosting !");

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParallelStrategy>(mortar, "PARALLEL_STRATEGY") ==
          INPAR::MORTAR::roundrobinevaluate or
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ParallelStrategy>(mortar, "PARALLEL_STRATEGY") ==
          INPAR::MORTAR::roundrobinghost)
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

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(mortar, "PARALLEL_REDIST") ==
          INPAR::MORTAR::parredist_dynamic and
      onlymeshtying)
    dserror("ERROR: Dynamic parallel redistribution not possible for meshtying");

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(mortar, "PARALLEL_REDIST") !=
          INPAR::MORTAR::parredist_none &&
      mortar.get<int>("MIN_ELEPROC") < 0)
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
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(mortar, "PARALLEL_REDIST") !=
          INPAR::MORTAR::parredist_none)
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
  mtparams.setParameters(mortar);
  mtparams.setParameters(meshtying);
  mtparams.setParameters(wearlist);

  // *********************************************************************
  // predefined params for meshtying and contact
  // *********************************************************************
  if (meshtyingandcontact and !DRT::INPUT::IntegralValue<int>(meshtying, "DISCR_SMOOTHING"))
  {
    // set options for mortar coupling
    mtparams.set<std::string>("SEARCH_ALGORITHM", "Binarytree");
    mtparams.set<double>("SEARCH_PARAM", 0.3);
    mtparams.set<std::string>("SEARCH_USE_AUX_POS", "no");
    mtparams.set<std::string>("PARALLEL_REDIST", "static");
    mtparams.set<std::string>("LM_SHAPEFCN", "dual");
    mtparams.set<std::string>("REDUNDANT_STORAGE", "Master");
    mtparams.set<std::string>("SYSTEM", "condensed");
    mtparams.set<bool>("NURBS", false);
    mtparams.set<int>("NUMGP_PER_DIM", -1);
    mtparams.set<std::string>("STRATEGY", "LagrangianMultipliers");
    mtparams.set<std::string>("INTTYPE", "segments");
  }
  else if (meshtyingandcontact and DRT::INPUT::IntegralValue<int>(meshtying, "DISCR_SMOOTHING"))
  {
    // set options for mortar coupling
    mtparams.set<std::string>("SEARCH_ALGORITHM", "Binarytree");
    mtparams.set<double>("SEARCH_PARAM", 0.3);
    mtparams.set<std::string>("SEARCH_USE_AUX_POS", "no");
    mtparams.set<std::string>("PARALLEL_REDIST", "static");
    mtparams.set<std::string>("LM_SHAPEFCN", "dual");
    mtparams.set<std::string>("REDUNDANT_STORAGE", "Master");
    mtparams.set<std::string>("SYSTEM", "condensed");
    mtparams.set<bool>("NURBS", false);
  }
  // *********************************************************************
  // smooth interfaces
  // *********************************************************************
  // NURBS PROBLEM?
  if (distype == "Nurbs")
    mtparams.set<bool>("NURBS", true);
  else
    mtparams.set<bool>("NURBS", false);

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
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(mortar, "PARALLEL_REDIST") !=
          INPAR::MORTAR::parredist_none)
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

  mtparams.setName("CONTACT DYNAMIC / MORTAR COUPLING");

  // no parallel redistribution in the serial case
  if (Comm().NumProc() == 1) mtparams.set<std::string>("PARALLEL_REDIST", "None");

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
        DRT::INPUT::IntegralValue<INPAR::MORTAR::RedundantStorage>(mtparams, "REDUNDANT_STORAGE");
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
    interface->FillComplete(maxdof);

  }  // for (int i=0; i<(int)contactconditions.size(); ++i)
  if (Comm().MyPID() == 0) std::cout << "done!" << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::MORTAR::MortarInterface> MORTAR::STRATEGY::FactoryMT::CreateInterface(const int id,
    const Epetra_Comm& comm, const int dim, Teuchos::ParameterList& icparams,
    const bool selfcontact, const enum INPAR::MORTAR::RedundantStorage redundant,
    const Teuchos::RCP<std::pair<enum XFEM::FieldName,
        Teuchos::RCP<const DRT::DiscretizationInterface>>>& parent_dis_pair,
    Teuchos::RCP<MORTAR::IDataContainer> idata_ptr)
{
  INPAR::CONTACT::SolvingStrategy stype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(icparams, "STRATEGY");

  return CreateInterface(
      stype, id, comm, dim, icparams, selfcontact, redundant, parent_dis_pair, idata_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::MORTAR::MortarInterface> MORTAR::STRATEGY::FactoryMT::CreateInterface(
    const enum INPAR::CONTACT::SolvingStrategy stype, const int id, const Epetra_Comm& comm,
    const int dim, Teuchos::ParameterList& icparams, const bool selfcontact,
    const enum INPAR::MORTAR::RedundantStorage redundant,
    const Teuchos::RCP<std::pair<enum XFEM::FieldName,
        Teuchos::RCP<const DRT::DiscretizationInterface>>>& parent_dis_pair,
    Teuchos::RCP<MORTAR::IDataContainer> idata_ptr)
{
  Teuchos::RCP<MORTAR::MortarInterface> newinterface =
      Teuchos::rcp(new MORTAR::MortarInterface(idata_ptr, id, comm, dim, icparams, redundant));

  return newinterface;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::MtAbstractStrategy> MORTAR::STRATEGY::FactoryMT::BuildStrategy(
    const Teuchos::ParameterList& cparams, const bool& poroslave, const bool& poromaster,
    const int& dof_offset, std::vector<Teuchos::RCP<MORTAR::MortarInterface>>& interfaces) const
{
  const INPAR::CONTACT::SolvingStrategy stype =
      DRT::INPUT::IntegralValue<enum INPAR::CONTACT::SolvingStrategy>(cparams, "STRATEGY");
  Teuchos::RCP<CONTACT::AbstractStratDataContainer> data_ptr = Teuchos::null;

  return BuildStrategy(stype, cparams, poroslave, poromaster, dof_offset, interfaces,
      Discret().DofRowMap(), Discret().NodeRowMap(), Dim(), CommPtr(), data_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::MtAbstractStrategy> MORTAR::STRATEGY::FactoryMT::BuildStrategy(
    const INPAR::CONTACT::SolvingStrategy stype, const Teuchos::ParameterList& mtparams,
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

  double alphaf = 0.0;
  const Teuchos::ParameterList& sdynparams = DRT::Problem::Instance()->StructuralDynamicParams();
  if (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams, "DYNAMICTYP") ==
      INPAR::STR::dyna_genalpha)
    alphaf = sdynparams.sublist("GENALPHA").get<double>("ALPHA_F");
  if (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams, "DYNAMICTYP") ==
      INPAR::STR::dyna_gemm)
    alphaf = sdynparams.sublist("GEMM").get<double>("ALPHA_F");
  if (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams, "DYNAMICTYP") ==
      INPAR::STR::dyna_onesteptheta)
    alphaf = 1.0 - sdynparams.sublist("ONESTEPTHETA").get<double>("THETA");

  if (stype == INPAR::CONTACT::solution_lagmult)
  {
    strategy_ptr = Teuchos::rcp(new CONTACT::MtLagrangeStrategy(
        dof_row_map, node_row_map, mtparams, interfaces, dim, comm_ptr, alphaf, dof_offset));
  }
  else if (stype == INPAR::CONTACT::solution_penalty or stype == INPAR::CONTACT::solution_uzawa)
    strategy_ptr = Teuchos::rcp(new CONTACT::MtPenaltyStrategy(
        dof_row_map, node_row_map, mtparams, interfaces, dim, comm_ptr, alphaf, dof_offset));
  else
    dserror("ERROR: Unrecognized strategy");

  if (comm_ptr->MyPID() == 0) std::cout << "done!" << std::endl;
  return strategy_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::FactoryMT::BuildSearchTree(
    const std::vector<Teuchos::RCP<MORTAR::MortarInterface>>& interfaces) const
{
  // create binary search tree
  for (unsigned i = 0; i < interfaces.size(); ++i) interfaces[i]->CreateSearchTree();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::FactoryMT::ExtractParentDiscretization(
    const DRT::DiscretizationInterface& full_discret,
    const std::vector<DRT::Condition*>& given_ccgroup,
    XSTR::MultiDiscretizationWrapper::cXDisPair& parent_dis_pair) const
{
  const XSTR::MultiDiscretizationWrapper* dis_wrapper =
      dynamic_cast<const XSTR::MultiDiscretizationWrapper*>(&full_discret);

  // default case: do nothing and just return the input object
  if (dis_wrapper == NULL)
  {
    parent_dis_pair = std::make_pair(
        XFEM::structure, Teuchos::rcp<const DRT::DiscretizationInterface>(&full_discret, false));
    return;
  }

  bool coincide = false;
  std::vector<std::vector<DRT::Condition*>> curr_ccgroup;
  XSTR::MultiDiscretizationWrapper::XDisMap::const_iterator cit;
  for (cit = dis_wrapper->DiscretMap().begin(); cit != dis_wrapper->DiscretMap().end(); ++cit)
  {
    // continue if no contact conditions could be found in the current discretization
    if (CONTACT::UTILS::GetContactConditionGroups(curr_ccgroup, *(cit->second), false)) continue;
    /* loop over the condition grps of the current discretization and try
     * to find the one with same interface ID's */
    std::vector<std::vector<DRT::Condition*>>::const_iterator vv_cit;
    for (vv_cit = curr_ccgroup.begin(); vv_cit != curr_ccgroup.end(); ++vv_cit)
    {
      // continue if the sizes do not fit of the contact condition groups
      if (vv_cit->size() != given_ccgroup.size()) continue;
      // loop over the entries of the contact condition groups
      for (unsigned j = 0; j < given_ccgroup.size(); ++j)
      {
        const std::vector<int>* ggroupv = given_ccgroup[j]->Get<std::vector<int>>("Interface ID");
        if (ggroupv == NULL)
          dserror("ERROR: Given Contact Conditions do not have value 'Interface ID'");
        const std::vector<int>* cgroupv = (*vv_cit)[j]->Get<std::vector<int>>("Interface ID");
        if (cgroupv == NULL or cgroupv->size() != ggroupv->size())
        {
          coincide = false;
          break;
        }
        // compare the Interface ID's of the contact condition grp's
        if ((*cgroupv)[0] != (*ggroupv)[0])
        {
          coincide = false;
          break;
        }
        // has to be true for all of them
        coincide = true;
      }
      if (coincide) break;
    }
    // if the search was successful, return the desired pair
    if (coincide)
    {
      parent_dis_pair = std::make_pair(cit->first, cit->second.getConst());
      return;
    }
  }
  // if the search failed, throw an error
  IO::cout << "\n:::: Given Contact Condition Group ::::\n";
  for (unsigned i = 0; i < given_ccgroup.size(); ++i) IO::cout << *(given_ccgroup[i]) << "\n";
  dserror(
      "There is no wrapped discretization which belongs to the given "
      "contact condition group!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::STRATEGY::FactoryMT::PrintStrategyBanner(
    const enum INPAR::CONTACT::SolvingStrategy soltype)
{
  // some parameters
  const Teuchos::ParameterList& smortar = DRT::Problem::Instance()->MortarCouplingParams();
  const Teuchos::ParameterList& scontact = DRT::Problem::Instance()->ContactDynamicParams();
  INPAR::MORTAR::ShapeFcn shapefcn =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(smortar, "LM_SHAPEFCN");
  INPAR::CONTACT::SystemType systype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(scontact, "SYSTEM");
  INPAR::MORTAR::AlgorithmType algorithm =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(smortar, "ALGORITHM");
  bool smoothing = DRT::INPUT::IntegralValue<int>(scontact, "DISCR_SMOOTHING");
  bool nonSmoothGeometries = DRT::INPUT::IntegralValue<int>(scontact, "NONSMOOTH_GEOMETRIES");

  if (nonSmoothGeometries)
  {
    if (soltype == INPAR::CONTACT::solution_lagmult)
    {
      IO::cout << "================================================================\n";
      IO::cout << "===== Lagrange Multiplier Strategy =============================\n";
      IO::cout << "===== NONSMOOTH - GEOMETRIES ===================================\n";
      IO::cout << "================================================================\n\n";
    }
    else
      dserror("ERROR: Invalid system type for contact/meshtying interface smoothing");
  }
  else if (smoothing)
  {
    if (soltype == INPAR::CONTACT::solution_lagmult)
    {
      IO::cout << "================================================================\n";
      IO::cout << "========= !!! EXPERIMENTAL VERSION  !!!     ====================\n";
      IO::cout << "================================================================\n\n";
      IO::cout << "================================================================\n";
      IO::cout << "===== Interface smoothing approach with     ====================\n";
      IO::cout << "===== Standard Lagrange multiplier strategy ====================\n";
      IO::cout << "===== (Saddle point formulation) ===============================\n";
      IO::cout << "================================================================\n\n";
    }
    else if (soltype == INPAR::CONTACT::solution_penalty)
    {
      IO::cout << "================================================================\n";
      IO::cout << "========= !!! EXPERIMENTAL VERSION  !!!     ====================\n";
      IO::cout << "================================================================\n\n";
      IO::cout << "================================================================\n";
      IO::cout << "===== Interface smoothing approach with     ====================\n";
      IO::cout << "===== Standard Penalty strategy             ====================\n";
      IO::cout << "===== (Pure displacement formulation)===========================\n";
      IO::cout << "================================================================\n\n";
    }
    else
      dserror("ERROR: Invalid system type for contact/meshtying interface smoothing");
  }
  else
  {
    if (algorithm == INPAR::MORTAR::algorithm_mortar)
    {
      // saddle point formulation
      if (systype == INPAR::CONTACT::system_saddlepoint)
      {
        if (soltype == INPAR::CONTACT::solution_lagmult &&
            shapefcn == INPAR::MORTAR::shape_standard)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Standard Lagrange multiplier strategy ====================\n";
          IO::cout << "===== (Saddle point formulation) ===============================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_lagmult &&
                 shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Lagrange multiplier strategy ========================\n";
          IO::cout << "===== (Saddle point formulation) ===============================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_lagmult &&
                 shapefcn == INPAR::MORTAR::shape_petrovgalerkin)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Petrov-Galerkin Lagrange multiplier strategy =============\n";
          IO::cout << "===== (Saddle point formulation) ===============================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_penalty &&
                 shapefcn == INPAR::MORTAR::shape_standard)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Standard Penalty strategy ================================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_penalty &&
                 shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Penalty strategy ====================================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_uzawa &&
                 shapefcn == INPAR::MORTAR::shape_standard)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Uzawa Augmented Lagrange strategy ========================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_uzawa && shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Uzawa Augmented Lagrange strategy ===================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_combo)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Combination of different Solving Strategies ==============\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_augmented)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Augmented Lagrange strategy ==============================\n";
          IO::cout << "===== (Saddle point formulation) ===============================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_std_lagrange)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Standard Lagrange strategy ===============================\n";
          IO::cout << "===== Derived from the Augmented formulation ===================\n";
          IO::cout << "===== (Saddle point formulation) ===============================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_steepest_ascent)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Steepest Ascent strategy =================================\n";
          IO::cout << "===== (Condensed formulation) ==================================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_xcontact &&
                 shapefcn == INPAR::MORTAR::shape_standard)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Extended contact strategy ================================\n";
          IO::cout << "===== (Saddle point formulation) ===============================\n";
          IO::cout << "================================================================\n\n";
        }
        else
          dserror("ERROR: Invalid strategy or shape function type for contact/meshtying");
      }

      // condensed formulation
      else if (systype == INPAR::CONTACT::system_condensed ||
               systype == INPAR::CONTACT::system_condensed_lagmult)
      {
        if (soltype == INPAR::CONTACT::solution_lagmult && shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Lagrange multiplier strategy ========================\n";
          IO::cout << "===== (Condensed formulation) ==================================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_lagmult &&
                 shapefcn == INPAR::MORTAR::shape_standard &&
                 DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(smortar, "LM_QUAD") ==
                     INPAR::MORTAR::lagmult_const)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== const Lagrange multiplier strategy =======================\n";
          IO::cout << "===== (Condensed formulation) ==================================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_lagmult &&
                 shapefcn == INPAR::MORTAR::shape_petrovgalerkin)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Petrov-Galerkin Lagrange multiplier strategy =============\n";
          IO::cout << "===== (Condensed formulation) ==================================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_penalty &&
                 shapefcn == INPAR::MORTAR::shape_standard)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Standard Penalty strategy ================================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_penalty &&
                 shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Penalty strategy ====================================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_uzawa &&
                 shapefcn == INPAR::MORTAR::shape_standard)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Uzawa Augmented Lagrange strategy ========================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_uzawa && shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Uzawa Augmented Lagrange strategy ===================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else
          dserror("ERROR: Invalid strategy or shape function type for contact/meshtying");
      }
    }
    else if (algorithm == INPAR::MORTAR::algorithm_nts)
    {
      IO::cout << "================================================================\n";
      IO::cout << "===== Node-To-Segment approach =================================\n";
      IO::cout << "================================================================\n\n";
    }
    else if (algorithm == INPAR::MORTAR::algorithm_lts)
    {
      IO::cout << "================================================================\n";
      IO::cout << "===== Line-To-Segment approach =================================\n";
      IO::cout << "================================================================\n\n";
    }
    else if (algorithm == INPAR::MORTAR::algorithm_stl)
    {
      IO::cout << "================================================================\n";
      IO::cout << "===== Segment-To-Line approach =================================\n";
      IO::cout << "================================================================\n\n";
    }
    else if (algorithm == INPAR::MORTAR::algorithm_gpts)
    {
      IO::cout << "================================================================\n";
      IO::cout << "===== Gauss-Point-To-Segment approach ==========================\n";
      IO::cout << "================================================================\n\n";
    }
    // invalid system type
    else
      dserror("ERROR: Invalid system type for contact/meshtying");
  }
  return;
}
