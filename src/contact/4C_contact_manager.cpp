/*----------------------------------------------------------------------*/
/*! \file
\brief 4C implementation of main class to control all contact

\level 1


*/
/*----------------------------------------------------------------------*/

#include "4C_contact_manager.hpp"

#include "4C_contact_aug_interface.hpp"
#include "4C_contact_aug_strategy.hpp"
#include "4C_contact_element.hpp"
#include "4C_contact_friction_node.hpp"
#include "4C_contact_interface.hpp"
#include "4C_contact_lagrange_strategy.hpp"
#include "4C_contact_lagrange_strategy_poro.hpp"
#include "4C_contact_lagrange_strategy_tsi.hpp"
#include "4C_contact_lagrange_strategy_wear.hpp"
#include "4C_contact_nitsche_strategy.hpp"
#include "4C_contact_nitsche_strategy_fpi.hpp"
#include "4C_contact_nitsche_strategy_fsi.hpp"
#include "4C_contact_nitsche_strategy_poro.hpp"
#include "4C_contact_node.hpp"
#include "4C_contact_penalty_strategy.hpp"
#include "4C_contact_strategy_factory.hpp"
#include "4C_contact_utils.hpp"
#include "4C_contact_utils_parallel.hpp"
#include "4C_contact_wear_interface.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_inpar_mortar.hpp"
#include "4C_inpar_wear.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_mortar_utils.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 03/08|
 *----------------------------------------------------------------------*/
CONTACT::Manager::Manager(Core::FE::Discretization& discret, double alphaf)
    : Mortar::ManagerBase(), discret_(discret)
{
  // overwrite base class communicator
  comm_ = Teuchos::rcp(Discret().Comm().Clone());

  // create some local variables (later to be stored in strategy)
  const int dim = Global::Problem::Instance()->NDim();
  if (dim != 2 && dim != 3) FOUR_C_THROW("Contact problem must be 2D or 3D");
  std::vector<Teuchos::RCP<CONTACT::Interface>> interfaces;
  Teuchos::ParameterList contactParams;

  // read and check contact input parameters
  if (Comm().MyPID() == 0) std::cout << "Checking contact input parameters..........." << std::endl;

  read_and_check_input(contactParams);
  if (Comm().MyPID() == 0) std::cout << "done!" << std::endl;

  // check for fill_complete of discretization
  if (!Discret().Filled()) FOUR_C_THROW("discretization is not fillcomplete");

  // let's check for contact boundary conditions in the discretization and and detect groups of
  // matching conditions. For each group, create a contact interface and store it.
  if (Comm().MyPID() == 0) std::cout << "Building contact interface(s)..............." << std::endl;

  // Vector that contains solid-to-solid and beam-to-solid contact pairs
  std::vector<Core::Conditions::Condition*> beamandsolidcontactconditions(0);
  Discret().GetCondition("Contact", beamandsolidcontactconditions);

  // Vector that solely contains solid-to-solid contact pairs
  std::vector<Core::Conditions::Condition*> contactconditions(0);

  // Sort out beam-to-solid contact pairs, since these are treated in the beam3contact framework
  for (const auto& beamSolidCondition : beamandsolidcontactconditions)
  {
    if (beamSolidCondition->parameters().get<std::string>("Application") != "Beamtosolidcontact")
      contactconditions.push_back(beamSolidCondition);
  }

  // there must be more than one contact condition
  // unless we have a self contact problem!
  if (contactconditions.size() < 1) FOUR_C_THROW("Not enough contact conditions in discretization");
  if (contactconditions.size() == 1)
  {
    const std::string side = contactconditions[0]->parameters().get<std::string>("Side");
    FOUR_C_THROW_UNLESS(side == "Selfcontact", "Not enough contact conditions in discretization");
  }

  // find all pairs of matching contact conditions
  // there is a maximum of (conditions / 2) groups
  std::vector<int> foundgroups(0);
  int numgroupsfound = 0;

  // maximum dof number in discretization
  // later we want to create NEW Lagrange multiplier degrees of
  // freedom, which of course must not overlap with existing displacement dofs
  int maxdof = Discret().dof_row_map()->MaxAllGID();

  // get input parameters
  Inpar::CONTACT::SolvingStrategy stype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contactParams, "STRATEGY");
  Inpar::Wear::WearLaw wearLaw =
      Core::UTILS::IntegralValue<Inpar::Wear::WearLaw>(contactParams, "WEARLAW");
  Inpar::Wear::WearType wearType =
      Core::UTILS::IntegralValue<Inpar::Wear::WearType>(contactParams, "WEARTYPE");
  Inpar::CONTACT::ConstraintDirection constr_direction =
      Core::UTILS::IntegralValue<Inpar::CONTACT::ConstraintDirection>(
          contactParams, "CONSTRAINT_DIRECTIONS");
  Inpar::CONTACT::FrictionType frictionType =
      Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(contactParams, "FRICTION");
  Inpar::CONTACT::AdhesionType adhesionType =
      Core::UTILS::IntegralValue<Inpar::CONTACT::AdhesionType>(contactParams, "ADHESION");
  const bool nurbs = contactParams.get<bool>("NURBS");
  Inpar::Mortar::AlgorithmType algo =
      Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(contactParams, "ALGORITHM");

  bool friplus = false;
  if ((wearLaw != Inpar::Wear::wear_none) ||
      (contactParams.get<int>("PROBTYPE") == Inpar::CONTACT::tsi))
    friplus = true;

  // only for poro
  bool poromaster = false;
  bool poroslave = false;
  bool structmaster = false;
  bool structslave = false;
  int slavetype = -1;
  int mastertype = -1;  // 1 poro, 0 struct, -1 default
  bool isanyselfcontact = false;

  for (unsigned i = 0; i < contactconditions.size(); ++i)
  {
    // initialize vector for current group of conditions and temp condition
    std::vector<Core::Conditions::Condition*> currentgroup(0);
    Core::Conditions::Condition* tempcond = nullptr;

    // try to build contact group around this condition
    currentgroup.push_back(contactconditions[i]);
    const int groupid1 = currentgroup[0]->parameters().get<int>("Interface ID");

    // In case of MultiScale contact this is the id of the interface's constitutive contact law
    int contactconstitutivelawid = currentgroup[0]->parameters().get<int>("ConstitutiveLawID");

    bool foundit = false;

    // only one surface per group is ok for self contact
    const std::string& side = contactconditions[i]->parameters().get<std::string>("Side");
    if (side == "Selfcontact") foundit = true;

    for (unsigned j = 0; j < contactconditions.size(); ++j)
    {
      if (j == i) continue;  // do not detect contactconditions[i] again
      tempcond = contactconditions[j];
      const int groupid2 = tempcond->parameters().get<int>("Interface ID");
      if (groupid1 != groupid2) continue;  // not in the group
      foundit = true;                      // found a group entry
      currentgroup.push_back(tempcond);    // store it in currentgroup
    }

    // now we should have found a group of conds
    if (!foundit) FOUR_C_THROW("Cannot find matching contact condition for id %d", groupid1);

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
    std::vector<bool> isslave(0);
    std::vector<bool> isself(0);
    CONTACT::UTILS::GetMasterSlaveSideInfo(isslave, isself, currentgroup);
    for (const bool is : isself)
    {
      if (is)
      {
        isanyselfcontact = true;
        break;
      }
    }

    // find out which sides are initialized as In/Active and other initalization data
    std::vector<bool> isactive(currentgroup.size());
    bool Two_half_pass(false);
    bool Check_nonsmooth_selfcontactsurface(false);
    bool Searchele_AllProc(false);

    CONTACT::UTILS::GetInitializationInfo(Two_half_pass, Check_nonsmooth_selfcontactsurface,
        Searchele_AllProc, isactive, isslave, isself, currentgroup);

    // create interface local parameter list (copy)
    Teuchos::ParameterList icparams = contactParams;

    // find out if interface-specific coefficients of friction are given
    if (frictionType == Inpar::CONTACT::friction_tresca ||
        frictionType == Inpar::CONTACT::friction_coulomb ||
        frictionType == Inpar::CONTACT::friction_stick)
    {
      // read interface COFs
      std::vector<double> frcoeff(currentgroup.size());

      for (unsigned j = 0; j < currentgroup.size(); ++j)
        frcoeff[j] = currentgroup[j]->parameters().get<double>("FrCoeffOrBound");

      // check consistency of interface COFs
      for (unsigned j = 1; j < currentgroup.size(); ++j)
      {
        if (frcoeff[j] != frcoeff[0])
          FOUR_C_THROW("Inconsistency in friction coefficients of interface %i", groupid1);
      }

      // check for infeasible value of COF
      if (frcoeff[0] < 0.0) FOUR_C_THROW("Negative FrCoeff / FrBound on interface %i", groupid1);

      // add COF locally to contact parameter list of this interface
      if (frictionType == Inpar::CONTACT::friction_tresca)
      {
        icparams.setEntry("FRBOUND", static_cast<Teuchos::ParameterEntry>(frcoeff[0]));
        icparams.setEntry("FRCOEFF", static_cast<Teuchos::ParameterEntry>(-1.0));
      }
      else if (frictionType == Inpar::CONTACT::friction_coulomb)
      {
        icparams.setEntry("FRCOEFF", static_cast<Teuchos::ParameterEntry>(frcoeff[0]));
        icparams.setEntry("FRBOUND", static_cast<Teuchos::ParameterEntry>(-1.0));
      }
      // dummy values for FRCOEFF and FRBOUND have to be set,
      // since entries are accessed regardless of the friction law
      else if (frictionType == Inpar::CONTACT::friction_stick)
      {
        icparams.setEntry("FRCOEFF", static_cast<Teuchos::ParameterEntry>(-1.0));
        icparams.setEntry("FRBOUND", static_cast<Teuchos::ParameterEntry>(-1.0));
      }
    }

    // find out if interface-specific coefficients of adhesion are given
    if (adhesionType == Inpar::CONTACT::adhesion_bound)
    {
      // read interface COFs
      std::vector<double> ad_bound(currentgroup.size());
      for (unsigned j = 0; j < currentgroup.size(); ++j)
        ad_bound[j] = currentgroup[j]->parameters().get<double>("AdhesionBound");

      // check consistency of interface COFs
      for (unsigned j = 1; j < currentgroup.size(); ++j)
      {
        if (ad_bound[j] != ad_bound[0])
          FOUR_C_THROW("Inconsistency in adhesion bounds of interface %i", groupid1);
      }

      // check for infeasible value of COF
      if (ad_bound[0] < 0.0) FOUR_C_THROW("Negative adhesion bound on interface %i", groupid1);

      // add COF locally to contact parameter list of this interface
      icparams.setEntry("ADHESION_BOUND", static_cast<Teuchos::ParameterEntry>(ad_bound[0]));
    }

    // add information to contact parameter list of this interface
    icparams.set<bool>("Two_half_pass", Two_half_pass);
    icparams.set<bool>("Check_nonsmooth_selfcontactsurface", Check_nonsmooth_selfcontactsurface);
    icparams.set<bool>("Searchele_AllProc", Searchele_AllProc);

    // Safety check for interface storage redundancy in case of self contact
    Inpar::Mortar::ExtendGhosting redundant =
        Teuchos::getIntegralValue<Inpar::Mortar::ExtendGhosting>(
            icparams.sublist("PARALLEL REDISTRIBUTION"), "GHOSTING_STRATEGY");
    if (isanyselfcontact == true && redundant != Inpar::Mortar::ExtendGhosting::redundant_all)
      FOUR_C_THROW("Manager: Self contact requires fully redundant slave and master storage");

    // Use factory to create an empty interface and store it in this Manager.
    Teuchos::RCP<CONTACT::Interface> newinterface = STRATEGY::Factory::CreateInterface(groupid1,
        Comm(), dim, icparams, isself[0], Teuchos::null, Teuchos::null, contactconstitutivelawid);
    interfaces.push_back(newinterface);

    // Get the RCP to the last created interface
    Teuchos::RCP<CONTACT::Interface> interface = interfaces.back();

    // note that the nodal ids are unique because they come from
    // one global problem discretization containing all nodes of the
    // contact interface.
    // We rely on this fact, therefore it is not possible to
    // do contact between two distinct discretizations here.

    // collect all initially active nodes
    std::vector<int> initialactive;

    //-------------------------------------------------- process nodes
    for (unsigned j = 0; j < currentgroup.size(); ++j)
    {
      // get all nodes and add them
      const std::vector<int>* nodeids = currentgroup[j]->GetNodes();
      if (!nodeids) FOUR_C_THROW("Condition does not have Node Ids");
      for (unsigned k = 0; k < (*nodeids).size(); ++k)
      {
        int gid = (*nodeids)[k];
        // do only nodes that I have in my discretization
        if (!Discret().NodeColMap()->MyGID(gid)) continue;
        Core::Nodes::Node* node = Discret().gNode(gid);
        if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);

        // store global IDs of initially active nodes
        if (isactive[j]) initialactive.push_back(gid);

        // find out if this node is initial active on another Condition
        // and do NOT overwrite this status then!
        bool foundinitialactive = false;
        if (!isactive[j])
        {
          for (unsigned k = 0; k < initialactive.size(); ++k)
          {
            if (gid == initialactive[k])
            {
              foundinitialactive = true;
              break;
            }
          }
        }

        // create Node object or FriNode object in the frictional case
        // for the boolean variable initactive we use isactive[j]+foundinitialactive,
        // as this is true for BOTH initial active nodes found for the first time
        // and found for the second, third, ... time!
        if (frictionType != Inpar::CONTACT::friction_none)
        {
          Teuchos::RCP<CONTACT::FriNode> cnode =
              Teuchos::rcp(new CONTACT::FriNode(node->Id(), node->X(), node->Owner(),
                  Discret().Dof(0, node), isslave[j], isactive[j] + foundinitialactive, friplus));
          //-------------------
          // get nurbs weight!
          if (nurbs) Mortar::UTILS::prepare_nurbs_node(node, cnode);

          // get edge and corner information:
          std::vector<Core::Conditions::Condition*> contactcornercond(0);
          Discret().GetCondition("mrtrcorner", contactcornercond);
          for (unsigned j = 0; j < contactcornercond.size(); j++)
          {
            if (contactcornercond.at(j)->ContainsNode(node->Id()))
            {
              cnode->SetOnCorner() = true;
            }
          }
          std::vector<Core::Conditions::Condition*> contactedgecond(0);
          Discret().GetCondition("mrtredge", contactedgecond);
          for (unsigned j = 0; j < contactedgecond.size(); j++)
          {
            if (contactedgecond.at(j)->ContainsNode(node->Id()))
            {
              cnode->SetOnEdge() = true;
            }
          }

          // Check, if this node (and, in case, which dofs) are in the contact symmetry condition
          std::vector<Core::Conditions::Condition*> contactSymconditions(0);
          Discret().GetCondition("mrtrsym", contactSymconditions);

          for (unsigned l = 0; l < contactSymconditions.size(); l++)
            if (contactSymconditions.at(l)->ContainsNode(node->Id()))
            {
              const std::vector<int>& onoff =
                  contactSymconditions.at(l)->parameters().get<std::vector<int>>("onoff");
              for (unsigned k = 0; k < onoff.size(); k++)
                if (onoff.at(k) == 1) cnode->DbcDofs()[k] = true;
              if (stype == Inpar::CONTACT::solution_lagmult &&
                  constr_direction != Inpar::CONTACT::constr_xyz)
                FOUR_C_THROW(
                    "Contact symmetry with Lagrange multiplier method"
                    " only with contact constraints in xyz direction.\n"
                    "Set CONSTRAINT_DIRECTIONS to xyz in CONTACT input section");
            }

          // note that we do not have to worry about double entries
          // as the AddNode function can deal with this case!
          // the only problem would have occurred for the initial active nodes,
          // as their status could have been overwritten, but is prevented
          // by the "foundinitialactive" block above!
          interface->AddNode(cnode);
        }
        else
        {
          Teuchos::RCP<CONTACT::Node> cnode = Teuchos::rcp(new CONTACT::Node(node->Id(), node->X(),
              node->Owner(), Discret().Dof(0, node), isslave[j], isactive[j] + foundinitialactive));
          //-------------------
          // get nurbs weight!
          if (nurbs)
          {
            Mortar::UTILS::prepare_nurbs_node(node, cnode);
          }

          // get edge and corner information:
          std::vector<Core::Conditions::Condition*> contactcornercond(0);
          Discret().GetCondition("mrtrcorner", contactcornercond);
          for (unsigned j = 0; j < contactcornercond.size(); j++)
          {
            if (contactcornercond.at(j)->ContainsNode(node->Id()))
            {
              cnode->SetOnCorner() = true;
            }
          }
          std::vector<Core::Conditions::Condition*> contactedgecond(0);
          Discret().GetCondition("mrtredge", contactedgecond);
          for (unsigned j = 0; j < contactedgecond.size(); j++)
          {
            if (contactedgecond.at(j)->ContainsNode(node->Id()))
            {
              cnode->SetOnEdge() = true;
            }
          }


          // Check, if this node (and, in case, which dofs) are in the contact symmetry condition
          std::vector<Core::Conditions::Condition*> contactSymconditions(0);
          Discret().GetCondition("mrtrsym", contactSymconditions);

          for (unsigned l = 0; l < contactSymconditions.size(); l++)
            if (contactSymconditions.at(l)->ContainsNode(node->Id()))
            {
              const std::vector<int>& onoff =
                  contactSymconditions.at(l)->parameters().get<std::vector<int>>("onoff");
              for (unsigned k = 0; k < onoff.size(); k++)
                if (onoff.at(k) == 1) cnode->DbcDofs()[k] = true;
              if (stype == Inpar::CONTACT::solution_lagmult &&
                  constr_direction != Inpar::CONTACT::constr_xyz)
                FOUR_C_THROW(
                    "Contact symmetry with Lagrange multiplier method"
                    " only with contact constraints in xyz direction.\n"
                    "Set CONSTRAINT_DIRECTIONS to xyz in CONTACT input section");
            }

          // note that we do not have to worry about double entries
          // as the AddNode function can deal with this case!
          // the only problem would have occured for the initial active nodes,
          // as their status could have been overwritten, but is prevented
          // by the "foundinitialactive" block above!
          interface->AddNode(cnode);
        }
      }
    }

    //----------------------------------------------- process elements
    int ggsize = 0;
    for (unsigned j = 0; j < currentgroup.size(); ++j)
    {
      // get elements from condition j of current group
      std::map<int, Teuchos::RCP<Core::Elements::Element>>& currele = currentgroup[j]->Geometry();

      // elements in a boundary condition have a unique id
      // but ids are not unique among 2 distinct conditions
      // due to the way elements in conditions are build.
      // We therefore have to give the second, third,... set of elements
      // different ids. ids do not have to be continuous, we just add a large
      // enough number ggsize to all elements of cond2, cond3,... so they are
      // different from those in cond1!!!
      // note that elements in ele1/ele2 already are in column (overlapping) map
      int lsize = (int)currele.size();
      int gsize = 0;
      Comm().SumAll(&lsize, &gsize, 1);

      std::map<int, Teuchos::RCP<Core::Elements::Element>>::iterator fool;
      for (fool = currele.begin(); fool != currele.end(); ++fool)
      {
        Teuchos::RCP<Core::Elements::Element> ele = fool->second;
        Teuchos::RCP<CONTACT::Element> cele = Teuchos::rcp(new CONTACT::Element(ele->Id() + ggsize,
            ele->Owner(), ele->Shape(), ele->num_node(), ele->NodeIds(), isslave[j], nurbs));

        if ((contactParams.get<int>("PROBTYPE") == Inpar::CONTACT::poroelast ||
                contactParams.get<int>("PROBTYPE") == Inpar::CONTACT::poroscatra) &&
            algo != Inpar::Mortar::algorithm_gpts)
          set_poro_parent_element(slavetype, mastertype, cele, ele);

        if (algo == Inpar::Mortar::algorithm_gpts)
        {
          Teuchos::RCP<Core::Elements::FaceElement> faceele =
              Teuchos::rcp_dynamic_cast<Core::Elements::FaceElement>(ele, true);
          if (faceele == Teuchos::null) FOUR_C_THROW("Cast to FaceElement failed!");
          if (faceele->parent_element() == nullptr) FOUR_C_THROW("face parent does not exist");
          if (Discret().ElementColMap()->LID(faceele->parent_element()->Id()) == -1)
            FOUR_C_THROW("vol dis does not have parent ele");
          cele->set_parent_master_element(faceele->parent_element(), faceele->FaceParentNumber());
        }

        //------------------------------------------------------------------
        // get knotvector, normal factor and zero-size information for nurbs
        if (nurbs)
        {
          Mortar::UTILS::prepare_nurbs_element(discret, ele, cele, dim);
        }

        interface->add_element(cele);
      }  // for (fool=ele1.start(); fool != ele1.end(); ++fool)

      ggsize += gsize;  // update global element counter
    }

    /* Finalize the contact interface construction
     *
     * Always assign degrees of freedom here, because we need a valid column map for further contact
     * setup. This is an initial one time cost, that does not matter compared to the repeated
     * fill_complete calls due to dynamic redistribution.
     */
    if (CONTACT::UTILS::UseSafeRedistributeAndGhosting(contactParams))
    {
      /* Finalize parallel layout of maps. Note: Do not redistribute here.
       *
       * Since this is the initial setup, we don't need redistribution here, just a proper extension
       * of the interface ghosting.
       */
      interface->update_parallel_layout_and_data_structures(false, true, maxdof, 0.0);
    }
    else
      interface->fill_complete(true, maxdof);

    if ((contactParams.get<int>("PROBTYPE") == Inpar::CONTACT::poroelast ||
            contactParams.get<int>("PROBTYPE") == Inpar::CONTACT::poroscatra) &&
        algo != Inpar::Mortar::algorithm_gpts)
      find_poro_interface_types(
          poromaster, poroslave, structmaster, structslave, slavetype, mastertype);
  }
  if (Comm().MyPID() == 0) std::cout << "done!" << std::endl;

  //**********************************************************************
  // create the solver strategy object and pass all necessary data to it
  if (Comm().MyPID() == 0)
  {
    std::cout << "Building contact strategy object............";
    fflush(stdout);
  }

  // build the correct data container
  Teuchos::RCP<CONTACT::AbstractStratDataContainer> data_ptr =
      Teuchos::rcp(new CONTACT::AbstractStratDataContainer());

  // create LagrangeStrategyWear for wear as non-distinct quantity
  if (stype == Inpar::CONTACT::solution_lagmult && wearLaw != Inpar::Wear::wear_none &&
      (wearType == Inpar::Wear::wear_intstate || wearType == Inpar::Wear::wear_primvar))
  {
    strategy_ = Teuchos::rcp(new Wear::LagrangeStrategyWear(data_ptr, Discret().dof_row_map(),
        Discret().NodeRowMap(), contactParams, interfaces, dim, comm_, alphaf, maxdof));
  }
  else if (stype == Inpar::CONTACT::solution_lagmult)
  {
    if (contactParams.get<int>("PROBTYPE") == Inpar::CONTACT::poroelast ||
        contactParams.get<int>("PROBTYPE") == Inpar::CONTACT::poroscatra)
    {
      strategy_ = Teuchos::rcp(
          new LagrangeStrategyPoro(data_ptr, Discret().dof_row_map(), Discret().NodeRowMap(),
              contactParams, interfaces, dim, comm_, alphaf, maxdof, poroslave, poromaster));
    }
    else if (contactParams.get<int>("PROBTYPE") == Inpar::CONTACT::tsi)
    {
      strategy_ = Teuchos::rcp(new LagrangeStrategyTsi(data_ptr, Discret().dof_row_map(),
          Discret().NodeRowMap(), contactParams, interfaces, dim, comm_, alphaf, maxdof));
    }
    else
    {
      strategy_ = Teuchos::rcp(new LagrangeStrategy(data_ptr, Discret().dof_row_map(),
          Discret().NodeRowMap(), contactParams, interfaces, dim, comm_, alphaf, maxdof));
    }
  }
  else if (((stype == Inpar::CONTACT::solution_penalty ||
                stype == Inpar::CONTACT::solution_multiscale) &&
               algo != Inpar::Mortar::algorithm_gpts) ||
           stype == Inpar::CONTACT::solution_uzawa)
  {
    strategy_ = Teuchos::rcp(new PenaltyStrategy(data_ptr, Discret().dof_row_map(),
        Discret().NodeRowMap(), contactParams, interfaces, dim, comm_, alphaf, maxdof));
  }
  else if (algo == Inpar::Mortar::algorithm_gpts &&
           (stype == Inpar::CONTACT::solution_nitsche || stype == Inpar::CONTACT::solution_penalty))
  {
    if ((contactParams.get<int>("PROBTYPE") == Inpar::CONTACT::poroelast ||
            contactParams.get<int>("PROBTYPE") == Inpar::CONTACT::poroscatra) &&
        stype == Inpar::CONTACT::solution_nitsche)
    {
      strategy_ = Teuchos::rcp(new NitscheStrategyPoro(data_ptr, Discret().dof_row_map(),
          Discret().NodeRowMap(), contactParams, interfaces, dim, comm_, alphaf, maxdof));
    }
    else if (contactParams.get<int>("PROBTYPE") == Inpar::CONTACT::fsi &&
             stype == Inpar::CONTACT::solution_nitsche)
    {
      strategy_ = Teuchos::rcp(new NitscheStrategyFsi(data_ptr, Discret().dof_row_map(),
          Discret().NodeRowMap(), contactParams, interfaces, dim, comm_, alphaf, maxdof));
    }
    else if (contactParams.get<int>("PROBTYPE") == Inpar::CONTACT::fpi &&
             stype == Inpar::CONTACT::solution_nitsche)
    {
      strategy_ = Teuchos::rcp(new NitscheStrategyFpi(data_ptr, Discret().dof_row_map(),
          Discret().NodeRowMap(), contactParams, interfaces, dim, comm_, alphaf, maxdof));
    }
    else
    {
      strategy_ = Teuchos::rcp(new NitscheStrategy(data_ptr, Discret().dof_row_map(),
          Discret().NodeRowMap(), contactParams, interfaces, dim, comm_, alphaf, maxdof));
    }
  }
  else if (stype == Inpar::CONTACT::solution_augmented)
  {
    FOUR_C_THROW(
        "The augmented contact formulation is no longer supported in the"
        " old structural time integrator!");
  }
  else
  {
    FOUR_C_THROW("Unrecognized contact strategy");
  }

  dynamic_cast<CONTACT::AbstractStrategy&>(*strategy_).setup(false, true);

  if (Comm().MyPID() == 0) std::cout << "done!" << std::endl;
  //**********************************************************************

  // print friction information of interfaces
  if (Comm().MyPID() == 0)
  {
    for (unsigned i = 0; i < interfaces.size(); ++i)
    {
      double checkfrcoeff = 0.0;
      if (frictionType == Inpar::CONTACT::friction_tresca)
      {
        checkfrcoeff = interfaces[i]->interface_params().get<double>("FRBOUND");
        std::cout << std::endl << "Interface         " << i + 1 << std::endl;
        std::cout << "FrBound (Tresca)  " << checkfrcoeff << std::endl;
      }
      else if (frictionType == Inpar::CONTACT::friction_coulomb)
      {
        checkfrcoeff = interfaces[i]->interface_params().get<double>("FRCOEFF");
        std::cout << std::endl << "Interface         " << i + 1 << std::endl;
        std::cout << "FrCoeff (Coulomb) " << checkfrcoeff << std::endl;
      }
    }
  }

  // print initial parallel redistribution
  if (Comm().MyPID() == 0 && Comm().NumProc() > 1)
    std::cout << "\nInitial parallel distribution of all contact interfaces:" << std::endl;
  for (auto& interface : interfaces) interface->print_parallel_distribution();

  // create binary search tree
  for (auto& interface : interfaces) interface->CreateSearchTree();

  return;
}


/*----------------------------------------------------------------------*
 |  read and check input parameters (public)                  popp 04/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Manager::read_and_check_input(Teuchos::ParameterList& cparams)
{
  // read parameter lists from Global::Problem
  const Teuchos::ParameterList& mortar = Global::Problem::Instance()->mortar_coupling_params();
  const Teuchos::ParameterList& contact = Global::Problem::Instance()->contact_dynamic_params();
  const Teuchos::ParameterList& wearlist = Global::Problem::Instance()->WearParams();
  const Teuchos::ParameterList& tsic = Global::Problem::Instance()->TSIContactParams();
  const Teuchos::ParameterList& stru = Global::Problem::Instance()->structural_dynamic_params();

  // read Problem Type and Problem Dimension from Global::Problem
  const Core::ProblemType problemtype = Global::Problem::Instance()->GetProblemType();
  Core::FE::ShapeFunctionType distype = Global::Problem::Instance()->spatial_approximation_type();
  const int dim = Global::Problem::Instance()->NDim();

  // in case just System type system_condensed_lagmult
  if (Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(contact, "SYSTEM") ==
      Inpar::CONTACT::system_condensed_lagmult)
    FOUR_C_THROW(
        "For Contact anyway just the lagrange multiplier can be condensed, choose SYSTEM = "
        "Condensed.");

  // *********************************************************************
  // invalid parallel strategies
  // *********************************************************************
  const Teuchos::ParameterList& mortarParallelRedistParams =
      mortar.sublist("PARALLEL REDISTRIBUTION");

  if (Teuchos::getIntegralValue<Inpar::Mortar::ParallelRedist>(mortarParallelRedistParams,
          "PARALLEL_REDIST") != Inpar::Mortar::ParallelRedist::redist_none &&
      mortarParallelRedistParams.get<int>("MIN_ELEPROC") < 0)
    FOUR_C_THROW(
        "Minimum number of elements per processor for parallel redistribution must be >= 0");

  if (Teuchos::getIntegralValue<Inpar::Mortar::ParallelRedist>(mortarParallelRedistParams,
          "PARALLEL_REDIST") == Inpar::Mortar::ParallelRedist::redist_dynamic &&
      mortarParallelRedistParams.get<double>("MAX_BALANCE_EVAL_TIME") < 1.0)
    FOUR_C_THROW(
        "Maximum allowed value of load balance for dynamic parallel redistribution must be "
        ">= 1.0");

  if (problemtype == Core::ProblemType::tsi &&
      Teuchos::getIntegralValue<Inpar::Mortar::ParallelRedist>(mortarParallelRedistParams,
          "PARALLEL_REDIST") != Inpar::Mortar::ParallelRedist::redist_none)
    FOUR_C_THROW("Parallel redistribution not yet implemented for TSI problems");

  // *********************************************************************
  // adhesive contact
  // *********************************************************************
  if (Core::UTILS::IntegralValue<Inpar::CONTACT::AdhesionType>(contact, "ADHESION") !=
          Inpar::CONTACT::adhesion_none and
      Core::UTILS::IntegralValue<Inpar::Wear::WearLaw>(wearlist, "WEARLAW") !=
          Inpar::Wear::wear_none)
    FOUR_C_THROW("Adhesion combined with wear not yet tested!");

  if (Core::UTILS::IntegralValue<Inpar::CONTACT::AdhesionType>(contact, "ADHESION") !=
          Inpar::CONTACT::adhesion_none and
      Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(contact, "FRICTION") !=
          Inpar::CONTACT::friction_none)
    FOUR_C_THROW("Adhesion combined with friction not yet tested!");

  // *********************************************************************
  // generally invalid combinations (nts/mortar)
  // *********************************************************************
  if ((Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
              Inpar::CONTACT::solution_penalty ||
          Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
              Inpar::CONTACT::solution_nitsche) &&
      contact.get<double>("PENALTYPARAM") <= 0.0)
    FOUR_C_THROW("Penalty parameter eps = 0, must be greater than 0");

  if ((Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
              Inpar::CONTACT::solution_penalty ||
          Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
              Inpar::CONTACT::solution_nitsche) &&
      Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(contact, "FRICTION") !=
          Inpar::CONTACT::friction_none &&
      contact.get<double>("PENALTYPARAMTAN") <= 0.0)
    FOUR_C_THROW("Tangential penalty parameter eps = 0, must be greater than 0");

  if (Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
          Inpar::CONTACT::solution_uzawa &&
      contact.get<double>("PENALTYPARAM") <= 0.0)
    FOUR_C_THROW("Penalty parameter eps = 0, must be greater than 0");

  if (Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
          Inpar::CONTACT::solution_uzawa &&
      Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(contact, "FRICTION") !=
          Inpar::CONTACT::friction_none &&
      contact.get<double>("PENALTYPARAMTAN") <= 0.0)
    FOUR_C_THROW("Tangential penalty parameter eps = 0, must be greater than 0");

  if (Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
          Inpar::CONTACT::solution_uzawa &&
      contact.get<int>("UZAWAMAXSTEPS") < 2)
    FOUR_C_THROW("Maximum number of Uzawa / Augmentation steps must be at least 2");

  if (Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
          Inpar::CONTACT::solution_uzawa &&
      contact.get<double>("UZAWACONSTRTOL") <= 0.0)
    FOUR_C_THROW("Constraint tolerance for Uzawa / Augmentation scheme must be greater than 0");

  if (Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(contact, "FRICTION") !=
          Inpar::CONTACT::friction_none &&
      contact.get<double>("SEMI_SMOOTH_CT") == 0.0)
    FOUR_C_THROW("Parameter ct = 0, must be greater than 0 for frictional contact");

  if (Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
          Inpar::CONTACT::solution_augmented &&
      contact.get<double>("SEMI_SMOOTH_CN") <= 0.0)
    FOUR_C_THROW("Regularization parameter cn, must be greater than 0 for contact problems");

  if (Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(contact, "FRICTION") ==
          Inpar::CONTACT::friction_tresca &&
      dim == 3 &&
      Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") !=
          Inpar::CONTACT::solution_nitsche)
    FOUR_C_THROW("3D frictional contact with Tresca's law not yet implemented");

  if (Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(contact, "FRICTION") !=
          Inpar::CONTACT::friction_none &&
      Core::UTILS::IntegralValue<int>(contact, "SEMI_SMOOTH_NEWTON") != 1 && dim == 3)
    FOUR_C_THROW("3D frictional contact only implemented with Semi-smooth Newton");

  if (Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
          Inpar::CONTACT::solution_augmented &&
      Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(contact, "FRICTION") !=
          Inpar::CONTACT::friction_none)
    FOUR_C_THROW(
        "Frictional contact is for the augmented Lagrange formulation not yet implemented!");

  if (Core::UTILS::IntegralValue<int>(mortar, "CROSSPOINTS") == true && dim == 3)
    FOUR_C_THROW("Crosspoints / edge node modification not yet implemented for 3D");

  if (Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(contact, "FRICTION") ==
          Inpar::CONTACT::friction_tresca &&
      Core::UTILS::IntegralValue<int>(contact, "FRLESS_FIRST") == true)
    // Hopefully coming soon, when Coulomb and Tresca are combined. Until then, throw error.
    FOUR_C_THROW("Frictionless first contact step with Tresca's law not yet implemented");

  if (Core::UTILS::IntegralValue<Inpar::CONTACT::Regularization>(
          contact, "CONTACT_REGULARIZATION") != Inpar::CONTACT::reg_none &&
      Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") !=
          Inpar::CONTACT::solution_lagmult)
    FOUR_C_THROW(
        "Regularized Contact just available for Dual Mortar Contact with Lagrangean "
        "Multiplier!");

  if (Core::UTILS::IntegralValue<Inpar::CONTACT::Regularization>(
          contact, "CONTACT_REGULARIZATION") != Inpar::CONTACT::reg_none &&
      Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(contact, "FRICTION") !=
          Inpar::CONTACT::friction_none)
    FOUR_C_THROW("Regularized Contact for contact with friction not implemented yet!");

  // *********************************************************************
  // warnings
  // *********************************************************************
  if (mortar.get<double>("SEARCH_PARAM") == 0.0 && Comm().MyPID() == 0)
    std::cout << ("Warning: Contact search called without inflation of bounding volumes\n")
              << std::endl;

  if (Core::UTILS::IntegralValue<Inpar::Wear::WearSide>(wearlist, "WEAR_SIDE") !=
      Inpar::Wear::wear_slave)
    std::cout << ("\n \n Warning: Contact with both-sided wear is still experimental !")
              << std::endl;


  // *********************************************************************
  //                       MORTAR-SPECIFIC CHECKS
  // *********************************************************************
  if (Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(mortar, "ALGORITHM") ==
      Inpar::Mortar::algorithm_mortar)
  {
    // *********************************************************************
    // invalid parameter combinations
    // *********************************************************************
    if (Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") !=
            Inpar::CONTACT::solution_lagmult &&
        Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
            Inpar::Mortar::shape_petrovgalerkin)
      FOUR_C_THROW("Petrov-Galerkin approach for LM only with Lagrange multiplier strategy");

    if (Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
            Inpar::CONTACT::solution_lagmult &&
        (Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
                Inpar::Mortar::shape_standard &&
            Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(mortar, "LM_QUAD") !=
                Inpar::Mortar::lagmult_const) &&
        Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(contact, "SYSTEM") ==
            Inpar::CONTACT::system_condensed)
      FOUR_C_THROW("Condensation of linear system only possible for dual Lagrange multipliers");

    if (Core::UTILS::IntegralValue<Inpar::Mortar::ConsistentDualType>(
            mortar, "LM_DUAL_CONSISTENT") != Inpar::Mortar::consistent_none &&
        Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") !=
            Inpar::CONTACT::solution_lagmult &&
        Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(mortar, "LM_SHAPEFCN") !=
            Inpar::Mortar::shape_standard)
      FOUR_C_THROW(
          "Consistent dual shape functions in boundary elements only for Lagrange "
          "multiplier strategy.");

    if (Core::UTILS::IntegralValue<Inpar::Mortar::ConsistentDualType>(
            mortar, "LM_DUAL_CONSISTENT") != Inpar::Mortar::consistent_none &&
        Core::UTILS::IntegralValue<Inpar::Mortar::IntType>(mortar, "INTTYPE") ==
            Inpar::Mortar::inttype_elements &&
        (Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
            Inpar::Mortar::shape_dual))
      FOUR_C_THROW(
          "Consistent dual shape functions in boundary elements not for purely "
          "element-based integration.");

    if (Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
            Inpar::CONTACT::solution_nitsche &&
        Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(mortar, "ALGORITHM") !=
            Inpar::Mortar::algorithm_gpts)
      FOUR_C_THROW("Nitsche contact only with GPTS algorithm.");


    // *********************************************************************
    // not (yet) implemented combinations
    // *********************************************************************

    if (Core::UTILS::IntegralValue<int>(mortar, "CROSSPOINTS") == true &&
        Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(mortar, "LM_QUAD") ==
            Inpar::Mortar::lagmult_lin)
      FOUR_C_THROW("Crosspoints and linear LM interpolation for quadratic FE not yet compatible");

    // check for self contact
    bool self = false;
    {
      std::vector<Core::Conditions::Condition*> contactCondition(0);
      Discret().GetCondition("Mortar", contactCondition);

      for (const auto& condition : contactCondition)
      {
        const std::string side = condition->parameters().get<std::string>("Side");
        if (side == "Selfcontact") self = true;
      }
    }

    if (self == true &&
        Teuchos::getIntegralValue<Inpar::Mortar::ParallelRedist>(mortarParallelRedistParams,
            "PARALLEL_REDIST") != Inpar::Mortar::ParallelRedist::redist_none)
      FOUR_C_THROW("Self contact and parallel redistribution not yet compatible");

    if (Core::UTILS::IntegralValue<int>(contact, "INITCONTACTBYGAP") == true &&
        contact.get<double>("INITCONTACTGAPVALUE") == 0.0)
      FOUR_C_THROW(
          "For initialization of init contact with gap, the INITCONTACTGAPVALUE is needed.");

    if (Core::UTILS::IntegralValue<Inpar::Wear::WearLaw>(wearlist, "WEARLAW") !=
            Inpar::Wear::wear_none &&
        Core::UTILS::IntegralValue<int>(contact, "FRLESS_FIRST") == true)
      FOUR_C_THROW("Frictionless first contact step with wear not yet implemented");

    if (problemtype != Core::ProblemType::ehl &&
        Core::UTILS::IntegralValue<int>(contact, "REGULARIZED_NORMAL_CONTACT") == true)
      FOUR_C_THROW("Regularized normal contact only implemented for EHL");

    // *********************************************************************
    // Augmented Lagrangian strategy
    // *********************************************************************
    if (Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
        Inpar::CONTACT::solution_augmented)
    {
      FOUR_C_THROW("No longer supported!");
    }

    // *********************************************************************
    // thermal-structure-interaction contact
    // *********************************************************************
    if (problemtype == Core::ProblemType::tsi &&
        Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
            Inpar::Mortar::shape_standard)
      FOUR_C_THROW("Thermal contact only for dual shape functions");

    if (problemtype == Core::ProblemType::tsi &&
        Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(contact, "SYSTEM") !=
            Inpar::CONTACT::system_condensed)
      FOUR_C_THROW("Thermal contact only for dual shape functions with condensed system");

    // no nodal scaling in for thermal-structure-interaction
    if (problemtype == Core::ProblemType::tsi &&
        tsic.get<double>("TEMP_DAMAGE") <= tsic.get<double>("TEMP_REF"))
      FOUR_C_THROW("damage temperature must be greater than reference temperature");

    // *********************************************************************
    // contact with wear
    // *********************************************************************
    if (Core::UTILS::IntegralValue<Inpar::Wear::WearLaw>(wearlist, "WEARLAW") ==
            Inpar::Wear::wear_none &&
        wearlist.get<double>("WEARCOEFF") != 0.0)
      FOUR_C_THROW("Wear coefficient only necessary in the context of wear.");

    if (problemtype == Core::ProblemType::structure and
        Core::UTILS::IntegralValue<Inpar::Wear::WearLaw>(wearlist, "WEARLAW") !=
            Inpar::Wear::wear_none and
        Core::UTILS::IntegralValue<Inpar::Wear::WearTimInt>(wearlist, "WEARTIMINT") !=
            Inpar::Wear::wear_expl)
      FOUR_C_THROW(
          "Wear calculation for pure structure problems only with explicit internal state "
          "variable approach reasonable!");

    if (Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(contact, "FRICTION") ==
            Inpar::CONTACT::friction_none &&
        Core::UTILS::IntegralValue<Inpar::Wear::WearLaw>(wearlist, "WEARLAW") !=
            Inpar::Wear::wear_none)
      FOUR_C_THROW("Wear models only applicable to frictional contact.");

    if (Core::UTILS::IntegralValue<Inpar::Wear::WearLaw>(wearlist, "WEARLAW") !=
            Inpar::Wear::wear_none &&
        wearlist.get<double>("WEARCOEFF") <= 0.0)
      FOUR_C_THROW("No valid wear coefficient provided, must be equal or greater 0.0");

    //    if (Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact,"STRATEGY") !=
    //    Inpar::CONTACT::solution_lagmult
    //        && Core::UTILS::IntegralValue<Inpar::Wear::WearLaw>(wearlist, "WEARLAW")     !=
    //        Inpar::Wear::wear_none)
    //      FOUR_C_THROW("Wear model only applicable in combination with Lagrange multiplier
    //      strategy.");

    if (Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(contact, "FRICTION") ==
            Inpar::CONTACT::friction_tresca &&
        Core::UTILS::IntegralValue<Inpar::Wear::WearLaw>(wearlist, "WEARLAW") !=
            Inpar::Wear::wear_none)
      FOUR_C_THROW("Wear only for Coulomb friction!");

    // *********************************************************************
    // 3D quadratic mortar (choice of interpolation and testing fcts.)
    // *********************************************************************
    if (Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(mortar, "LM_QUAD") ==
            Inpar::Mortar::lagmult_pwlin &&
        Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
            Inpar::Mortar::shape_dual)
      FOUR_C_THROW(
          "No piecewise linear approach (for LM) implemented for quadratic contact with "
          "DUAL shape fct.");

    // *********************************************************************
    // poroelastic contact
    // *********************************************************************
    if (problemtype == Core::ProblemType::poroelast ||
        problemtype == Core::ProblemType::poroscatra || problemtype == Core::ProblemType::fpsi ||
        problemtype == Core::ProblemType::fpsi_xfem)
    {
      const Teuchos::ParameterList& porodyn =
          Global::Problem::Instance()->poroelast_dynamic_params();
      if ((Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(mortar, "LM_SHAPEFCN") !=
                  Inpar::Mortar::shape_dual &&
              Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(mortar, "LM_SHAPEFCN") !=
                  Inpar::Mortar::shape_petrovgalerkin) &&
          Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
              Inpar::CONTACT::solution_lagmult)
        FOUR_C_THROW("POROCONTACT: Only dual and petrovgalerkin shape functions implemented yet!");

      if (Teuchos::getIntegralValue<Inpar::Mortar::ParallelRedist>(mortarParallelRedistParams,
              "PARALLEL_REDIST") != Inpar::Mortar::ParallelRedist::redist_none &&
          Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
              Inpar::CONTACT::solution_lagmult)
        FOUR_C_THROW(
            "POROCONTACT: Parallel Redistribution not implemented yet!");  // Since we use Pointers
                                                                           // to Parent Elements,
                                                                           // which are not copied
                                                                           // to other procs!

      if (Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") !=
              Inpar::CONTACT::solution_lagmult &&
          Core::UTILS::IntegralValue<int>(porodyn, "CONTACTNOPEN"))
        FOUR_C_THROW("POROCONTACT: Use Lagrangean Strategy for poro contact!");

      if (Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(contact, "FRICTION") !=
              Inpar::CONTACT::friction_none &&
          Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
              Inpar::CONTACT::solution_lagmult)
        FOUR_C_THROW("POROCONTACT: is_friction for poro contact not implemented!");

      if (Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(contact, "SYSTEM") !=
              Inpar::CONTACT::system_condensed &&
          Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
              Inpar::CONTACT::solution_lagmult)
        FOUR_C_THROW("POROCONTACT: System has to be condensed for poro contact!");

      if ((dim != 3) && (dim != 2))
      {
        const Teuchos::ParameterList& porodyn =
            Global::Problem::Instance()->poroelast_dynamic_params();
        if (Core::UTILS::IntegralValue<int>(porodyn, "CONTACTNOPEN"))
          FOUR_C_THROW("POROCONTACT: PoroContact with no penetration just tested for 3d (and 2d)!");
      }
    }

    // *********************************************************************
    // element-based vs. segment-based mortar integration
    // *********************************************************************
    Inpar::Mortar::IntType inttype =
        Core::UTILS::IntegralValue<Inpar::Mortar::IntType>(mortar, "INTTYPE");

    if (inttype == Inpar::Mortar::inttype_elements && mortar.get<int>("NUMGP_PER_DIM") <= 0)
      FOUR_C_THROW("Invalid Gauss point number NUMGP_PER_DIM for element-based integration.");

    if (inttype == Inpar::Mortar::inttype_elements_BS && mortar.get<int>("NUMGP_PER_DIM") <= 0)
      FOUR_C_THROW(
          "Invalid Gauss point number NUMGP_PER_DIM for element-based integration with "
          "boundary segmentation."
          "\nPlease note that the value you have to provide only applies to the element-based "
          "integration"
          "\ndomain, while pre-defined default values will be used in the segment-based boundary "
          "domain.");

    if ((inttype == Inpar::Mortar::inttype_elements ||
            inttype == Inpar::Mortar::inttype_elements_BS) &&
        mortar.get<int>("NUMGP_PER_DIM") <= 1)
      FOUR_C_THROW("Invalid Gauss point number NUMGP_PER_DIM for element-based integration.");
  }  // END MORTAR CHECKS

  // *********************************************************************
  //                       NTS-SPECIFIC CHECKS
  // *********************************************************************
  else if (Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(mortar, "ALGORITHM") ==
           Inpar::Mortar::algorithm_nts)
  {
    if (problemtype == Core::ProblemType::poroelast or problemtype == Core::ProblemType::fpsi or
        problemtype == Core::ProblemType::tsi)
      FOUR_C_THROW("NTS only for problem type: structure");
  }  // END NTS CHECKS

  // *********************************************************************
  //                       GPTS-SPECIFIC CHECKS
  // *********************************************************************
  else if (Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(mortar, "ALGORITHM") ==
           Inpar::Mortar::algorithm_gpts)
  {
    const_cast<Teuchos::ParameterList&>(Global::Problem::Instance()->contact_dynamic_params())
        .set("SYSTEM", "none");

    if (contact.get<double>("PENALTYPARAM") <= 0.0)
      FOUR_C_THROW("Penalty parameter eps = 0, must be greater than 0");

    if (problemtype != Core::ProblemType::structure &&
        problemtype != Core::ProblemType::poroelast && problemtype != Core::ProblemType::fsi_xfem &&
        problemtype != Core::ProblemType::fpsi_xfem)
      FOUR_C_THROW(
          "GPTS algorithm only tested for structural, FSI-CutFEM, FPSI-CutFEM, and "
          "poroelastic problems");

    if (Core::UTILS::IntegralValue<Inpar::Wear::WearLaw>(wearlist, "WEARLAW") !=
        Inpar::Wear::wear_none)
      FOUR_C_THROW("GPTS algorithm not implemented for wear");

  }  // END GPTS CHECKS

  // *********************************************************************
  // store contents of BOTH ParameterLists in local parameter list
  // *********************************************************************
  cparams.setParameters(mortar);
  cparams.setParameters(contact);
  cparams.setParameters(wearlist);
  cparams.setParameters(tsic);
  if (problemtype == Core::ProblemType::tsi)
    cparams.set<double>(
        "TIMESTEP", Global::Problem::Instance()->TSIDynamicParams().get<double>("TIMESTEP"));
  else if (problemtype != Core::ProblemType::structure)
  {
    // rauch 01/16
    if (Comm().MyPID() == 0)
      std::cout << "\n \n  Warning: CONTACT::Manager::read_and_check_input() reads TIMESTEP = "
                << stru.get<double>("TIMESTEP") << " from --STRUCTURAL DYNAMIC \n"
                << std::endl;
    cparams.set<double>("TIMESTEP", stru.get<double>("TIMESTEP"));
  }
  else
    cparams.set<double>("TIMESTEP", stru.get<double>("TIMESTEP"));

  // *********************************************************************
  // NURBS contact
  // *********************************************************************
  switch (distype)
  {
    case Core::FE::ShapeFunctionType::nurbs:
    {
      cparams.set<bool>("NURBS", true);
      break;
    }
    default:
    {
      cparams.set<bool>("NURBS", false);
      break;
    }
  }

  // *********************************************************************
  cparams.setName("CONTACT DYNAMIC / MORTAR COUPLING");

  // store relevant problem types
  if (problemtype == Core::ProblemType::structure)
  {
    cparams.set<int>("PROBTYPE", Inpar::CONTACT::structure);
  }
  else if (problemtype == Core::ProblemType::tsi)
  {
    cparams.set<int>("PROBTYPE", Inpar::CONTACT::tsi);
  }
  else if (problemtype == Core::ProblemType::struct_ale)
  {
    cparams.set<int>("PROBTYPE", Inpar::CONTACT::structalewear);
  }
  else if (problemtype == Core::ProblemType::poroelast or problemtype == Core::ProblemType::fpsi or
           problemtype == Core::ProblemType::poroscatra)
  {
    const Teuchos::ParameterList& porodyn = Global::Problem::Instance()->poroelast_dynamic_params();
    if (problemtype == Core::ProblemType::poroelast or problemtype == Core::ProblemType::fpsi)
      cparams.set<int>("PROBTYPE", Inpar::CONTACT::poroelast);
    else if (problemtype == Core::ProblemType::poroscatra)
      cparams.set<int>("PROBTYPE", Inpar::CONTACT::poroscatra);
    // porotimefac = 1/(theta*dt) --- required for derivation of structural displacements!
    double porotimefac =
        1 / (stru.sublist("ONESTEPTHETA").get<double>("THETA") * stru.get<double>("TIMESTEP"));
    cparams.set<double>("porotimefac", porotimefac);
    cparams.set<bool>("CONTACTNOPEN",
        Core::UTILS::IntegralValue<int>(porodyn, "CONTACTNOPEN"));  // used in the integrator
  }
  else if (problemtype == Core::ProblemType::fsi_xfem)
  {
    cparams.set<int>("PROBTYPE", Inpar::CONTACT::fsi);
  }
  else if (problemtype == Core::ProblemType::fpsi_xfem)
  {
    const Teuchos::ParameterList& porodyn = Global::Problem::Instance()->poroelast_dynamic_params();
    cparams.set<int>("PROBTYPE", Inpar::CONTACT::fpi);
    // porotimefac = 1/(theta*dt) --- required for derivation of structural displacements!
    double porotimefac =
        1 / (stru.sublist("ONESTEPTHETA").get<double>("THETA") * stru.get<double>("TIMESTEP"));
    cparams.set<double>("porotimefac", porotimefac);
    cparams.set<bool>("CONTACTNOPEN",
        Core::UTILS::IntegralValue<int>(porodyn, "CONTACTNOPEN"));  // used in the integrator
  }
  else
  {
    cparams.set<int>("PROBTYPE", Inpar::CONTACT::other);
  }

  // no parallel redistribution in the serial case
  if (Comm().NumProc() == 1)
    cparams.sublist("PARALLEL REDISTRIBUTION").set<std::string>("PARALLEL_REDIST", "None");

  // set dimension
  cparams.set<int>("DIMENSION", dim);
  return true;
}

/*----------------------------------------------------------------------*
 |  write restart information for contact (public)            popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::write_restart(Core::IO::DiscretizationWriter& output, bool forcedrestart)
{
  // clear cache of maps due to varying vector size
  output.clear_map_cache();

  // quantities to be written for restart
  std::map<std::string, Teuchos::RCP<Epetra_Vector>> restart_vectors;

  // quantities to be written for restart
  GetStrategy().DoWriteRestart(restart_vectors, forcedrestart);

  if (GetStrategy().lagrange_multiplier_old() != Teuchos::null)
    output.write_vector("lagrmultold", GetStrategy().lagrange_multiplier_old());

  // write all vectors specified by used strategy
  for (std::map<std::string, Teuchos::RCP<Epetra_Vector>>::const_iterator p =
           restart_vectors.begin();
       p != restart_vectors.end(); ++p)
    output.write_vector(p->first, p->second);

  return;
}

/*----------------------------------------------------------------------*
 |  read restart information for contact (public)             popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::read_restart(Core::IO::DiscretizationReader& reader,
    Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<Epetra_Vector> zero)
{
  // If Parent Elements are required, we need to reconnect them before contact restart!
  Inpar::Mortar::AlgorithmType atype =
      Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(GetStrategy().Params(), "ALGORITHM");
  if (atype == Inpar::Mortar::algorithm_gpts)
  {
    for (unsigned i = 0;
         i < dynamic_cast<CONTACT::AbstractStrategy&>(GetStrategy()).contact_interfaces().size();
         ++i)
      dynamic_cast<CONTACT::AbstractStrategy&>(GetStrategy())
          .contact_interfaces()[i]
          ->create_volume_ghosting();
  }

  // If Parent Elements are required, we need to reconnect them before contact restart!
  if ((GetStrategy().Params().get<int>("PROBTYPE") == Inpar::CONTACT::poroelast ||
          GetStrategy().Params().get<int>("PROBTYPE") == Inpar::CONTACT::poroscatra) ||
      GetStrategy().Params().get<int>("PROBTYPE") == Inpar::CONTACT::fpi)
    reconnect_parent_elements();

  // this is contact, thus we need the displacement state for restart
  // let strategy object do all the work
  GetStrategy().DoReadRestart(reader, dis);

  return;
}

/*----------------------------------------------------------------------*
 |  write interface tractions for postprocessing (public)     popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::postprocess_quantities(Core::IO::DiscretizationWriter& output)
{
  if (GetStrategy().IsNitsche()) return;

  // *********************************************************************
  // active contact set and slip set
  // *********************************************************************

  // evaluate active set and slip set
  Teuchos::RCP<Epetra_Vector> activeset =
      Teuchos::rcp(new Epetra_Vector(*GetStrategy().active_row_nodes()));
  activeset->PutScalar(1.0);
  if (GetStrategy().is_friction())
  {
    Teuchos::RCP<Epetra_Vector> slipset =
        Teuchos::rcp(new Epetra_Vector(*GetStrategy().slip_row_nodes()));
    slipset->PutScalar(1.0);
    Teuchos::RCP<Epetra_Vector> slipsetexp =
        Teuchos::rcp(new Epetra_Vector(*GetStrategy().active_row_nodes()));
    Core::LinAlg::Export(*slipset, *slipsetexp);
    activeset->Update(1.0, *slipsetexp, 1.0);
  }

  // export to problem node row map
  Teuchos::RCP<Epetra_Map> problemnodes = GetStrategy().ProblemNodes();
  Teuchos::RCP<Epetra_Vector> activesetexp = Teuchos::rcp(new Epetra_Vector(*problemnodes));
  Core::LinAlg::Export(*activeset, *activesetexp);

  if (GetStrategy().WearBothDiscrete())
  {
    Teuchos::RCP<Epetra_Vector> mactiveset =
        Teuchos::rcp(new Epetra_Vector(*GetStrategy().MasterActiveNodes()));
    mactiveset->PutScalar(1.0);
    Teuchos::RCP<Epetra_Vector> slipset =
        Teuchos::rcp(new Epetra_Vector(*GetStrategy().MasterSlipNodes()));
    slipset->PutScalar(1.0);
    Teuchos::RCP<Epetra_Vector> slipsetexp =
        Teuchos::rcp(new Epetra_Vector(*GetStrategy().MasterActiveNodes()));
    Core::LinAlg::Export(*slipset, *slipsetexp);
    mactiveset->Update(1.0, *slipsetexp, 1.0);

    Teuchos::RCP<Epetra_Vector> mactivesetexp = Teuchos::rcp(new Epetra_Vector(*problemnodes));
    Core::LinAlg::Export(*mactiveset, *mactivesetexp);
    activesetexp->Update(1.0, *mactivesetexp, 1.0);
  }

  output.write_vector("activeset", activesetexp);

  // *********************************************************************
  //  weighted gap
  // *********************************************************************
  // export to problem dof row map
  Teuchos::RCP<Epetra_Map> gapnodes = GetStrategy().ProblemNodes();
  Teuchos::RCP<Epetra_Vector> gaps =
      Teuchos::rcp_dynamic_cast<CONTACT::AbstractStrategy>(strategy_)->contact_wgap();
  if (gaps != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> gapsexp = Teuchos::rcp(new Epetra_Vector(*gapnodes));
    Core::LinAlg::Export(*gaps, *gapsexp);

    output.write_vector("gap", gapsexp);
  }

  // *********************************************************************
  // contact tractions
  // *********************************************************************

  // evaluate contact tractions
  GetStrategy().compute_contact_stresses();

  // export to problem dof row map
  Teuchos::RCP<Epetra_Map> problemdofs = GetStrategy().ProblemDofs();

  // normal direction
  Teuchos::RCP<Epetra_Vector> normalstresses = GetStrategy().contact_normal_stress();
  Teuchos::RCP<Epetra_Vector> normalstressesexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  Core::LinAlg::Export(*normalstresses, *normalstressesexp);

  // tangential plane
  Teuchos::RCP<Epetra_Vector> tangentialstresses = GetStrategy().contact_tangential_stress();
  Teuchos::RCP<Epetra_Vector> tangentialstressesexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  Core::LinAlg::Export(*tangentialstresses, *tangentialstressesexp);

  // write to output
  // contact tractions in normal and tangential direction
  output.write_vector("norcontactstress", normalstressesexp);
  output.write_vector("tancontactstress", tangentialstressesexp);

  if (GetStrategy().contact_normal_force() != Teuchos::null)
  {
    // normal direction
    Teuchos::RCP<Epetra_Vector> normalforce = GetStrategy().contact_normal_force();
    Teuchos::RCP<Epetra_Vector> normalforceexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
    Core::LinAlg::Export(*normalforce, *normalforceexp);

    // tangential plane
    Teuchos::RCP<Epetra_Vector> tangentialforce = GetStrategy().contact_tangential_force();
    Teuchos::RCP<Epetra_Vector> tangentialforceexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
    Core::LinAlg::Export(*tangentialforce, *tangentialforceexp);

    // write to output
    // contact tractions in normal and tangential direction
    output.write_vector("norslaveforce", normalforceexp);
    output.write_vector("tanslaveforce", tangentialforceexp);
  }


#ifdef CONTACTFORCEOUTPUT

  // *********************************************************************
  // contact forces on slave non master side,
  // in normal and tangential direction
  // *********************************************************************
  // vectors for contact forces
  Teuchos::RCP<Epetra_Vector> fcslavenor =
      Teuchos::rcp(new Epetra_Vector(GetStrategy().DMatrix()->RowMap()));
  Teuchos::RCP<Epetra_Vector> fcslavetan =
      Teuchos::rcp(new Epetra_Vector(GetStrategy().DMatrix()->RowMap()));
  Teuchos::RCP<Epetra_Vector> fcmasternor =
      Teuchos::rcp(new Epetra_Vector(GetStrategy().MMatrix()->DomainMap()));
  Teuchos::RCP<Epetra_Vector> fcmastertan =
      Teuchos::rcp(new Epetra_Vector(GetStrategy().MMatrix()->DomainMap()));

  // vectors with problem dof row map
  Teuchos::RCP<Epetra_Vector> fcslavenorexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  Teuchos::RCP<Epetra_Vector> fcslavetanexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  Teuchos::RCP<Epetra_Vector> fcmasternorexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  Teuchos::RCP<Epetra_Vector> fcmastertanexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));

  // multiplication
  GetStrategy().DMatrix()->Multiply(true, *normalstresses, *fcslavenor);
  GetStrategy().DMatrix()->Multiply(true, *tangentialstresses, *fcslavetan);
  GetStrategy().MMatrix()->Multiply(true, *normalstresses, *fcmasternor);
  GetStrategy().MMatrix()->Multiply(true, *tangentialstresses, *fcmastertan);

#ifdef MASTERNODESINCONTACT
  // BEGIN: to output the global ID's of the master nodes in contact - devaal 02.2011

  int dim = Global::Problem::Instance()->NDim();

  if (dim == 2) FOUR_C_THROW("Only working for 3D");

  std::vector<int> lnid, gnid;

  // std::cout << "MasterNor" << fcmasternor->MyLength() << std::endl;

  for (int i = 0; i < fcmasternor->MyLength(); i = i + 3)
  {
    // check if master node in contact
    if (sqrt(((*fcmasternor)[i]) * ((*fcmasternor)[i]) +
             ((*fcmasternor)[i + 1]) * ((*fcmasternor)[i + 1]) +
             ((*fcmasternor)[i + 2]) * ((*fcmasternor)[i] + 2)) > 0.00001)
    {
      lnid.push_back((fcmasternor->Map()).GID(i) / 3);
    }
  }

  // we want to gather data from on all procs
  std::vector<int> allproc(Comm().NumProc());
  for (int i = 0; i < Comm().NumProc(); ++i) allproc[i] = i;

  // communicate all data to proc 0
  Core::LinAlg::Gather<int>(lnid, gnid, static_cast<int>(llproc.size()), allproc.data(), Comm());

  // std::cout << " size of gnid:" << gnid.size() << std::endl;

  ////////////////
  ///// attempt at obtaining the nid and relative displacement u of master nodes in contact - devaal
  // define my own interface
  Mortar::StrategyBase& myStrategy = GetStrategy();
  AbstractStrategy& myContactStrategy = dynamic_cast<AbstractStrategy&>(myStrategy);

  std::vector<Teuchos::RCP<CONTACT::Interface>> myInterface = myContactStrategy.ContactInterfaces();

  // check interface size - just doing this now for a single interface

  if (myInterface.size() != 1) FOUR_C_THROW("Interface size should be 1");

  std::cout << "OUTPUT OF MASTER NODE IN CONTACT" << std::endl;
  for (const auto& globalNodeId : gnid) std::cout << globalNodeId << std::endl;

#endif

  //  // when we do a boundary modification we shift slave entries to the M matrix with
  //  // negative sign. Therefore, we have to extract the right force entries from the
  //  // master force which correcpond to the slave force!
  //  Teuchos::RCP<Epetra_Vector> slavedummy =
  //      Teuchos::rcp(new Epetra_Vector(GetStrategy().d_matrix()->RowMap(),true));
  //  Core::LinAlg::Export(*fcmasternor,*slavedummy);
  //  int err = fcslavenor->Update(-1.0,*slavedummy,1.0);
  //  if(err!=0)
  //    FOUR_C_THROW("ERROR");
  //
  //  Teuchos::RCP<Epetra_Vector> masterdummy =
  //      Teuchos::rcp(new Epetra_Vector(GetStrategy().m_matrix()->DomainMap(),true));
  //  Core::LinAlg::Export(*slavedummy,*masterdummy);
  //  err = fcmasternor->Update(-1.0,*masterdummy,1.0);
  //  if(err!=0)
  //    FOUR_C_THROW("ERROR");

  // export
  Core::LinAlg::Export(*fcslavenor, *fcslavenorexp);
  Core::LinAlg::Export(*fcslavetan, *fcslavetanexp);
  Core::LinAlg::Export(*fcmasternor, *fcmasternorexp);
  Core::LinAlg::Export(*fcmastertan, *fcmastertanexp);

  // contact forces on slave and master side
  output.write_vector("norslaveforce", fcslavenorexp);
  output.write_vector("tanslaveforce", fcslavetanexp);
  output.write_vector("normasterforce", fcmasternorexp);
  output.write_vector("tanmasterforce", fcmastertanexp);

#ifdef CONTACTEXPORT
  // export averaged node forces to xxx.force
  double resultnor[fcslavenor->NumVectors()];
  double resulttan[fcslavetan->NumVectors()];
  fcslavenor->Norm2(resultnor);
  fcslavetan->Norm2(resulttan);

  if (Comm().MyPID() == 0)
  {
    std::cout << "resultnor= " << resultnor[0] << std::endl;
    std::cout << "resulttan= " << resulttan[0] << std::endl;

    FILE* MyFile = nullptr;
    std::ostringstream filename;
    const std::string filebase =
        Global::Problem::Instance()->OutputControlFile()->file_name_only_prefix();
    filename << filebase << ".force";
    MyFile = fopen(filename.str().c_str(), "at+");
    if (MyFile)
    {
      // fprintf(MyFile,valuename.c_str());
      fprintf(MyFile, "%g\t", resultnor[0]);
      fprintf(MyFile, "%g\n", resulttan[0]);
      fclose(MyFile);
    }
    else
      FOUR_C_THROW("File for Output could not be opened.");
  }
#endif  // CONTACTEXPORT
#endif  // CONTACTFORCEOUTPUT

  // *********************************************************************
  // wear with internal state variable approach
  // *********************************************************************
  bool wwear = GetStrategy().WeightedWear();
  if (wwear)
  {
    // ***************************************************************************
    // we do not compute the non-weighted wear here. we just write    farah 06/13
    // the output. the non-weighted wear will be used as dirichlet-b.
    // for the ale problem. n.w.wear will be called in stru_ale_algorithm.cpp
    // and computed in GetStrategy().OutputWear();
    // ***************************************************************************

    // evaluate wear (not weighted)
    GetStrategy().OutputWear();

    // write output
    Teuchos::RCP<Epetra_Vector> wearoutput = GetStrategy().ContactWear();
    Teuchos::RCP<Epetra_Vector> wearoutputexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
    Core::LinAlg::Export(*wearoutput, *wearoutputexp);
    output.write_vector("wear", wearoutputexp);
    GetStrategy().ContactWear()->PutScalar(0.0);
  }

  // *********************************************************************
  // poro contact
  // *********************************************************************
  bool poro = GetStrategy().has_poro_no_penetration();
  if (poro)
  {
    // output of poro no penetration lagrange multiplier!
    CONTACT::LagrangeStrategyPoro& costrategy =
        dynamic_cast<CONTACT::LagrangeStrategyPoro&>(GetStrategy());
    Teuchos::RCP<Epetra_Vector> lambdaout = costrategy.LambdaNoPen();
    Teuchos::RCP<Epetra_Vector> lambdaoutexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
    Core::LinAlg::Export(*lambdaout, *lambdaoutexp);
    output.write_vector("poronopen_lambda", lambdaoutexp);
  }
  return;
}

/*-----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Manager::postprocess_quantities_per_interface(
    Teuchos::RCP<Teuchos::ParameterList> outputParams)
{
  GetStrategy().postprocess_quantities_per_interface(outputParams);
}

/*----------------------------------------------------------------------------------------------*
 |  Reconnect Contact Element -- Parent Element Pointers (required for restart)       ager 04/16|
 *---------------------------------------------------------------------------------------------*/
void CONTACT::Manager::reconnect_parent_elements()
{
  {
    const Epetra_Map* elecolmap = discret_.ElementColMap();

    CONTACT::AbstractStrategy& strategy = dynamic_cast<CONTACT::AbstractStrategy&>(GetStrategy());

    for (auto& interface : strategy.contact_interfaces())
    {
      const Epetra_Map* ielecolmap = interface->Discret().ElementColMap();

      for (int i = 0; i < ielecolmap->NumMyElements(); ++i)
      {
        int gid = ielecolmap->GID(i);

        Core::Elements::Element* ele = interface->Discret().gElement(gid);
        if (!ele) FOUR_C_THROW("Cannot find element with gid %", gid);
        Core::Elements::FaceElement* faceele = dynamic_cast<Core::Elements::FaceElement*>(ele);

        int volgid = faceele->ParentElementId();
        if (elecolmap->LID(volgid) == -1)  // Volume discretization has not Element
          FOUR_C_THROW(
              "Manager::reconnect_parent_elements: Element %d does not exist on this Proc!",
              volgid);

        Core::Elements::Element* vele = discret_.gElement(volgid);
        if (!vele) FOUR_C_THROW("Cannot find element with gid %", volgid);

        faceele->set_parent_master_element(vele, faceele->FaceParentNumber());
      }
    }
  }
}

/*----------------------------------------------------------------------*
 |  Set Parent Elements for Poro Face Elements                ager 11/15|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::set_poro_parent_element(int& slavetype, int& mastertype,
    Teuchos::RCP<CONTACT::Element>& cele, Teuchos::RCP<Core::Elements::Element>& ele)
{
  // ints to communicate decision over poro bools between processors on every interface
  // safety check - because there may not be mixed interfaces and structural slave elements
  // slavetype ... 1 poro, 0 struct, -1 default
  // mastertype ... 1 poro, 0 struct, -1 default
  Teuchos::RCP<Core::Elements::FaceElement> faceele =
      Teuchos::rcp_dynamic_cast<Core::Elements::FaceElement>(ele, true);
  if (faceele == Teuchos::null) FOUR_C_THROW("Cast to FaceElement failed!");
  cele->PhysType() = Mortar::Element::other;
  std::vector<Teuchos::RCP<Core::Conditions::Condition>> porocondvec;
  discret_.GetCondition("PoroCoupling", porocondvec);
  if (!cele->IsSlave())  // treat an element as a master element if it is no slave element
  {
    for (unsigned int i = 0; i < porocondvec.size(); ++i)
    {
      std::map<int, Teuchos::RCP<Core::Elements::Element>>::const_iterator eleitergeometry;
      for (eleitergeometry = porocondvec[i]->Geometry().begin();
           eleitergeometry != porocondvec[i]->Geometry().end(); ++eleitergeometry)
      {
        if (faceele->parent_element()->Id() == eleitergeometry->second->Id())
        {
          if (mastertype == 0)
            FOUR_C_THROW(
                "struct and poro master elements on the same processor - no mixed interface "
                "supported");
          cele->PhysType() = Mortar::Element::poro;
          mastertype = 1;
          break;
        }
      }
    }
    if (cele->PhysType() == Mortar::Element::other)
    {
      if (mastertype == 1)
        FOUR_C_THROW(
            "struct and poro master elements on the same processor - no mixed interface supported");
      cele->PhysType() = Mortar::Element::structure;
      mastertype = 0;
    }
  }
  else if (cele->IsSlave())  // treat an element as slave element if it is one
  {
    for (unsigned int i = 0; i < porocondvec.size(); ++i)
    {
      std::map<int, Teuchos::RCP<Core::Elements::Element>>::const_iterator eleitergeometry;
      for (eleitergeometry = porocondvec[i]->Geometry().begin();
           eleitergeometry != porocondvec[i]->Geometry().end(); ++eleitergeometry)
      {
        if (faceele->parent_element()->Id() == eleitergeometry->second->Id())
        {
          if (slavetype == 0)
            FOUR_C_THROW(
                "struct and poro master elements on the same processor - no mixed interface "
                "supported");
          cele->PhysType() = Mortar::Element::poro;
          slavetype = 1;
          break;
        }
      }
    }
    if (cele->PhysType() == Mortar::Element::other)
    {
      if (slavetype == 1)
        FOUR_C_THROW(
            "struct and poro master elements on the same processor - no mixed interface supported");
      cele->PhysType() = Mortar::Element::structure;
      slavetype = 0;
    }
  }
  // store information about parent for porous contact (required for calculation of deformation
  // gradient!) in every contact element although only really needed for phystype poro
  cele->set_parent_master_element(faceele->parent_element(), faceele->FaceParentNumber());
  return;
}

/*----------------------------------------------------------------------*
 |  Find Physical Type (Poro or Structure) of Poro Interface  ager 11/15|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::find_poro_interface_types(bool& poromaster, bool& poroslave,
    bool& structmaster, bool& structslave, int& slavetype, int& mastertype)
{
  // find poro and structure elements when a poro coupling condition is applied on an element
  // and restrict to pure poroelastic or pure structural interfaces' sides.
  //(only poro slave elements AND (only poro master elements or only structure master elements)
  // Tell the contact element which physical type it is to extract PhysType in contact integrator
  // bools to decide which side is structural and which side is poroelastic to manage all 4
  // constellations
  // s-s, p-s, s-p, p-p
  // wait for all processors to determine if they have poro or structural master or slave elements
  comm_->Barrier();
  std::vector<int> slaveTypeList(comm_->NumProc());
  std::vector<int> masterTypeList(comm_->NumProc());
  comm_->GatherAll(&slavetype, slaveTypeList.data(), 1);
  comm_->GatherAll(&mastertype, masterTypeList.data(), 1);
  comm_->Barrier();

  for (int i = 0; i < comm_->NumProc(); ++i)
  {
    switch (slaveTypeList[i])
    {
      case -1:
        break;
      case 1:
        if (structslave)
          FOUR_C_THROW(
              "struct and poro slave elements in the same problem - no mixed interface "
              "constellations supported");
        // adjust FOUR_C_THROW text, when more than one interface is supported
        poroslave = true;
        break;
      case 0:
        if (poroslave)
          FOUR_C_THROW(
              "struct and poro slave elements in the same problem - no mixed interface "
              "constellations supported");
        structslave = true;
        break;
      default:
        FOUR_C_THROW("this cannot happen");
        break;
    }
  }

  for (int i = 0; i < comm_->NumProc(); ++i)
  {
    switch (masterTypeList[i])
    {
      case -1:
        break;
      case 1:
        if (structmaster)
          FOUR_C_THROW(
              "struct and poro master elements in the same problem - no mixed interface "
              "constellations supported");
        // adjust FOUR_C_THROW text, when more than one interface is supported
        poromaster = true;
        break;
      case 0:
        if (poromaster)
          FOUR_C_THROW(
              "struct and poro master elements in the same problem - no mixed interface "
              "constellations supported");
        structmaster = true;
        break;
      default:
        FOUR_C_THROW("this cannot happen");
        break;
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
