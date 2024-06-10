/*----------------------------------------------------------------------*/
/*! \file

\brief Basis of all WEAR algorithms that perform a coupling between the
       structural field equation and ALE field equations

\level 2


*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  farah 11/13 |
 *----------------------------------------------------------------------*/
#include "4C_wear_algorithm.hpp"

#include "4C_adapter_ale.hpp"
#include "4C_adapter_ale_wear.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_adapter_str_structure.hpp"
#include "4C_contact_aug_interface.hpp"
#include "4C_contact_element.hpp"
#include "4C_contact_friction_node.hpp"
#include "4C_contact_lagrange_strategy_wear.hpp"
#include "4C_contact_meshtying_contact_bridge.hpp"
#include "4C_contact_strategy_factory.hpp"
#include "4C_contact_utils.hpp"
#include "4C_contact_wear_interface.hpp"
#include "4C_fem_nurbs_discretization.hpp"
#include "4C_fem_nurbs_discretization_control_point.hpp"
#include "4C_fem_nurbs_discretization_knotvector.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_ale.hpp"
#include "4C_inpar_wear.hpp"
#include "4C_mortar_manager_base.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | Constructor                                              farah 11/13 |
 *----------------------------------------------------------------------*/
Wear::Algorithm::Algorithm(const Epetra_Comm& comm)
    : AlgorithmBase(comm, Global::Problem::Instance()->structural_dynamic_params())

{
  /*--------------------------------------------------------------------*
   | first create structure then ale --> important for discretization   |
   | numbering and therefore for the post_ensight.cpp                   |
   *--------------------------------------------------------------------*/

  // create structure
  Teuchos::RCP<Adapter::StructureBaseAlgorithm> structure = Teuchos::rcp(
      new Adapter::StructureBaseAlgorithm(Global::Problem::Instance()->structural_dynamic_params(),
          const_cast<Teuchos::ParameterList&>(
              Global::Problem::Instance()->structural_dynamic_params()),
          Global::Problem::Instance()->GetDis("structure")));
  structure_ =
      Teuchos::rcp_dynamic_cast<Adapter::FSIStructureWrapper>(structure->structure_field());
  structure_->Setup();

  if (structure_ == Teuchos::null)
    FOUR_C_THROW("cast from Adapter::Structure to Adapter::FSIStructureWrapper failed");

  // ask base algorithm for the ale time integrator
  Teuchos::RCP<Adapter::AleBaseAlgorithm> ale = Teuchos::rcp(
      new Adapter::AleBaseAlgorithm(Global::Problem::Instance()->structural_dynamic_params(),
          Global::Problem::Instance()->GetDis("ale")));
  ale_ = Teuchos::rcp_dynamic_cast<Adapter::AleWearWrapper>(ale->ale_field());
  if (ale_ == Teuchos::null)
    FOUR_C_THROW("cast from Adapter::Ale to Adapter::AleFsiWrapper failed");

  // create empty operator
  ale_->create_system_matrix();

  // contact/meshtying manager
  cmtman_ = structure_field()->meshtying_contact_bridge()->ContactManager();

  // copy interfaces for material configuration
  // stactic cast of mortar strategy to contact strategy
  Mortar::StrategyBase& strategy = cmtman_->GetStrategy();
  Wear::LagrangeStrategyWear& cstrategy = static_cast<Wear::LagrangeStrategyWear&>(strategy);

  // get dimension
  dim_ = strategy.Dim();

  // get vector of contact interfaces
  interfaces_ = cstrategy.ContactInterfaces();

  // create contact interfaces for material conf.
  create_material_interface();

  // input
  check_input();
}



/*----------------------------------------------------------------------*
 | Check compatibility of input parameters                  farah 09/14 |
 *----------------------------------------------------------------------*/
void Wear::Algorithm::check_input()
{
  //  Teuchos::ParameterList apara = Global::Problem::Instance()->AleDynamicParams();
  //
  //  Inpar::ALE::AleDynamic aletype =
  //      Core::UTILS::IntegralValue<Inpar::ALE::AleDynamic>(apara, "ALE_TYPE");

  return;
}


/*----------------------------------------------------------------------*
 | Create interfaces for material conf.                     farah 09/14 |
 *----------------------------------------------------------------------*/
void Wear::Algorithm::create_material_interface()
{
  Mortar::StrategyBase& strategy = cmtman_->GetStrategy();
  Wear::LagrangeStrategyWear& cstrategy = static_cast<Wear::LagrangeStrategyWear&>(strategy);

  // create some local variables (later to be stored in strategy)
  int dim = Global::Problem::Instance()->NDim();
  if (dim != 2 && dim != 3) FOUR_C_THROW("Contact problem must be 2D or 3D");
  Teuchos::ParameterList cparams = cstrategy.Params();

  // check for fill_complete of discretization
  if (!structure_->discretization()->Filled()) FOUR_C_THROW("discretization is not fillcomplete");

  // let's check for contact boundary conditions in discret
  // and detect groups of matching conditions
  // for each group, create a contact interface and store it
  if (Comm().MyPID() == 0)
  {
    std::cout << "Building contact interface(s) for Mat. conf. ...............";
    fflush(stdout);
  }

  std::vector<Core::Conditions::Condition*> contactconditions(0);
  structure_->discretization()->GetCondition("Contact", contactconditions);

  // there must be more than one contact condition
  // unless we have a self contact problem!
  if ((int)contactconditions.size() < 1)
    FOUR_C_THROW("Not enough contact conditions in discretization");
  if ((int)contactconditions.size() == 1)
  {
    const std::string& side = contactconditions[0]->parameters().Get<std::string>("Side");
    if (side != "Selfcontact") FOUR_C_THROW("Not enough contact conditions in discretization");
  }

  // find all pairs of matching contact conditions
  // there is a maximum of (conditions / 2) groups
  std::vector<int> foundgroups(0);
  int numgroupsfound = 0;

  // maximum dof number in discretization
  // later we want to create NEW Lagrange multiplier degrees of
  // freedom, which of course must not overlap with displacement dofs
  int maxdof = structure_->discretization()->dof_row_map()->MaxAllGID();

  // get input par.
  Inpar::CONTACT::SolvingStrategy stype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(cparams, "STRATEGY");
  Inpar::Wear::WearLaw wlaw = Core::UTILS::IntegralValue<Inpar::Wear::WearLaw>(cparams, "WEARLAW");
  Inpar::CONTACT::ConstraintDirection constr_direction =
      Core::UTILS::IntegralValue<Inpar::CONTACT::ConstraintDirection>(
          cparams, "CONSTRAINT_DIRECTIONS");

  bool friplus = false;
  if ((wlaw != Inpar::Wear::wear_none) || (cparams.get<int>("PROBTYPE") == Inpar::CONTACT::tsi))
    friplus = true;

  bool isanyselfcontact = false;

  for (int i = 0; i < (int)contactconditions.size(); ++i)
  {
    // initialize vector for current group of conditions and temp condition
    std::vector<Core::Conditions::Condition*> currentgroup(0);
    Core::Conditions::Condition* tempcond = nullptr;

    // try to build contact group around this condition
    currentgroup.push_back(contactconditions[i]);
    int groupid1 = currentgroup[0]->parameters().Get<int>("Interface ID");
    bool foundit = false;

    // only one surface per group is ok for self contact
    const std::string& side = contactconditions[i]->parameters().Get<std::string>("Side");
    if (side == "Selfcontact") foundit = true;

    for (int j = 0; j < (int)contactconditions.size(); ++j)
    {
      if (j == i) continue;  // do not detect contactconditions[i] again
      tempcond = contactconditions[j];
      int groupid2 = currentgroup[0]->parameters().Get<int>("Interface ID");
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
    std::vector<bool> isslave(0);
    std::vector<bool> isself(0);
    CONTACT::UTILS::GetMasterSlaveSideInfo(isslave, isself, currentgroup);
    for (const bool is : isself)
      if (is)
      {
        isanyselfcontact = true;
        break;
      }

    // find out which sides are initialized as Active
    std::vector<bool> isactive(currentgroup.size());
    bool Two_half_pass(false);
    bool Check_nonsmooth_selfcontactsurface(false);
    bool Searchele_AllProc(false);

    CONTACT::UTILS::GetInitializationInfo(Two_half_pass, Check_nonsmooth_selfcontactsurface,
        Searchele_AllProc, isactive, isslave, isself, currentgroup);

    // create interface local parameter list (copy)
    Teuchos::ParameterList icparams = cparams;

    // find out if interface-specific coefficients of friction are given
    Inpar::CONTACT::FrictionType fric =
        Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(cparams, "FRICTION");
    if (fric == Inpar::CONTACT::friction_tresca || fric == Inpar::CONTACT::friction_coulomb)
    {
      // read interface COFs
      std::vector<double> frcoeff((int)currentgroup.size());
      for (int j = 0; j < (int)currentgroup.size(); ++j)
        frcoeff[j] = currentgroup[j]->parameters().Get<double>("FrCoeffOrBound");

      // check consistency of interface COFs
      for (int j = 1; j < (int)currentgroup.size(); ++j)
        if (frcoeff[j] != frcoeff[0])
          FOUR_C_THROW("Inconsistency in friction coefficients of interface %i", groupid1);

      // check for infeasible value of COF
      if (frcoeff[0] < 0.0) FOUR_C_THROW("Negative FrCoeff / FrBound on interface %i", groupid1);

      // add COF locally to contact parameter list of this interface
      if (fric == Inpar::CONTACT::friction_tresca)
      {
        icparams.setEntry("FRBOUND", static_cast<Teuchos::ParameterEntry>(frcoeff[0]));
        icparams.setEntry("FRCOEFF", static_cast<Teuchos::ParameterEntry>(-1.0));
      }
      else if (fric == Inpar::CONTACT::friction_coulomb)
      {
        icparams.setEntry("FRCOEFF", static_cast<Teuchos::ParameterEntry>(frcoeff[0]));
        icparams.setEntry("FRBOUND", static_cast<Teuchos::ParameterEntry>(-1.0));
      }
    }

    // find out if interface-specific coefficients of friction are given
    Inpar::CONTACT::AdhesionType ad =
        Core::UTILS::IntegralValue<Inpar::CONTACT::AdhesionType>(cparams, "ADHESION");
    if (ad == Inpar::CONTACT::adhesion_bound)
    {
      // read interface COFs
      std::vector<double> ad_bound((int)currentgroup.size());
      for (int j = 0; j < (int)currentgroup.size(); ++j)
        ad_bound[j] = currentgroup[j]->parameters().Get<double>("AdhesionBound");

      // check consistency of interface COFs
      for (int j = 1; j < (int)currentgroup.size(); ++j)
        if (ad_bound[j] != ad_bound[0])
          FOUR_C_THROW("Inconsistency in adhesion bounds of interface %i", groupid1);

      // check for infeasible value of COF
      if (ad_bound[0] < 0.0) FOUR_C_THROW("Negative adhesion bound on interface %i", groupid1);

      // add COF locally to contact parameter list of this interface
      icparams.setEntry("ADHESION_BOUND", static_cast<Teuchos::ParameterEntry>(ad_bound[0]));
    }

    // add information to parameter list of this interface
    icparams.set<bool>("Two_half_pass", Two_half_pass);
    icparams.set<bool>("Check_nonsmooth_selfcontactsurface", Check_nonsmooth_selfcontactsurface);
    icparams.set<bool>("Searchele_AllProc", Searchele_AllProc);

    // for structural contact we currently choose redundant master storage
    // the only exception is self contact where a redundant slave is needed, too
    Inpar::Mortar::ExtendGhosting redundant =
        Teuchos::getIntegralValue<Inpar::Mortar::ExtendGhosting>(
            icparams.sublist("PARALLEL REDISTRIBUTION"), "GHOSTING_STRATEGY");
    if (isanyselfcontact == true && redundant != Inpar::Mortar::ExtendGhosting::redundant_all)
      FOUR_C_THROW("Self contact requires fully redundant slave and master storage");

    // decide between contactinterface, augmented interface and wearinterface
    Teuchos::RCP<CONTACT::Interface> newinterface = CONTACT::STRATEGY::Factory::CreateInterface(
        groupid1, Comm(), dim, icparams, isself[0], Teuchos::null);
    interfacesMat_.push_back(newinterface);

    // get it again
    Teuchos::RCP<CONTACT::Interface> interface = interfacesMat_[(int)interfacesMat_.size() - 1];

    // note that the nodal ids are unique because they come from
    // one global problem discretization containing all nodes of the
    // contact interface.
    // We rely on this fact, therefore it is not possible to
    // do contact between two distinct discretizations here.

    // collect all intial active nodes
    std::vector<int> initialactive;

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
        if (!structure_->discretization()->NodeColMap()->MyGID(gid)) continue;
        Core::Nodes::Node* node = structure_->discretization()->gNode(gid);
        if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);

        // store initial active node gids
        if (isactive[j]) initialactive.push_back(gid);

        // find out if this node is initial active on another Condition
        // and do NOT overwrite this status then!
        bool foundinitialactive = false;
        if (!isactive[j])
        {
          for (int k = 0; k < (int)initialactive.size(); ++k)
            if (gid == initialactive[k])
            {
              foundinitialactive = true;
              break;
            }
        }

        // create Node object or FriNode object in the frictional case
        Inpar::CONTACT::FrictionType ftype =
            Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(cparams, "FRICTION");

        // for the boolean variable initactive we use isactive[j]+foundinitialactive,
        // as this is true for BOTH initial active nodes found for the first time
        // and found for the second, third, ... time!
        if (ftype != Inpar::CONTACT::friction_none)
        {
          Teuchos::RCP<CONTACT::FriNode> cnode = Teuchos::rcp(new CONTACT::FriNode(node->Id(),
              node->X(), node->Owner(), structure_->discretization()->Dof(0, node), isslave[j],
              isactive[j] + foundinitialactive, friplus));
          //-------------------
          // get nurbs weight!
          if (cparams.get<bool>("NURBS") == true)
          {
            Discret::Nurbs::ControlPoint* cp = dynamic_cast<Discret::Nurbs::ControlPoint*>(node);

            cnode->NurbsW() = cp->W();
          }

          // Check, if this node (and, in case, which dofs) are in the contact symmetry condition
          std::vector<Core::Conditions::Condition*> contactSymconditions(0);
          structure_->discretization()->GetCondition("mrtrsym", contactSymconditions);

          for (unsigned j = 0; j < contactSymconditions.size(); j++)
            if (contactSymconditions.at(j)->ContainsNode(node->Id()))
            {
              const std::vector<int>& onoff =
                  contactSymconditions.at(j)->parameters().Get<std::vector<int>>("onoff");
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
        else
        {
          Teuchos::RCP<CONTACT::Node> cnode = Teuchos::rcp(new CONTACT::Node(node->Id(), node->X(),
              node->Owner(), structure_->discretization()->Dof(0, node), isslave[j],
              isactive[j] + foundinitialactive));
          //-------------------
          // get nurbs weight!
          if (cparams.get<bool>("NURBS") == true)
          {
            Discret::Nurbs::ControlPoint* cp = dynamic_cast<Discret::Nurbs::ControlPoint*>(node);

            cnode->NurbsW() = cp->W();
          }

          // Check, if this node (and, in case, which dofs) are in the contact symmetry condition
          std::vector<Core::Conditions::Condition*> contactSymconditions(0);
          structure_->discretization()->GetCondition("mrtrsym", contactSymconditions);

          for (unsigned j = 0; j < contactSymconditions.size(); j++)
            if (contactSymconditions.at(j)->ContainsNode(node->Id()))
            {
              const std::vector<int>& onoff =
                  contactSymconditions.at(j)->parameters().Get<std::vector<int>>("onoff");
              for (unsigned k = 0; k < onoff.size(); k++)
                if (onoff.at(k) == 1) cnode->DbcDofs()[k] = true;
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
    for (int j = 0; j < (int)currentgroup.size(); ++j)
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
        Teuchos::RCP<CONTACT::Element> cele =
            Teuchos::rcp(new CONTACT::Element(ele->Id() + ggsize, ele->Owner(), ele->Shape(),
                ele->num_node(), ele->NodeIds(), isslave[j], cparams.get<bool>("NURBS")));

        //------------------------------------------------------------------
        // get knotvector, normal factor and zero-size information for nurbs
        if (cparams.get<bool>("NURBS") == true)
        {
          Discret::Nurbs::NurbsDiscretization* nurbsdis =
              dynamic_cast<Discret::Nurbs::NurbsDiscretization*>(
                  &(*(structure_->discretization())));

          Teuchos::RCP<Discret::Nurbs::Knotvector> knots = (*nurbsdis).GetKnotVector();
          std::vector<Core::LinAlg::SerialDenseVector> parentknots(dim);
          std::vector<Core::LinAlg::SerialDenseVector> mortarknots(dim - 1);

          Teuchos::RCP<Core::Elements::FaceElement> faceele =
              Teuchos::rcp_dynamic_cast<Core::Elements::FaceElement>(ele, true);
          double normalfac = 0.0;
          bool zero_size = knots->get_boundary_ele_and_parent_knots(parentknots, mortarknots,
              normalfac, faceele->ParentMasterElement()->Id(), faceele->FaceMasterNumber());

          // store nurbs specific data to node
          cele->ZeroSized() = zero_size;
          cele->Knots() = mortarknots;
          cele->NormalFac() = normalfac;
        }

        interface->add_element(cele);
      }  // for (fool=ele1.start(); fool != ele1.end(); ++fool)

      ggsize += gsize;  // update global element counter
    }

    //-------------------- finalize the contact interface construction
    interface->fill_complete(maxdof);

  }  // for (int i=0; i<(int)contactconditions.size(); ++i)
  if (Comm().MyPID() == 0) std::cout << "done!" << std::endl;

  return;
}

FOUR_C_NAMESPACE_CLOSE
