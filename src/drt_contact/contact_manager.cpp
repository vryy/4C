/*!----------------------------------------------------------------------
\file contact_manager.cpp
\brief BACI implementation of main class to control all contact

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "contact_manager.H"
#include "contact_interface.H"
#include "contact_node.H"
#include "contact_element.H"
#include "contact_lagrange_strategy.H"
#include "contact_penalty_strategy.H"
#include "contact_defines.H"
#include "friction_node.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_mortar.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 03/08|
 *----------------------------------------------------------------------*/
CONTACT::CoManager::CoManager(DRT::Discretization& discret, double alphaf) :
MORTAR::ManagerBase(),
discret_(discret)
{
  // overwrite base class communicator
  comm_ = rcp(Discret().Comm().Clone());

  // welcome message

  // create some local variables (later to be stored in strategy)
  RCP<Epetra_Map> problemrowmap = rcp(new Epetra_Map(*(Discret().DofRowMap())));
  const Teuchos::ParameterList& psize = DRT::Problem::Instance()->ProblemSizeParams();
  int dim = psize.get<int>("DIM");
  if (dim!= 2 && dim!=3) dserror("ERROR: Contact problem must be 2D or 3D");
  vector<RCP<CONTACT::CoInterface> > interfaces;
  Teuchos::ParameterList cparams;

  // read and check contact input parameters
  if(Comm().MyPID()==0)
  {
    cout << "Checking contact input parameters...........";
    fflush(stdout);
  }
  ReadAndCheckInput(cparams);
  if(Comm().MyPID()==0) cout << "done!" << endl;

  // check for FillComplete of discretization
  if (!Discret().Filled()) dserror("Discretization is not fillcomplete");

  // let's check for contact boundary conditions in discret
  // and detect groups of matching conditions
  // for each group, create a contact interface and store it
  if(Comm().MyPID()==0)
  {
    cout << "Building contact interface(s)...............";
    fflush(stdout);
  }

  vector<DRT::Condition*> contactconditions(0);
  Discret().GetCondition("Contact",contactconditions);

  // there must be more than one contact condition
  // unless we have a self contact problem!
  if ((int)contactconditions.size()<1)
    dserror("Not enough contact conditions in discretization");
  if ((int)contactconditions.size()==1)
  {
    const string* side = contactconditions[0]->Get<string>("Side");
    if (*side != "Selfcontact")
      dserror("Not enough contact conditions in discretization");
  }

  // find all pairs of matching contact conditions
  // there is a maximum of (conditions / 2) groups
  vector<int> foundgroups(0);
  int numgroupsfound = 0;

  // maximum dof number in discretization
  // later we want to create NEW Lagrange multiplier degrees of
  // freedom, which of course must not overlap with displacement dofs
  int maxdof = Discret().DofRowMap()->MaxAllGID();

  for (int i=0; i<(int)contactconditions.size(); ++i)
  {
    // initialize vector for current group of conditions and temp condition
    vector<DRT::Condition*> currentgroup(0);
    DRT::Condition* tempcond = NULL;

    // try to build contact group around this condition
    currentgroup.push_back(contactconditions[i]);
    const vector<int>* group1v = currentgroup[0]->Get<vector<int> >("contact id");
    if (!group1v) dserror("Contact Conditions does not have value 'contact id'");
    int groupid1 = (*group1v)[0];
    bool foundit = false;

    // only one surface per group is ok for self contact
    const string* side = contactconditions[i]->Get<string>("Side");
    if (*side == "Selfcontact") foundit = true;

    for (int j=0; j<(int)contactconditions.size(); ++j)
    {
      if (j==i) continue; // do not detect contactconditions[i] again
      tempcond = contactconditions[j];
      const vector<int>* group2v = tempcond->Get<vector<int> >("contact id");
      if (!group2v) dserror("Contact Conditions does not have value 'contact id'");
      int groupid2 = (*group2v)[0];
      if (groupid1 != groupid2) continue; // not in the group
      foundit = true; // found a group entry
      currentgroup.push_back(tempcond); // store it in currentgroup
    }

    // now we should have found a group of conds
    if (!foundit) dserror("Cannot find matching contact condition for id %d",groupid1);

    // see whether we found this group before
    bool foundbefore = false;
    for (int j=0; j<numgroupsfound; ++j)
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
    vector<const string*> sides((int)currentgroup.size());
    vector<bool> isslave((int)currentgroup.size());
    vector<bool> isself((int)currentgroup.size());

    for (int j=0;j<(int)sides.size();++j)
    {
      sides[j] = currentgroup[j]->Get<string>("Side");
      if (*sides[j] == "Slave")
      {
        hasslave = true;
        isslave[j] = true;
        isself[j] = false;
      }
      else if (*sides[j] == "Master")
      {
        hasmaster = true;
        isslave[j] = false;
        isself[j] = false;
      }
      else if (*sides[j] == "Selfcontact")
      {
        hasmaster = true;
        hasslave = true;
        isslave[j] = false;
        isself[j] = true;
      }
      else
        dserror("ERROR: CoManager: Unknown contact side qualifier!");
    }

    if (!hasslave) dserror("Slave side missing in contact condition group!");
    if (!hasmaster) dserror("Master side missing in contact condition group!");

    // check for self contact group
    if (isself[0])
    {
      for (int j=1;j<(int)isself.size();++j)
        if (!isself[j]) dserror("Inconsistent definition of self contact condition group!");
    }

    // find out which sides are initialized as Active
    vector<const string*> active((int)currentgroup.size());
    vector<bool> isactive((int)currentgroup.size());

    for (int j=0;j<(int)sides.size();++j)
    {
      active[j] = currentgroup[j]->Get<string>("Initialization");
      if (*sides[j] == "Slave")
      {
        // slave sides may be initialized as "Active" or as "Inactive"
        if (*active[j] == "Active")        isactive[j] = true;
        else if (*active[j] == "Inactive") isactive[j] = false;
        else                               dserror("ERROR: Unknown contact init qualifier!");
      }
      else if (*sides[j] == "Master")
      {
        // master sides must NOT be initialized as "Active" as this makes no sense
        if (*active[j] == "Active")        dserror("ERROR: Master side cannot be active!");
        else if (*active[j] == "Inactive") isactive[j] = false;
        else                               dserror("ERROR: Unknown contact init qualifier!");
      }
      else if (*sides[j] == "Selfcontact")
      {
        // Selfcontact surfs must NOT be initialized as "Active" as this makes no sense
        if (*active[j] == "Active")        dserror("ERROR: Selfcontact surface cannot be active!");
        else if (*active[j] == "Inactive") isactive[j] = false;
        else                               dserror("ERROR: Unknown contact init qualifier!");
      }
      else
        dserror("ERROR: CoManager: Unknown contact side qualifier!");
    }

    // create an empty interface and store it in this Manager
    // create an empty contact interface and store it in this Manager
    // (for structural contact we do NOT want redundant storage)
    // (the only exception is self contact where a redundant slave is needed)
    bool redundant = false;
    if (isself[0]) redundant = true;
    interfaces.push_back(rcp(new CONTACT::CoInterface(groupid1,Comm(),dim,cparams,isself[0],redundant)));

    // get it again
    RCP<CONTACT::CoInterface> interface = interfaces[(int)interfaces.size()-1];

    // note that the nodal ids are unique because they come from
    // one global problem discretization conatining all nodes of the
    // contact interface
    // We rely on this fact, therefore it is not possible to
    // do contact between two distinct discretizations here

    // collect all intial active nodes
    std::vector<int> initialactive;

    //-------------------------------------------------- process nodes
    for (int j=0;j<(int)currentgroup.size();++j)
    {
      // get all nodes and add them
      const vector<int>* nodeids = currentgroup[j]->Nodes();
      if (!nodeids) dserror("Condition does not have Node Ids");
      for (int k=0; k<(int)(*nodeids).size(); ++k)
      {
        int gid = (*nodeids)[k];
        // do only nodes that I have in my discretization
        if (!Discret().NodeColMap()->MyGID(gid)) continue;
        DRT::Node* node = Discret().gNode(gid);
        if (!node) dserror("Cannot find node with gid %",gid);

        // store initial active node gids
        if (isactive[j]) initialactive.push_back(gid);

        // find out if this node is initial active on another Condition
        // and do NOT overwrite this status then!
        bool foundinitialactive = false;
        if (!isactive[j])
        {
          for (int k=0;k<(int)initialactive.size();++k)
            if (gid == initialactive[k])
            {
              foundinitialactive = true;
              break;
            }
        }

        // create CoNode object or FriNode object in the frictional case
        INPAR::CONTACT::FrictionType ftype = Teuchos::getIntegralValue<INPAR::CONTACT::FrictionType>(cparams,"FRICTION");

        // for the boolean variable initactive we use isactive[j]+foundinitialactive,
        // as this is true for BOTH initial active nodes found for the first time
        // and found for the second, third, ... time!
        if (ftype != INPAR::CONTACT::friction_none)
        {
           RCP<CONTACT::FriNode> cnode = rcp(new CONTACT::FriNode(node->Id(),node->X(),
                                                             node->Owner(),
                                                             Discret().NumDof(node),
                                                             Discret().Dof(node),
                                                             isslave[j],isactive[j]+foundinitialactive));

          // note that we do not have to worry about double entries
          // as the AddNode function can deal with this case!
          // the only problem would have occured for the initial active nodes,
          // as their status could have been overwritten, but is prevented
          // by the "foundinitialactive" block above!
          interface->AddCoNode(cnode);
        }
        else
        {
          RCP<CONTACT::CoNode> cnode = rcp(new CONTACT::CoNode(node->Id(),node->X(),
                                                           node->Owner(),
                                                           Discret().NumDof(node),
                                                           Discret().Dof(node),
                                                           isslave[j],isactive[j]+foundinitialactive));
          // note that we do not have to worry about double entries
          // as the AddNode function can deal with this case!
          // the only problem would have occured for the initial active nodes,
          // as their status could have been overwritten, but is prevented
          // by the "foundinitialactive" block above!
          interface->AddCoNode(cnode);
        }
      }
    }

    //----------------------------------------------- process elements
    int ggsize = 0;
    for (int j=0;j<(int)currentgroup.size();++j)
    {
      // get elements from condition j of current group
      map<int,RCP<DRT::Element> >& currele = currentgroup[j]->Geometry();

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
      Comm().SumAll(&lsize,&gsize,1);


      map<int,RCP<DRT::Element> >::iterator fool;
      for (fool=currele.begin(); fool != currele.end(); ++fool)
      {
        RCP<DRT::Element> ele = fool->second;
        RCP<CONTACT::CoElement> cele = rcp(new CONTACT::CoElement(ele->Id()+ggsize,
                                                                ele->Owner(),
                                                                ele->Shape(),
                                                                ele->NumNode(),
                                                                ele->NodeIds(),
                                                                isslave[j]));
        interface->AddCoElement(cele);
      } // for (fool=ele1.start(); fool != ele1.end(); ++fool)

      ggsize += gsize; // update global element counter
    }

    //-------------------- finalize the contact interface construction
    interface->FillComplete(maxdof);

  } // for (int i=0; i<(int)contactconditions.size(); ++i)
  if(Comm().MyPID()==0) cout << "done!" << endl;

  //**********************************************************************
  // create the solver strategy object
  // and pass all necessary data to it
  if(Comm().MyPID()==0)
  {
    cout << "Building contact strategy object............";
    fflush(stdout);
  }
  INPAR::CONTACT::SolvingStrategy stype =
      Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(cparams,"STRATEGY");
  if (stype == INPAR::CONTACT::solution_lagmult)
    strategy_ = rcp(new CoLagrangeStrategy(problemrowmap,cparams,interfaces,dim,comm_,alphaf,maxdof));
  else if (stype == INPAR::CONTACT::solution_penalty)
    strategy_ = rcp(new CoPenaltyStrategy(problemrowmap,cparams,interfaces,dim,comm_,alphaf,maxdof));
  else if (stype == INPAR::CONTACT::solution_auglag)
    strategy_ = rcp(new CoPenaltyStrategy(problemrowmap,cparams,interfaces,dim,comm_,alphaf,maxdof));
  else
    dserror("Unrecognized strategy");
  if(Comm().MyPID()==0) cout << "done!" << endl;
  //**********************************************************************

  // print initial parallel redistribution
  for (int i=0; i<(int)interfaces.size();++i)
    interfaces[i]->PrintParallelDistribution(i+1);

  // create binary search tree
  for (int i=0; i<(int)interfaces.size();++i)
    interfaces[i]->CreateSearchTree();

  // print parameter list to screen
  if (Comm().MyPID()==0)
  {
    cout << "\ngiven parameters in list '" << GetStrategy().Params().name() << "':\n";
    cout << GetStrategy().Params() << endl;
  }

  return;
}


/*----------------------------------------------------------------------*
 |  read and check input parameters (public)                  popp 04/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoManager::ReadAndCheckInput(Teuchos::ParameterList& cparams)
{
  // read parameter list and problemtype from DRT::Problem
  const Teuchos::ParameterList& input = DRT::Problem::Instance()->MeshtyingAndContactParams();
  const Teuchos::ParameterList& psize = DRT::Problem::Instance()->ProblemSizeParams();
  std::string problemtype = DRT::Problem::Instance()->ProblemType();
  int dim = psize.get<int>("DIM");

  // *********************************************************************
  // this is mortar contact
  // *********************************************************************
  if (Teuchos::getIntegralValue<INPAR::CONTACT::ApplicationType>(input,"APPLICATION") != INPAR::CONTACT::app_mortarcontact)
    dserror("You should not be here...");

  // *********************************************************************
  // invalid parameter combinations
  // *********************************************************************
  if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(input,"STRATEGY") == INPAR::CONTACT::solution_penalty &&
                                                 input.get<double>("PENALTYPARAM") <= 0.0)
    dserror("Penalty parameter eps = 0, must be greater than 0");

  if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(input,"STRATEGY") == INPAR::CONTACT::solution_auglag &&
                                                 input.get<double>("PENALTYPARAM") <= 0.0)
    dserror("Penalty parameter eps = 0, must be greater than 0");

  if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(input,"STRATEGY") == INPAR::CONTACT::solution_auglag &&
                                                   input.get<int>("UZAWAMAXSTEPS") < 2)
    dserror("Maximum number of Uzawa / Augmentation steps must be at least 2");

  if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(input,"STRATEGY") == INPAR::CONTACT::solution_auglag &&
                                               input.get<double>("UZAWACONSTRTOL") <= 0.0)
    dserror("Constraint tolerance for Uzawa / Augmentation scheme must be greater than 0");

  if (Teuchos::getIntegralValue<INPAR::CONTACT::FrictionType>(input,"FRICTION") != INPAR::CONTACT::friction_none &&
                                            input.get<double>("SEMI_SMOOTH_CT") == 0.0)
    dserror("Parameter ct = 0, must be greater than 0 for frictional contact");

  if (Teuchos::getIntegralValue<INPAR::CONTACT::FrictionType>(input,"FRICTION") == INPAR::CONTACT::friction_tresca &&
                                                   input.get<double>("FRBOUND") <= 0.0)
    dserror("No valid Tresca friction bound provided, must be greater than 0");

  if (Teuchos::getIntegralValue<INPAR::CONTACT::FrictionType>(input,"FRICTION") == INPAR::CONTACT::friction_coulomb &&
                                                   input.get<double>("FRCOEFF") <= 0.0)
    dserror("No valid Coulomb friction coefficient provided, must be greater than 0");

  if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(input,"STRATEGY") == INPAR::CONTACT::solution_lagmult &&
      Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(input,"SHAPEFCN") == INPAR::MORTAR::shape_standard &&
      Teuchos::getIntegralValue<INPAR::CONTACT::SystemType>(input,"SYSTEM") == INPAR::CONTACT::system_condensed)
    dserror("Condensation of linear system only possible for dual Lagrange multipliers");

  if (Teuchos::getIntegralValue<INPAR::MORTAR::ParRedist>(input,"PARALLEL_REDIST") != INPAR::MORTAR::parredist_none &&
                                                     input.get<int>("MIN_ELEPROC") <  0)
    dserror("Minimum number of elements per processor for parallel redistribution must be >= 0");
  
  if (Teuchos::getIntegralValue<INPAR::MORTAR::ParRedist>(input,"PARALLEL_REDIST") == INPAR::MORTAR::parredist_dynamic &&
                                                     input.get<double>("MAX_BALANCE") <  1.0)
    dserror("Maximum allowed value of load balance for dynamic parallel redistribution must be >= 1.0");
  
  // *********************************************************************
  // not (yet) implemented combinations
  // *********************************************************************
  if (Teuchos::getIntegralValue<INPAR::CONTACT::FrictionType>(input,"FRICTION") == INPAR::CONTACT::friction_tresca &&
                                                                            dim == 3)
    dserror("3D frictional contact with Tresca's law not yet implemented");

  if (Teuchos::getIntegralValue<INPAR::CONTACT::FrictionType>(input,"FRICTION") != INPAR::CONTACT::friction_none &&
                     Teuchos::getIntegralValue<int>(input,"SEMI_SMOOTH_NEWTON") != 1 &&
                                                                            dim == 3)
    dserror("3D frictional contact only implemented with Semi-smooth Newton");

#ifndef CONTACTCOMPHUEBER
  if (Teuchos::getIntegralValue<INPAR::CONTACT::FrictionType>(input,"FRICTION") != INPAR::CONTACT::friction_none &&
                                                                            dim == 3)
    dserror("3D frictional contact without flag CONTACTCOMPHUEBER not yet implemented");
#endif

  if (Teuchos::getIntegralValue<int>(input,"CROSSPOINTS") == true && dim == 3)
    dserror("ERROR: Crosspoints / edge node modification not yet implemented for 3D");

  // check for self contact
  vector<DRT::Condition*> coco(0);
  Discret().GetCondition("Contact",coco);
  bool self = false;

  for (int k=0;k<(int)coco.size();++k)
  {
    const string* side = coco[k]->Get<string>("Side");
    if (*side == "Selfcontact") self = true;
  }

  if (self == true &&
      Teuchos::getIntegralValue<INPAR::MORTAR::ParRedist>(input,"PARALLEL_REDIST") != INPAR::MORTAR::parredist_none)
    dserror("ERROR: Self contact and parallel redistribution not yet compatible");
  
  // *********************************************************************
  // thermal-structure-interaction contact
  // *********************************************************************
  
  if (problemtype=="tsi" &&
      Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(input,"SHAPEFCN") != INPAR::MORTAR::shape_standard &&
      Teuchos::getIntegralValue<int>(input,"THERMOLAGMULT")==false)
    dserror("Thermal contact without Lagrange Multipliers only for standard shape functions");

  if (problemtype=="tsi" &&
      Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(input,"SHAPEFCN") == INPAR::MORTAR::shape_standard &&
      Teuchos::getIntegralValue<int>(input,"THERMOLAGMULT")==true)
    dserror("Thermal contact with Lagrange Multipliers only for dual shape functions");
  
  // no parallel redistribution in for thermal-structure-interaction
  if (problemtype=="tsi" && 
      Teuchos::getIntegralValue<INPAR::MORTAR::ParRedist>(input,"PARALLEL_REDIST") != INPAR::MORTAR::parredist_none)
    dserror("ERROR: Parallel redistribution not yet implemented for TSI problems");  

  // *********************************************************************
  // contact with wear
  // *********************************************************************
  
  if (Teuchos::getIntegralValue<INPAR::CONTACT::WearType>(input,"WEAR") == INPAR::CONTACT::wear_none &&
      input.get<double>("WEARCOEFF") != 0.0)
    dserror("ERROR: Wear coefficient only necessary in the context of wear.");
  
  if (Teuchos::getIntegralValue<INPAR::CONTACT::FrictionType>(input,"FRICTION") == INPAR::CONTACT::friction_none &&
      Teuchos::getIntegralValue<INPAR::CONTACT::WearType>(input,"WEAR") != INPAR::CONTACT::wear_none)
    dserror("ERROR: Wear models only applicable to frictional contact.");

  if (Teuchos::getIntegralValue<INPAR::CONTACT::WearType>(input,"WEAR") != INPAR::CONTACT::wear_none &&
      input.get<double>("WEARCOEFF") <= 0.0)
    dserror("ERROR: No valid wear coefficient provided, must be equal or greater 0.");

  // *********************************************************************
  // 3D quadratic mortar (choice of interpolation and testing fcts.)
  // *********************************************************************
  if ((Teuchos::getIntegralValue<INPAR::MORTAR::LagMultQuad3D>(input,"LAGMULT_QUAD3D") == INPAR::MORTAR::lagmult_quad_pwlin ||
       Teuchos::getIntegralValue<INPAR::MORTAR::LagMultQuad3D>(input,"LAGMULT_QUAD3D") == INPAR::MORTAR::lagmult_quad_lin) &&
       Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(input,"SHAPEFCN") == INPAR::MORTAR::shape_standard)
    dserror("No Petrov-Galerkin approach (for LM) implemented for 3D quadratic contact with STANDARD shape fct.");

  if ((Teuchos::getIntegralValue<INPAR::MORTAR::LagMultQuad3D>(input,"LAGMULT_QUAD3D") == INPAR::MORTAR::lagmult_pwlin_pwlin ||
      Teuchos::getIntegralValue<INPAR::MORTAR::LagMultQuad3D>(input,"LAGMULT_QUAD3D") == INPAR::MORTAR::lagmult_lin_lin ||
      Teuchos::getIntegralValue<INPAR::MORTAR::LagMultQuad3D>(input,"LAGMULT_QUAD3D") == INPAR::MORTAR::lagmult_quad_pwlin ||
      Teuchos::getIntegralValue<INPAR::MORTAR::LagMultQuad3D>(input,"LAGMULT_QUAD3D") == INPAR::MORTAR::lagmult_quad_lin) &&
      Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(input,"SHAPEFCN") == INPAR::MORTAR::shape_dual)
    dserror("Only quadratic/quadratic approach (for LM) implemented for 3D quadratic contact with DUAL shape fct.");

  // *********************************************************************
  // warnings
  // *********************************************************************
  if (input.get<double>("SEARCH_PARAM") == 0.0)
    cout << ("Warning: Contact search called without inflation of bounding volumes\n") << endl;

  // store ParameterList in local parameter list
  cparams = input;

  // no parallel redistribution in the serial case
  if (Comm().NumProc()==1) cparams.set<string>("PARALLEL_REDIST","None");

  return true;
}

/*----------------------------------------------------------------------*
 |  write restart information for contact (public)            popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoManager::WriteRestart(IO::DiscretizationWriter& output)
{
  // quantities to be written for restart
  RCP<Epetra_Vector> activetoggle;
  RCP<Epetra_Vector> sliptoggle;
  RCP<Epetra_Vector> weightedwear;

  // quantities to be written for restart
  GetStrategy().DoWriteRestart(activetoggle,sliptoggle,weightedwear);

  // write restart information for contact
  output.WriteVector("lagrmultold",GetStrategy().LagrMultOld());
  output.WriteVector("activetoggle",activetoggle);

  // friction
  if(GetStrategy().Friction())
    output.WriteVector("sliptoggle",sliptoggle);
  
  if (weightedwear != Teuchos::null)
    output.WriteVector("weightedwear", weightedwear);

  return;
}

/*----------------------------------------------------------------------*
 |  read restart information for contact (public)             popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoManager::ReadRestart(IO::DiscretizationReader& reader,
                                     RCP<Epetra_Vector> dis, RCP<Epetra_Vector> zero)
{
  // this is contact, thus we need the displacement state for restart
  // let strategy object do all the work
  GetStrategy().DoReadRestart(reader, dis);

  return;
}

/*----------------------------------------------------------------------*
 |  write interface tractions for postprocessing (public)     popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoManager::PostprocessTractions(IO::DiscretizationWriter& output)
{

  // *********************************************************************
  // contact tractions
  // *********************************************************************
  
  // evaluate contact tractions
  GetStrategy().OutputStresses();
  
  // problem row map  
  RCP<Epetra_Map> problem = GetStrategy().ProblemRowMap();

  // normal direction
  RCP<Epetra_Vector> normalstresses = GetStrategy().ContactNorStress();
  RCP<Epetra_Vector> normalstressesexp = rcp(new Epetra_Vector(*problem));
  LINALG::Export(*normalstresses, *normalstressesexp);

  // tangential plane
  RCP<Epetra_Vector> tangentialstresses = GetStrategy().ContactTanStress();
  RCP<Epetra_Vector> tangentialstressesexp = rcp(new Epetra_Vector(*problem));
  LINALG::Export(*tangentialstresses, *tangentialstressesexp);
  
  // write to output
  // contact tractions in normal and tangential direction
  output.WriteVector("norcontactstress",normalstressesexp);
  output.WriteVector("tancontactstress",tangentialstressesexp);

#ifdef CONTACTFORCEOUTPUT  

  // *********************************************************************
  // contact forces on slave non master side,
  // in normal and tangential direction
  // *********************************************************************
 
  // vectors for contact forces
  RCP<Epetra_Vector> fcslavenor = rcp(new Epetra_Vector(GetStrategy().DMatrix()->RowMap()));
  RCP<Epetra_Vector> fcslavetan = rcp(new Epetra_Vector(GetStrategy().DMatrix()->RowMap()));
  RCP<Epetra_Vector> fcmasternor = rcp(new Epetra_Vector(GetStrategy().MMatrix()->DomainMap()));
  RCP<Epetra_Vector> fcmastertan = rcp(new Epetra_Vector(GetStrategy().MMatrix()->DomainMap()));
  
  // vectors with problemrowmap
  RCP<Epetra_Vector> fcslavenorexp = rcp(new Epetra_Vector(*problem));
  RCP<Epetra_Vector> fcslavetanexp = rcp(new Epetra_Vector(*problem));
  RCP<Epetra_Vector> fcmasternorexp = rcp(new Epetra_Vector(*problem));
  RCP<Epetra_Vector> fcmastertanexp = rcp(new Epetra_Vector(*problem));

  // multiplication
  GetStrategy().DMatrix()->Multiply(true, *normalstresses, *fcslavenor);
  GetStrategy().DMatrix()->Multiply(true, *tangentialstresses, *fcslavetan);
  GetStrategy().MMatrix()->Multiply(true, *normalstresses, *fcmasternor);
  GetStrategy().MMatrix()->Multiply(true, *tangentialstresses, *fcmastertan);

  // export  
  LINALG::Export(*fcslavenor,*fcslavenorexp);
  LINALG::Export(*fcslavetan,*fcslavetanexp);
  LINALG::Export(*fcmasternor,*fcmasternorexp);
  LINALG::Export(*fcmastertan,*fcmastertanexp);

  // contact forces on slave and master side
  output.WriteVector("norslaveforce",fcslavenorexp);
  output.WriteVector("tanslaveforce",fcslavetanexp);
  output.WriteVector("normasterforce",fcmasternorexp);
  output.WriteVector("tanmasterforce",fcmastertanexp);
#endif
 
  // *********************************************************************
  // wear
  // *********************************************************************
 
  bool wear = GetStrategy().Wear();
  if (wear)
  {
    // evaluate wear (not weighted)
    GetStrategy().OutputWear();

    // write output
    RCP<Epetra_Vector> wearoutput = GetStrategy().ContactWear();
    RCP<Epetra_Vector> wearoutputexp = rcp(new Epetra_Vector(*problem));
    LINALG::Export(*wearoutput, *wearoutputexp);
    output.WriteVector("wear",wearoutputexp);
  }
 return;
}

#endif  // #ifdef CCADISCRET
