/*!----------------------------------------------------------------------
\file drt_contact_manager.cpp
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
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "drt_contact_manager.H"
#include "drt_contact_lagrange_strategy.H"
#include "drt_contact_penalty_strategy.H"
#include "drt_cnode.H"
#include "drt_celement.H"
#include "contactdefines.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/linalg_utils.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 03/08|
 *----------------------------------------------------------------------*/
CONTACT::Manager::Manager(DRT::Discretization& discret, double alphaf) :
CONTACT::ManagerBase(),
discret_(discret)
{
  // overwrite base class communicator
  comm_ = rcp(Discret().Comm().Clone());
  
  // create some local variables (later to be stored in strategy)
  RCP<Epetra_Map> problemrowmap = rcp(new Epetra_Map(*(Discret().DofRowMap())));
  const Teuchos::ParameterList& psize = DRT::Problem::Instance()->ProblemSizeParams();
  int dim = psize.get<int>("DIM");
  if (dim!= 2 && dim!=3) dserror("ERROR: Contact problem must be 2D or 3D");
  vector<RCP<CONTACT::Interface> > interfaces;
  Teuchos::ParameterList cparams;

  // read and check contact input parameters
  ReadAndCheckInput(cparams);

  // check for FillComplete of discretization
  if (!Discret().Filled()) dserror("Discretization is not fillcomplete");

  // let's check for contact boundary conditions in discret
  // and detect groups of matching conditions
  // for each group, create a contact interface and store it
  vector<DRT::Condition*> contactconditions(0);
  Discret().GetCondition("Contact",contactconditions);
  if ((int)contactconditions.size()<=1)  dserror("Not enough contact conditions in discretization");

  // find all pairs of matching contact conditions
  // there is a maximum of (conditions / 2) groups
  vector<int> foundgroups(0);
  int numgroupsfound = 0;

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

    // create an empty interface and store it local vector
    interfaces.push_back(rcp(new CONTACT::Interface(groupid1,Comm(),dim,cparams)));
        
    // get it again
    RCP<CONTACT::Interface> interface = interfaces[(int)interfaces.size()-1];

    // find out which sides are Master and Slave
    bool hasslave = false;
    bool hasmaster = false;
    vector<const string*> sides((int)currentgroup.size());
    vector<bool> isslave((int)currentgroup.size());

    for (int j=0;j<(int)sides.size();++j)
    {
      sides[j] = currentgroup[j]->Get<string>("Side");
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
        dserror("ERROR: ContactManager: Unknown contact side qualifier!");
    }

    if (!hasslave) dserror("Slave side missing in contact condition group!");
    if (!hasmaster) dserror("Master side missing in contact condition group!");

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
      else
        dserror("ERROR: ContactManager: Unknown contact side qualifier!");
    }

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

        // create CNode object
        // for the boolean variable initactive we use isactive[j]+foundinitialactive,
        // as this is true for BOTH initial active nodes found for the first time
        // and found for the second, third, ... time!
        RCP<CONTACT::CNode> cnode = rcp(new CONTACT::CNode(node->Id(),node->X(),
                                                           node->Owner(),
                                                           Discret().NumDof(node),
                                                           Discret().Dof(node),
                                                           isslave[j],isactive[j]+foundinitialactive));

        // note that we do not have to worry about double entries
        // as the AddNode function can deal with this case!
        // the only problem would have occured for the initial active nodes,
        // as their status could have been overwritten, but is prevented
        // by the "foundinitialactive" block above!
        interface->AddCNode(cnode);
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
        RCP<CONTACT::CElement> cele = rcp(new CONTACT::CElement(ele->Id()+ggsize,
                                                                DRT::Element::element_contact,
                                                                ele->Owner(),
                                                                ele->Shape(),
                                                                ele->NumNode(),
                                                                ele->NodeIds(),
                                                                isslave[j]));
        interface->AddCElement(cele);
      } // for (fool=ele1.start(); fool != ele1.end(); ++fool)

      ggsize += gsize; // update global element counter
    }

    //-------------------- finalize the contact interface construction
    interface->FillComplete();

  } // for (int i=0; i<(int)contactconditions.size(); ++i)

  // create the solver strategy object
  // and pass all necessary data to it
  INPAR::CONTACT::SolvingStrategy stype =
      Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(cparams,"STRATEGY");
  if (stype == INPAR::CONTACT::solution_lagmult)
    strategy_ = rcp(new LagrangeStrategy(problemrowmap,cparams,interfaces,dim,comm_,alphaf)); 
  else if (stype == INPAR::CONTACT::solution_penalty)
    strategy_ = rcp(new PenaltyStrategy(problemrowmap,cparams,interfaces,dim,comm_,alphaf));
  else if (stype == INPAR::CONTACT::solution_auglag)
    dserror("Cannot cope with augmented lagrange strategy yet");
  else
    dserror("Unrecognized strategy");
    
  // **** initialization of row/column maps moved to AbstractStrategy **** //
  // since the manager does not operate over nodes, elements, dofs anymore 
  // ********************************************************************* //

  return;
}


/*----------------------------------------------------------------------*
 |  read and check input parameters (public)                  popp 04/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Manager::ReadAndCheckInput(Teuchos::ParameterList& cparams)
{
  // read parameter list from DRT::Problem
  const Teuchos::ParameterList& input = DRT::Problem::Instance()->StructuralContactParams();

  // invalid parameter combinations
  if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(input,"STRATEGY") == INPAR::CONTACT::solution_auglag)
    dserror("Augmented Lagrange strategy not yet implemented");
  
  if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(input,"STRATEGY") == INPAR::CONTACT::solution_lagmult &&
      Teuchos::getIntegralValue<INPAR::CONTACT::ShapeFcn>(input,"SHAPEFCN") != INPAR::CONTACT::shape_dual )
      dserror("Lagrange multiplier strategy only implemented for dual shape fct.");
  
  if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(input,"STRATEGY") == INPAR::CONTACT::solution_penalty &&
      Teuchos::getIntegralValue<INPAR::CONTACT::ShapeFcn>(input,"SHAPEFCN") != INPAR::CONTACT::shape_standard )
      dserror("Penalty strategy only implemented for standard shape fct.");
  
  if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(input,"STRATEGY") == INPAR::CONTACT::solution_penalty &&
                                                 input.get<double>("PENALTYPARAM") <= 0.0)
      dserror("Penalty parameter eps = 0, must be greater than 0");
  
  if (Teuchos::getIntegralValue<INPAR::CONTACT::ContactType>(input,"CONTACT")   == INPAR::CONTACT::contact_normal &&
      Teuchos::getIntegralValue<INPAR::CONTACT::ContactFrictionType>(input,"FRICTION") != INPAR::CONTACT::friction_none)
    dserror("Friction law supplied for normal contact");

  if (Teuchos::getIntegralValue<INPAR::CONTACT::ContactType>(input,"CONTACT")   == INPAR::CONTACT::contact_frictional &&
      Teuchos::getIntegralValue<INPAR::CONTACT::ContactFrictionType>(input,"FRICTION") == INPAR::CONTACT::friction_none)
    dserror("No friction law supplied for frictional contact");

  if (Teuchos::getIntegralValue<INPAR::CONTACT::ContactType>(input,"CONTACT") == INPAR::CONTACT::contact_frictional &&
                                          input.get<double>("SEMI_SMOOTH_CT") == 0.0)
  	dserror("Friction Parameter ct = 0, must be greater than 0");

  if (Teuchos::getIntegralValue<INPAR::CONTACT::ContactFrictionType>(input,"FRICTION") == INPAR::CONTACT::friction_tresca &&
                                                   input.get<double>("FRBOUND") <= 0.0)
    dserror("No valid Tresca friction bound provided, must be greater than 0");

  if (Teuchos::getIntegralValue<INPAR::CONTACT::ContactFrictionType>(input,"FRICTION") == INPAR::CONTACT::friction_coulomb &&
                                                   input.get<double>("FRCOEFF") <= 0.0)
    dserror("No valid Coulomb friction coefficient provided, must be greater than 0");

  if (Teuchos::getIntegralValue<INPAR::CONTACT::ContactType>(input,"CONTACT")   == INPAR::CONTACT::contact_meshtying &&
      Teuchos::getIntegralValue<INPAR::CONTACT::ContactFrictionType>(input,"FRICTION") != INPAR::CONTACT::friction_stick)
    dserror("Friction law other than stick supplied for mesh tying");

  if (Teuchos::getIntegralValue<INPAR::CONTACT::ContactSearchAlgorithm>(input,"SEARCH_ALGORITHM") == INPAR::CONTACT::search_bfnode &&
                                                                input.get<double>("SEARCH_PARAM") == 0.0)
    dserror("Search radius sp = 0, must be greater than 0 for node-based search");

#ifdef CONTACTRELVELMATERIAL
  // check full linearization
  bool fulllin   = Teuchos::getIntegralValue<int>(input,"FULL_LINEARIZATION");
  if (fulllin) dserror ("Full linearization only running for evaluating\n"
  		                  "the relative velocity with change of projection");
#endif

  // warnings
  if ((Teuchos::getIntegralValue<INPAR::CONTACT::ContactSearchAlgorithm>(input,"SEARCH_ALGORITHM") == INPAR::CONTACT::search_bfele ||
      Teuchos::getIntegralValue<INPAR::CONTACT::ContactSearchAlgorithm>(input,"SEARCH_ALGORITHM")  == INPAR::CONTACT::search_binarytree) &&
                                                                 input.get<double>("SEARCH_PARAM") == 0.0)
    cout << ("Warning: Ele-based / binary tree search called without inflation of bounding volumes\n") << endl;

  // store ParameterList in local parameter list
  cparams = input;

  return true;
}

/*----------------------------------------------------------------------*
 |  write restart information for contact (public)            popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::WriteRestart(IO::DiscretizationWriter& output)
{
  // quantities to be written for restart
  RCP<Epetra_Vector> activetoggle;
  RCP<Epetra_Vector> sliptoggle;

  // quantities to be written for restart
  GetStrategy().DoWriteRestart(activetoggle, sliptoggle);
    
  // write restart information for contact
  output.WriteVector("lagrmultold",GetStrategy().LagrMultOld());
  output.WriteVector("activetoggle",activetoggle);
  output.WriteVector("sliptoggle",sliptoggle);

  return;
}

/*----------------------------------------------------------------------*
 |  read restart information for contact (public)             popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::ReadRestart(IO::DiscretizationReader& reader,
                                   RCP<Epetra_Vector> dis)
{
  // read restart information for contact
  reader.ReadVector(GetStrategy().LagrMultOld(),"lagrmultold");
  reader.ReadVector(GetStrategy().LagrMult(),"lagrmultold");
  
  GetStrategy().StoreNodalQuantities(AbstractStrategy::lmold);
  GetStrategy().StoreNodalQuantities(AbstractStrategy::lmcurrent);

  RCP<Epetra_Vector> activetoggle =rcp(new Epetra_Vector(*(GetStrategy().SlaveRowNodes())));
  reader.ReadVector(activetoggle,"activetoggle");

  RCP<Epetra_Vector> sliptoggle =rcp(new Epetra_Vector(*(GetStrategy().SlaveRowNodes())));
  reader.ReadVector(sliptoggle,"sliptoggle");

  GetStrategy().DoReadRestart(activetoggle, sliptoggle, dis);
  
  return;
}


#endif  // #ifdef CCADISCRET
