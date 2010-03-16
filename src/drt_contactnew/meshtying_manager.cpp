/*!----------------------------------------------------------------------
\file meshtying_manager.cpp
\brief BACI implementation of main class to control all meshtying

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
#include "meshtying_manager.H"
#include "meshtying_lagrange_strategy.H"
#include "meshtying_penalty_strategy.H"
#include "meshtying_defines.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_mortar.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 03/08|
 *----------------------------------------------------------------------*/
CONTACT::MtManager::MtManager(DRT::Discretization& discret, double alphaf) :
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
  if (dim!= 2 && dim!=3) dserror("ERROR: Meshtying problem must be 2D or 3D");
  vector<RCP<MORTAR::MortarInterface> > interfaces;
  Teuchos::ParameterList mtparams;

  // read and check meshtying input parameters
  if(Comm().MyPID()==0)
  {
    cout << "Checking meshtying input parameters...........";
    fflush(stdout);
  }
  ReadAndCheckInput(mtparams);
  if(Comm().MyPID()==0) cout << "done!" << endl;

  // check for FillComplete of discretization
  if (!Discret().Filled()) dserror("Discretization is not fillcomplete");

  // let's check for meshtying boundary conditions in discret
  // and detect groups of matching conditions
  // for each group, create a contact interface and store it
  if(Comm().MyPID()==0)
  {
    cout << "Building meshtying interface(s)...............";
    fflush(stdout);
  }
  
  vector<DRT::Condition*> contactconditions(0);
  Discret().GetCondition("Contact",contactconditions);
  
  // there must be more than one meshtying condition
  if ((int)contactconditions.size()<2)
    dserror("Not enough contact conditions in discretization");

  // find all pairs of matching meshtying conditions
  // there is a maximum of (conditions / 2) groups
  vector<int> foundgroups(0);
  int numgroupsfound = 0;

  for (int i=0; i<(int)contactconditions.size(); ++i)
  {
    // initialize vector for current group of conditions and temp condition
    vector<DRT::Condition*> currentgroup(0);
    DRT::Condition* tempcond = NULL;

    // try to build meshtying group around this condition
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
        dserror("ERROR: MtManager: Unknown contact side qualifier!");
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
        // slave sides must be initialized as "Active"
        if (*active[j] == "Active")        isactive[j] = true;
        else if (*active[j] == "Inactive") dserror("ERROR: Slave side must be active for meshtying!");
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
        dserror("ERROR: MtManager: Unknown contact side qualifier!");
    }

    // create an empty meshtying interface and store it in this Manager
    interfaces.push_back(rcp(new MORTAR::MortarInterface(groupid1,Comm(),dim,mtparams)));
    
    // get it again
    RCP<MORTAR::MortarInterface> interface = interfaces[(int)interfaces.size()-1];

    // note that the nodal ids are unique because they come from
    // one global problem discretization conatining all nodes of the
    // contact interface
    // We rely on this fact, therefore it is not possible to
    // do meshtying between two distinct discretizations here

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

        // create MortarNode object
        RCP<MORTAR::MortarNode> mtnode = rcp(new MORTAR::MortarNode(node->Id(),node->X(),
                                                              node->Owner(),
                                                              Discret().NumDof(node),
                                                              Discret().Dof(node),
                                                              isslave[j]));

        // note that we do not have to worry about double entries
        // as the AddNode function can deal with this case!
        interface->AddMortarNode(mtnode);
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
        RCP<MORTAR::MortarElement> mtele = rcp(new MORTAR::MortarElement(ele->Id()+ggsize,
                                                                   DRT::Element::element_mortar,
                                                                   ele->Owner(),
                                                                   ele->Shape(),
                                                                   ele->NumNode(),
                                                                   ele->NodeIds(),
                                                                   isslave[j]));
        interface->AddMortarElement(mtele);
      } // for (fool=ele1.start(); fool != ele1.end(); ++fool)

      ggsize += gsize; // update global element counter
    }

    //-------------------- finalize the meshtying interface construction
    interface->FillComplete();
    
    //---------------------------------------- create binary search tree
    interface->CreateSearchTree();

  } // for (int i=0; i<(int)contactconditions.size(); ++i)
  if(Comm().MyPID()==0) cout << "done!" << endl;
  
  // create the solver strategy object
  // and pass all necessary data to it
  if(Comm().MyPID()==0)
  {
    cout << "Building meshtying strategy object............";
    fflush(stdout);
  }
  INPAR::CONTACT::SolvingStrategy stype =
      Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(mtparams,"STRATEGY");
  if (stype == INPAR::CONTACT::solution_lagmult)
    strategy_ = rcp(new MtLagrangeStrategy(Discret(),problemrowmap,mtparams,interfaces,dim,comm_,alphaf)); 
  else if (stype == INPAR::CONTACT::solution_penalty)
    strategy_ = rcp(new MtPenaltyStrategy(Discret(),problemrowmap,mtparams,interfaces,dim,comm_,alphaf));
  else if (stype == INPAR::CONTACT::solution_auglag)
    strategy_ = rcp(new MtPenaltyStrategy(Discret(),problemrowmap,mtparams,interfaces,dim,comm_,alphaf));
  else
    dserror("Unrecognized strategy");
  if(Comm().MyPID()==0) cout << "done!" << endl << endl;
  // **** initialization of row/column maps moved to AbstractStrategy **** //
  // since the manager does not operate over nodes, elements, dofs anymore 
  // ********************************************************************* //

  // print parameter list to screen
  //if (Comm().MyPID()==0) cout << GetStrategy().Params() << endl;
  
  return;
}


/*----------------------------------------------------------------------*
 |  read and check input parameters (public)                  popp 04/08|
 *----------------------------------------------------------------------*/
bool CONTACT::MtManager::ReadAndCheckInput(Teuchos::ParameterList& mtparams)
{
  // read parameter list from DRT::Problem
  const Teuchos::ParameterList& input = DRT::Problem::Instance()->MeshtyingAndContactParams();
  const Teuchos::ParameterList& psize = DRT::Problem::Instance()->ProblemSizeParams();
  int dim = psize.get<int>("DIM");
  
  // *********************************************************************
  // this is mortar meshtying
  // *********************************************************************
  if (Teuchos::getIntegralValue<INPAR::CONTACT::ApplicationType>(input,"APPLICATION") != INPAR::CONTACT::app_mortarmeshtying)
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
                                                   input.get<int>("UZAWAMAXSTEPS") <  2)
    dserror("Maximum number of Uzawa / Augmentation steps must be at least 2");
  
  if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(input,"STRATEGY") == INPAR::CONTACT::solution_auglag &&
                                               input.get<double>("UZAWACONSTRTOL") <= 0.0)
    dserror("Constraint tolerance for Uzawa / Augmentation scheme must be greater than 0");
  
  if (Teuchos::getIntegralValue<INPAR::CONTACT::ApplicationType>(input,"APPLICATION") == INPAR::CONTACT::app_mortarmeshtying &&
            Teuchos::getIntegralValue<INPAR::CONTACT::FrictionType>(input,"FRICTION") != INPAR::CONTACT::friction_none)
    dserror("Friction law supplied for mortar meshtying");

  if (Teuchos::getIntegralValue<INPAR::MORTAR::SearchAlgorithm>(input,"SEARCH_ALGORITHM") == INPAR::MORTAR::search_bfnode &&
                                                        input.get<double>("SEARCH_PARAM") == 0.0)
    dserror("Search radius sp = 0, must be greater than 0 for node-based search");
  
  // *********************************************************************
  // not (yet) implemented combinations
  // *********************************************************************
  if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(input,"STRATEGY") == INPAR::CONTACT::solution_lagmult &&
      Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(input,"SHAPEFCN") != INPAR::MORTAR::shape_dual )
    dserror("Lagrange multiplier strategy only implemented for dual shape fct.");
  
  if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(input,"STRATEGY") == INPAR::CONTACT::solution_penalty &&
      Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(input,"SHAPEFCN") != INPAR::MORTAR::shape_standard )
    dserror("Penalty strategy only implemented for standard shape fct.");
  
  if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(input,"STRATEGY") == INPAR::CONTACT::solution_auglag &&
      Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(input,"SHAPEFCN") != INPAR::MORTAR::shape_standard )
    dserror("Augmented Lagrange strategy only implemented for standard shape fct.");
  
  if (Teuchos::getIntegralValue<int>(input,"CROSSPOINTS") == true && dim == 3)
    dserror("ERROR: Crosspoints / edge node modification not yet implemented for 3D");
  
  // *********************************************************************
  // 3D quadratic mortar (choice of interpolation and testing fcts.)
  // *********************************************************************
  if ((Teuchos::getIntegralValue<INPAR::MORTAR::LagMultQuad3D>(input,"LAGMULT_QUAD3D") == INPAR::MORTAR::lagmult_quad_pwlin ||
       Teuchos::getIntegralValue<INPAR::MORTAR::LagMultQuad3D>(input,"LAGMULT_QUAD3D") == INPAR::MORTAR::lagmult_quad_lin) &&
       Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(input,"SHAPEFCN") == INPAR::MORTAR::shape_standard)
    dserror("No Petrov-Galerkin approach (for LM) implemented for 3D quadratic meshtying with STANDARD shape fct.");
  
  if ((Teuchos::getIntegralValue<INPAR::MORTAR::LagMultQuad3D>(input,"LAGMULT_QUAD3D") == INPAR::MORTAR::lagmult_pwlin_pwlin ||
       Teuchos::getIntegralValue<INPAR::MORTAR::LagMultQuad3D>(input,"LAGMULT_QUAD3D") == INPAR::MORTAR::lagmult_lin_lin ||
       Teuchos::getIntegralValue<INPAR::MORTAR::LagMultQuad3D>(input,"LAGMULT_QUAD3D") == INPAR::MORTAR::lagmult_quad_pwlin ||
       Teuchos::getIntegralValue<INPAR::MORTAR::LagMultQuad3D>(input,"LAGMULT_QUAD3D") == INPAR::MORTAR::lagmult_quad_lin) &&
       Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(input,"SHAPEFCN") == INPAR::MORTAR::shape_dual)
    dserror("Only quadratic/quadratic approach (for LM) implemented for 3D quadratic meshtying with DUAL shape fct.");
  
  // *********************************************************************
  // warnings
  // *********************************************************************
  if ((Teuchos::getIntegralValue<INPAR::MORTAR::SearchAlgorithm>(input,"SEARCH_ALGORITHM") == INPAR::MORTAR::search_bfele ||
       Teuchos::getIntegralValue<INPAR::MORTAR::SearchAlgorithm>(input,"SEARCH_ALGORITHM") == INPAR::MORTAR::search_binarytree) &&
                                                         input.get<double>("SEARCH_PARAM") == 0.0)
    cout << ("Warning: Ele-based / binary tree search called without inflation of bounding volumes\n") << endl;

  // store ParameterList in local parameter list
  mtparams = input;

  return true;
}

/*----------------------------------------------------------------------*
 |  write restart information for meshtying (public)          popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtManager::WriteRestart(IO::DiscretizationWriter& output)
{
  // write restart information for meshtying
  output.WriteVector("lagrmultold",GetStrategy().LagrMultOld());

  return;
}

/*----------------------------------------------------------------------*
 |  read restart information for meshtying (public)           popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtManager::ReadRestart(IO::DiscretizationReader& reader,
                                     RCP<Epetra_Vector> dis)
{
  // let strategy object do all the work
  GetStrategy().DoReadRestart(reader, dis);
  
  return;
}

/*----------------------------------------------------------------------*
 |  write interface tractions for postprocessing (public)     popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtManager::PostprocessTractions(IO::DiscretizationWriter& output)
{
  // evaluate interface tractions
  RCP<Epetra_Map> problem = GetStrategy().ProblemRowMap();
  RCP<Epetra_Vector> traction = rcp(new Epetra_Vector(*(GetStrategy().LagrMultOld())));
  RCP<Epetra_Vector> tractionexp = rcp(new Epetra_Vector(*problem));
  LINALG::Export(*traction, *tractionexp);
  
  // write to output
  output.WriteVector("interfacetraction",tractionexp);
  
  return;
}

#endif  // #ifdef CCADISCRET
