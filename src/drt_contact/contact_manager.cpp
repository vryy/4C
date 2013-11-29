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
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "contact_manager.H"
#include "contact_interface.H"
#include "contact_node.H"
#include "contact_element.H"
#include "contact_lagrange_strategy.H"
#include "contact_wear_lagrange_strategy.H"
#include "contact_wear_interface.H"
#include "contact_penalty_strategy.H"
#include "contact_defines.H"
#include "friction_node.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_inpar/inpar_wear.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_io/io_control.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 03/08|
 *----------------------------------------------------------------------*/
CONTACT::CoManager::CoManager(DRT::Discretization& discret, double alphaf) :
MORTAR::ManagerBase(),
discret_(discret)
{
  // overwrite base class communicator
  comm_ = Teuchos::rcp(Discret().Comm().Clone());

  // welcome message

  // create some local variables (later to be stored in strategy)
  int dim = DRT::Problem::Instance()->NDim();
  if (dim!= 2 && dim!=3) dserror("ERROR: Contact problem must be 2D or 3D");
  std::vector<Teuchos::RCP<CONTACT::CoInterface> > interfaces;
  Teuchos::ParameterList cparams;

  // read and check contact input parameters
  if(Comm().MyPID()==0)
  {
    std::cout << "Checking contact input parameters...........";
    fflush(stdout);
  }
  ReadAndCheckInput(cparams);
  if(Comm().MyPID()==0) std::cout << "done!" << std::endl;

  // check for FillComplete of discretization
  if (!Discret().Filled()) dserror("Discretization is not fillcomplete");

  // let's check for contact boundary conditions in discret
  // and detect groups of matching conditions
  // for each group, create a contact interface and store it
  if(Comm().MyPID()==0)
  {
    std::cout << "Building contact interface(s)...............";
    fflush(stdout);
  }

  std::vector<DRT::Condition*> contactconditions(0);
  Discret().GetCondition("Mortar",contactconditions);

  // there must be more than one contact condition
  // unless we have a self contact problem!
  if ((int)contactconditions.size()<1)
    dserror("Not enough contact conditions in discretization");
  if ((int)contactconditions.size()==1)
  {
    const std::string* side = contactconditions[0]->Get<std::string>("Side");
    if (*side != "Selfcontact")
      dserror("Not enough contact conditions in discretization");
  }

  // find all pairs of matching contact conditions
  // there is a maximum of (conditions / 2) groups
  std::vector<int> foundgroups(0);
  int numgroupsfound = 0;

  // maximum dof number in discretization
  // later we want to create NEW Lagrange multiplier degrees of
  // freedom, which of course must not overlap with displacement dofs
  int maxdof = Discret().DofRowMap()->MaxAllGID();

  // get input par.
  INPAR::CONTACT::SolvingStrategy stype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(cparams,"STRATEGY");
  INPAR::CONTACT::WearLaw wlaw =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::WearLaw>(cparams,"WEARLAW");
  INPAR::CONTACT::WearType wtype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::WearType>(cparams,"WEARTYPE");

  bool friplus = false;
  if ((wlaw!=INPAR::CONTACT::wear_none) || (cparams.get<int>("PROBTYPE")==INPAR::CONTACT::tsi))
    friplus=true;

  for (int i=0; i<(int)contactconditions.size(); ++i)
  {
    // initialize vector for current group of conditions and temp condition
    std::vector<DRT::Condition*> currentgroup(0);
    DRT::Condition* tempcond = NULL;

    // try to build contact group around this condition
    currentgroup.push_back(contactconditions[i]);
    const std::vector<int>* group1v = currentgroup[0]->Get<std::vector<int> >("Interface ID");
    if (!group1v) dserror("Contact Conditions does not have value 'Interface ID'");
    int groupid1 = (*group1v)[0];
    bool foundit = false;

    // only one surface per group is ok for self contact
    const std::string* side = contactconditions[i]->Get<std::string>("Side");
    if (*side == "Selfcontact") foundit = true;

    for (int j=0; j<(int)contactconditions.size(); ++j)
    {
      if (j==i) continue; // do not detect contactconditions[i] again
      tempcond = contactconditions[j];
      const std::vector<int>* group2v = tempcond->Get<std::vector<int> >("Interface ID");
      if (!group2v) dserror("Contact Conditions does not have value 'Interface ID'");
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
    std::vector<const std::string*> sides((int)currentgroup.size());
    std::vector<bool> isslave((int)currentgroup.size());
    std::vector<bool> isself((int)currentgroup.size());

    for (int j=0;j<(int)sides.size();++j)
    {
      sides[j] = currentgroup[j]->Get<std::string>("Side");
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
    std::vector<const std::string*> active((int)currentgroup.size());
    std::vector<bool> isactive((int)currentgroup.size());

    for (int j=0;j<(int)sides.size();++j)
    {
      active[j] = currentgroup[j]->Get<std::string>("Initialization");
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

    // create interface local parameter list (copy)
    Teuchos::ParameterList icparams = cparams;

    // find out if interface-specific coefficients of friction are given
    INPAR::CONTACT::FrictionType fric = DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(cparams,"FRICTION");
    if (fric == INPAR::CONTACT::friction_tresca || fric == INPAR::CONTACT::friction_coulomb)
    {
      // read interface COFs
      std::vector<double> frcoeff((int)currentgroup.size());
      for (int j=0;j<(int)currentgroup.size();++j)
        frcoeff[j] = currentgroup[j]->GetDouble("FrCoeffOrBound");

      // check consistency of interface COFs
      for (int j=1;j<(int)currentgroup.size();++j)
        if (frcoeff[j] != frcoeff[0])
          dserror("ERROR: Inconsistency in friction coefficients of interface %i",groupid1);

      // check for infeasible value of COF
      if (frcoeff[0]<0.0)
        dserror("ERROR: Negative FrCoeff / FrBound on interface %i",groupid1);

      // add COF locally to contact parameter list of this interface
      if (fric == INPAR::CONTACT::friction_tresca)
      {
        icparams.setEntry("FRBOUND",static_cast<Teuchos::ParameterEntry>(frcoeff[0]));
        icparams.setEntry("FRCOEFF",static_cast<Teuchos::ParameterEntry>(-1.0));
      }
      else if (fric == INPAR::CONTACT::friction_coulomb)
      {
        icparams.setEntry("FRCOEFF",static_cast<Teuchos::ParameterEntry>(frcoeff[0]));
        icparams.setEntry("FRBOUND",static_cast<Teuchos::ParameterEntry>(-1.0));
      }
    }

    // create an empty interface and store it in this Manager
    // create an empty contact interface and store it in this Manager
    // (for structural contact we currently choose redundant master storage)
    // (the only exception is self contact where a redundant slave is needed, too)
    INPAR::MORTAR::RedundantStorage redundant = DRT::INPUT::IntegralValue<INPAR::MORTAR::RedundantStorage>(icparams,"REDUNDANT_STORAGE");
//    if (isself[0]==false && redundant != INPAR::MORTAR::redundant_master)
//      dserror("ERROR: CoManager: Contact requires redundant master storage");
    if (isself[0]==true && redundant != INPAR::MORTAR::redundant_all)
      dserror("ERROR: CoManager: Self contact requires redundant slave and master storage");

    // decide between contactinterface and wearinterface
    Teuchos::RCP<CONTACT::CoInterface> newinterface=Teuchos::null;
    if(wlaw!=INPAR::CONTACT::wear_none)
      newinterface=Teuchos::rcp(new CONTACT::WearInterface(groupid1,Comm(),dim,icparams,isself[0],redundant));
    else
      newinterface=Teuchos::rcp(new CONTACT::CoInterface(groupid1,Comm(),dim,icparams,isself[0],redundant));
    interfaces.push_back(newinterface);

    // get it again
    Teuchos::RCP<CONTACT::CoInterface> interface = interfaces[(int)interfaces.size()-1];

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
      const std::vector<int>* nodeids = currentgroup[j]->Nodes();
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
        INPAR::CONTACT::FrictionType ftype = DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(cparams,"FRICTION");

        // for the boolean variable initactive we use isactive[j]+foundinitialactive,
        // as this is true for BOTH initial active nodes found for the first time
        // and found for the second, third, ... time!
        if (ftype != INPAR::CONTACT::friction_none)
        {
           Teuchos::RCP<CONTACT::FriNode> cnode = Teuchos::rcp(new CONTACT::FriNode(node->Id(),node->X(),
                                                             node->Owner(),
                                                             Discret().NumDof(node),
                                                             Discret().Dof(node),
                                                             isslave[j],isactive[j]+foundinitialactive,
                                                             friplus));

#ifdef CONTACTCONSTRAINTXYZ
           // Check, if this node (and, in case, which dofs) are in the contact symmetry condition
           std::vector<DRT::Condition*> contactSymconditions(0);
           Discret().GetCondition("mrtrsym",contactSymconditions);

           for (unsigned j=0; j<contactSymconditions.size(); j++)
             if (contactSymconditions.at(j)->ContainsNode(node->Id()))
             {
               const std::vector<int>*    onoff  = contactSymconditions.at(j)->Get<std::vector<int> >("onoff");
               for (unsigned k=0; k<onoff->size(); k++)
                 if (onoff->at(k)==1)
                   cnode->DbcDofs()[k]=true;
             }
#else
           std::vector<DRT::Condition*> contactSymconditions(0);
           Discret().GetCondition("mrtrsym",contactSymconditions);
           if (contactSymconditions.size()!=0)
             dserror("Contact symmetry condition only with CONTACTCONSTRAINTXYZ flag");
#endif

          // note that we do not have to worry about double entries
          // as the AddNode function can deal with this case!
          // the only problem would have occured for the initial active nodes,
          // as their status could have been overwritten, but is prevented
          // by the "foundinitialactive" block above!
          interface->AddCoNode(cnode);
        }
        else
        {
          Teuchos::RCP<CONTACT::CoNode> cnode = Teuchos::rcp(new CONTACT::CoNode(node->Id(),node->X(),
                                                           node->Owner(),
                                                           Discret().NumDof(node),
                                                           Discret().Dof(node),
                                                           isslave[j],isactive[j]+foundinitialactive));

          // Check, if this node (and, in case, which dofs) are in the contact symmetry condition
          std::vector<DRT::Condition*> contactSymconditions(0);
          Discret().GetCondition("mrtrsym",contactSymconditions);

          for (unsigned j=0; j<contactSymconditions.size(); j++)
            if (contactSymconditions.at(j)->ContainsNode(node->Id()))
            {
              const std::vector<int>*    onoff  = contactSymconditions.at(j)->Get<std::vector<int> >("onoff");
              for (unsigned k=0; k<onoff->size(); k++)
                if (onoff->at(k)==1)
                  cnode->DbcDofs()[k]=true;
            }

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
      std::map<int,Teuchos::RCP<DRT::Element> >& currele = currentgroup[j]->Geometry();

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


      std::map<int,Teuchos::RCP<DRT::Element> >::iterator fool;
      for (fool=currele.begin(); fool != currele.end(); ++fool)
      {
        Teuchos::RCP<DRT::Element> ele = fool->second;
        Teuchos::RCP<CONTACT::CoElement> cele = Teuchos::rcp(new CONTACT::CoElement(ele->Id()+ggsize,
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
  if(Comm().MyPID()==0) std::cout << "done!" << std::endl;

  //**********************************************************************
  // create the solver strategy object
  // and pass all necessary data to it
  if(Comm().MyPID()==0)
  {
    std::cout << "Building contact strategy object............";
    fflush(stdout);
  }

  // create WearLagrangeStrategy for wear as non-distinct quantity
  if (stype == INPAR::CONTACT::solution_lagmult && wlaw!=INPAR::CONTACT::wear_none &&
      (wtype==INPAR::CONTACT::wear_expl || wtype==INPAR::CONTACT::wear_impl || wtype==INPAR::CONTACT::wear_discr))
    strategy_ = Teuchos::rcp(new WearLagrangeStrategy(Discret(),cparams,interfaces,dim,comm_,alphaf,maxdof));
  else if (stype == INPAR::CONTACT::solution_lagmult)
    strategy_ = Teuchos::rcp(new CoLagrangeStrategy(Discret(),cparams,interfaces,dim,comm_,alphaf,maxdof));
  else if (stype == INPAR::CONTACT::solution_penalty)
    strategy_ = Teuchos::rcp(new CoPenaltyStrategy(Discret(),cparams,interfaces,dim,comm_,alphaf,maxdof));
  else if (stype == INPAR::CONTACT::solution_auglag)
    strategy_ = Teuchos::rcp(new CoPenaltyStrategy(Discret(),cparams,interfaces,dim,comm_,alphaf,maxdof));
  else
    dserror("Unrecognized strategy");
  if(Comm().MyPID()==0) std::cout << "done!" << std::endl;
  //**********************************************************************

  // print friction information of interfaces
  INPAR::CONTACT::FrictionType fric = DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(cparams,"FRICTION");
  if (Comm().MyPID()==0)
  {
    for (int i=0; i<(int)interfaces.size();++i)
    {
      double checkfrcoeff = 0.0;
      if (fric == INPAR::CONTACT::friction_tresca)
      {
        checkfrcoeff = interfaces[i]->IParams().get<double>("FRBOUND");
        std::cout << std::endl << "Interface         " << i+1 << std::endl;
        std::cout <<         "FrBound (Tresca)  " << checkfrcoeff << std::endl;
      }
      else if (fric == INPAR::CONTACT::friction_coulomb)
      {
        checkfrcoeff = interfaces[i]->IParams().get<double>("FRCOEFF");
        std::cout << std::endl << "Interface         " << i+1 << std::endl;
        std::cout <<         "FrCoeff (Coulomb) " << checkfrcoeff << std::endl;
      }
    }
  }

  // print initial parallel redistribution
  for (int i=0; i<(int)interfaces.size();++i)
    interfaces[i]->PrintParallelDistribution(i+1);

  // create binary search tree
  for (int i=0; i<(int)interfaces.size();++i)
    interfaces[i]->CreateSearchTree();

  // show default parameters
  if (Comm().MyPID()==0)
  {
    std::cout << std::endl;
    DRT::INPUT::PrintDefaultParameters(IO::cout,GetStrategy().Params());
  }

  return;
}


/*----------------------------------------------------------------------*
 |  read and check input parameters (public)                  popp 04/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoManager::ReadAndCheckInput(Teuchos::ParameterList& cparams)
{
  // read parameter lists from DRT::Problem
  const Teuchos::ParameterList& mortar   = DRT::Problem::Instance()->MortarCouplingParams();
  const Teuchos::ParameterList& contact  = DRT::Problem::Instance()->ContactDynamicParams();
  const Teuchos::ParameterList& wearlist = DRT::Problem::Instance()->WearParams();
  const Teuchos::ParameterList& tsic     = DRT::Problem::Instance()->TSIContactParams();

  // structure params only for delta_t
  const Teuchos::ParameterList& stru     = DRT::Problem::Instance()->StructuralDynamicParams();

  // read Problem Type and Problem Dimension from DRT::Problem
  const PROBLEM_TYP problemtype = DRT::Problem::Instance()->ProblemType();
  int dim = DRT::Problem::Instance()->NDim();

  // *********************************************************************
  // this is mortar contact
  // *********************************************************************
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::ApplicationType>(contact,"APPLICATION") != INPAR::CONTACT::app_mortarcontact)
    dserror("You should not be here...");

  // *********************************************************************
  // invalid parameter combinations
  // *********************************************************************
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") != INPAR::CONTACT::solution_lagmult &&
              DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"SHAPEFCN") == INPAR::MORTAR::shape_petrovgalerkin)
    dserror("Petrov-Galerkin approach for LM only with Lagrange multiplier strategy");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") == INPAR::CONTACT::solution_penalty &&
                                                 contact.get<double>("PENALTYPARAM") <= 0.0)
    dserror("Penalty parameter eps = 0, must be greater than 0");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") == INPAR::CONTACT::solution_penalty &&
         DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") != INPAR::CONTACT::friction_none &&
                                              contact.get<double>("PENALTYPARAMTAN") <= 0.0)
    dserror("Tangential penalty parameter eps = 0, must be greater than 0");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") == INPAR::CONTACT::solution_auglag &&
                                                 contact.get<double>("PENALTYPARAM") <= 0.0)
    dserror("Penalty parameter eps = 0, must be greater than 0");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") == INPAR::CONTACT::solution_auglag &&
         DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") != INPAR::CONTACT::friction_none &&
                                              contact.get<double>("PENALTYPARAMTAN") <= 0.0)
    dserror("Tangential penalty parameter eps = 0, must be greater than 0");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") == INPAR::CONTACT::solution_auglag &&
                                                   contact.get<int>("UZAWAMAXSTEPS") < 2)
    dserror("Maximum number of Uzawa / Augmentation steps must be at least 2");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") == INPAR::CONTACT::solution_auglag &&
                                               contact.get<double>("UZAWACONSTRTOL") <= 0.0)
    dserror("Constraint tolerance for Uzawa / Augmentation scheme must be greater than 0");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") != INPAR::CONTACT::friction_none &&
                                            contact.get<double>("SEMI_SMOOTH_CT") == 0.0)
    dserror("Parameter ct = 0, must be greater than 0 for frictional contact");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") == INPAR::CONTACT::solution_lagmult &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"SHAPEFCN") == INPAR::MORTAR::shape_standard &&
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(contact,"SYSTEM") == INPAR::CONTACT::system_condensed)
    dserror("Condensation of linear system only possible for dual Lagrange multipliers");

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(mortar,"PARALLEL_REDIST") != INPAR::MORTAR::parredist_none &&
                                                     mortar.get<int>("MIN_ELEPROC") <  0)
    dserror("Minimum number of elements per processor for parallel redistribution must be >= 0");
  
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(mortar,"PARALLEL_REDIST") == INPAR::MORTAR::parredist_dynamic &&
                                                     mortar.get<double>("MAX_BALANCE") <  1.0)
    dserror("Maximum allowed value of load balance for dynamic parallel redistribution must be >= 1.0");

  if(DRT::INPUT::IntegralValue<int>(mortar,"LM_DUAL_CONSISTENT")==true &&
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") != INPAR::CONTACT::solution_lagmult &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"SHAPEFCN") != INPAR::MORTAR::shape_standard )
    dserror("ERROR: Consistent dual shape functions in boundary elements only for Lagrange multiplier strategy.");

  if(DRT::INPUT::IntegralValue<int>(mortar,"LM_DUAL_CONSISTENT")==true &&
       DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(mortar,"INTTYPE") == INPAR::MORTAR::inttype_elements &&
       (DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"SHAPEFCN") == INPAR::MORTAR::shape_dual))
    dserror("ERROR: Consistent dual shape functions in boundary elements not for purely element-based integration.");

  if(DRT::INPUT::IntegralValue<int>(mortar,"LM_NODAL_SCALE")==true &&
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") != INPAR::CONTACT::solution_lagmult)
    dserror("ERROR: Nodal scaling of Lagrange multipliers only for Lagrange multiplier strategy.");

  if(DRT::INPUT::IntegralValue<int>(mortar,"LM_NODAL_SCALE")==true &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(mortar,"INTTYPE") == INPAR::MORTAR::inttype_elements)
   dserror("ERROR: Nodal scaling of Lagrange multipliers not for purely element-based integration.");

  // *********************************************************************
  // not (yet) implemented combinations
  // *********************************************************************
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") == INPAR::CONTACT::friction_tresca &&
                                                                            dim == 3)
    dserror("3D frictional contact with Tresca's law not yet implemented");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") != INPAR::CONTACT::friction_none &&
                     DRT::INPUT::IntegralValue<int>(contact,"SEMI_SMOOTH_NEWTON") != 1 &&
                                                                            dim == 3)
    dserror("3D frictional contact only implemented with Semi-smooth Newton");

  if (DRT::INPUT::IntegralValue<int>(mortar,"CROSSPOINTS") == true && dim == 3)
    dserror("ERROR: Crosspoints / edge node modification not yet implemented for 3D");

  if (DRT::INPUT::IntegralValue<int>(mortar,"CROSSPOINTS") == true &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar,"LAGMULT_QUAD") == INPAR::MORTAR::lagmult_lin_lin)
    dserror("ERROR: Crosspoints and linear LM interpolation for quadratic FE not yet compatible");

  // check for self contact
  std::vector<DRT::Condition*> coco(0);
  Discret().GetCondition("Mortar",coco);
  bool self = false;

  for (int k=0;k<(int)coco.size();++k)
  {
    const std::string* side = coco[k]->Get<std::string>("Side");
    if (*side == "Selfcontact") self = true;
  }

  if (self == true &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(mortar,"PARALLEL_REDIST") != INPAR::MORTAR::parredist_none)
    dserror("ERROR: Self contact and parallel redistribution not yet compatible");
  
  if(DRT::INPUT::IntegralValue<int>(contact,"INITCONTACTBYGAP")==true &&
     contact.get<double>("INITCONTACTGAPVALUE") == 0.0)
    dserror("ERROR: For initialization of init contact with gap, the INITCONTACTGAPVALUE is needed."); 

  if(DRT::INPUT::IntegralValue<int>(mortar,"LM_DUAL_CONSISTENT")==true &&
     DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar,"LAGMULT_QUAD") != INPAR::MORTAR::lagmult_undefined)
    dserror("ERROR: Consistent dual shape functions in boundary elements only for linear shape functions.");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") == INPAR::CONTACT::friction_tresca &&
      DRT::INPUT::IntegralValue<int>(contact,"FRLESS_FIRST")==true)
    dserror("Frictionless first contact step with Tresca's law not yet implemented"); // hopefully coming soon, when Coulomb and Tresca are combined

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearLaw>(wearlist,"WEARLAW") != INPAR::CONTACT::wear_none &&
      DRT::INPUT::IntegralValue<int>(contact,"FRLESS_FIRST")==true)
    dserror("Frictionless first contact step with wear not yet implemented");


  // *********************************************************************
  // thermal-structure-interaction contact
  // *********************************************************************
  
  if (problemtype==prb_tsi &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"SHAPEFCN") != INPAR::MORTAR::shape_standard &&
      DRT::INPUT::IntegralValue<int>(tsic,"THERMOLAGMULT")==false)
    dserror("Thermal contact without Lagrange Multipliers only for standard shape functions");

  if (problemtype==prb_tsi &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"SHAPEFCN") == INPAR::MORTAR::shape_standard &&
      DRT::INPUT::IntegralValue<int>(tsic,"THERMOLAGMULT")==true)
    dserror("Thermal contact with Lagrange Multipliers only for dual shape functions");
  
  // no parallel redistribution in for thermal-structure-interaction
  if (problemtype==prb_tsi &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(mortar,"PARALLEL_REDIST") != INPAR::MORTAR::parredist_none)
    dserror("ERROR: Parallel redistribution not yet implemented for TSI problems");  

  // no nodal scaling in for thermal-structure-interaction
  if (problemtype==prb_tsi &&
      DRT::INPUT::IntegralValue<int>(mortar,"LM_NODAL_SCALE")==true)
    dserror("ERROR: Nodal scaling not yet implemented for TSI problems");

  // *********************************************************************
  // contact with wear
  // *********************************************************************
  
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearLaw>(wearlist,"WEARLAW") == INPAR::CONTACT::wear_none &&
      wearlist.get<double>("WEARCOEFF") != 0.0)
    dserror("ERROR: Wear coefficient only necessary in the context of wear.");
  
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") == INPAR::CONTACT::friction_none &&
      DRT::INPUT::IntegralValue<INPAR::CONTACT::WearLaw>(wearlist,"WEARLAW") != INPAR::CONTACT::wear_none)
    dserror("ERROR: Wear models only applicable to frictional contact.");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearLaw>(wearlist,"WEARLAW") != INPAR::CONTACT::wear_none &&
      wearlist.get<double>("WEARCOEFF") <= 0.0)
    dserror("ERROR: No valid wear coefficient provided, must be equal or greater 0.");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearLaw>(wearlist,"WEARLAW") != INPAR::CONTACT::wear_none &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(mortar,"INTTYPE") != INPAR::MORTAR::inttype_segments)
    dserror("ERROR: Calculation of wear only possible by employing segment-based integration!");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearLaw>(wearlist,"WEARLAW") != INPAR::CONTACT::wear_none &&
      DRT::INPUT::IntegralValue<int>(mortar,"LM_NODAL_SCALE") ==true)
    dserror("ERROR: Combination of LM_NODAL_SCALE and WEAR not (yet) implemented.");

  if(DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") != INPAR::CONTACT::solution_lagmult &&
     DRT::INPUT::IntegralValue<INPAR::CONTACT::WearLaw>(wearlist,"WEARLAW") != INPAR::CONTACT::wear_none)
    dserror("ERROR: Wear model only applicable in combination with Lagrange multiplier strategy.");

  if(DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") == INPAR::CONTACT::friction_tresca &&
      DRT::INPUT::IntegralValue<INPAR::CONTACT::WearLaw>(wearlist,"WEARLAW") != INPAR::CONTACT::wear_none)
    dserror("ERROR: Wear only for Coulomb friction!");

  // *********************************************************************
  // 3D quadratic mortar (choice of interpolation and testing fcts.)
  // *********************************************************************
  if ((DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar,"LAGMULT_QUAD") == INPAR::MORTAR::lagmult_quad_pwlin ||
       DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar,"LAGMULT_QUAD") == INPAR::MORTAR::lagmult_quad_lin) &&
       DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"SHAPEFCN") == INPAR::MORTAR::shape_standard)
    dserror("Only quad/quad, lin/lin, pwlin/pwlin (for LM) implemented for quadratic contact with STANDARD shape fct.");

  if ((DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar,"LAGMULT_QUAD") == INPAR::MORTAR::lagmult_pwlin_pwlin ||
      DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar,"LAGMULT_QUAD") == INPAR::MORTAR::lagmult_lin_lin ||
      DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar,"LAGMULT_QUAD") == INPAR::MORTAR::lagmult_quad_pwlin ||
      DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar,"LAGMULT_QUAD") == INPAR::MORTAR::lagmult_quad_lin) &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"SHAPEFCN") == INPAR::MORTAR::shape_dual)
    dserror("Only quadratic/quadratic approach (for LM) implemented for quadratic contact with DUAL shape fct.");

#ifdef MORTARTRAFO
    dserror("MORTARTRAFO not yet implemented for contact, only for meshtying");
#endif // #ifndef MORTARTRAFO

  // *********************************************************************
  // element-based vs. segment-based mortar integration
  // *********************************************************************
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(mortar,"INTTYPE") == INPAR::MORTAR::inttype_segments &&
                                         mortar.get<int>("NUMGP_PER_DIM") != 0)
    dserror("It is not possible to choose a Gauss rule with NUMGP_PER_DIM for segment-based integration."
            "\nSegment-based integration always uses pre-defined default values (5 GP per segment in 2D "
            "\nand 7 GP per triangular cell in 3D). Ask a 'contact person' if you want to change this.");

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(mortar,"INTTYPE") == INPAR::MORTAR::inttype_elements &&
                                         mortar.get<int>("NUMGP_PER_DIM") <= 0)
    dserror("Invalid Gauss point number NUMGP_PER_DIM for element-based integration.");

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(mortar,"INTTYPE") == INPAR::MORTAR::inttype_elements_BS &&
                                         mortar.get<int>("NUMGP_PER_DIM") <= 0)
    dserror("Invalid Gauss point number NUMGP_PER_DIM for element-based integration with boundary segmentation."
            "\nPlease note that the value you have to provide only applies to the element-based integration"
            "\ndomain, while pre-defined default values will be used in the segment-based boundary domain.");

  // *********************************************************************
  // warnings
  // *********************************************************************
  if (mortar.get<double>("SEARCH_PARAM") == 0.0 && Comm().MyPID()==0)
    std::cout << ("Warning: Contact search called without inflation of bounding volumes\n") << std::endl;

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearSide>(wearlist,"BOTH_SIDED_WEAR") !=  INPAR::CONTACT::wear_slave)
    std::cout << ("\n \n Warning: Contact with both-sided wear is still experimental !") << std::endl;

  // store contents of BOTH ParameterLists in local parameter list
  cparams.setParameters(mortar);
  cparams.setParameters(contact);
  cparams.setParameters(wearlist);
  cparams.setParameters(tsic);
  cparams.set<double>("TIMESTEP",stru.get<double>("TIMESTEP"));
  cparams.setName("CONTACT DYNAMIC / MORTAR COUPLING");

  // store relevant problem types
  if (problemtype==prb_structure)
      cparams.set<int> ("PROBTYPE",INPAR::CONTACT::structure);
  else if (problemtype==prb_tsi)
    cparams.set<int> ("PROBTYPE",INPAR::CONTACT::tsi);
  else if (problemtype==prb_struct_ale)
    cparams.set<int> ("PROBTYPE",INPAR::CONTACT::structalewear);
  else
    cparams.set<int> ("PROBTYPE",INPAR::CONTACT::other);

  // no parallel redistribution in the serial case
  if (Comm().NumProc()==1) cparams.set<std::string>("PARALLEL_REDIST","None");

  return true;
}

/*----------------------------------------------------------------------*
 |  write restart information for contact (public)            popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoManager::WriteRestart(IO::DiscretizationWriter& output, bool forcedrestart)
{
  // quantities to be written for restart
  Teuchos::RCP<Epetra_Vector> activetoggle;
  Teuchos::RCP<Epetra_Vector> sliptoggle;
  Teuchos::RCP<Epetra_Vector> weightedwear;
  Teuchos::RCP<Epetra_Vector> realwear;


  // quantities to be written for restart
  GetStrategy().DoWriteRestart(activetoggle,sliptoggle,weightedwear,realwear, forcedrestart);

  // export restart information for contact to problem dof row map
  Teuchos::RCP<Epetra_Map> problemdofs = GetStrategy().ProblemDofs();
  Teuchos::RCP<Epetra_Vector> lagrmultoldexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  LINALG::Export(*(GetStrategy().LagrMultOld()),*lagrmultoldexp);
  Teuchos::RCP<Epetra_Vector> activetoggleexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  LINALG::Export(*activetoggle,*activetoggleexp);

  // write restart information for contact
  output.WriteVector("lagrmultold",lagrmultoldexp);
  output.WriteVector("activetoggle",activetoggleexp);

  // friction
  if(GetStrategy().Friction())
  {
    Teuchos::RCP<Epetra_Vector> sliptoggleexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
    LINALG::Export(*sliptoggle,*sliptoggleexp);
    output.WriteVector("sliptoggle",sliptoggleexp);
  }
  
  // weighted wear
  if (weightedwear != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> weightedwearexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
    LINALG::Export(*weightedwear,*weightedwearexp);
    output.WriteVector("weightedwear", weightedwearexp);
  }

  // unweighted  wear
  if (realwear != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> realwearexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
    LINALG::Export(*realwear,*realwearexp);
    output.WriteVector("realwear", realwearexp);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  read restart information for contact (public)             popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoManager::ReadRestart(IO::DiscretizationReader& reader,
                                     Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<Epetra_Vector> zero)
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
  // active contact set and slip set
  // *********************************************************************

  // evaluate active set and slip set
  Teuchos::RCP<Epetra_Vector> activeset = Teuchos::rcp(new Epetra_Vector(*GetStrategy().ActiveRowNodes()));
  activeset->PutScalar(1.0);
  if (GetStrategy().Friction())
  {
    Teuchos::RCP<Epetra_Vector> slipset = Teuchos::rcp(new Epetra_Vector(*GetStrategy().SlipRowNodes()));
    slipset->PutScalar(1.0);
    Teuchos::RCP<Epetra_Vector> slipsetexp = Teuchos::rcp(new Epetra_Vector(*GetStrategy().ActiveRowNodes()));
    LINALG::Export(*slipset, *slipsetexp);
    activeset->Update(1.0,*slipsetexp,1.0);
  }

  // export to problem node row map
  Teuchos::RCP<Epetra_Map> problemnodes = GetStrategy().ProblemNodes();
  Teuchos::RCP<Epetra_Vector> activesetexp = Teuchos::rcp(new Epetra_Vector(*problemnodes));
  LINALG::Export(*activeset,*activesetexp);

  output.WriteVector("activeset",activesetexp);

  // *********************************************************************
  // contact tractions
  // *********************************************************************
  
  // evaluate contact tractions
  GetStrategy().OutputStresses();
  
  // export to problem dof row map
  Teuchos::RCP<Epetra_Map> problemdofs = GetStrategy().ProblemDofs();

  // normal direction
  Teuchos::RCP<Epetra_Vector> normalstresses = GetStrategy().ContactNorStress();
  Teuchos::RCP<Epetra_Vector> normalstressesexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  LINALG::Export(*normalstresses, *normalstressesexp);

  // tangential plane
  Teuchos::RCP<Epetra_Vector> tangentialstresses = GetStrategy().ContactTanStress();
  Teuchos::RCP<Epetra_Vector> tangentialstressesexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
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
  Teuchos::RCP<Epetra_Vector> fcslavenor = Teuchos::rcp(new Epetra_Vector(GetStrategy().DMatrix()->RowMap()));
  Teuchos::RCP<Epetra_Vector> fcslavetan = Teuchos::rcp(new Epetra_Vector(GetStrategy().DMatrix()->RowMap()));
  Teuchos::RCP<Epetra_Vector> fcmasternor = Teuchos::rcp(new Epetra_Vector(GetStrategy().MMatrix()->DomainMap()));
  Teuchos::RCP<Epetra_Vector> fcmastertan = Teuchos::rcp(new Epetra_Vector(GetStrategy().MMatrix()->DomainMap()));
  
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
  //BEGIN: to output the global ID's of the master nodes in contact - devaal 02.2011
  
  int dim = DRT::Problem::Instance()->NDim();

  if (dim == 2)
    dserror("Only working for 3D");

  std::vector<int>  lnid, gnid;

  //std::cout << "MasterNor" << fcmasternor->MyLength() << std::endl;

  for (int i=0; i<fcmasternor->MyLength(); i=i+3)
  {
   
    //check if master node in contact
    if (sqrt(((*fcmasternor)[i])*((*fcmasternor)[i])+((*fcmasternor)[i+1])*((*fcmasternor)[i+1])+((*fcmasternor)[i+2])*((*fcmasternor)[i]+2)) > 0.00001)
    {
      lnid.push_back((fcmasternor->Map()).GID(i)/3);
     }
  }

  // we want to gather data from on all procs
  std::vector<int> allproc(Comm().NumProc());
  for (int i=0; i<Comm().NumProc(); ++i) allproc[i] = i;

  // communicate all data to proc 0
  LINALG::Gather<int>(lnid,gnid,(int)allproc.size(),&allproc[0],Comm());
  
  //std::cout << " size of gnid:" << gnid.size() << std::endl;
  
  ////////////////
  ///// attempt at obtaining the nid and relative displacement u of master nodes in contact - devaal
  // define my own interface
  MORTAR::StrategyBase& myStrategy = GetStrategy();
  CoAbstractStrategy& myContactStrategy = static_cast<CoAbstractStrategy&>(myStrategy);
  
  std::vector<Teuchos::RCP<CONTACT::CoInterface> > myInterface = myContactStrategy.ContactInterfaces();
  
  //check interface size - just doing this now for a single interface
  
  if (myInterface.size() != 1)
    dserror("Interface size should be 1");
  
  std::cout << "OUTPUT OF MASTER NODE IN CONTACT" << std::endl;
  //std::cout << "Master_node_in_contact x_dis y_dis z_dis" << std::endl;
  for (int i=0; i<(int)gnid.size(); ++i)
  {
      int myGid = gnid[i];
      std::cout << gnid[i] << std::endl; // << " " << myUx << " " << myUy << " " << myUz << std::endl;
  }
  
#endif  //MASTERNODESINCONTACT: to output the global ID's of the master nodes in contact 
    
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

#ifdef CONTACTEXPORT
  // export averaged node forces to xxx.force
  double resultnor[fcslavenor->NumVectors()];
  double resulttan[fcslavetan->NumVectors()];
  fcslavenor->Norm2(resultnor);
  fcslavetan->Norm2(resulttan);

  if(Comm().MyPID()==0)
  {
    std::cout << "resultnor= " << resultnor[0] << std::endl;
    std::cout << "resulttan= " << resulttan[0] << std::endl;

    FILE* MyFile = NULL;
    std::ostringstream filename;
    const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
    filename << filebase << ".force";
    MyFile = fopen(filename.str().c_str(), "at+");
    if (MyFile)
    {
      //fprintf(MyFile,valuename.c_str());
      fprintf(MyFile, "%g\t",resultnor[0]);
      fprintf(MyFile, "%g\n",resulttan[0]);
      fclose(MyFile);
    }
    else
      dserror("ERROR: File for Output could not be opened.");
  }
#endif //CONTACTEXPORT

#endif //CONTACTFORCEOUTPUT
 
 
  // *********************************************************************
  // wear
  // *********************************************************************
  bool wear = GetStrategy().Wear();
  if (wear)
  {

    // ***************************************************************************
    // we do not compute the non-weighted wear here. we just write    farah 06/13
    // the output. the non-weighted wear will be used as dirichlet-b.
    // for the ale problem. n.w.wear will be called in stru_ale_algorithm.cpp
    // and computed in GetStrategy().OutputWear();
    // ***************************************************************************

    // evaluate wear (not weighted)
    //GetStrategy().OutputWear();

    // write output
    Teuchos::RCP<Epetra_Vector> wearoutput = GetStrategy().ContactWear();
    Teuchos::RCP<Epetra_Vector> wearoutputexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
    LINALG::Export(*wearoutput, *wearoutputexp);
    output.WriteVector("wear",wearoutputexp);
  }
 return;
}

