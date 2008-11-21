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
#include "drt_cnode.H"
#include "drt_celement.H"
#include "contactdefines.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_lib/linalg_utils.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 03/08|
 *----------------------------------------------------------------------*/
CONTACT::Manager::Manager(DRT::Discretization& discret, double alphaf) :
CONTACT::ManagerBase(),
discret_(discret)
{
  // overwrite some base class variables
  alphaf_ = alphaf;
  comm_ = rcp(Discret().Comm().Clone());
  problemrowmap_ = rcp(new Epetra_Map(*(Discret().DofRowMap())));
  
  // get problem dimension (2D or 3D) and store into dim_
  const Teuchos::ParameterList& psize = DRT::Problem::Instance()->ProblemSizeParams();
  dim_=psize.get<int>("DIM");
  if (Dim()!= 2 && Dim()!=3) dserror("ERROR: Contact problem must be 2D or 3D");
  
  // read and check contact input parameters
  ReadAndCheckInput();
  
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

    // create an empty interface and store it in this Manager
    interface_.push_back(rcp(new CONTACT::Interface(groupid1,Comm(),Dim())));
    
    // get it again
    RCP<CONTACT::Interface> interface = interface_[(int)interface_.size()-1];

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
  
#ifdef DEBUG
  //for (int i=0;i<(int)interface_.size();++i)
  //  cout << *interface_[i];
#endif // #ifdef DEBUG
  
  // setup global slave / master Epetra_Maps
  // (this is done by looping over all interfaces and merging)
  for (int i=0;i<(int)interface_.size();++i)
  {
  	gsnoderowmap_ = LINALG::MergeMap(gsnoderowmap_,interface_[i]->SlaveRowNodes());
  	gsdofrowmap_ = LINALG::MergeMap(gsdofrowmap_,interface_[i]->SlaveRowDofs());
  	gmdofrowmap_ = LINALG::MergeMap(gmdofrowmap_,interface_[i]->MasterRowDofs());
  }
  
  // setup global dof map of non-slave-ormaster dofs
  // (this is done by splitting from the dicretization dof map) 
  gndofrowmap_ = LINALG::SplitMap(*(discret_.DofRowMap()),*gsdofrowmap_);
  gndofrowmap_ = LINALG::SplitMap(*gndofrowmap_,*gmdofrowmap_);
  
  // initialize active sets and slip sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i=0;i<(int)interface_.size();++i)
  {
    interface_[i]->InitializeActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_,interface_[i]->ActiveNodes(),false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_,interface_[i]->ActiveDofs(),false);
    gactiven_ = LINALG::MergeMap(gactiven_,interface_[i]->ActiveNDofs(),false);
    gactivet_ = LINALG::MergeMap(gactivet_,interface_[i]->ActiveTDofs(),false);
    gslipnodes_ = LINALG::MergeMap(gslipnodes_,interface_[i]->SlipNodes(),false);
    gslipdofs_ = LINALG::MergeMap(gslipdofs_,interface_[i]->SlipDofs(),false);
    gslipt_ = LINALG::MergeMap(gslipt_,interface_[i]->SlipTDofs(),false);
  }
    
  // setup Lagrange muliplier vectors
  z_          = rcp(new Epetra_Vector(*gsdofrowmap_));
  zold_       = rcp(new Epetra_Vector(*gsdofrowmap_));
  
  // setup global Mortar matrices Dold and Mold
  dold_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_));
  dold_->Zero();
  dold_->Complete();
  mold_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_));
  mold_->Zero();
  mold_->Complete(*gmdofrowmap_,*gsdofrowmap_);
  
  // friction
  // setup vector of displacement jumps (slave dof) 
  jump_       = rcp(new Epetra_Vector(*gsdofrowmap_));
  
  return;
}


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CONTACT::Manager& manager)
{
  manager.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print manager (public)                                   mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::Print(ostream& os) const
{
  if (Comm().MyPID()==0)
  {
    os << "---------------------------------------------CONTACT::Manager\n"
       << "Contact interfaces: " << (int)interface_.size() << endl
       << "-------------------------------------------------------------\n";
  }
  Comm().Barrier();
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    cout << *(interface_[i]);
  }
  Comm().Barrier();

  return;
}

/*----------------------------------------------------------------------*
 |  read and check input parameters (public)                  popp 04/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Manager::ReadAndCheckInput()
{
  // read parameter list from DRT::Problem
  const Teuchos::ParameterList& input = DRT::Problem::Instance()->StructuralContactParams();
  
  // read contact type
  switch (Teuchos::getIntegralValue<int>(input,"CONTACT"))
  {
    case INPUTPARAMS::contact_none:
      scontact_.set<string>("contact type","none");
      break;
    case INPUTPARAMS::contact_normal:
      scontact_.set<string>("contact type","normal");
      break;
    case INPUTPARAMS::contact_frictional:
      scontact_.set<string>("contact type","frictional");
      break;
    case INPUTPARAMS::contact_meshtying:
      scontact_.set<string>("contact type","meshtying");
      break;
    default:
      dserror("Cannot cope with choice of contact type");
      break;
  }
  
  // read basis trafo flag
  scontact_.set<bool> ("basis transformation",Teuchos::getIntegralValue<int>(input,"BASISTRAFO"));
  
  // read friction type
  switch (Teuchos::getIntegralValue<int>(input,"FRICTION"))
  {
    case INPUTPARAMS::friction_none:
      scontact_.set<string>("friction type","none");
      break;
    case INPUTPARAMS::friction_stick:
      scontact_.set<string>("friction type","stick");
      break;
    case INPUTPARAMS::friction_tresca:
      scontact_.set<string>("friction type","tresca");
      break;
    case INPUTPARAMS::friction_coulomb:
      scontact_.set<string>("friction type","coulomb");
      break;
    default:
      dserror("Cannot cope with choice of friction type");
      break;
  }
  
  // read friction parameters
  scontact_.set<double> ("friction bound",input.get<double>("FRBOUND"));
  scontact_.set<double> ("friction coefficient",input.get<double>("FRCOEFF"));
  
  // read full linearization flag
  scontact_.set<bool> ("full linearization",Teuchos::getIntegralValue<int>(input,"FULL_LINEARIZATION"));
  
  // read semi-smooth Newton flag and parameters
  scontact_.set<bool> ("semismooth newton",Teuchos::getIntegralValue<int>(input,"SEMI_SMOOTH_NEWTON"));
  scontact_.set<double>("semismooth cn",input.get<double>("SEMI_SMOOTH_CN"));
  scontact_.set<double>("semismooth ct",input.get<double>("SEMI_SMOOTH_CT"));
    
  // check contact input parameters
  string ctype   = scontact_.get<string>("contact type","none");
  bool btrafo    = scontact_.get<bool>("basis transformation",false);
  string ftype   = scontact_.get<string>("friction type","none");
  //double frbound = scontact_.get<double>("friction bound",0.0);
  //double frcoeff = scontact_.get<double>("friction coeffiecient",0.0);
  bool fulllin   = scontact_.get<bool>("full linearization",false);
  double ct      = scontact_.get<double>("semismooth ct",0.0);
  
  // invalid parameter combinations
  if (btrafo)
    dserror("Basis transformed versions are not up to date");
  if (ctype=="normal" && ftype !="none")
    dserror("Friction law supplied for normal contact");
  if (ctype=="frictional" && ftype=="none")
    dserror("No friction law supplied for frictional contact");
  if (ctype=="frictional" && ftype=="coulomb")
    dserror("Coulomb friction law not yet implemented");
  if (ctype=="frictional" && fulllin)
    dserror("Full linearization not yet implemented for friction");
  if (btrafo && fulllin)
    dserror("Full linearization not yet implemented for basis trafo case");
  if (ctype=="frictional" && ct == 0)
  	dserror("Friction Parameter ct = 0, must be greater than 0");
  
  // overrule input in certain cases
  if (ctype=="meshtying" && ftype!="stick")
    scontact_.set<string>("friction type","stick");
    
  return true;
}
  
/*----------------------------------------------------------------------*
 |  write restart information for contact (public)            popp 03/08|
 *----------------------------------------------------------------------*/
RCP<Epetra_Vector> CONTACT::Manager::WriteRestart()
{
  RCP<Epetra_Vector> activetoggle = rcp(new Epetra_Vector(*gsnoderowmap_));
  
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* cnode = static_cast<CNode*>(node);
      
      // set value active / inactive in toggle vector
      if (cnode->Active()) (*activetoggle)[j]=1;
    }
  }
  
  return activetoggle;
}

/*----------------------------------------------------------------------*
 |  read restart information for contact (public)             popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::ReadRestart(const RCP<Epetra_Vector> activetoggle)
{
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      if ((*activetoggle)[j]==1)
      {
        int gid = interface_[i]->SlaveRowNodes()->GID(j);
        DRT::Node* node = interface_[i]->Discret().gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %",gid);
        CNode* cnode = static_cast<CNode*>(node);
      
        // set value active / inactive in cnode
        cnode->Active()=true;
      }
    }
  }
  
  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i=0;i<(int)interface_.size();++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_,interface_[i]->ActiveNodes(),false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_,interface_[i]->ActiveDofs(),false);
    gactiven_ = LINALG::MergeMap(gactiven_,interface_[i]->ActiveNDofs(),false);
    gactivet_ = LINALG::MergeMap(gactivet_,interface_[i]->ActiveTDofs(),false);
    gslipnodes_ = LINALG::MergeMap(gslipnodes_,interface_[i]->SlipNodes(),false);
    gslipdofs_ = LINALG::MergeMap(gslipdofs_,interface_[i]->SlipDofs(),false);
    gslipt_ = LINALG::MergeMap(gslipt_,interface_[i]->SlipTDofs(),false);
  }
  
  // update flag for global contact status
  if (gactivenodes_->NumGlobalElements())
    IsInContact()=true;
  
  return;
}


#endif  // #ifdef CCADISCRET
