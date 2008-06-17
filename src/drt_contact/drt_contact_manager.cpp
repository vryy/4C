/*!----------------------------------------------------------------------
\file drt_contact_manager.cpp
\brief Main class to control all contact

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"

#include "drt_contact_manager.H"
#include "drt_cnode.H"
#include "drt_celement.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_solver.H"
#include "contactdefines.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 03/08|
 *----------------------------------------------------------------------*/
CONTACT::Manager::Manager(DRT::Discretization& discret, double alphaf) :
discret_(discret),
alphaf_(alphaf),
activesetconv_(false),
activesetsteps_(0),
isincontact_(false)
{
#ifdef CONTACTMESHTYING
#ifndef CONTACTSTICKING
  dserror("ERROR: Contact_Manager: Mesh Tying option only in combination with sticking!");
#endif
#endif
  
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
    interface_.push_back(rcp(new CONTACT::Interface(groupid1,Comm())));
    
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
    
    // note that the nodal ids are unique because they come from
    // one global problem discretization conatining all nodes of the
    // contact interface
    // We rely on this fact, therefore it is not possible to
    // do contact between two distinct discretizations here
    
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
        RCP<CONTACT::CNode> cnode = rcp(new CONTACT::CNode(node->Id(),node->X(),
                                                           node->Owner(),
                                                           Discret().NumDof(node),
                                                           Discret().Dof(node),isslave[j]));
        
        // note that we do not have to worry about double entries
        // as the AddNode function can deal with this case!
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
  bool initialcontact = scontact_.get("initial contact",false);
  for (int i=0;i<(int)interface_.size();++i)
  {
    interface_[i]->InitializeActiveSet(initialcontact);
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
  
  // read initial contact and basis trafo flags
  scontact_.set<bool> ("initial contact",Teuchos::getIntegralValue<int>(input,"INIT_CONTACT"));
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
  bool init      = scontact_.get<bool>("initial contact",false);
  bool btrafo    = scontact_.get<bool>("basis transformation",false);
  string ftype   = scontact_.get<string>("friction type","none");
  //double frbound = scontact_.get<double>("friction bound",0.0);
  //double frcoeff = scontact_.get<double>("friction coeffiecient",0.0);
  bool fulllin   = scontact_.get<bool>("full linearization",false);
  bool semismooth= scontact_.get<bool>("semismooth newton",false);
  
  // invalid parameter combinations
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
  if (ctype=="frictional" && semismooth)
    dserror("Semi-smooth Newton approach not yet implemented for friction");
  
  // overrule input in certain cases
  if (ctype=="meshtying" && !init)
    scontact_.set<bool>("initial contact",true);
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
  }
  
  // update flag for global contact status
  if (gactivenodes_->NumGlobalElements())
    IsInContact()=true;
  
  return;
}

/*----------------------------------------------------------------------*
 |  initialize contact for next Newton step (public)          popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::Initialize(int numiter)
{
#ifdef CONTACTCHECKHUEEBER
  if (numiter==0)
  {
#endif // #ifdef CONTACTCHECKHUEEBER
    
  // initialize / reset interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->Initialize();
  }

  // intitialize Dold and Mold if not done already
   if (dold_==null)
   {
     dold_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
     dold_->Zero();
     dold_->Complete();
   }
   if (mold_==null)
   {
     mold_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
     mold_->Zero();
     mold_->Complete(*gmdofrowmap_,*gsdofrowmap_);
   }
   
  // (re)setup global Mortar LINALG::SparseMatrices and Epetra_Vectors
  dmatrix_    = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
  mmatrix_    = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
  mhatmatrix_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
  g_          = LINALG::CreateVector(*gsnoderowmap_,true);
    
  // (re)setup global normal and tangent matrices
  nmatrix_ = rcp(new LINALG::SparseMatrix(*gactiven_,3));
  tmatrix_ = rcp(new LINALG::SparseMatrix(*gactivet_,3));
  
  // (re)setup global Tresca friction matrix L and vector R
  string ftype   = scontact_.get<string>("friction type","none");
  if (ftype=="tresca")
  {
	// FIXME: at the moment we do NOT seperate stick and slip mode
	// for Tresca friction. This is pure slip! Later we have to split
	// the active set into stick and split set!
	lmatrix_ = rcp(new LINALG::SparseMatrix(*gactivet_,10));
	r_       = LINALG::CreateVector(*gactivet_,true);
  }
  
  // (re)setup global matrices containing derivatives
  smatrix_ = rcp(new LINALG::SparseMatrix(*gactiven_,3));
  pmatrix_ = rcp(new LINALG::SparseMatrix(*gactivet_,3));
  
  // (re)setup global matrices containing fc derivatives
  lindmatrix_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
  linmmatrix_ = rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
  
  
#ifdef CONTACTCHECKHUEEBER
  }
#endif // #ifdef CONTACTCHECKHUEEBER
  return;
}

/*----------------------------------------------------------------------*
 |  set current deformation state (public)                    popp 11/07|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::SetState(const string& statename,
                                const RCP<Epetra_Vector> vec)
{
  // set state on interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->SetState(statename,vec);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate contact Mortar matrices D,M only (public)        popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::EvaluateMortar()
{ 
  // call interfaces to evaluate Mortar coupling
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->Evaluate();
    interface_[i]->AssembleDMG(*dmatrix_,*mmatrix_,*g_);
  }
    
  // FillComplete() global Mortar matrices
  dmatrix_->Complete();
  mmatrix_->Complete(*gmdofrowmap_,*gsdofrowmap_);
    
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate contact (public)                                 popp 04/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::Evaluate(RCP<LINALG::SparseMatrix> kteff,
                                RCP<Epetra_Vector> feff, int numiter)
{ 
  // check if Tresca friction and/or basis transformation should be applied
  string ftype   = scontact_.get<string>("friction type","none");
  bool btrafo = scontact_.get<bool>("basis transformation",false);
  
  // Tresca friction cases
  if (ftype=="tresca")
  {
	  if (btrafo)
		  EvaluateTrescaBasisTrafo(kteff,feff,numiter);
	  else
		  dserror("ERROR: Tresca friction not yet implemented without basis trafo!");
  }
  
  // Other cases (Frictionless, Stick, MeshTying)
  else
  {
  if (btrafo)
    EvaluateBasisTrafo(kteff,feff,numiter);
  else
    EvaluateNoBasisTrafo(kteff,feff,numiter);
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate Tresca friction with basis trafo (public)        mgit 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::EvaluateTrescaBasisTrafo(RCP<LINALG::SparseMatrix> kteff,
                                                RCP<Epetra_Vector> feff, int numiter)
{ 
  // FIXME: Currently only the old LINALG::Multiply method is used,
  // because there are still problems with the transposed version of
  // MLMultiply if a row has no entries! One day we should use ML...
  
  // input parameters
  string ctype   = scontact_.get<string>("contact type","none");
  string ftype   = scontact_.get<string>("friction type","none");
  
  /**********************************************************************/
  /* evaluate interfaces                                                */
  /* (nodal normals, projections, Mortar integration, Mortar assembly)  */
  /**********************************************************************/
#ifdef CONTACTCHECKHUEEBER
  if (numiter==0)
  {
#endif // #ifdef CONTACTCHECKHUEEBER
  for (int i=0; i<(int)interface_.size(); ++i)
  { 
	  interface_[i]->Evaluate();
    interface_[i]->AssembleDMG(*dmatrix_,*mmatrix_,*g_);
  }

  // FillComplete() global Mortar matrices
  dmatrix_->Complete();
  mmatrix_->Complete(*gmdofrowmap_,*gsdofrowmap_);
#ifdef CONTACTCHECKHUEEBER
  }
#endif // #ifdef CONTACTCHECKHUEEBER
  // export weighted gap vector to gactiveN-map
  RCP<Epetra_Vector> gact = LINALG::CreateVector(*gactivenodes_,true);
  if (gact->GlobalLength())
  {
    LINALG::Export(*g_,*gact);
    gact->ReplaceMap(*gactiven_);
  }
  
  /**********************************************************************/
  /* build global matrix n with normal vectors of active nodes,         */
  /* global matrix t with tangent vectors of active nodes               */
  /* and global matrix l and vector r for frictional contact            */
  /**********************************************************************/
  
  // read tresca friction bound
  double frbound = scontact_.get<double>("friction bound",0.0);
  
  // read weighting factor ct
  // (this is necessary in semi-smooth Newton case, as the search for the
  // active set is now part of the Newton iteration. Thus, we do not know
  // the active / inactive status in advance and we can have a state in
  // which both firctional conditions are violated. Here we have to weigh
  // the two violations via ct!
  double ct = scontact_.get<double>("semismooth ct",0.0);

#ifdef CONTACTCHECKHUEEBER
  if (numiter==0)
  {
#endif // #ifdef CONTACTCHECKHUEEBER
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->AssembleNT(*nmatrix_,*tmatrix_);
    interface_[i]->AssembleTresca(*lmatrix_,*r_,frbound,ct); 
  }
  
  // FillComplete() global matrices N, T and L
  nmatrix_->Complete(*gactivedofs_,*gactiven_);
  tmatrix_->Complete(*gactivedofs_,*gactivet_);
  lmatrix_->Complete(*gactivet_,*gactivet_);
#ifdef CONTACTCHECKHUEEBER
  }
#endif // #ifdef CONTACTCHECKHUEEBER

  /**********************************************************************/
  /* Multiply Mortar matrices: m^ = inv(d) * m                          */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> invd = rcp(new LINALG::SparseMatrix(*dmatrix_));
  RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_,true);
  int err = 0;
  
  // extract diagonal of invd into diag
  invd->ExtractDiagonalCopy(*diag);
  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
    if ((*diag)[i]==0.0) (*diag)[i]=1.0;
  
  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");
  
  // re-insert inverted diagonal into invd
  err = invd->ReplaceDiagonalValues(*diag);
  // we cannot use this check, as we deliberately replaced zero entries
  //if (err>0) dserror("ERROR: ReplaceDiagonalValues: Missing diagonal entry!");
  
  // do the multiplication M^ = inv(D) * M
  mhatmatrix_ = LINALG::Multiply(*invd,false,*mmatrix_,false);
  
  /**********************************************************************/
  /* Split kteff into 3x3 block matrix                                  */
  /**********************************************************************/
  // we want to split k into 3 groups s,m,n = 9 blocks
  RCP<LINALG::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;
  
  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!) 
  RCP<LINALG::SparseMatrix> ksmsm, ksmn, knsm;
  
  // we also need the combined sm rowmap
  // (this map is NOT allowed to have an overlap !!!)
  RCP<Epetra_Map> gsmdofs = LINALG::MergeMap(gsdofrowmap_,gmdofrowmap_,false);
  
  // some temporary RCPs
  RCP<Epetra_Map> tempmap;
  RCP<LINALG::SparseMatrix> tempmtx1;
  RCP<LINALG::SparseMatrix> tempmtx2;
  RCP<LINALG::SparseMatrix> tempmtx3;
  
  // split into slave/master part + structure part
  LINALG::SplitMatrix2x2(kteff,gsmdofs,gndofrowmap_,gsmdofs,gndofrowmap_,ksmsm,ksmn,knsm,knn);
  
  // further splits into slave part + master part
  LINALG::SplitMatrix2x2(ksmsm,gsdofrowmap_,gmdofrowmap_,gsdofrowmap_,gmdofrowmap_,kss,ksm,kms,kmm);
  LINALG::SplitMatrix2x2(ksmn,gsdofrowmap_,gmdofrowmap_,gndofrowmap_,tempmap,ksn,tempmtx1,kmn,tempmtx2);
  LINALG::SplitMatrix2x2(knsm,gndofrowmap_,tempmap,gsdofrowmap_,gmdofrowmap_,kns,knm,tempmtx1,tempmtx2);
  
  // output for checking everything
#ifdef CONTACTDIMOUTPUT
  if(Comm().MyPID()==0)
  {
    cout << endl << "*****************" << endl;
    cout << "k:   " << kteff->EpetraMatrix()->NumGlobalRows() << " x " << kteff->EpetraMatrix()->NumGlobalCols() << endl;
    if (kss!=null) cout << "kss: " << kss->EpetraMatrix()->NumGlobalRows() << " x " << kss->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kss: null" << endl;
    if (ksm!=null) cout << "ksm: " << ksm->EpetraMatrix()->NumGlobalRows() << " x " << ksm->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "ksm: null" << endl;
    if (ksn!=null) cout << "ksn: " << ksn->EpetraMatrix()->NumGlobalRows() << " x " << ksn->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "ksn: null" << endl;
    if (kms!=null) cout << "kms: " << kms->EpetraMatrix()->NumGlobalRows() << " x " << kms->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kms: null" << endl;
    if (kmm!=null) cout << "kmm: " << kmm->EpetraMatrix()->NumGlobalRows() << " x " << kmm->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kmm: null" << endl;
    if (kmn!=null) cout << "kmn: " << kmn->EpetraMatrix()->NumGlobalRows() << " x " << kmn->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kmn: null" << endl;
    if (kns!=null) cout << "kns: " << kns->EpetraMatrix()->NumGlobalRows() << " x " << kns->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kns: null" << endl;
    if (knm!=null) cout << "knm: " << knm->EpetraMatrix()->NumGlobalRows() << " x " << knm->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "knm: null" << endl;
    if (knn!=null) cout << "knn: " << knn->EpetraMatrix()->NumGlobalRows() << " x " << knn->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "knn: null" << endl;
    cout << "*****************" << endl;
  }
#endif // #ifdef CONTACTDIMOUTPUT
  
  /**********************************************************************/
  /* Split feff into 3 subvectors                                       */
  /**********************************************************************/
  // we want to split f into 3 groups s.m,n
  RCP<Epetra_Vector> fs, fm, fn;
  
  // temporarily we need the group sm
  RCP<Epetra_Vector> fsm;
  
  // do the vector splitting smn -> sm+n -> s+m+n
  LINALG::SplitVector(*feff,*gsmdofs,fsm,*gndofrowmap_,fn);
  LINALG::SplitVector(*fsm,*gsdofrowmap_,fs,*gmdofrowmap_,fm);
  
  // output for checking everything
#ifdef CONTACTDIMOUTPUT
  if(Comm().MyPID()==0)
  {
    cout << endl << "**********" << endl;
    cout << "f:  " << feff->GlobalLength() << endl;
    cout << "fs: " << fs->GlobalLength() << endl;
    cout << "fm: " << fm->GlobalLength() << endl;
    cout << "fn: " << fn->GlobalLength() << endl;
    cout << "**********" << endl;
  }
#endif // #ifdef CONTACTDIMOUTPUT
  
  // store some stuff for static condensation of LM
  fs_   = fs;
  invd_ = invd;
  ksn_  = ksn;
  ksm_  = ksm;
  kss_  = kss;
  
  /**********************************************************************/
  /* Apply basis transformation to k                                    */
  /**********************************************************************/
  // define temporary RCP mod
  RCP<LINALG::SparseMatrix> mod;
  
  // kss: nothing to do
  RCP<LINALG::SparseMatrix> kssmod = kss;
  
  // ksm: add kss*T(mbar)
  RCP<LINALG::SparseMatrix> ksmmod = LINALG::Multiply(*kss,false,*mhatmatrix_,false,false);
  ksmmod->Add(*ksm,false,1.0,1.0);
  ksmmod->Complete(ksm->DomainMap(),ksm->RowMap());
  
  // ksn: nothing to do
  RCP<LINALG::SparseMatrix> ksnmod = ksn;
  
  // kms: add T(mbar)*kss
  RCP<LINALG::SparseMatrix> kmsmod = LINALG::Multiply(*mhatmatrix_,true,*kss,false,false);
  kmsmod->Add(*kms,false,1.0,1.0);
  kmsmod->Complete(kms->DomainMap(),kms->RowMap());
  
  // kmm: add kms*T(mbar) + T(mbar)*ksm + T(mbar)*kss*mbar
  RCP<LINALG::SparseMatrix> kmmmod = LINALG::Multiply(*kms,false,*mhatmatrix_,false,false);
  mod = LINALG::Multiply(*mhatmatrix_,true,*ksm,false);
  kmmmod->Add(*mod,false,1.0,1.0);
  mod = LINALG::Multiply(*mhatmatrix_,true,*kss,false);
  mod = LINALG::Multiply(*mod,false,*mhatmatrix_,false);
  kmmmod->Add(*mod,false,1.0,1.0);
  kmmmod->Add(*kmm,false,1.0,1.0);
  kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());
  
  // kmn: add T(mbar)*ksn
  RCP<LINALG::SparseMatrix> kmnmod = LINALG::Multiply(*mhatmatrix_,true,*ksn,false,false);
  kmnmod->Add(*kmn,false,1.0,1.0);
  kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());
  
  // kns: nothing to do
  RCP<LINALG::SparseMatrix> knsmod = kns;
  
  // knm: add kns*mbar
  RCP<LINALG::SparseMatrix> knmmod = LINALG::Multiply(*kns,false,*mhatmatrix_,false,false);
  knmmod->Add(*knm,false,1.0,1.0);
  knmmod->Complete(knm->DomainMap(),knm->RowMap());
  
  // knn: nothing to do
  RCP<LINALG::SparseMatrix> knnmod = knn;
  
  // output for checking everything
#ifdef CONTACTDIMOUTPUT
  if(Comm().MyPID()==0)
  {
    cout << endl << "*****************" << endl;
    cout << "kmod:   " << kteff->EpetraMatrix()->NumGlobalRows() << " x " << kteff->EpetraMatrix()->NumGlobalCols() << endl;
    if (kssmod!=null) cout << "kssmod: " << kssmod->EpetraMatrix()->NumGlobalRows() << " x " << kssmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kssmod: null" << endl;
    if (ksmmod!=null) cout << "ksmmod: " << ksmmod->EpetraMatrix()->NumGlobalRows() << " x " << ksmmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "ksmmod: null" << endl;
    if (ksnmod!=null) cout << "ksnmod: " << ksnmod->EpetraMatrix()->NumGlobalRows() << " x " << ksnmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "ksnmod: null" << endl;
    if (kmsmod!=null) cout << "kmsmod: " << kmsmod->EpetraMatrix()->NumGlobalRows() << " x " << kmsmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kmsmod: null" << endl;
    if (kmmmod!=null) cout << "kmmmod: " << kmmmod->EpetraMatrix()->NumGlobalRows() << " x " << kmmmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kmmmod: null" << endl;
    if (kmnmod!=null) cout << "kmnmod: " << kmnmod->EpetraMatrix()->NumGlobalRows() << " x " << kmnmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kmnmod: null" << endl;
    if (knsmod!=null) cout << "knsmod: " << knsmod->EpetraMatrix()->NumGlobalRows() << " x " << knsmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "knsmod: null" << endl;
    if (knmmod!=null) cout << "knmmod: " << knmmod->EpetraMatrix()->NumGlobalRows() << " x " << knmmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "knmmod: null" << endl;
    if (knnmod!=null) cout << "knnmod: " << knnmod->EpetraMatrix()->NumGlobalRows() << " x " << knnmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "knnmod: null" << endl;
    cout << "*****************" << endl;
  }
#endif // #ifdef CONTACTDIMOUTPUT
  
  /**********************************************************************/
  /* Apply basis transformation to f                                    */
  /**********************************************************************/
  // fs: nothing to be done
  RCP<Epetra_Vector> fsmod = fs;
  
  // fm: add T(mbar)*fs
  RCP<Epetra_Vector> fmmod = rcp(new Epetra_Vector(*gmdofrowmap_));
  mhatmatrix_->Multiply(true,*fs,*fmmod);
  fmmod->Update(1.0,*fm,1.0);
  
  // fn: nothing to be done
  RCP<Epetra_Vector> fnmod = fn;
  
  // output for checking everything
#ifdef CONTACTDIMOUTPUT
  if(Comm().MyPID()==0)
  {
    cout << endl << "**********" << endl;
    cout << "fmod:  " << feff->GlobalLength() << endl;
    cout << "fsmod: " << fsmod->GlobalLength() << endl;
    cout << "fmmod: " << fmmod->GlobalLength() << endl;
    cout << "fnmod: " << fnmod->GlobalLength() << endl;
    cout << "**********" << endl;
  }
#endif // #ifdef CONTACTDIMOUTPUT
  
  /**********************************************************************/
  /* Split slave quantities into active / inactive                      */
  /**********************************************************************/
  // we want to split kssmod into 2 groups a,i = 4 blocks
  RCP<LINALG::SparseMatrix> kaamod, kaimod, kiamod, kiimod;
  
  // we want to split ksnmod / ksmmod into 2 groups a,i = 2 blocks
  RCP<LINALG::SparseMatrix> kanmod, kinmod, kammod, kimmod;
    
  // we will get the i rowmap as a by-product
  RCP<Epetra_Map> gidofs;
    
  // do the splitting
  LINALG::SplitMatrix2x2(kssmod,gactivedofs_,gidofs,gactivedofs_,gidofs,kaamod,kaimod,kiamod,kiimod);
  LINALG::SplitMatrix2x2(ksnmod,gactivedofs_,gidofs,gndofrowmap_,tempmap,kanmod,tempmtx1,kinmod,tempmtx2);
  LINALG::SplitMatrix2x2(ksmmod,gactivedofs_,gidofs,gmdofrowmap_,tempmap,kammod,tempmtx1,kimmod,tempmtx2);
  
  // output for checking everything
#ifdef CONTACTDIMOUTPUT
  if(Comm().MyPID()==0)
  {
    cout << endl << "*****************" << endl;
    cout << "kssmod: " << kssmod->EpetraMatrix()->NumGlobalRows() << " x " << kssmod->EpetraMatrix()->NumGlobalCols() << endl;
    if (kaamod!=null) cout << "kaamod: " << kaamod->EpetraMatrix()->NumGlobalRows() << " x " << kaamod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kaamod: null" << endl;
    if (kaimod!=null) cout << "kaimod: " << kaimod->EpetraMatrix()->NumGlobalRows() << " x " << kaimod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kaimod: null" << endl;
    if (kiamod!=null) cout << "kiamod: " << kiamod->EpetraMatrix()->NumGlobalRows() << " x " << kiamod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kiamod: null" << endl;
    if (kiimod!=null) cout << "kiimod: " << kiimod->EpetraMatrix()->NumGlobalRows() << " x " << kiimod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kiimod: null" << endl;
    cout << "*****************" << endl;
    
    cout << endl << "*****************" << endl;
    cout << "ksnmod: " << ksnmod->EpetraMatrix()->NumGlobalRows() << " x " << ksnmod->EpetraMatrix()->NumGlobalCols() << endl;
    if (kanmod!=null) cout << "kanmod: " << kanmod->EpetraMatrix()->NumGlobalRows() << " x " << kanmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kanmod: null" << endl;
    if (kinmod!=null) cout << "kinmod: " << kinmod->EpetraMatrix()->NumGlobalRows() << " x " << kinmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kinmod: null" << endl;
    cout << "*****************" << endl;
    
    cout << endl << "*****************" << endl;
    cout << "ksmmod: " << ksmmod->EpetraMatrix()->NumGlobalRows() << " x " << ksmmod->EpetraMatrix()->NumGlobalCols() << endl;
    if (kammod!=null) cout << "kammod: " << kammod->EpetraMatrix()->NumGlobalRows() << " x " << kammod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kammod: null" << endl;
    if (kimmod!=null) cout << "kimmod: " << kimmod->EpetraMatrix()->NumGlobalRows() << " x " << kimmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kimmod: null" << endl;
    cout << "*****************" << endl;
  }
#endif // #ifdef CONTACTDIMOUTPUT
  
  // we want to split fsmod into 2 groups a,i
  RCP<Epetra_Vector> famod, fimod;
  
  // do the vector splitting s -> a+i
  if (!gidofs->NumGlobalElements())
    famod = rcp(new Epetra_Vector(*fsmod));
  else if (!gactivedofs_->NumGlobalElements())
    fimod = rcp(new Epetra_Vector(*fsmod));
  else
  {
    LINALG::SplitVector(*fsmod,*gactivedofs_,famod,*gidofs,fimod);
  }
  
  // output for checking everything
#ifdef CONTACTDIMOUTPUT
  if(Comm().MyPID()==0)
  {
    cout << endl << "**********" << endl;
    if (fsmod!=null) cout << "fsmod: " << fsmod->GlobalLength() << endl;
    else cout << "fsmod: null" << endl;
    if (famod!=null) cout << "famod: " << famod->GlobalLength() << endl;
    else cout << "famod: null" << endl;
    if (fimod!=null) cout << "fimod: " << fimod->GlobalLength() << endl;
    else cout << "fimod: null" << endl;
    cout << "**********" << endl;
  }
#endif // #ifdef CONTACTDIMOUTPUT
  
  /**********************************************************************/
  /* Isolate active part from invd and dold                             */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> invda;
  LINALG::SplitMatrix2x2(invd_,gactivedofs_,gidofs,gactivedofs_,gidofs,invda,tempmtx1,tempmtx2,tempmtx3);
  invda->Scale(1/(1-alphaf_));
  
  RCP<LINALG::SparseMatrix> dolda, doldi;
  LINALG::SplitMatrix2x2(dold_,gactivedofs_,gidofs,gactivedofs_,gidofs,dolda,tempmtx1,tempmtx2,doldi);
    
  /**********************************************************************/
  /* Gen-alpha modifications                                            */
  /**********************************************************************/
  // fi: subtract alphaf * old contact forces (t_n)
  if (gidofs->NumGlobalElements())
  {
    RCP<Epetra_Vector> modi = rcp(new Epetra_Vector(*gidofs));
    LINALG::Export(*zold_,*modi);
    RCP<Epetra_Vector> tempveci = rcp(new Epetra_Vector(*gidofs));
    doldi->Multiply(false,*modi,*tempveci);
    fimod->Update(-alphaf_,*tempveci,1.0);
  }
   
  // fa: subtract alphaf * old contact forces (t_n)
  if (gactivedofs_->NumGlobalElements())
  {
    RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(*gactivedofs_));
    LINALG::Export(*zold_,*mod);
    RCP<Epetra_Vector> tempvec = rcp(new Epetra_Vector(*gactivedofs_));
    dolda->Multiply(false,*mod,*tempvec);
    famod->Update(-alphaf_,*tempvec,1.0);
  }
     
  // do the multiplications with t matrix
  RCP<LINALG::SparseMatrix> tkanmod, tkammod, tkaimod, tkaamod, tlmatrix;
  RCP<Epetra_Vector> tfamod;
  
  if(gactivedofs_->NumGlobalElements())
  {
  	// kanmod: multiply with tmatrix
  	tkanmod = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
  	tkanmod = LINALG::Multiply(*tkanmod,false,*kanmod,false,true);
  	
  	// kammod: multiply with tmatrix
  	tkammod = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
  	tkammod = LINALG::Multiply(*tkammod,false,*kammod,false,true);
  	  	
  	// friction
  	// lmatrix: multiply with tmatrix 
  	tlmatrix = LINALG::Multiply(*lmatrix_,false,*tmatrix_,false,true);
  	 	
  	// kaamod: multiply with tmatrix
  	tkaamod = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
  	tkaamod = LINALG::Multiply(*tkaamod,false,*kaamod,false,false);
  	
  	// add tlmatirx to tkaamod
  	tkaamod->Add(*tlmatrix,false,1.0,1.0);
  	tkaamod->Complete(kaamod->DomainMap(),kaamod->RowMap());
  	
  	if (gidofs->NumGlobalElements())
    {
      //kaimod: multiply with tmatrix
    	tkaimod = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
    	tkaimod = LINALG::Multiply(*tkaimod,false,*kaimod,false,true);
    }	  		
    
    // famod: multiply with tmatrix    
    tfamod = rcp(new Epetra_Vector(*gactivet_));
    RCP<LINALG::SparseMatrix> temp = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
    temp->Multiply(false,*famod,*tfamod); 
    
    // friction
    // add r to famod
    tfamod->Update(1.0,*r_,1.0);
  }
  
  // output for checking everything
#ifdef CONTACTDIMOUTPUT
  if(Comm().MyPID()==0)
  {
    cout << endl << "*****************" << endl;
    if (tkanmod!=null) cout << "tkanmod: " << tkanmod->EpetraMatrix()->NumGlobalRows() << " x " << tkanmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "tkanmod: null" << endl;
    if (tkammod!=null) cout << "tkammod: " << tkammod->EpetraMatrix()->NumGlobalRows() << " x " << tkammod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "tkammod: null" << endl;
    if (tkaimod!=null) cout << "tkaimod: " << tkaimod->EpetraMatrix()->NumGlobalRows() << " x " << tkaimod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "tkaimod: null" << endl;
    if (tkaamod!=null) cout << "tkaamod: " << tkaamod->EpetraMatrix()->NumGlobalRows() << " x " << tkaamod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "tkaamod: null" << endl;
    cout << "*****************" << endl;
    
    cout << endl << "**********" << endl;
    if (tfamod!=null) cout << "tfamod: " << tfamod->GlobalLength() << endl;
    else cout << "tfamod: null" << endl;
    cout << "**********" << endl;
  }
#endif // #ifdef CONTACTDIMOUTPUT
  
  /**********************************************************************/
  /* Global setup of kteffnew, feffnew (including contact)              */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> kteffnew = rcp(new LINALG::SparseMatrix(*(discret_.DofRowMap()),81));
  RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*(discret_.DofRowMap()));
  
  // add n / m submatrices to kteffnew
  kteffnew->Add(*knnmod,false,1.0,1.0);
  kteffnew->Add(*knmmod,false,1.0,1.0);
  kteffnew->Add(*kmnmod,false,1.0,1.0);
  kteffnew->Add(*kmmmod,false,1.0,1.0);
  
  // add a / i submatrices to kteffnew, if existing
  if (knsmod!=null) kteffnew->Add(*knsmod,false,1.0,1.0);
  if (kmsmod!=null) kteffnew->Add(*kmsmod,false,1.0,1.0);
  if (kinmod!=null) kteffnew->Add(*kinmod,false,1.0,1.0);
  if (kimmod!=null) kteffnew->Add(*kimmod,false,1.0,1.0);
  if (kiimod!=null) kteffnew->Add(*kiimod,false,1.0,1.0);
  if (kiamod!=null) kteffnew->Add(*kiamod,false,1.0,1.0);
  
  // add matrix of normals to kteffnew
  kteffnew->Add(*nmatrix_,false,1.0,1.0);

 
  // add submatrices with tangents to kteffnew, if existing
  if (tkanmod!=null) kteffnew->Add(*tkanmod,false,1.0,1.0);
  if (tkammod!=null) kteffnew->Add(*tkammod,false,1.0,1.0);
  if (tkaimod!=null) kteffnew->Add(*tkaimod,false,1.0,1.0);
  if (tkaamod!=null) kteffnew->Add(*tkaamod,false,1.0,1.0);
  
  // FillComplete kteffnew (square)
  kteffnew->Complete();
  
  // add n / m subvectors to feffnew
  RCP<Epetra_Vector> fnmodexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  RCP<Epetra_Vector> fmmodexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  LINALG::Export(*fnmod,*fnmodexp);
  LINALG::Export(*fmmod,*fmmodexp);
  feffnew->Update(1.0,*fnmodexp,1.0,*fmmodexp,1.0);
  
  // add i / ta subvectors to feffnew, if existing
  RCP<Epetra_Vector> fimodexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  RCP<Epetra_Vector> tfamodexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  if (fimod!=null) LINALG::Export(*fimod,*fimodexp);
  if (tfamod!=null) LINALG::Export(*tfamod,*tfamodexp);
  feffnew->Update(1.0,*fimodexp,1.0,*tfamodexp,1.0);
  
  // add weighted gap vector to feffnew, if existing
  RCP<Epetra_Vector> gexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  if (gact->GlobalLength()) LINALG::Export(*gact,*gexp);
  feffnew->Update(1.0,*gexp,1.0);
  
  /**********************************************************************/
  /* Replace kteff and feff by kteffnew and feffnew                   */
  /**********************************************************************/
  *kteff = *kteffnew;
  *feff = *feffnew;
  //LINALG::PrintSparsityToPostscript(*(kteff->EpetraMatrix()));
  
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate contact with basis transformation (public)       popp 11/07|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::EvaluateBasisTrafo(RCP<LINALG::SparseMatrix> kteff,
                                          RCP<Epetra_Vector> feff, int numiter)
{ 
  // FIXME: Currently only the old LINALG::Multiply method is used,
  // because there are still problems with the transposed version of
  // MLMultiply if a row has no entries! One day we should use ML...
  
  // input parameters
  string ctype   = scontact_.get<string>("contact type","none");
  string ftype   = scontact_.get<string>("friction type","none");
  
  /**********************************************************************/
  /* evaluate interfaces                                                */
  /* (nodal normals, projections, Mortar integration, Mortar assembly)  */
  /**********************************************************************/
#ifdef CONTACTCHECKHUEEBER
  if (numiter==0)
  {
#endif // #ifdef CONTACTCHECKHUEEBER
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->Evaluate();
    interface_[i]->AssembleDMG(*dmatrix_,*mmatrix_,*g_);
  }
  
  // FillComplete() global Mortar matrices
  dmatrix_->Complete();
  mmatrix_->Complete(*gmdofrowmap_,*gsdofrowmap_);
#ifdef CONTACTCHECKHUEEBER
  }
#endif // #ifdef CONTACTCHECKHUEEBER
  // export weighted gap vector to gactiveN-map
  RCP<Epetra_Vector> gact = LINALG::CreateVector(*gactivenodes_,true);
  if (gact->GlobalLength())
  {
    LINALG::Export(*g_,*gact);
    gact->ReplaceMap(*gactiven_);
  }
  
  /**********************************************************************/
  /* build global matrix n with normal vectors of active nodes          */
  /* and global matrix t with tangent vectors of active nodes           */
  /**********************************************************************/
#ifdef CONTACTCHECKHUEEBER
  if (numiter==0)
  {
#endif // #ifdef CONTACTCHECKHUEEBER
  for (int i=0; i<(int)interface_.size(); ++i)
    interface_[i]->AssembleNT(*nmatrix_,*tmatrix_);
    
  // FillComplete() global matrices N and T
  nmatrix_->Complete(*gactivedofs_,*gactiven_);
  tmatrix_->Complete(*gactivedofs_,*gactivet_);
#ifdef CONTACTCHECKHUEEBER
  }
#endif // #ifdef CONTACTCHECKHUEEBER
  /**********************************************************************/
  /* Multiply Mortar matrices: m^ = inv(d) * m                          */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> invd = rcp(new LINALG::SparseMatrix(*dmatrix_));
  RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_,true);
  int err = 0;
  
  // extract diagonal of invd into diag
  invd->ExtractDiagonalCopy(*diag);
  
  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
    if ((*diag)[i]==0.0) (*diag)[i]=1.0;
  
  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");
  
  // re-insert inverted diagonal into invd
  err = invd->ReplaceDiagonalValues(*diag);
  // we cannot use this check, as we deliberately replaced zero entries
  //if (err>0) dserror("ERROR: ReplaceDiagonalValues: Missing diagonal entry!");
  
  // do the multiplication M^ = inv(D) * M
  mhatmatrix_ = LINALG::Multiply(*invd,false,*mmatrix_,false);
  
  /**********************************************************************/
  /* Split kteff into 3x3 block matrix                                  */
  /**********************************************************************/
  // we want to split k into 3 groups s,m,n = 9 blocks
  RCP<LINALG::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;
  
  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!) 
  RCP<LINALG::SparseMatrix> ksmsm, ksmn, knsm;
  
  // we also need the combined sm rowmap
  // (this map is NOT allowed to have an overlap !!!)
  RCP<Epetra_Map> gsmdofs = LINALG::MergeMap(gsdofrowmap_,gmdofrowmap_,false);
  
  // some temporary RCPs
  RCP<Epetra_Map> tempmap;
  RCP<LINALG::SparseMatrix> tempmtx1;
  RCP<LINALG::SparseMatrix> tempmtx2;
  
  // split into slave/master part + structure part
  LINALG::SplitMatrix2x2(kteff,gsmdofs,gndofrowmap_,gsmdofs,gndofrowmap_,ksmsm,ksmn,knsm,knn);
  
  // further splits into slave part + master part
  LINALG::SplitMatrix2x2(ksmsm,gsdofrowmap_,gmdofrowmap_,gsdofrowmap_,gmdofrowmap_,kss,ksm,kms,kmm);
  LINALG::SplitMatrix2x2(ksmn,gsdofrowmap_,gmdofrowmap_,gndofrowmap_,tempmap,ksn,tempmtx1,kmn,tempmtx2);
  LINALG::SplitMatrix2x2(knsm,gndofrowmap_,tempmap,gsdofrowmap_,gmdofrowmap_,kns,knm,tempmtx1,tempmtx2);
  
  // output for checking everything
#ifdef CONTACTDIMOUTPUT
  if(Comm().MyPID()==0)
  {
    cout << endl << "*****************" << endl;
    cout << "k:   " << kteff->EpetraMatrix()->NumGlobalRows() << " x " << kteff->EpetraMatrix()->NumGlobalCols() << endl;
    if (kss!=null) cout << "kss: " << kss->EpetraMatrix()->NumGlobalRows() << " x " << kss->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kss: null" << endl;
    if (ksm!=null) cout << "ksm: " << ksm->EpetraMatrix()->NumGlobalRows() << " x " << ksm->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "ksm: null" << endl;
    if (ksn!=null) cout << "ksn: " << ksn->EpetraMatrix()->NumGlobalRows() << " x " << ksn->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "ksn: null" << endl;
    if (kms!=null) cout << "kms: " << kms->EpetraMatrix()->NumGlobalRows() << " x " << kms->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kms: null" << endl;
    if (kmm!=null) cout << "kmm: " << kmm->EpetraMatrix()->NumGlobalRows() << " x " << kmm->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kmm: null" << endl;
    if (kmn!=null) cout << "kmn: " << kmn->EpetraMatrix()->NumGlobalRows() << " x " << kmn->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kmn: null" << endl;
    if (kns!=null) cout << "kns: " << kns->EpetraMatrix()->NumGlobalRows() << " x " << kns->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kns: null" << endl;
    if (knm!=null) cout << "knm: " << knm->EpetraMatrix()->NumGlobalRows() << " x " << knm->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "knm: null" << endl;
    if (knn!=null) cout << "knn: " << knn->EpetraMatrix()->NumGlobalRows() << " x " << knn->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "knn: null" << endl;
    cout << "*****************" << endl;
  }
#endif // #ifdef CONTACTDIMOUTPUT
  
  /**********************************************************************/
  /* Split feff into 3 subvectors                                       */
  /**********************************************************************/
  // we want to split f into 3 groups s.m,n
  RCP<Epetra_Vector> fs, fm, fn;
  
  // temporarily we need the group sm
  RCP<Epetra_Vector> fsm;
  
  // do the vector splitting smn -> sm+n -> s+m+n
  LINALG::SplitVector(*feff,*gsmdofs,fsm,*gndofrowmap_,fn);
  LINALG::SplitVector(*fsm,*gsdofrowmap_,fs,*gmdofrowmap_,fm);
  
  // output for checking everything
#ifdef CONTACTDIMOUTPUT
  if(Comm().MyPID()==0)
  {
    cout << endl << "**********" << endl;
    cout << "f:  " << feff->GlobalLength() << endl;
    cout << "fs: " << fs->GlobalLength() << endl;
    cout << "fm: " << fm->GlobalLength() << endl;
    cout << "fn: " << fn->GlobalLength() << endl;
    cout << "**********" << endl;
  }
#endif // #ifdef CONTACTDIMOUTPUT
  
  // store some stuff for static condensation of LM
  fs_   = fs;
  invd_ = invd;
  ksn_  = ksn;
  ksm_  = ksm;
  kss_  = kss;
    
  /**********************************************************************/
  /* Apply basis transformation to k                                    */
  /**********************************************************************/
  // define temporary RCP mod
  RCP<LINALG::SparseMatrix> mod;
  
  // kss: nothing to do
  RCP<LINALG::SparseMatrix> kssmod = kss;
  
  // ksm: add kss*T(mbar)
  RCP<LINALG::SparseMatrix> ksmmod = LINALG::Multiply(*kss,false,*mhatmatrix_,false,false);
  ksmmod->Add(*ksm,false,1.0,1.0);
  ksmmod->Complete(ksm->DomainMap(),ksm->RowMap());
  
  // ksn: nothing to do
  RCP<LINALG::SparseMatrix> ksnmod = ksn;
  
  // kms: add T(mbar)*kss
  RCP<LINALG::SparseMatrix> kmsmod = LINALG::Multiply(*mhatmatrix_,true,*kss,false,false);
  kmsmod->Add(*kms,false,1.0,1.0);
  kmsmod->Complete(kms->DomainMap(),kms->RowMap());
  
  // kmm: add kms*T(mbar) + T(mbar)*ksm + T(mbar)*kss*mbar
  RCP<LINALG::SparseMatrix> kmmmod = LINALG::Multiply(*kms,false,*mhatmatrix_,false,false);
  mod = LINALG::Multiply(*mhatmatrix_,true,*ksm,false);
  kmmmod->Add(*mod,false,1.0,1.0);
  mod = LINALG::Multiply(*mhatmatrix_,true,*kss,false);
  mod = LINALG::Multiply(*mod,false,*mhatmatrix_,false);
  kmmmod->Add(*mod,false,1.0,1.0);
  kmmmod->Add(*kmm,false,1.0,1.0);
  kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());
  
  // kmn: add T(mbar)*ksn
  RCP<LINALG::SparseMatrix> kmnmod = LINALG::Multiply(*mhatmatrix_,true,*ksn,false,false);
  kmnmod->Add(*kmn,false,1.0,1.0);
  kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());
  
  // kns: nothing to do
  RCP<LINALG::SparseMatrix> knsmod = kns;
  
  // knm: add kns*mbar
  RCP<LINALG::SparseMatrix> knmmod = LINALG::Multiply(*kns,false,*mhatmatrix_,false,false);
  knmmod->Add(*knm,false,1.0,1.0);
  knmmod->Complete(knm->DomainMap(),knm->RowMap());
  
  // knn: nothing to do
  RCP<LINALG::SparseMatrix> knnmod = knn;
  
  // output for checking everything
#ifdef CONTACTDIMOUTPUT
  if(Comm().MyPID()==0)
  {
    cout << endl << "*****************" << endl;
    cout << "kmod:   " << kteff->EpetraMatrix()->NumGlobalRows() << " x " << kteff->EpetraMatrix()->NumGlobalCols() << endl;
    if (kssmod!=null) cout << "kssmod: " << kssmod->EpetraMatrix()->NumGlobalRows() << " x " << kssmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kssmod: null" << endl;
    if (ksmmod!=null) cout << "ksmmod: " << ksmmod->EpetraMatrix()->NumGlobalRows() << " x " << ksmmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "ksmmod: null" << endl;
    if (ksnmod!=null) cout << "ksnmod: " << ksnmod->EpetraMatrix()->NumGlobalRows() << " x " << ksnmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "ksnmod: null" << endl;
    if (kmsmod!=null) cout << "kmsmod: " << kmsmod->EpetraMatrix()->NumGlobalRows() << " x " << kmsmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kmsmod: null" << endl;
    if (kmmmod!=null) cout << "kmmmod: " << kmmmod->EpetraMatrix()->NumGlobalRows() << " x " << kmmmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kmmmod: null" << endl;
    if (kmnmod!=null) cout << "kmnmod: " << kmnmod->EpetraMatrix()->NumGlobalRows() << " x " << kmnmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kmnmod: null" << endl;
    if (knsmod!=null) cout << "knsmod: " << knsmod->EpetraMatrix()->NumGlobalRows() << " x " << knsmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "knsmod: null" << endl;
    if (knmmod!=null) cout << "knmmod: " << knmmod->EpetraMatrix()->NumGlobalRows() << " x " << knmmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "knmmod: null" << endl;
    if (knnmod!=null) cout << "knnmod: " << knnmod->EpetraMatrix()->NumGlobalRows() << " x " << knnmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "knnmod: null" << endl;
    cout << "*****************" << endl;
  }
#endif // #ifdef CONTACTDIMOUTPUT
  
  /**********************************************************************/
  /* Apply basis transformation to f                                    */
  /**********************************************************************/
  // fs: nothing to be done
  RCP<Epetra_Vector> fsmod = fs;
  
  // fm: add T(mbar)*fs
  RCP<Epetra_Vector> fmmod = rcp(new Epetra_Vector(*gmdofrowmap_));
  mhatmatrix_->Multiply(true,*fs,*fmmod);
  fmmod->Update(1.0,*fm,1.0);
  
  // fn: nothing to be done
  RCP<Epetra_Vector> fnmod = fn;
  
  // output for checking everything
#ifdef CONTACTDIMOUTPUT
  if(Comm().MyPID()==0)
  {
    cout << endl << "**********" << endl;
    cout << "fmod:  " << feff->GlobalLength() << endl;
    cout << "fsmod: " << fsmod->GlobalLength() << endl;
    cout << "fmmod: " << fmmod->GlobalLength() << endl;
    cout << "fnmod: " << fnmod->GlobalLength() << endl;
    cout << "**********" << endl;
  }
#endif // #ifdef CONTACTDIMOUTPUT
  
  /**********************************************************************/
  /* Split slave quantities into active / inactive                      */
  /**********************************************************************/
  // we want to split kssmod into 2 groups a,i = 4 blocks
  RCP<LINALG::SparseMatrix> kaamod, kaimod, kiamod, kiimod;
  
  // we want to split ksnmod / ksmmod into 2 groups a,i = 2 blocks
  RCP<LINALG::SparseMatrix> kanmod, kinmod, kammod, kimmod;
    
  // we will get the i rowmap as a by-product
  RCP<Epetra_Map> gidofs;
    
  // do the splitting
  LINALG::SplitMatrix2x2(kssmod,gactivedofs_,gidofs,gactivedofs_,gidofs,kaamod,kaimod,kiamod,kiimod);
  LINALG::SplitMatrix2x2(ksnmod,gactivedofs_,gidofs,gndofrowmap_,tempmap,kanmod,tempmtx1,kinmod,tempmtx2);
  LINALG::SplitMatrix2x2(ksmmod,gactivedofs_,gidofs,gmdofrowmap_,tempmap,kammod,tempmtx1,kimmod,tempmtx2);
  
  // output for checking everything
#ifdef CONTACTDIMOUTPUT
  if(Comm().MyPID()==0)
  {
    cout << endl << "*****************" << endl;
    cout << "kssmod: " << kssmod->EpetraMatrix()->NumGlobalRows() << " x " << kssmod->EpetraMatrix()->NumGlobalCols() << endl;
    if (kaamod!=null) cout << "kaamod: " << kaamod->EpetraMatrix()->NumGlobalRows() << " x " << kaamod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kaamod: null" << endl;
    if (kaimod!=null) cout << "kaimod: " << kaimod->EpetraMatrix()->NumGlobalRows() << " x " << kaimod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kaimod: null" << endl;
    if (kiamod!=null) cout << "kiamod: " << kiamod->EpetraMatrix()->NumGlobalRows() << " x " << kiamod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kiamod: null" << endl;
    if (kiimod!=null) cout << "kiimod: " << kiimod->EpetraMatrix()->NumGlobalRows() << " x " << kiimod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kiimod: null" << endl;
    cout << "*****************" << endl;
    
    cout << endl << "*****************" << endl;
    cout << "ksnmod: " << ksnmod->EpetraMatrix()->NumGlobalRows() << " x " << ksnmod->EpetraMatrix()->NumGlobalCols() << endl;
    if (kanmod!=null) cout << "kanmod: " << kanmod->EpetraMatrix()->NumGlobalRows() << " x " << kanmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kanmod: null" << endl;
    if (kinmod!=null) cout << "kinmod: " << kinmod->EpetraMatrix()->NumGlobalRows() << " x " << kinmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kinmod: null" << endl;
    cout << "*****************" << endl;
    
    cout << endl << "*****************" << endl;
    cout << "ksmmod: " << ksmmod->EpetraMatrix()->NumGlobalRows() << " x " << ksmmod->EpetraMatrix()->NumGlobalCols() << endl;
    if (kammod!=null) cout << "kammod: " << kammod->EpetraMatrix()->NumGlobalRows() << " x " << kammod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kammod: null" << endl;
    if (kimmod!=null) cout << "kimmod: " << kimmod->EpetraMatrix()->NumGlobalRows() << " x " << kimmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kimmod: null" << endl;
    cout << "*****************" << endl;
  }
#endif // #ifdef CONTACTDIMOUTPUT
  
  // we want to split fsmod into 2 groups a,i
  RCP<Epetra_Vector> famod, fimod;
  
  // do the vector splitting s -> a+i
  if (!gidofs->NumGlobalElements())
    famod = rcp(new Epetra_Vector(*fsmod));
  else if (!gactivedofs_->NumGlobalElements())
    fimod = rcp(new Epetra_Vector(*fsmod));
  else
  {
    LINALG::SplitVector(*fsmod,*gactivedofs_,famod,*gidofs,fimod);
  }
  
  // output for checking everything
#ifdef CONTACTDIMOUTPUT
  if(Comm().MyPID()==0)
  {
    cout << endl << "**********" << endl;
    if (fsmod!=null) cout << "fsmod: " << fsmod->GlobalLength() << endl;
    else cout << "fsmod: null" << endl;
    if (famod!=null) cout << "famod: " << famod->GlobalLength() << endl;
    else cout << "famod: null" << endl;
    if (fimod!=null) cout << "fimod: " << fimod->GlobalLength() << endl;
    else cout << "fimod: null" << endl;
    cout << "**********" << endl;
  }
#endif // #ifdef CONTACTDIMOUTPUT
  
  /**********************************************************************/
  /* Isolate active / inactive part from dold                           */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> dolda, doldi;
  LINALG::SplitMatrix2x2(dold_,gactivedofs_,gidofs,gactivedofs_,gidofs,dolda,tempmtx1,tempmtx2,doldi);
  
  /**********************************************************************/
  /* Gen-alpha modifications                                            */
  /**********************************************************************/
  // fi: subtract alphaf * old contact forces (t_n)
  if (gidofs->NumGlobalElements())
  {
    RCP<Epetra_Vector> modi = rcp(new Epetra_Vector(*gidofs));
    LINALG::Export(*zold_,*modi);
    RCP<Epetra_Vector> tempveci = rcp(new Epetra_Vector(*gidofs));
    doldi->Multiply(false,*modi,*tempveci);
    fimod->Update(-alphaf_,*tempveci,1.0);
  }
   
  // fa: subtract alphaf * old contact forces (t_n)
  if (gactivedofs_->NumGlobalElements())
  {
    RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(*gactivedofs_));
    LINALG::Export(*zold_,*mod);
    RCP<Epetra_Vector> tempvec = rcp(new Epetra_Vector(*gactivedofs_));
    dolda->Multiply(false,*mod,*tempvec);
    famod->Update(-alphaf_,*tempvec,1.0);
  }
   
  // do the multiplications with t matrix
  RCP<LINALG::SparseMatrix> tkanmod, tkammod, tkaimod, tkaamod;
  RCP<Epetra_Vector> tfamod;
  
  if(gactivedofs_->NumGlobalElements())
  {
    tkanmod = LINALG::Multiply(*tmatrix_,false,*kanmod,false);
    tkammod = LINALG::Multiply(*tmatrix_,false,*kammod,false);
    tkaamod = LINALG::Multiply(*tmatrix_,false,*kaamod,false);
    
    if (gidofs->NumGlobalElements())
      tkaimod = LINALG::Multiply(*tmatrix_,false,*kaimod,false);
    
    tfamod = rcp(new Epetra_Vector(tmatrix_->RowMap()));
    tmatrix_->Multiply(false,*famod,*tfamod);
  }
  
  // output for checking everything
#ifdef CONTACTDIMOUTPUT
  if(Comm().MyPID()==0)
  {
    cout << endl << "*****************" << endl;
    if (tkanmod!=null) cout << "tkanmod: " << tkanmod->EpetraMatrix()->NumGlobalRows() << " x " << tkanmod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "tkanmod: null" << endl;
    if (tkammod!=null) cout << "tkammod: " << tkammod->EpetraMatrix()->NumGlobalRows() << " x " << tkammod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "tkammod: null" << endl;
    if (tkaimod!=null) cout << "tkaimod: " << tkaimod->EpetraMatrix()->NumGlobalRows() << " x " << tkaimod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "tkaimod: null" << endl;
    if (tkaamod!=null) cout << "tkaamod: " << tkaamod->EpetraMatrix()->NumGlobalRows() << " x " << tkaamod->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "tkaamod: null" << endl;
    cout << "*****************" << endl;
    
    cout << endl << "**********" << endl;
    if (tfamod!=null) cout << "tfamod: " << tfamod->GlobalLength() << endl;
    else cout << "tfamod: null" << endl;
    cout << "**********" << endl;
  }
#endif // #ifdef CONTACTDIMOUTPUT
  
  /**********************************************************************/
  /* Global setup of kteffnew, feffnew (including contact)              */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> kteffnew = rcp(new LINALG::SparseMatrix(*(discret_.DofRowMap()),81));
  RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*(discret_.DofRowMap()));
  
  // add n / m submatrices to kteffnew
  kteffnew->Add(*knnmod,false,1.0,1.0);
  kteffnew->Add(*knmmod,false,1.0,1.0);
  kteffnew->Add(*kmnmod,false,1.0,1.0);
  kteffnew->Add(*kmmmod,false,1.0,1.0);
  
  // add a / i submatrices to kteffnew, if existing
  if (knsmod!=null) kteffnew->Add(*knsmod,false,1.0,1.0);
  if (kmsmod!=null) kteffnew->Add(*kmsmod,false,1.0,1.0);
  if (kinmod!=null) kteffnew->Add(*kinmod,false,1.0,1.0);
  if (kimmod!=null) kteffnew->Add(*kimmod,false,1.0,1.0);
  if (kiimod!=null) kteffnew->Add(*kiimod,false,1.0,1.0);
  if (kiamod!=null) kteffnew->Add(*kiamod,false,1.0,1.0);
  
  // add matrix of normals to kteffnew
  kteffnew->Add(*nmatrix_,false,1.0,1.0);

  if (ftype=="none")
  {
    // add submatrices with tangents to kteffnew, if existing
    if (tkanmod!=null) kteffnew->Add(*tkanmod,false,1.0,1.0);
    if (tkammod!=null) kteffnew->Add(*tkammod,false,1.0,1.0);
    if (tkaimod!=null) kteffnew->Add(*tkaimod,false,1.0,1.0);
    if (tkaamod!=null) kteffnew->Add(*tkaamod,false,1.0,1.0);
  }
  else if (ftype=="stick")
  {
    // add matrix of tangents to kteffnew
    if (tmatrix_!=null) kteffnew->Add(*tmatrix_,false,1.0,1.0);
  }
  else
    dserror("ERROR: Evaluate: Invalid type of friction law");
  
  // FillComplete kteffnew (square)
  kteffnew->Complete();
  
  // add n / m subvectors to feffnew
  RCP<Epetra_Vector> fnmodexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  RCP<Epetra_Vector> fmmodexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  LINALG::Export(*fnmod,*fnmodexp);
  LINALG::Export(*fmmod,*fmmodexp);
  feffnew->Update(1.0,*fnmodexp,1.0,*fmmodexp,1.0);
  
  // add i / ta subvectors to feffnew, if existing
  RCP<Epetra_Vector> fimodexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  RCP<Epetra_Vector> tfamodexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  if (fimod!=null) LINALG::Export(*fimod,*fimodexp);
  if (tfamod!=null) LINALG::Export(*tfamod,*tfamodexp);
  
  if (ftype=="none")
   feffnew->Update(1.0,*fimodexp,1.0,*tfamodexp,1.0);
  else if (ftype=="stick")
    feffnew->Update(1.0,*fimodexp,0.0,*tfamodexp,1.0);
  else
    dserror("ERROR: Evaluate: Invalid type of friction law");
  
  if (ctype!="meshtying")
  {
    // add weighted gap vector to feffnew, if existing
    RCP<Epetra_Vector> gexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
    if (gact->GlobalLength()) LINALG::Export(*gact,*gexp);
    feffnew->Update(1.0,*gexp,1.0);
  }
  
  /**********************************************************************/
  /* Replace kteff and feff by kteffnew and feffnew                   */
  /**********************************************************************/
  *kteff = *kteffnew;
  *feff = *feffnew;
  //LINALG::PrintSparsityToPostscript(*(kteff->EpetraMatrix()));
  
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate contact without basis transformation (public)    popp 04/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::EvaluateNoBasisTrafo(RCP<LINALG::SparseMatrix> kteff,
                                            RCP<Epetra_Vector> feff, int numiter)
{ 
  // FIXME: Currently only the old LINALG::Multiply method is used,
  // because there are still problems with the transposed version of
  // MLMultiply if a row has no entries! One day we should use ML...
  
  // input parameters
  string ctype   = scontact_.get<string>("contact type","none");
  string ftype   = scontact_.get<string>("friction type","none");
  bool fulllin   = scontact_.get<bool>("full linearization",false);
  
  /**********************************************************************/
  /* evaluate interfaces                                                */
  /* (nodal normals, projections, Mortar integration, Mortar assembly)  */
  /**********************************************************************/
#ifdef CONTACTCHECKHUEEBER
  if (numiter==0)
  {
#endif // #ifdef CONTACTCHECKHUEEBER
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->Evaluate();
    interface_[i]->AssembleDMG(*dmatrix_,*mmatrix_,*g_);
  }
  
  // FillComplete() global Mortar matrices
  dmatrix_->Complete();
  mmatrix_->Complete(*gmdofrowmap_,*gsdofrowmap_);
  
  
#ifdef CONTACTCHECKHUEEBER
  }
#endif // #ifdef CONTACTCHECKHUEEBER
  // export weighted gap vector to gactiveN-map
  RCP<Epetra_Vector> gact = LINALG::CreateVector(*gactivenodes_,true);
  if (gact->GlobalLength())
  {
    LINALG::Export(*g_,*gact);
    gact->ReplaceMap(*gactiven_);
  }
  
  /**********************************************************************/
  /* build global matrix n with normal vectors of active nodes          */
  /* and global matrix t with tangent vectors of active nodes           */
  /* and global matrix s with normal derivatives of active nodes        */
  /**********************************************************************/
  // here and for the splitting later, we need the combined sm rowmap
  // (this map is NOT allowed to have an overlap !!!)
  RCP<Epetra_Map> gsmdofs = LINALG::MergeMap(gsdofrowmap_,gmdofrowmap_,false);
  
#ifdef CONTACTCHECKHUEEBER
  if (numiter==0)
  {
#endif // #ifdef CONTACTCHECKHUEEBER
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->AssembleNT(*nmatrix_,*tmatrix_);
    interface_[i]->AssembleS(*smatrix_);
    interface_[i]->AssembleP(*pmatrix_);
    interface_[i]->AssembleLinDM(*lindmatrix_,*linmmatrix_);
  }
    
  // FillComplete() global matrices N and T
  nmatrix_->Complete(*gactivedofs_,*gactiven_);
  tmatrix_->Complete(*gactivedofs_,*gactivet_);
  
  // FillComplete() global matrix S
  smatrix_->Complete(*gsmdofs,*gactiven_);
  
  // FillComplete() global matrix P
  pmatrix_->Complete(*gsdofrowmap_,*gactivet_);
  
  // FillComplete() global matrices LinD, LinM
  lindmatrix_->Complete();
  linmmatrix_->Complete(*gsmdofs,*gmdofrowmap_);
  
#ifdef CONTACTCHECKHUEEBER
  }
#endif // #ifdef CONTACTCHECKHUEEBER
  /**********************************************************************/
  /* Multiply Mortar matrices: m^ = inv(d) * m                          */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> invd = rcp(new LINALG::SparseMatrix(*dmatrix_));
  RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_,true);
  int err = 0;
  
  // extract diagonal of invd into diag
  invd->ExtractDiagonalCopy(*diag);
  
  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
    if ((*diag)[i]==0.0) (*diag)[i]=1.0;
  
  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");
  
  // re-insert inverted diagonal into invd
  err = invd->ReplaceDiagonalValues(*diag);
  // we cannot use this check, as we deliberately replaced zero entries
  //if (err>0) dserror("ERROR: ReplaceDiagonalValues: Missing diagonal entry!");
  
  // do the multiplication M^ = inv(D) * M
  mhatmatrix_ = LINALG::Multiply(*invd,false,*mmatrix_,false);
  
  /**********************************************************************/
  /* Split kteff into 3x3 block matrix                                  */
  /**********************************************************************/
  // we want to split k into 3 groups s,m,n = 9 blocks
  RCP<LINALG::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;
  
  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!) 
  RCP<LINALG::SparseMatrix> ksmsm, ksmn, knsm;
  
  // some temporary RCPs
  RCP<Epetra_Map> tempmap;
  RCP<LINALG::SparseMatrix> tempmtx1;
  RCP<LINALG::SparseMatrix> tempmtx2;
  RCP<LINALG::SparseMatrix> tempmtx3;
  
  // split into slave/master part + structure part
  LINALG::SplitMatrix2x2(kteff,gsmdofs,gndofrowmap_,gsmdofs,gndofrowmap_,ksmsm,ksmn,knsm,knn);
  
  // further splits into slave part + master part
  LINALG::SplitMatrix2x2(ksmsm,gsdofrowmap_,gmdofrowmap_,gsdofrowmap_,gmdofrowmap_,kss,ksm,kms,kmm);
  LINALG::SplitMatrix2x2(ksmn,gsdofrowmap_,gmdofrowmap_,gndofrowmap_,tempmap,ksn,tempmtx1,kmn,tempmtx2);
  LINALG::SplitMatrix2x2(knsm,gndofrowmap_,tempmap,gsdofrowmap_,gmdofrowmap_,kns,knm,tempmtx1,tempmtx2);
  
  // output for checking everything
#ifdef CONTACTDIMOUTPUT
  if(Comm().MyPID()==0)
  {
    cout << endl << "*****************" << endl;
    cout << "k:   " << kteff->EpetraMatrix()->NumGlobalRows() << " x " << kteff->EpetraMatrix()->NumGlobalCols() << endl;
    if (kss!=null) cout << "kss: " << kss->EpetraMatrix()->NumGlobalRows() << " x " << kss->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kss: null" << endl;
    if (ksm!=null) cout << "ksm: " << ksm->EpetraMatrix()->NumGlobalRows() << " x " << ksm->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "ksm: null" << endl;
    if (ksn!=null) cout << "ksn: " << ksn->EpetraMatrix()->NumGlobalRows() << " x " << ksn->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "ksn: null" << endl;
    if (kms!=null) cout << "kms: " << kms->EpetraMatrix()->NumGlobalRows() << " x " << kms->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kms: null" << endl;
    if (kmm!=null) cout << "kmm: " << kmm->EpetraMatrix()->NumGlobalRows() << " x " << kmm->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kmm: null" << endl;
    if (kmn!=null) cout << "kmn: " << kmn->EpetraMatrix()->NumGlobalRows() << " x " << kmn->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kmn: null" << endl;
    if (kns!=null) cout << "kns: " << kns->EpetraMatrix()->NumGlobalRows() << " x " << kns->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kns: null" << endl;
    if (knm!=null) cout << "knm: " << knm->EpetraMatrix()->NumGlobalRows() << " x " << knm->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "knm: null" << endl;
    if (knn!=null) cout << "knn: " << knn->EpetraMatrix()->NumGlobalRows() << " x " << knn->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "knn: null" << endl;
    cout << "*****************" << endl;
  }
#endif // #ifdef CONTACTDIMOUTPUT
  
  /**********************************************************************/
  /* Split feff into 3 subvectors                                       */
  /**********************************************************************/
  // we want to split f into 3 groups s.m,n
  RCP<Epetra_Vector> fs, fm, fn;
  
  // temporarily we need the group sm
  RCP<Epetra_Vector> fsm;
  
  // do the vector splitting smn -> sm+n -> s+m+n
  LINALG::SplitVector(*feff,*gsmdofs,fsm,*gndofrowmap_,fn);
  LINALG::SplitVector(*fsm,*gsdofrowmap_,fs,*gmdofrowmap_,fm);
  
  // output for checking everything
#ifdef CONTACTDIMOUTPUT
  if(Comm().MyPID()==0)
  {
    cout << endl << "**********" << endl;
    cout << "f:  " << feff->GlobalLength() << endl;
    cout << "fs: " << fs->GlobalLength() << endl;
    cout << "fm: " << fm->GlobalLength() << endl;
    cout << "fn: " << fn->GlobalLength() << endl;
    cout << "**********" << endl;
  }
#endif // #ifdef CONTACTDIMOUTPUT
  
  // store some stuff for static condensation of LM
  fs_   = fs;
  invd_ = invd;
  ksn_  = ksn;
  ksm_  = ksm;
  kss_  = kss;
    
  /**********************************************************************/
  /* Split slave quantities into active / inactive                      */
  /**********************************************************************/
  // we want to split kssmod into 2 groups a,i = 4 blocks
  RCP<LINALG::SparseMatrix> kaa, kai, kia, kii;
  
  // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
  RCP<LINALG::SparseMatrix> kan, kin, kam, kim, kma, kmi;
    
  // we will get the i rowmap as a by-product
  RCP<Epetra_Map> gidofs;
    
  // do the splitting
  LINALG::SplitMatrix2x2(kss,gactivedofs_,gidofs,gactivedofs_,gidofs,kaa,kai,kia,kii);
  LINALG::SplitMatrix2x2(ksn,gactivedofs_,gidofs,gndofrowmap_,tempmap,kan,tempmtx1,kin,tempmtx2);
  LINALG::SplitMatrix2x2(ksm,gactivedofs_,gidofs,gmdofrowmap_,tempmap,kam,tempmtx1,kim,tempmtx2);
  LINALG::SplitMatrix2x2(kms,gmdofrowmap_,tempmap,gactivedofs_,gidofs,kma,kmi,tempmtx1,tempmtx2);
  
  // output for checking everything
#ifdef CONTACTDIMOUTPUT
  if(Comm().MyPID()==0)
  {
    cout << endl << "*****************" << endl;
    cout << "kss: " << kss->EpetraMatrix()->NumGlobalRows() << " x " << kss->EpetraMatrix()->NumGlobalCols() << endl;
    if (kaa!=null) cout << "kaa: " << kaa->EpetraMatrix()->NumGlobalRows() << " x " << kaa->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kaa: null" << endl;
    if (kai!=null) cout << "kai: " << kai->EpetraMatrix()->NumGlobalRows() << " x " << kai->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kai: null" << endl;
    if (kia!=null) cout << "kia: " << kia->EpetraMatrix()->NumGlobalRows() << " x " << kia->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kia: null" << endl;
    if (kii!=null) cout << "kii: " << kii->EpetraMatrix()->NumGlobalRows() << " x " << kii->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kii: null" << endl;
    cout << "*****************" << endl;
    
    cout << endl << "*****************" << endl;
    cout << "ksn: " << ksn->EpetraMatrix()->NumGlobalRows() << " x " << ksn->EpetraMatrix()->NumGlobalCols() << endl;
    if (kan!=null) cout << "kan: " << kan->EpetraMatrix()->NumGlobalRows() << " x " << kan->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kan: null" << endl;
    if (kin!=null) cout << "kin: " << kin->EpetraMatrix()->NumGlobalRows() << " x " << kin->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kin: null" << endl;
    cout << "*****************" << endl;
    
    cout << endl << "*****************" << endl;
    cout << "ksm: " << ksm->EpetraMatrix()->NumGlobalRows() << " x " << ksm->EpetraMatrix()->NumGlobalCols() << endl;
    if (kam!=null) cout << "kam: " << kam->EpetraMatrix()->NumGlobalRows() << " x " << kam->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kam: null" << endl;
    if (kim!=null) cout << "kim: " << kim->EpetraMatrix()->NumGlobalRows() << " x " << kim->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kim: null" << endl;
    cout << "*****************" << endl;
    
    cout << endl << "*****************" << endl;
    cout << "kms: " << kms->EpetraMatrix()->NumGlobalRows() << " x " << kms->EpetraMatrix()->NumGlobalCols() << endl;
    if (kma!=null) cout << "kma: " << kma->EpetraMatrix()->NumGlobalRows() << " x " << kma->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kma: null" << endl;
    if (kmi!=null) cout << "kmi: " << kmi->EpetraMatrix()->NumGlobalRows() << " x " << kmi->EpetraMatrix()->NumGlobalCols() << endl;
    else cout << "kmi: null" << endl;
    cout << "*****************" << endl;
  }
#endif // #ifdef CONTACTDIMOUTPUT
  
  // we want to split fsmod into 2 groups a,i
  RCP<Epetra_Vector> fa = rcp(new Epetra_Vector(*gactivedofs_));
  RCP<Epetra_Vector> fi = rcp(new Epetra_Vector(*gidofs));
  
  // do the vector splitting s -> a+i
  if (!gidofs->NumGlobalElements())
    *fa = *fs;
  else if (!gactivedofs_->NumGlobalElements())
    *fi = *fs;
  else
  {
    LINALG::SplitVector(*fs,*gactivedofs_,fa,*gidofs,fi);
  }
  
  // output for checking everything
#ifdef CONTACTDIMOUTPUT
  if(Comm().MyPID()==0)
  {
    cout << endl << "**********" << endl;
    if (fs!=null) cout << "fs: " << fs->GlobalLength() << endl;
    else cout << "fs: null" << endl;
    if (fa!=null) cout << "fa: " << fa->GlobalLength() << endl;
    else cout << "fa: null" << endl;
    if (fi!=null) cout << "fi: " << fi->GlobalLength() << endl;
    else cout << "fi: null" << endl;
    cout << "**********" << endl;
  }
#endif // #ifdef CONTACTDIMOUTPUT
  
  /**********************************************************************/
  /* Isolate active part from mhat, invd and dold                       */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> mhata;
  LINALG::SplitMatrix2x2(mhatmatrix_,gactivedofs_,gidofs,gmdofrowmap_,tempmap,mhata,tempmtx1,tempmtx2,tempmtx3);
  mhata_=mhata;
  
  RCP<LINALG::SparseMatrix> invda;
  LINALG::SplitMatrix2x2(invd_,gactivedofs_,gidofs,gactivedofs_,gidofs,invda,tempmtx1,tempmtx2,tempmtx3);
  invda->Scale(1/(1-alphaf_));
  
  RCP<LINALG::SparseMatrix> dolda, doldi;
  LINALG::SplitMatrix2x2(dold_,gactivedofs_,gidofs,gactivedofs_,gidofs,dolda,tempmtx1,tempmtx2,doldi);
  
  /**********************************************************************/
  /* Split LinD and LinM into blocks                                    */
  /**********************************************************************/
  // we want to split lindmatrix_ into 2 groups a,i = 4 blocks
  RCP<LINALG::SparseMatrix> lindai, lindaa;
  
  // we want to split linmmatrix_ into 3 groups a,i,m = 3 blocks
  RCP<LINALG::SparseMatrix> linmmi, linmma, linmmm, linmms;
  
  if (fulllin)
  {
    // do the splitting
    LINALG::SplitMatrix2x2(lindmatrix_,gactivedofs_,gidofs,gactivedofs_,gidofs,lindaa,lindai,tempmtx1,tempmtx2);
    LINALG::SplitMatrix2x2(linmmatrix_,gmdofrowmap_,tempmap,gmdofrowmap_,gsdofrowmap_,linmmm,linmms,tempmtx1,tempmtx2);
    LINALG::SplitMatrix2x2(linmms,gmdofrowmap_,tempmap,gactivedofs_,gidofs,linmma,linmmi,tempmtx1,tempmtx2);
  
    // modification of kai, kaa
    // (this has to be done first as they are needed below)
    kai->Add(*lindai,false,1.0-alphaf_,1.0);
    kaa->Add(*lindaa,false,1.0-alphaf_,1.0);
  }  
  
  /**********************************************************************/
  /* Build the final K and f blocks                                     */
  /**********************************************************************/
  // knn: nothing to do
  
  // knm: nothing to do
  
  // kns: nothing to do
  
  // kmn: add T(mbaractive)*kan
  RCP<LINALG::SparseMatrix> kmnmod = LINALG::Multiply(*mhata,true,*kan,false,false);
  kmnmod->Add(*kmn,false,1.0,1.0);
  kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());
  
  // kmm: add T(mbaractive)*kam
  RCP<LINALG::SparseMatrix> kmmmod = LINALG::Multiply(*mhata,true,*kam,false,false);
  kmmmod->Add(*kmm,false,1.0,1.0);
  if (fulllin) kmmmod->Add(*linmmm,false,1.0-alphaf_,1.0);
  kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());
  
  // kmi: add T(mbaractive)*kai
  RCP<LINALG::SparseMatrix> kmimod = LINALG::Multiply(*mhata,true,*kai,false,false);
  kmimod->Add(*kmi,false,1.0,1.0);
  if (fulllin) kmimod->Add(*linmmi,false,1.0-alphaf_,1.0);
  kmimod->Complete(kmi->DomainMap(),kmi->RowMap());
  
  // kma: add T(mbaractive)*kaa
  RCP<LINALG::SparseMatrix> kmamod = LINALG::Multiply(*mhata,true,*kaa,false,false);
  kmamod->Add(*kma,false,1.0,1.0);
  if (fulllin) kmamod->Add(*linmma,false,1.0-alphaf_,1.0);
  kmamod->Complete(kma->DomainMap(),kma->RowMap());
  
  // kin: nothing to do
  
  // kim: nothing to do
  
  // kii: nothing to do
  
  // kia: nothing to do
  
  // n*mbaractive: do the multiplication
  RCP<LINALG::SparseMatrix> nmhata = LINALG::Multiply(*nmatrix_,false,*mhata,false,true);
  
  // nmatrix: nothing to do
  
  // kan: multiply with tmatrix
  RCP<LINALG::SparseMatrix> kanmod = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
  kanmod = LINALG::Multiply(*kanmod,false,*kan,false,true);
  
  // kam: multiply with tmatrix
  RCP<LINALG::SparseMatrix> kammod = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
  kammod = LINALG::Multiply(*kammod,false,*kam,false,true);
  
  // kai: multiply with tmatrix
  RCP<LINALG::SparseMatrix> kaimod = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
  kaimod = LINALG::Multiply(*kaimod,false,*kai,false,true);
  
  // kaa: multiply with tmatrix
  RCP<LINALG::SparseMatrix> kaamod = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
  kaamod = LINALG::Multiply(*kaamod,false,*kaa,false,true);
  
  // t*mbaractive: do the multiplication
  RCP<LINALG::SparseMatrix> tmhata = LINALG::Multiply(*tmatrix_,false,*mhata,false,true);
  
  // fn: nothing to do
  
  // fi: subtract alphaf * old contact forces (t_n)
  if (gidofs->NumGlobalElements())
  {
    RCP<Epetra_Vector> modi = rcp(new Epetra_Vector(*gidofs));
    LINALG::Export(*zold_,*modi);
    RCP<Epetra_Vector> tempveci = rcp(new Epetra_Vector(*gidofs));
    doldi->Multiply(false,*modi,*tempveci);
    fi->Update(-alphaf_,*tempveci,1.0);
  }
  
  // fa: subtract alphaf * old contact forces (t_n)
  if (gactivedofs_->NumGlobalElements())
  {
    RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(*gactivedofs_));
    LINALG::Export(*zold_,*mod);
    RCP<Epetra_Vector> tempvec = rcp(new Epetra_Vector(*gactivedofs_));
    dolda->Multiply(false,*mod,*tempvec);
    fa->Update(-alphaf_,*tempvec,1.0);
  }
    
  // fm: add alphaf * old contact forces (t_n)
  RCP<Epetra_Vector> tempvecm = rcp(new Epetra_Vector(*gmdofrowmap_));
  mold_->Multiply(true,*zold_,*tempvecm);
  fm->Update(alphaf_,*tempvecm,1.0);
    
  // fm: add T(mbaractive)*fa
  RCP<Epetra_Vector> fmmod = rcp(new Epetra_Vector(*gmdofrowmap_));
  mhata->Multiply(true,*fa,*fmmod);
  fmmod->Update(1.0,*fm,1.0);
    
  // fa: mutliply with tmatrix
  // (this had to wait as we had to modify fm first)
  RCP<Epetra_Vector> famod = rcp(new Epetra_Vector(*gactivet_));
  RCP<LINALG::SparseMatrix> tinvda = LINALG::Multiply(*tmatrix_,false,*invda,false,true);
  tinvda->Multiply(false,*fa,*famod);
  
  // gactive: nothing to do
  
#ifdef CONTACTFDGAP
  // FD check of weighted gap g derivatives
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    RCP<LINALG::SparseMatrix> deriv = rcp(new LINALG::SparseMatrix(*gactiven_,81));
    deriv->Add(*nmatrix_,false,1.0,1.0);
    deriv->Add(*smatrix_,false,1.0,1.0);
    deriv->Add(*nmhata,false,-1.0,1.0);
    deriv->Complete(*gsmdofs,*gactiven_);
    cout << *deriv << endl;
    interface_[i]->FDCheckGapDeriv();
  }
#endif // #ifdef CONTACTFDGAP
  
#ifdef CONTACTFDTANGLM
  // FD check of tangential LM derivatives (frictionless condition)
  for (int i=0; i<(int)interface_.size();++i)
  {
    cout << *pmatrix_ << endl;
    interface_[i]->FDCheckTangLMDeriv();
  }
#endif // #ifdef CONTACTFDTANGLM
    
  /**********************************************************************/
  /* Global setup of kteffnew, feffnew (including contact)              */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> kteffnew = rcp(new LINALG::SparseMatrix(*(discret_.DofRowMap()),81));
  RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*(discret_.DofRowMap()));
  
  // add n submatrices to kteffnew
  kteffnew->Add(*knn,false,1.0,1.0);
  kteffnew->Add(*knm,false,1.0,1.0);
  kteffnew->Add(*kns,false,1.0,1.0);

  // add m submatrices to kteffnew
  kteffnew->Add(*kmnmod,false,1.0,1.0);
  kteffnew->Add(*kmmmod,false,1.0,1.0);
  kteffnew->Add(*kmimod,false,1.0,1.0);
  kteffnew->Add(*kmamod,false,1.0,1.0);
  
  // add i submatrices to kteffnew
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kin,false,1.0,1.0);
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kim,false,1.0,1.0);
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kii,false,1.0,1.0);
  if (gidofs->NumGlobalElements()) kteffnew->Add(*kia,false,1.0,1.0);
  
  // add matrix nmhata to keteffnew
  if (gactiven_->NumGlobalElements()) kteffnew->Add(*nmhata,false,-1.0,1.0);

  // add matrix n to kteffnew
  if (gactiven_->NumGlobalElements()) kteffnew->Add(*nmatrix_,false,1.0,1.0);
  
  // add full linearization terms to kteffnew
  if (fulllin)
  {
   if (gactiven_->NumGlobalElements()) kteffnew->Add(*smatrix_,false,1.0,1.0);
   if (gactivet_->NumGlobalElements()) kteffnew->Add(*pmatrix_,false,-1.0,1.0);
  }
  
  if (ftype=="none")
  {
    // add a submatrices to kteffnew
    if (gactivet_->NumGlobalElements()) kteffnew->Add(*kanmod,false,1.0,1.0);
    if (gactivet_->NumGlobalElements()) kteffnew->Add(*kammod,false,1.0,1.0);
    if (gactivet_->NumGlobalElements()) kteffnew->Add(*kaimod,false,1.0,1.0);
    if (gactivet_->NumGlobalElements()) kteffnew->Add(*kaamod,false,1.0,1.0);
  }
  else if (ftype=="stick")
  {
    // add matrices t and tmhata to kteffnew
    if (gactivet_->NumGlobalElements()) kteffnew->Add(*tmatrix_,false,1.0,1.0);
    if (gactivet_->NumGlobalElements()) kteffnew->Add(*tmhata,false,-1.0,1.0);
  }
  else
    dserror("ERROR: Evaluate: Invalid type of friction law");
  
  // FillComplete kteffnew (square)
  kteffnew->Complete();
  
  // add n subvector to feffnew
  RCP<Epetra_Vector> fnexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  LINALG::Export(*fn,*fnexp);
  feffnew->Update(1.0,*fnexp,1.0);
  
  // add m subvector to feffnew
  RCP<Epetra_Vector> fmmodexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  LINALG::Export(*fmmod,*fmmodexp);
  feffnew->Update(1.0,*fmmodexp,1.0);
  
  // add i subvector to feffnew
  RCP<Epetra_Vector> fiexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  LINALG::Export(*fi,*fiexp);
  if (gidofs->NumGlobalElements()) feffnew->Update(1.0,*fiexp,1.0);
  
  if (ctype!="meshtying")
  {
    // add weighted gap vector to feffnew, if existing
    RCP<Epetra_Vector> gexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
    LINALG::Export(*gact,*gexp);
    if (gact->GlobalLength()) feffnew->Update(1.0,*gexp,1.0);
  }
  
  if (ftype=="none")
  {
    // add a subvector to feffnew
    RCP<Epetra_Vector> famodexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
    LINALG::Export(*famod,*famodexp);
    if (gactivenodes_->NumGlobalElements())feffnew->Update(1.0,*famodexp,1.0);
  }
  else if (ftype=="stick")
  {
    // do nothing here
  }
  else
    dserror("ERROR: Invalid type of friction law");  
  
  /**********************************************************************/
  /* Replace kteff and feff by kteffnew and feffnew                     */
  /**********************************************************************/
  *kteff = *kteffnew;
  *feff = *feffnew;
  
  return;
}

/*----------------------------------------------------------------------*
 |  Recovery method for displacements and LM (public)         popp 04/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::Recover(RCP<Epetra_Vector> disi)
{ 
  // check if basis transformation should be applied
  bool btrafo = scontact_.get<bool>("basis transformation",false);
  
  if (btrafo)
    RecoverBasisTrafo(disi);
  else
    RecoverNoBasisTrafo(disi);
    
  return;
}

/*----------------------------------------------------------------------*
 |  Recovery method (basis trafo case)                        popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::RecoverBasisTrafo(RCP<Epetra_Vector> disi)
{ 
#ifdef CONTACTCHECKHUEEBER
  // debugging (check S. Hüeber)
  if (IsInContact())
  {
    RCP<Epetra_Vector> gupdate = rcp(new Epetra_Vector(*gactiven_));
    RCP<Epetra_Vector> activedisi =rcp(new Epetra_Vector(*gactivedofs_));
    LINALG::Export(*disi,*activedisi);
    nmatrix_->Multiply(false,*activedisi,*gupdate);
    gupdate->ReplaceMap(*gactivenodes_);
    RCP<Epetra_Vector> gupdateexp = rcp(new Epetra_Vector(*gsnoderowmap_));
    LINALG::Export(*gupdate,*gupdateexp);
    gupdateexp->Update(1.0,*g_,-1.0);
    g_=gupdateexp;
    //cout << *g_;
    //cout << *gupdateexp;
  } 
#endif // #ifdef CONTACTCHECKHUEEBER
  
  // extract incremental jump from disi (for active set)
  incrjump_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*disi,*incrjump_);
  
  // friction
  // sum up incremental jumps from active set nodes
  jump_->Update(1.0,*incrjump_,1.0);
  
  // friction
  // store updaded jumps to nodes
  StoreNodalQuantities("jump");
  
  // extract master displacements from disi
  RCP<Epetra_Vector> disim = rcp(new Epetra_Vector(*gmdofrowmap_));
  LINALG::Export(*disi,*disim);
  
  // recover slave displacement increments
  RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(mhatmatrix_->RowMap()));
  mhatmatrix_->Multiply(false,*disim,*mod);
  
  RCP<Epetra_Vector> modexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  LINALG::Export(*mod,*modexp);
  disi->Update(1.0,*modexp,1.0);
  
  /**********************************************************************/
  /* Update Lagrange multipliers                                        */
  /**********************************************************************/
  // approximate update
  //invd_->Multiply(false,*fs_,*z_);
  
  // full update
  z_->Update(1.0,*fs_,0.0);
  RCP<Epetra_Vector> mod2 = rcp(new Epetra_Vector(*gsdofrowmap_));
  RCP<Epetra_Vector> slavedisp = rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*disi,*slavedisp);
  kss_->Multiply(false,*slavedisp,*mod2);
  z_->Update(-1.0,*mod2,1.0);
  RCP<Epetra_Vector> masterdisp = rcp(new Epetra_Vector(*gmdofrowmap_));
  LINALG::Export(*disi,*masterdisp);
  ksm_->Multiply(false,*masterdisp,*mod2);
  z_->Update(-1.0,*mod2,1.0);
  RCP<Epetra_Vector> innerdisp = rcp(new Epetra_Vector(*gndofrowmap_));
  LINALG::Export(*disi,*innerdisp);
  ksn_->Multiply(false,*innerdisp,*mod2);
  z_->Update(-1.0,*mod2,1.0);
  dold_->Multiply(false,*zold_,*mod);
  z_->Update(-alphaf_,*mod,1.0);
  RCP<Epetra_Vector> zcopy = rcp(new Epetra_Vector(*z_));
  invd_->Multiply(false,*zcopy,*z_);
  z_->Scale(1/(1-alphaf_));
  
  // store updated LM into nodes
  StoreNodalQuantities("lmcurrent");
  
  /*
  // CHECK OF CONTACT COUNDARY CONDITIONS---------------------------------
#ifdef DEBUG
  //debugging (check for z_i = 0)
  RCP<Epetra_Map> gidofs = LINALG::SplitMap(*gsdofrowmap_,*gactivedofs_);
  if (gidofs->NumGlobalElements())
  {
    RCP<Epetra_Vector> zinactive = rcp(new Epetra_Vector(*gidofs));
    LINALG::Export(*z_,*zinactive);
    cout << *zinactive << endl;
  }
  
  // debugging (check for N*[d_a] = g_a and T*z_a = 0)
  if (gactivedofs_->NumGlobalElements())
  { 
    RCP<Epetra_Vector> activejump = rcp(new Epetra_Vector(*gactivedofs_));
    RCP<Epetra_Vector> gtest = rcp(new Epetra_Vector(*gactiven_));
    LINALG::Export(*incrjump_,*activejump);
    nmatrix_->Multiply(false,*activejump,*gtest);
    cout << *gtest << endl << *g_ << endl;
    
    RCP<Epetra_Vector> zactive = rcp(new Epetra_Vector(*gactivedofs_));
    RCP<Epetra_Vector> zerotest = rcp(new Epetra_Vector(*gactivet_));
    LINALG::Export(*z_,*zactive);
    tmatrix_->Multiply(false,*zactive,*zerotest);
    cout << *zerotest << endl;
    cout << *zactive << endl;
  }
#endif // #ifdef DEBUG
  // CHECK OF CONTACT BOUNDARY CONDITIONS---------------------------------
  */
  
  return;
}

/*----------------------------------------------------------------------*
 | Recovery method (no basis trafo case) (public)             popp 04/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::RecoverNoBasisTrafo(RCP<Epetra_Vector> disi)
{
#ifdef CONTACTCHECKHUEEBER
  // debugging (check S. Hüeber)
  if (IsInContact())
  {
    RCP<LINALG::SparseMatrix> nmhata = LINALG::Multiply(*nmatrix_,false,*mhata_,false,true);
    
    RCP<Epetra_Vector> gupdate = rcp(new Epetra_Vector(*gactiven_));
    RCP<Epetra_Vector> gupdate2 = rcp(new Epetra_Vector(*gactiven_));
    RCP<Epetra_Vector> activedisi =rcp(new Epetra_Vector(*gactivedofs_));
    RCP<Epetra_Vector> masterdisi =rcp(new Epetra_Vector(*gmdofrowmap_));
    LINALG::Export(*disi,*activedisi);
    LINALG::Export(*disi,*masterdisi);
    nmatrix_->Multiply(false,*activedisi,*gupdate);
    nmhata->Multiply(false,*masterdisi,*gupdate2);
    gupdate->ReplaceMap(*gactivenodes_);
    gupdate2->ReplaceMap(*gactivenodes_);
    RCP<Epetra_Vector> gupdateexp = rcp(new Epetra_Vector(*gsnoderowmap_));
    RCP<Epetra_Vector> gupdateexp2 = rcp(new Epetra_Vector(*gsnoderowmap_));
    LINALG::Export(*gupdate,*gupdateexp);
    LINALG::Export(*gupdate2,*gupdateexp2);
    gupdateexp->Update(1.0,*g_,-1.0);
    gupdateexp->Update(1.0,*gupdateexp2,1.0);
    g_=gupdateexp;
  } 
#endif // #ifdef CONTACTCHECKHUEEBER
  
  // extract slave displacements from disi
  RCP<Epetra_Vector> disis = rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*disi,*disis);
  
  // extract master displacements from disi
  RCP<Epetra_Vector> disim = rcp(new Epetra_Vector(*gmdofrowmap_));
  LINALG::Export(*disi,*disim);
  
  // recover incremental jump (for active set)
  incrjump_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  mhatmatrix_->Multiply(false,*disim,*incrjump_);
  incrjump_->Update(1.0,*disis,-1.0);
  
  /**********************************************************************/
  /* Update Lagrange multipliers z_n+1                                  */
  /**********************************************************************/
  // approximate update
  //invd_->Multiply(false,*fs_,*z_);
    
  // full update
  z_->Update(1.0,*fs_,0.0);
  RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(*gsdofrowmap_));
  kss_->Multiply(false,*disis,*mod);
  z_->Update(-1.0,*mod,1.0);
  ksm_->Multiply(false,*disim,*mod);
  z_->Update(-1.0,*mod,1.0);
  RCP<Epetra_Vector> disin = rcp(new Epetra_Vector(*gndofrowmap_));
  LINALG::Export(*disi,*disin);
  ksn_->Multiply(false,*disin,*mod);
  z_->Update(-1.0,*mod,1.0);
  dold_->Multiply(false,*zold_,*mod);
  z_->Update(-alphaf_,*mod,1.0);
  RCP<Epetra_Vector> zcopy = rcp(new Epetra_Vector(*z_));
  invd_->Multiply(false,*zcopy,*z_);
  z_->Scale(1/(1-alphaf_));
  
  // store updated LM into nodes
  StoreNodalQuantities("lmcurrent");
    
  /* 
  // CHECK OF CONTACT COUNDARY CONDITIONS---------------------------------
#ifdef DEBUG
  //debugging (check for z_i = 0)
  RCP<Epetra_Map> gidofs = LINALG::SplitMap(*gsdofrowmap_,*gactivedofs_);
  if (gidofs->NumGlobalElements())
  {
    RCP<Epetra_Vector> zinactive = rcp(new Epetra_Vector(*gidofs));
    LINALG::Export(*z_,*zinactive);
    cout << *zinactive << endl;
  }
  
  bool fulllin   = scontact_.get<bool>("full linearization",false);
  
  // debugging (check for N*[d_a] = g_a and T*z_a = 0)
  if (gactivedofs_->NumGlobalElements())
  { 
    RCP<Epetra_Vector> activejump = rcp(new Epetra_Vector(*gactivedofs_));
    RCP<Epetra_Vector> gtest = rcp(new Epetra_Vector(*gactiven_));
    RCP<Epetra_Vector> gtest2 = rcp(new Epetra_Vector(*gactiven_));
    LINALG::Export(*incrjump_,*activejump);
    nmatrix_->Multiply(false,*activejump,*gtest);
    
    RCP<Epetra_Map> gsmdofs = LINALG::MergeMap(gsdofrowmap_,gmdofrowmap_,false);
    RCP<Epetra_Vector> disism = rcp(new Epetra_Vector(*gsmdofs));
    LINALG::Export(*disi,*disism);
    if (fulllin)
    {
      smatrix_->Multiply(false,*disism,*gtest2);
      gtest->Update(1.0,*gtest2,1.0);
    }
    cout << *gtest << endl << *g_ << endl;
    
    RCP<Epetra_Vector> zactive = rcp(new Epetra_Vector(*gactivedofs_));
    RCP<Epetra_Vector> zerotest = rcp(new Epetra_Vector(*gactivet_));
    RCP<Epetra_Vector> zerotest2 = rcp(new Epetra_Vector(*gactivet_));
    LINALG::Export(*z_,*zactive);
    tmatrix_->Multiply(false,*zactive,*zerotest);
    if (fulllin)
    {
    pmatrix_->Multiply(false,*disis,*zerotest2);
    zerotest->Update(1.0,*zerotest2,1.0);
    }
    cout << *zerotest << endl;
  }
#endif // #ifdef DEBUG
  // CHECK OF CONTACT BOUNDARY CONDITIONS---------------------------------
  */
  
  return;
}

/*----------------------------------------------------------------------*
 |  Update active set and check for convergence (public)      popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::UpdateActiveSet(RCP<Epetra_Vector> disn)
{
  // get input parameter ctype
  string ctype   = scontact_.get<string>("contact type","none");
  
  // assume that active set has converged and check for opposite
  activesetconv_=true;
  
#ifdef CONTACTCHECKHUEEBER
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    RCP<LINALG::SparseMatrix> temp1 = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));;
    RCP<LINALG::SparseMatrix> temp2 = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));;
    interface_[i]->SetState("displacement",disn);
    interface_[i]->Evaluate();
    interface_[i]->AssembleDMG(*temp1,*temp2,*g_);
  }
#endif
  
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    if (i>0) dserror("ERROR: UpdateActiveSet: Double active node check needed for n interfaces!");
    
    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* cnode = static_cast<CNode*>(node);
      
      // get weighting factor from nodal D-map
      double wii;
      if ((int)((cnode->GetD()).size())==0) wii = 0.0;
      else wii = (cnode->GetD()[0])[cnode->Dofs()[0]];
      
      // compute incr. normal displacement and weighted gap
      double nincr = 0.0;
      for (int k=0;k<3;++k)
        nincr += wii * cnode->n()[k] * (*incrjump_)[incrjump_->Map().LID(2*gid)+k];
      double wgap = (*g_)[g_->Map().LID(gid)];
      
      if (cnode->n()[2] != 0.0) dserror("ERROR: UpdateActiveSet: Not yet implemented for 3D!");
      
      // compute normal part of Lagrange multiplier
      double nz = 0.0;
      double nzold = 0.0;
      for (int k=0;k<2;++k)
      {
        nz += cnode->n()[k] * cnode->lm()[k];
        nzold += cnode->n()[k] * cnode->lmold()[k];
      }
      
      // check nodes of inactive set *************************************
      // (by definition they fulfill the condition z_j = 0)
      // (thus we only have to check ncr.disp. jump and weighted gap)
      if (cnode->Active()==false)
      {
        // check for fulfilment of contact condition
        // if (abs(nz) > 1e-8)
        //  dserror("ERROR: UpdateActiveSet: Exact inactive node condition violated "
        //          "for node ID: %i ", cnode->Id());
        
        // check for penetration
        if (nincr-wgap > 0)
        {
          cnode->Active() = true;
          activesetconv_ = false;
        }
      }
      
      // check nodes of active set ***************************************
      // (by definition they fulfill the non-penetration condition)
      // (thus we only have to check for positive Lagrange multipliers)
      else
      {
        // check for fulfilment of contact condition
        // if (abs(nincr-wgap) > 1e-8)
        // dserror("ERROR: UpdateActiveSet: Exact active node condition violated "
        //         "for node ID: %i ", cnode->Id());
        
        // check for tensile contact forces
        //if (nz <= 0) // no averaging of Lagrange multipliers
        if (0.5*nz+0.5*nzold <= 0) // averaging of Lagrange multipliers
        {
          if (ctype!="meshtying")
          {
            cnode->Active() = false;
            activesetconv_ = false;
          }
          else
          {
            cnode->Active() = true;   // set all nodes active for mesh tying
            activesetconv_ = true;    // no active set loop for mesh tying
          }
        }
      }
    }
  }

  // broadcast convergence status among processors
  int convcheck = 0;
  int localcheck = activesetconv_;
  Comm().SumAll(&localcheck,&convcheck,1);
  
  // active set is only converged, if converged on all procs
  // if not, increase no. of active set steps too
  if (convcheck!=Comm().NumProc())
  {
    activesetconv_=false;
    ActiveSetSteps() += 1;
  }
  
  // update zig-zagging history (shift by one)
  if (zigzagtwo_!=null) zigzagthree_  = rcp(new Epetra_Map(*zigzagtwo_));
  if (zigzagone_!=null) zigzagtwo_    = rcp(new Epetra_Map(*zigzagone_));
  if (gactivenodes_!=null) zigzagone_ = rcp(new Epetra_Map(*gactivenodes_));
    
  // (re)setup active global Epetra_Maps
  gactivenodes_ = null;
  gactivedofs_ = null;
  gactiven_ = null;
  gactivet_ = null;

  // friction  
  // reset displacement jumps (slave dofs)
   jump_->Scale(0.0); 
   StoreNodalQuantities("jump");
      
  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i=0;i<(int)interface_.size();++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_,interface_[i]->ActiveNodes(),false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_,interface_[i]->ActiveDofs(),false);
    gactiven_ = LINALG::MergeMap(gactiven_,interface_[i]->ActiveNDofs(),false);
    gactivet_ = LINALG::MergeMap(gactivet_,interface_[i]->ActiveTDofs(),false);
  }
  
  // CHECK FOR ZIG-ZAGGING / JAMMING OF THE ACTIVE SET
  // *********************************************************************
  // A problem of the active set strategy which sometimes arises is known
  // from optimization literature as jamming or zig-zagging. This means
  // that within a load/time-step the algorithm can have more than one
  // solution due to the fact that the active set is not unique. Hence the
  // algorithm jumps between the solutions of the active set. The non-
  // uniquenesss results either from hoghly curved contact surfaces or
  // from the FE discretization, Thus the uniqueness of the closest-point-
  // projection cannot be guaranteed.
  // *********************************************************************
  // To overcome this problem we monitor the development of the active
  // set scheme in our contact algorithms. We can identify zig-zagging by
  // comparing the current active set with the active set of the second-
  // and third-last iteration. If an identity occurs, we consider the
  // active set strategy as converged instantly, accepting the current
  // version of the active set and proceeding with the next time/load step.
  // This very simple approach helps stabilizing the contact algorithm!
  // *********************************************************************
  if (ActiveSetSteps()>2)
  {
    if (zigzagtwo_!=null)
    {
      if (zigzagtwo_->SameAs(*gactivenodes_))
      {
        // set active set converged
        activesetconv_=true;
        
        // output to screen
        if (Comm().MyPID()==0)
          cout << "DETECTED 1-2 ZIG-ZAGGING OF ACTIVE SET................." << endl;
      }
    }
    
    if (zigzagthree_!=null)
    {
      if (zigzagthree_->SameAs(*gactivenodes_))
      {
        // set active set converged
        activesetconv_=true;
        
        // output to screen
        if (Comm().MyPID()==0)
          cout << "DETECTED 1-2-3 ZIG-ZAGGING OF ACTIVE SET................" << endl;
      }
    }
  }
  
  // reset zig-zagging history
  if (activesetconv_==true)
  {
    zigzagone_  = null;
    zigzagtwo_  = null;
    zigzagthree_= null;
  }
  
  // output of active set status to screen
  if (Comm().MyPID()==0 && activesetconv_==false)
    cout << "ACTIVE SET ITERATION " << ActiveSetSteps()-1
         << " NOT CONVERGED - REPEAT TIME STEP................." << endl;
  else if (Comm().MyPID()==0 && activesetconv_==true)
    cout << "ACTIVE SET CONVERGED IN " << ActiveSetSteps()
         << " STEP(S)................." << endl;
    
  // update flag for global contact status
  if (gactivenodes_->NumGlobalElements())
    IsInContact()=true;
  /*
#ifdef DEBUG
  // visualization with gmsh
  if (activesetconv_)
    for (int i=0;i<(int)interface_.size();++i)
      interface_[i]->VisualizeGmsh(interface_[i]->CSegs());
#endif // #ifdef DEBUG
  */
  return;
}

/*----------------------------------------------------------------------*
 |  Update active set and check for convergence (public)      popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::UpdateActiveSetSemiSmooth(RCP<Epetra_Vector> disn)
{
  // FIXME: Here we do not consider zig-zagging yet!
  
  // get input parameter ctype
  string ctype   = scontact_.get<string>("contact type","none");
  
  // read weighting factor cn
  // (this is necessary in semi-smooth Newton case, as the search for the
  // active set is now part of the Newton iteration. Thus, we do not know
  // the active / inactive status in advance and we can have a state in
  // which both the condition znormal = 0 and wgap = 0 are violated. Here
  // we have to weigh the two violations via cn!
  double cn = scontact_.get<double>("semismooth cn",0.0);
        
  // assume that active set has converged and check for opposite
  activesetconv_=true;
  
#ifdef CONTACTCHECKHUEEBER
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    RCP<LINALG::SparseMatrix> temp1 = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));;
    RCP<LINALG::SparseMatrix> temp2 = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));;
    interface_[i]->SetState("displacement",disn);
    interface_[i]->Evaluate();
    interface_[i]->AssembleDMG(*temp1,*temp2,*g_);
  }
#endif
  
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    if (i>0) dserror("ERROR: UpdateActiveSet: Double active node check needed for n interfaces!");
    
    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* cnode = static_cast<CNode*>(node);
      
      // get weighting factor from nodal D-map
      double wii;
      if ((int)((cnode->GetD()).size())==0) wii = 0.0;
      else wii = (cnode->GetD()[0])[cnode->Dofs()[0]];
      
      // compute incr. normal displacement and weighted gap
      double nincr = 0.0;
      for (int k=0;k<3;++k)
        nincr += wii * cnode->n()[k] * (*incrjump_)[incrjump_->Map().LID(2*gid)+k];
      double wgap = (*g_)[g_->Map().LID(gid)];
      
      if (cnode->n()[2] != 0.0) dserror("ERROR: UpdateActiveSet: Not yet implemented for 3D!");
      
      // compute normal part of Lagrange multiplier
      double nz = 0.0;
      double nzold = 0.0;
      for (int k=0;k<2;++k)
      {
        nz += cnode->n()[k] * cnode->lm()[k];
        nzold += cnode->n()[k] * cnode->lmold()[k];
      }
      
      // check nodes of inactive set *************************************
      if (cnode->Active()==false)
      {
        // check for penetration and/or tensile contact forces
        if (nz + cn*(nincr-wgap) > 0)
        {
          cnode->Active() = true;
          activesetconv_ = false;
        }
      }
      
      // check nodes of active set ***************************************
      else
      {
        // check for tensile contact forces and/or penetration
        //if (nz + cn*(nincr-wgap) <= 0) // no averaging of Lagrange multipliers
        if ((0.5*nz+0.5*nzold) + cn*(nincr-wgap) <= 0) // averaging of Lagrange multipliers
        {
          if (ctype!="meshtying")
          {
            cnode->Active() = false;
            activesetconv_ = false;
          }
          else
          {
            cnode->Active() = true;   // set all nodes active for mesh tying
            activesetconv_ = true;    // no active set loop for mesh tying
          }
        } 
      } 
    }
  }

  // broadcast convergence status among processors
  int convcheck = 0;
  int localcheck = activesetconv_;
  Comm().SumAll(&localcheck,&convcheck,1);
  
  // active set is only converged, if converged on all procs
  // if not, increase no. of active set steps too
  if (convcheck!=Comm().NumProc())
  {
    activesetconv_=false;
    ActiveSetSteps() += 1;
  }
  
  // (re)setup active global Epetra_Maps
  gactivenodes_ = null;
  gactivedofs_ = null;
  gactiven_ = null;
  gactivet_ = null;
  
  // friction  
  // reset displacement jumps (slave dofs)
  jump_->Scale(0.0); 
  StoreNodalQuantities("jump");
   
  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i=0;i<(int)interface_.size();++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_,interface_[i]->ActiveNodes(),false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_,interface_[i]->ActiveDofs(),false);
    gactiven_ = LINALG::MergeMap(gactiven_,interface_[i]->ActiveNDofs(),false);
    gactivet_ = LINALG::MergeMap(gactivet_,interface_[i]->ActiveTDofs(),false);
  }
  
  // output of active set status to screen
  if (Comm().MyPID()==0 && activesetconv_==false)
    cout << "ACTIVE SET HAS CHANGED... CHANGE No. " << ActiveSetSteps()-1 << endl;
  
  // update flag for global contact status
  if (gactivenodes_->NumGlobalElements())
    IsInContact()=true;
  /*
#ifdef DEBUG
  // visualization with gmsh
  if (activesetconv_)
    for (int i=0;i<(int)interface_.size();++i)
      interface_[i]->VisualizeGmsh(interface_[i]->CSegs());
#endif // #ifdef DEBUG
  */
  return;
}
/*----------------------------------------------------------------------*
 |  Compute contact forces (public)                           popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::ContactForces(RCP<Epetra_Vector> fresm)
{
  // Note that we ALWAYS use a TR-like approach to compute the contact
  // forces. This means we never explicitly compute fc at the generalized
  // mid-point n+1-alphaf, but use a linear combination of the old end-
  // point n and the new end-point n+1 instead:
  // F_{c;n+1-alpha_f} := (1-alphaf) * F_{c;n+1} +  alpha_f * F_{c;n}
  
  // FIXME: fresm is only here for debugging purposes!
  // compute two subvectors of fc each via Lagrange multipliers z_n+1, z_n
  RCP<Epetra_Vector> fcslavetemp  = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  RCP<Epetra_Vector> fcmastertemp = rcp(new Epetra_Vector(mmatrix_->DomainMap()));
  RCP<Epetra_Vector> fcslavetempend  = rcp(new Epetra_Vector(dold_->RowMap()));
  RCP<Epetra_Vector> fcmastertempend = rcp(new Epetra_Vector(mold_->DomainMap()));
  dmatrix_->Multiply(false,*z_,*fcslavetemp);
  mmatrix_->Multiply(true,*z_,*fcmastertemp);
  dold_->Multiply(false,*zold_,*fcslavetempend);
  mold_->Multiply(true,*zold_,*fcmastertempend);
  
  // export the contact forces to full dof layout
  RCP<Epetra_Vector> fcslave  = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  RCP<Epetra_Vector> fcmaster = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  RCP<Epetra_Vector> fcslaveend  = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  RCP<Epetra_Vector> fcmasterend = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  LINALG::Export(*fcslavetemp,*fcslave);
  LINALG::Export(*fcmastertemp,*fcmaster);
  LINALG::Export(*fcslavetempend,*fcslaveend);
  LINALG::Export(*fcmastertempend,*fcmasterend);
  
  // build total contact force vector (TR-like!!!)
  fc_=fcslave;
  fc_->Update(-(1.0-alphaf_),*fcmaster,1.0-alphaf_);
  fc_->Update(alphaf_,*fcslaveend,1.0);
  fc_->Update(-alphaf_,*fcmasterend,1.0);
  
  /*
  // CHECK OF CONTACT FORCE EQUILIBRIUM ----------------------------------
#ifdef DEBUG
  RCP<Epetra_Vector> fresmslave  = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  RCP<Epetra_Vector> fresmmaster = rcp(new Epetra_Vector(mmatrix_->DomainMap()));
  LINALG::Export(*fresm,*fresmslave);
  LINALG::Export(*fresm,*fresmmaster);
  
  vector<double> gfcs(3);
  vector<double> ggfcs(3);
  vector<double> gfcm(3);
  vector<double> ggfcm(3);
  int dimcheck = (gsdofrowmap_->NumGlobalElements())/(gsnoderowmap_->NumGlobalElements());
  if (dimcheck!=2 && dimcheck!=3) dserror("ERROR: ContactForces: Debugging for 3D not implemented yet");
  
  for (int i=0;i<fcslavetemp->MyLength();++i)
  {
    if ((i+dimcheck)%dimcheck == 0) gfcs[0]+=(*fcslavetemp)[i];
    else if ((i+dimcheck)%dimcheck == 1) gfcs[1]+=(*fcslavetemp)[i];
    else if ((i+dimcheck)%dimcheck == 2) gfcs[2]+=(*fcslavetemp)[i];
    else dserror("ERROR: Contact Forces: Dim. error in debugging part!");
  }
  
  for (int i=0;i<fcmastertemp->MyLength();++i)
  {
    if ((i+dimcheck)%dimcheck == 0) gfcm[0]-=(*fcmastertemp)[i];
    else if ((i+dimcheck)%dimcheck == 1) gfcm[1]-=(*fcmastertemp)[i];
    else if ((i+dimcheck)%dimcheck == 2) gfcm[2]-=(*fcmastertemp)[i];
    else dserror("ERROR: Contact Forces: Dim. error in debugging part!");
  }
  
  for (int i=0;i<3;++i)
  {
    Comm().SumAll(&gfcs[i],&ggfcs[i],1);
    Comm().SumAll(&gfcm[i],&ggfcm[i],1);
  }
  
  double slavenorm = 0.0;
  fcslavetemp->Norm2(&slavenorm);
  double slavenormend = 0.0;
  fcslavetempend->Norm2(&slavenormend);
  double fresmslavenorm = 0.0;
  fresmslave->Norm2(&fresmslavenorm);
  if (Comm().MyPID()==0)
  {
    cout << "Slave Contact Force Norm (n+1):  " << slavenorm << endl;
    cout << "Slave Contact Force Norm (n):  " << slavenormend << endl;
    cout << "Slave Residual Force Norm: " << fresmslavenorm << endl;
    cout << "Slave Contact Force Vector: " << ggfcs[0] << " " << ggfcs[1] << " " << ggfcs[2] << endl;
  }
  double masternorm = 0.0;
  fcmastertemp->Norm2(&masternorm);
  double masternormend = 0.0;
  fcmastertempend->Norm2(&masternormend);
  double fresmmasternorm = 0.0;
  fresmmaster->Norm2(&fresmmasternorm);
  if (Comm().MyPID()==0)
  {
    cout << "Master Contact Force Norm (n+1): " << masternorm << endl;
    cout << "Master Contact Force Norm (n): " << masternormend << endl;
    cout << "Master Residual Force Norm " << fresmmasternorm << endl;
    cout << "Master Contact Force Vector: " << ggfcm[0] << " " << ggfcm[1] << " " << ggfcm[2] << endl;
  }
#endif // #ifdef DEBUG
  // CHECK OF CONTACT FORCE EQUILIBRIUM ----------------------------------
  */
  
  return;
}

/*----------------------------------------------------------------------*
 |  Store Lagrange mulitpliers and displacent jumps into CNode popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::StoreNodalQuantities(const string& state)
{
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    if (i>0) dserror("ERROR: StoreNodalQuantities: Double active node check needed for n interfaces!");
    
    // export global LM vector or Jump vector to current interface slave dof row map
    RCP<Epetra_Vector> vectorglobal = null;
    if (state=="lmcurrent") vectorglobal = LagrMult();
    else if (state=="lmold") vectorglobal = LagrMultOld();
    else if (state=="jump") vectorglobal = Jump(); 
    else dserror("ERROR: StoreNodalQuantities: Unknown state string variable!");
    RCP<Epetra_Map> sdofrowmap = interface_[i]->SlaveRowDofs();
    RCP<Epetra_Vector> vectorinterface = rcp(new Epetra_Vector(*sdofrowmap));
    LINALG::Export(*vectorglobal,*vectorinterface);
    
    // loop over all slave row nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* cnode = static_cast<CNode*>(node);
      
      if (cnode->NumDof() > 2) dserror("ERROR: StoreNodalQuantities: Not yet implemented for 3D!");
      
      // extract this node's LM from zcurr
      for (int k=0;k<2;++k)
      {
        if (state=="lmcurrent")
          cnode->lm()[k] = (*vectorinterface)[vectorinterface->Map().LID(2*gid)+k];
        else if (state=="lmold")
          cnode->lmold()[k] = (*vectorinterface)[vectorinterface->Map().LID(2*gid)+k];
        else if (state=="jump")
        	cnode->jump()[k] = (*vectorinterface)[vectorinterface->Map().LID(2*gid)+k];
        else
          dserror("ERROR: StoreNodalQuantities: Unknown state string variable!");
      }
    }
  }
    
  return;
}

/*----------------------------------------------------------------------*
 |  Store D and M last coverged step <-> current step        popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::StoreDM(const string& state)
{
  //store Dold and Mold matrix in D and M
  if (state=="current")
  {
    dmatrix_ = dold_;
    mmatrix_ = mold_;
  } 
    
  // store D and M matrix in Dold and Mold
  else if (state=="old")
  {
    dold_ = dmatrix_;
    mold_ = mmatrix_;
  }
  
  // unknown conversion
  else
    dserror("ERROR: StoreDM: Unknown conversion requested!");
  
  return;
}

/*----------------------------------------------------------------------*
 |  Print current active set to screen                        popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::PrintActiveSet()
{
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    if (i>0) dserror("ERROR: UpdateActiveSet: Double active node check needed for n interfaces!");
    
    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* cnode = static_cast<CNode*>(node);
      
      // get weighting factor from nodal D-map
      double wii;
      if ((int)((cnode->GetD()).size())==0) wii = 0.0;
      else wii = (cnode->GetD()[0])[cnode->Dofs()[0]];
      
      // compute incr. normal displacement and weighted gap
      double nincr = 0.0;
      for (int k=0;k<3;++k)
        nincr += wii * cnode->n()[k] * (*incrjump_)[incrjump_->Map().LID(2*gid)+k];
      double wgap = (*g_)[g_->Map().LID(gid)];
      
      if (cnode->n()[2] != 0.0) dserror("ERROR: UpdateActiveSet: Not yet implemented for 3D!");
      
      // compute normal part of Lagrange multiplier
      double nz = 0.0;
      double nzold = 0.0;
      for (int k=0;k<2;++k)
      {
        nz += cnode->n()[k] * cnode->lm()[k];
        nzold += cnode->n()[k] * cnode->lmold()[k];
      }
      
      // print nodes of inactive set *************************************
      if (cnode->Active()==false)
      {
        cout << "INACTIVE: " << gid << " " << wgap << " " << nz << endl;
      }
      
      // print nodes of active set ***************************************
      else
      {
        cout << "ACTIVE:   " << gid << " " << nz << " " << nzold << " "
             << 0.5*nz+0.5*nzold << " " << wgap << endl;
      }
    }
  }

  return;
}

#endif  // #ifdef CCADISCRET
