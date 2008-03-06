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

#include "drt_contact_manager.H"
#include "drt_cnode.H"
#include "drt_celement.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_solver.H"
#include "contactdefines.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::Manager::Manager(DRT::Discretization& discret) :
discret_(discret),
activesetconv_(false)
{
#if 0 // OLD VERSION
  if (!Discret().Filled()) dserror("Discretization is not fillcomplete");

  // let's check for contact boundary conditions in discret
  // and detect pairs of matching conditions
  // for each pair, create a contact interface and store it
  vector<DRT::Condition*> contactconditions(0);
  Discret().GetCondition("Contact",contactconditions);
  if ((int)contactconditions.size()<=1)  dserror("Not enough contact conditions in discretization");
  if ((int)contactconditions.size() %2) dserror("Odd number of contact conditions is impossible");

  // find all pairs of matching contact conditions
  // there are num conditions / 2 pairs
  vector<int> foundpairs(contactconditions.size()/2);
  int numpairsfound = 0;
  for (int i=0; i<(int)contactconditions.size(); ++i)
  {
    DRT::Condition* cond1 = contactconditions[i];
    DRT::Condition* cond2 = NULL;
    const vector<int>* pair1v = cond1->Get<vector<int> >("contact id");
    if (!pair1v) dserror("Contact Conditions does not have value 'contact id'");
    int pairid1 = (*pair1v)[0];
    bool foundit = false;
    const vector<int>* pair2v = NULL;
    int pairid2 = -1;
    for (int j=0; j<(int)contactconditions.size(); ++j)
    {
      if (j==i) continue; // do not detect contactconditions[i] again
      cond2 = contactconditions[j];
      pair2v = cond2->Get<vector<int> >("contact id");
      if (!pair2v) dserror("Contact Conditions does not have value 'contact id'");
      pairid2 = (*pair2v)[0];
      if (pairid1 != pairid2) continue; // not a pair
      foundit = true; // found a pair
      break;
    }

    // now we should have found a pair cond1/cond2
    if (!foundit) dserror("Cannot find matching contact condition for id %d",pairid1);

    // see whether we found this pair before
    bool foundbefore = false;
    for (int j=0; j<numpairsfound; ++j)
      if (pairid1 == foundpairs[j])
      {
        foundbefore = true;
        break;
      }

    // if we have processed this pair before, do nothing
    if (foundbefore) continue;

    // we have not found this pair pairid1/pairid2 before, process it
    foundpairs[numpairsfound] = pairid1;
    ++numpairsfound;

    // create an empty interface and store it in this Manager
    interface_.push_back(rcp(new CONTACT::Interface(pairid1,Comm())));

    // get it again
    RCP<CONTACT::Interface> interface = interface_[(int)interface_.size()-1];

    // find out which side is Master and Slave
    const string* side1 = cond1->Get<string>("Side");
    const string* side2 = cond2->Get<string>("Side");
    if (!side1 || !side2) dserror("Data for Master/Slave side missing in condition");
    if (*side1 == *side2) dserror("2 Slave sides or 2 Master sides not allowed");
    bool is1slave = false;
    if (*side1 == "Slave") is1slave = true;

    // note that the nodal ids are unique because they come from
    // one global problem discretization conatining all nodes of the
    // contact interface
    // We rely on this fact, therefore it is not possible to
    // do contact between two distinct discretizations here
    
    //--------------------------------------------- process side 1 nodes
    // get all nodes and add them
    const vector<int>* nodeids1 = cond1->Nodes();
    if (!nodeids1) dserror("Condition does not have Node Ids");
    for (int j=0; j<(int)(*nodeids1).size(); ++j)
    {
      int gid = (*nodeids1)[j];
      // do only nodes that I have in my discretization
      if (!Discret().NodeColMap()->MyGID(gid)) continue;
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("Cannot find node with gid %",gid);
      RCP<CONTACT::CNode> cnode = rcp(new CONTACT::CNode(node->Id(),node->X(),
                                                         node->Owner(),
                                                         Discret().NumDof(node),
                                                         Discret().Dof(node),is1slave));
      interface->AddCNode(cnode);
    } // for (int j=0; j<(int)(*nodeids1).size(); ++j)

    //----------------------------------------------- process side 2 nodes
    // get all nodes and add them
    const vector<int>* nodeids2 = cond2->Nodes();
    if (!nodeids2) dserror("Condition does not have Node Ids");
    for (int j=0; j<(int)(*nodeids2).size(); ++j)
    {
      int gid = (*nodeids2)[j];
      if (!Discret().NodeColMap()->MyGID(gid)) continue;
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("Cannot find node with gid %",gid);
      RCP<CONTACT::CNode> cnode = rcp(new CONTACT::CNode(node->Id(),node->X(),
                                                         node->Owner(),
                                                         Discret().NumDof(node),
                                                         Discret().Dof(node),
                                                         !is1slave));
      interface->AddCNode(cnode);
    } // for (int j=0; j<(int)(*nodeids2).size(); ++j)

    //-------------------------------------------- process elements
    // get elements from condition cond1/cond2
    map<int,RCP<DRT::Element> >& ele1 = cond1->Geometry();
    map<int,RCP<DRT::Element> >& ele2 = cond2->Geometry();
    // elements in a boundary condition have a unique id
    // but ids are not unique among 2 distinct conditions
    // due to the way elements in conditions are build.
    // We therefore have to give the second set of elements different ids
    // ids do not have to be continous, we just add a large enough number
    // gsize to all elements of cond2 so they are different from those in cond1.
    // note that elements in ele1/ele2 already are in column (overlapping) map
    int lsize = (int)ele1.size();
    int gsize = 0;
    Comm().SumAll(&lsize,&gsize,1);

    //---------------------------------------- process elements in ele1
    map<int,RCP<DRT::Element> >::iterator fool;
    for (fool=ele1.begin(); fool != ele1.end(); ++fool)
    {
      RCP<DRT::Element> ele = fool->second;
      RCP<CONTACT::CElement> cele = rcp(new CONTACT::CElement(ele->Id(),
                                                              DRT::Element::element_contact,
                                                              ele->Owner(),
                                                              ele->Shape(),
                                                              ele->NumNode(),
                                                              ele->NodeIds(),
                                                              is1slave));
      interface->AddCElement(cele);
    } // for (fool=ele1.start(); fool != ele1.end(); ++fool)

    //----------------------------------------- process element in ele2
    // note the change in the id of the elements
    for (fool=ele2.begin(); fool != ele2.end(); ++fool)
    {
      RCP<DRT::Element> ele = fool->second;
      RCP<CONTACT::CElement> cele = rcp(new CONTACT::CElement(ele->Id()+gsize,
                                                              DRT::Element::element_contact,
                                                              ele->Owner(),
                                                              ele->Shape(),
                                                              ele->NumNode(),
                                                              ele->NodeIds(),
                                                              !is1slave));
      interface->AddCElement(cele);
    } // for (fool=ele1.start(); fool != ele1.end(); ++fool)


    //-------------------- finalize the contact interface construction
    interface->FillComplete();

  } // for (int i=0; i<(int)contactconditions.size(); ++i)
#endif // #if 1 (OLD VERSION)
  
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
    
    //--------------------------------------------- process nodes
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

    //-------------------------------------------- process elements
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
      // enough number ggsize to all elements of cond2, cnod3,... so they are
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
  
  // setup global active node/dof maps
  gactivenodes_ = rcp(new Epetra_Map(0,0,Comm()));
  gactivedofs_ = rcp(new Epetra_Map(0,0,Comm()));
  gactiven_ = rcp(new Epetra_Map(0,0,Comm()));
  gactivet_ = rcp(new Epetra_Map(0,0,Comm()));
  
  // setup Lagrange muliplier vectors
  zold_       = rcp(new Epetra_Vector(*gsdofrowmap_));
  z_          = rcp(new Epetra_Vector(*gsdofrowmap_));
  		
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
 |  initialize contact for next Newton step (public)          popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::Initialize()
{
  // initialize / reset interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->Initialize();
  }

  // (re)setup global Mortar LINALG::SparseMatrices and Epetra_Vectors
  dmatrix_    = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
  mmatrix_    = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
  mhatmatrix_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
  g_          = LINALG::CreateVector(*gsnoderowmap_,true);
    
  // (re)setup global normal and tangent matrices
  nmatrix_ = rcp(new LINALG::SparseMatrix(*gactiven_,3));
  tmatrix_ = rcp(new LINALG::SparseMatrix(*gactivet_,3));
  
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
 |  evaluate contact (public)                                 popp 11/07|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::Evaluate(RCP<LINALG::SparseMatrix>& kteff,
                                RCP<Epetra_Vector>& feff)
{  
  /**********************************************************************/
  /* evaluate interfaces                                                */
  /* (nodal normals, projections, Mortar integration, Mortar assembly)  */
  /**********************************************************************/
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->Evaluate();
    interface_[i]->AssembleDMG(*dmatrix_,*mmatrix_,*g_);
  }
  
  // FillComplete() global Mortar matrices
  dmatrix_->Complete();
  mmatrix_->Complete(*gmdofrowmap_,*gsdofrowmap_);
  
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
  for (int i=0; i<(int)interface_.size(); ++i)
    interface_[i]->AssembleNT(*nmatrix_,*tmatrix_);
    
  // FillComplete() global matrices N and T
  nmatrix_->Complete(*gactivedofs_,*gactiven_);
  tmatrix_->Complete(*gactivedofs_,*gactivet_);
  
  /**********************************************************************/
  /* Multiply Mortar matrices: m^ = inv(d) * m                          */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> invd = rcp(new LINALG::SparseMatrix(*dmatrix_));
  RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_,true);
  int err = 0;
  
  // extract diagonal of invd into diag
  invd->ExtractDiagonalCopy(*diag);
  
  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");
  
  // re-insert inverted diagonal into invd
  err = invd->ReplaceDiagonalValues(*diag);
  if (err>0) dserror("ERROR: ReplaceDiagonalValues: Missing diagonal entry!");
  
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
  RCP<Epetra_Map> gsmdofs = LINALG::MergeMap(gsdofrowmap_,gmdofrowmap_);
  
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
  
  // store some stuff for static condensation of LM
  ksn_  = ksn;
  ksm_  = ksm;
  kss_  = kss;
  invd_ = invd;
  
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
  kteff = kteffnew;
  feff = feffnew;
  //LINALG::PrintSparsityToPostscript(*(kteff->EpetraMatrix()));
  
  /**********************************************************************/
  /* Update Lagrange multipliers                                        */
  /**********************************************************************/
  //invd->Multiply(false,*fsmod,*z_);
  z_->Update(1.0,*fsmod,0.0);

  return;
}

/*----------------------------------------------------------------------*
 |  Transform displacement increment vector (public)          popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::RecoverDisp(RCP<Epetra_Vector>& disi)
{  
  // extract incremental jump from disi (for active set)
  incrjump_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*disi,*incrjump_);
  
  // extract master displacements from disi
  RCP<Epetra_Vector> disim = rcp(new Epetra_Vector(*gmdofrowmap_));
  LINALG::Export(*disi,*disim);
  
  // recover slave displacement increments
  RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(mhatmatrix_->RowMap()));
  mhatmatrix_->Multiply(false,*disim,*mod);
  
  RCP<Epetra_Vector> modexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  LINALG::Export(*mod,*modexp);
  disi->Update(1.0,*modexp,1.0);
  
  // update Lagrange multipliers
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
  RCP<Epetra_Vector> zcopy = rcp(new Epetra_Vector(*z_));
  invd_->Multiply(false,*zcopy,*z_);
  
  /*
  // CHECK OF CONTACT COUNDARY CONDITIONS---------------------------------
#ifdef DEBUG
  //debugging (check for z_i = 0)
  RCP<Epetra_Map> gidofs = LINALG::SplitMap(*gsdofrowmap_,*gactivedofs_);
  RCP<Epetra_Vector> zinactive = rcp(new Epetra_Vector(*gidofs));
  LINALG::Export(*z_,*zinactive);
  cout << *zinactive << endl;
  
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
  }
#endif // #ifdef DEBUG
  // CHECK OF CONTACT BOUNDARY CONDITIONS---------------------------------
  */
  
  return;
}

/*----------------------------------------------------------------------*
 |  Update active set and check for convergence (public)      popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::UpdateActiveSet()
{
  // assume that active set has converged and check for opposite
  activesetconv_=true;
  
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
      
      // check nodes of inactive set *************************************
      // (by definition they fulfill the condition z_j = 0)
      // (thus we only have to check ncr.disp. jump and weighted gap)
      if (cnode->Active()==false)
      {
        double wii = (cnode->GetD()[0])[cnode->Dofs()[0]];
        double nincr = 0.0;
        for (int k=0;k<3;++k)
          nincr += wii * cnode->n()[k] * (*incrjump_)[incrjump_->Map().LID(2*gid)+k];
        double wgap = (*g_)[g_->Map().LID(gid)];
        
        if (cnode->n()[2] != 0.0) dserror("ERROR: UpdateActiveSet: Not yet implemented for 3D!");
        
        // check for penetration
        if (nincr>wgap)
        {
          cnode->Active() = true;
          activesetconv_ = false;
        }
        
        //cout << "INACTIVE: " << i << " " << j << " " << gid << " "
        //     << nincr << " " << (*g_)[g_->Map().LID(gid)] << " "
        //     << cnode->HasProj() << endl;
      }
      
      // check nodes of active set ***************************************
      // (by definition they fulfill the non-penetration condition)
      // (thus we only have to check for positive Lagrange multipliers)
      else
      {
        double nz = 0.0;
        double nzold = 0.0;
        for (int k=0;k<3;++k)
        {
          nz += cnode->n()[k] * (*z_)[z_->Map().LID(2*gid)+k];
          nzold += cnode->n()[k] * (*zold_)[zold_->Map().LID(2*gid)+k];
        }
        
        // check for tensile contact forces
        //if (nz<0) // no averaging of Lagrange multipliers
        if ((0.5*nz+0.5*nzold)<0) // averaging of Lagrange multipliers
        {
          cnode->Active() = false;
          activesetconv_ = false;
        }
        
        cout << "ACTIVE: " << i << " " << j << " " << gid << " "
             << nz << " " << nzold << " " << 0.5*nz+0.5*nzold
             << " " << cnode->Getg() << endl;  
      }
      
    }
  }

  // broadcast convergence status among processors
  int convcheck = 0;
  int localcheck = activesetconv_;
  Comm().SumAll(&localcheck,&convcheck,1);
  
  // active set is only converged, if converged on all procs
  if (convcheck!=Comm().NumProc())
    activesetconv_=false;
  
  // output of active set status to screen
  if (Comm().MyPID()==0 && activesetconv_==false)
    cout << "ACTIVE SET NOT CONVERGED - REPEAT TIME STEP................." << endl;
  else if (Comm().MyPID()==0 && activesetconv_==true)
    cout << "ACTIVE SET CONVERGED................." << endl;
  
  // (re)setup active global Epetra_Maps
  gactivenodes_ = null;
  gactivedofs_ = null;
  gactiven_ = null;
  gactivet_ = null;
  
  // update active sets of all interfaces
  for (int i=0;i<(int)interface_.size();++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_,interface_[i]->ActiveNodes());
    gactivedofs_ = LINALG::MergeMap(gactivedofs_,interface_[i]->ActiveDofs());
    gactiven_ = LINALG::MergeMap(gactiven_,interface_[i]->ActiveNDofs());;
    gactivet_ = LINALG::MergeMap(gactivet_,interface_[i]->ActiveTDofs());;
  }
  
  // create empty maps, if active set = null
  if (gactivenodes_==null)
  {
    gactivenodes_ = rcp(new Epetra_Map(0,0,Comm()));
    gactivedofs_ = rcp(new Epetra_Map(0,0,Comm()));
    gactiven_ = rcp(new Epetra_Map(0,0,Comm()));
    gactivet_ = rcp(new Epetra_Map(0,0,Comm()));
  }
  
  // store Lagrange multipliers zold_ if active set converged
  if (activesetconv_) zold_=rcp(new Epetra_Vector(*z_));
  
#ifdef DEBUG
  // visualization with gmsh
  if (activesetconv_)
    for (int i=0;i<(int)interface_.size();++i)
      interface_[i]->VisualizeGmsh(interface_[i]->CSegs());
#endif // #ifdef DEBUG
  
  return;
}

/*----------------------------------------------------------------------*
 |  Compute contact forces (public)                           popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::ContactForces(RCP<Epetra_Vector>& fc,
                                     RCP<Epetra_Vector>& fresm)
{
  // FIXME: fresm is only here for debugging purposes!
  // compute two subvectors of fc via Lagrange multipliers z
  RCP<Epetra_Vector> fcslavetemp  = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  RCP<Epetra_Vector> fcmastertemp = rcp(new Epetra_Vector(mmatrix_->DomainMap()));
  dmatrix_->Multiply(false,*z_,*fcslavetemp);
  mmatrix_->Multiply(true,*z_,*fcmastertemp);
  
  // export the contact forces to full dof layout
  RCP<Epetra_Vector> fcslave  = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  RCP<Epetra_Vector> fcmaster = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
  LINALG::Export(*fcslavetemp,*fcslave);
  LINALG::Export(*fcmastertemp,*fcmaster);
  
  // build total contact force vector
  fc->Update(1.0,*fcslave,0.0);
  fc->Update(-1.0,*fcmaster,1.0);
  
  // store into member variable
  fc_=fc;
  
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
  double fresmslavenorm = 0.0;
  fresmslave->Norm2(&fresmslavenorm);
  if (Comm().MyPID()==0)
  {
    cout << "Slave Contact Force Norm:  " << slavenorm << endl;
    cout << "Slave Residual Force Norm: " << fresmslavenorm << endl;
    cout << "Slave Contact Force Vector: " << ggfcs[0] << " " << ggfcs[1] << " " << ggfcs[2] << endl;
  }
  double masternorm = 0.0;
  fcmastertemp->Norm2(&masternorm);
  double fresmmasternorm = 0.0;
  fresmmaster->Norm2(&fresmmasternorm);
  if (Comm().MyPID()==0)
  {
    cout << "Master Contact Force Norm: " << masternorm << endl;
    cout << "Master Residual Force Norm " << fresmmasternorm << endl;
    cout << "Master Contact Force Vector: " << ggfcm[0] << " " << ggfcm[1] << " " << ggfcm[2] << endl;
  }
#endif // #ifdef DEBUG
  // CHECK OF CONTACT FORCE EQUILIBRIUM ----------------------------------
  */
  
  return;
}

#endif  // #ifdef CCADISCRET
