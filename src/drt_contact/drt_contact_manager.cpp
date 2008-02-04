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

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::Manager::Manager(DRT::Discretization& discret) :
discret_(discret)
{
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

  // setup global Epetra_Maps
  // (this is done by looping over all interfaces and merging)
  for (int i=0;i<(int)interface_.size();++i)
  {
  	gsnoderowmap_ = LINALG::MergeMap(gsnoderowmap_,interface_[i]->SlaveRowNodes());
  	gsdofrowmap_ = LINALG::MergeMap(gsdofrowmap_,interface_[i]->SlaveRowDofs());
  	gmdofrowmap_ = LINALG::MergeMap(gmdofrowmap_,interface_[i]->MasterRowDofs());
  }
  
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
	
	// (re)setup global Epetra_CrsMatrices and Epetra_Vectors
	D_ = LINALG::CreateMatrix(*gsdofrowmap_,10);
	M_ = LINALG::CreateMatrix(*gsdofrowmap_,100);
	g_ = LINALG::CreateVector(*gsnoderowmap_,true);
		
	// update active global Epetra_Maps
	for (int i=0;i<(int)interface_.size();++i)
	{
		gactivenodes_ = LINALG::MergeMap(gactivenodes_,interface_[i]->ActiveNodes());
		gactivedofs_ = LINALG::MergeMap(gactivedofs_,interface_[i]->ActiveDofs());
	}
	
  return;
}

/*----------------------------------------------------------------------*
 |  set current deformation state (public)                    popp 11/07|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::SetState(const string& statename, const RCP<Epetra_Vector> vec)
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
void CONTACT::Manager::Evaluate(RCP<Epetra_CrsMatrix>& Kteff,
																RCP<Epetra_Vector>& feff)
{	
	/**********************************************************************/
	/* evaluate interfaces                                                */
	/* (nodal normals, projections, Mortar integration, Mortar assembly)  */
	/**********************************************************************/
	for (int i=0; i<(int)interface_.size(); ++i)
  {
	  interface_[i]->Evaluate();
	  interface_[i]->Assemble_DMG(*D_,*M_,*g_);
  }
	
	// FillComplete() global Mortar matrices
	LINALG::Complete(*D_);
	LINALG::Complete(*M_,*gmdofrowmap_,*gsdofrowmap_);
	
	/**********************************************************************/
	/* build global matrix N with normal vectors of active nodes          */
	/* and global matrix T with tangent vectors of active nodes           */
	/**********************************************************************/
	if (gactivenodes_!=null)
	{ 
		N_ = LINALG::CreateMatrix(*gactivenodes_,3);
	  T_ = LINALG::CreateMatrix(*gactivenodes_,3);
	  
	  for (int i=0; i<(int)interface_.size(); ++i)
	    interface_[i]->Assemble_NT(*N_,*T_);
	  
	  // FillComplete() global matrices N and T
	  LINALG::Complete(*N_,*gactivedofs_,*gactivenodes_);
	  LINALG::Complete(*T_,*gactivedofs_,*gactivenodes_);
	}
	
	/**********************************************************************/
	/* Multiply Mortar matrices: M^ = inv(D) * M                          */
	/**********************************************************************/
	RCP<Epetra_CrsMatrix> invD = rcp(new Epetra_CrsMatrix(*D_));
	RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_,true);
	int err = 0;
	
	invD->ExtractDiagonalCopy(*diag);
	err = diag->Reciprocal(*diag);
	if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");
	err = invD->ReplaceDiagonalValues(*diag);
	if (err>0) dserror("ERROR: ReplaceDiagonalValues: Missing diagonal entry!");
	
	RCP<Epetra_CrsMatrix> Mbar = LINALG::MatMatMult(*invD,false,*M_,false);
	
	/**********************************************************************/
	/* Split Kteff into 3x3 block matrix                                  */
	/**********************************************************************/
	// we want to split K into 3 groups s,m,n = 9 blocks
	RCP<Epetra_CrsMatrix> Kss, Ksm, Ksn, Kms, Kmm, Kmn, Kns, Knm, Knn;
	
	// temporarily we need the blocks K_sm_sm, K_sm_n, K_n_sm
	RCP<Epetra_CrsMatrix> Ksmsm, Ksmn, Knsm;
	
	// we also need the combined sm rowmap
	RCP<Epetra_Map> gsmdofs = LINALG::MergeMap(gsdofrowmap_,gmdofrowmap_);
	
	// we will get the n rowmap as a by-product
	RCP<Epetra_Map> gndofs;
	
	// some temporary RCPs
	RCP<Epetra_Map> tempmap = null;
	RCP<Epetra_CrsMatrix> tempmtx1 = null;
	RCP<Epetra_CrsMatrix> tempmtx2 = null;
	
	// split into slave/master part + structure part
	LINALG::SplitMatrix2x2(Kteff,gsmdofs,gndofs,Ksmsm,Ksmn,Knsm,Knn);
	
	// further splits into slave part + master part
	LINALG::SplitMatrix2x2(Ksmsm,gsdofrowmap_,gmdofrowmap_,Kss,Ksm,Kms,Kmm);
	LINALG::SplitMatrix2x2(Ksmn,gsdofrowmap_,gmdofrowmap_,gndofs,tempmap,Ksn,tempmtx1,Kmn,tempmtx2);
	LINALG::SplitMatrix2x2(Knsm,gndofs,tempmap,gsdofrowmap_,gmdofrowmap_,Kns,Knm,tempmtx1,tempmtx2);
	
	// output for checking everything
	if(Kteff->Comm().MyPID()==0)
	{
		cout << endl << "*****************" << endl;
		cout << "K:   " << Kteff->NumGlobalRows() << " x " << Kteff->NumGlobalCols() << endl;
		if (Kss!=null) cout << "Kss: " << Kss->NumGlobalRows() << " x " << Kss->NumGlobalCols() << endl;
		else cout << "Kss: null" << endl;
		if (Ksm!=null) cout << "Ksm: " << Ksm->NumGlobalRows() << " x " << Ksm->NumGlobalCols() << endl;
		else cout << "Ksm: null" << endl;
		if (Ksn!=null) cout << "Ksn: " << Ksn->NumGlobalRows() << " x " << Ksn->NumGlobalCols() << endl;
		else cout << "Ksn: null" << endl;
		if (Kms!=null) cout << "Kms: " << Kms->NumGlobalRows() << " x " << Kms->NumGlobalCols() << endl;
		else cout << "Kms: null" << endl;
		if (Kmm!=null) cout << "Kmm: " << Kmm->NumGlobalRows() << " x " << Kmm->NumGlobalCols() << endl;
		else cout << "Kmm: null" << endl;
		if (Kmn!=null) cout << "Kmn: " << Kmn->NumGlobalRows() << " x " << Kmn->NumGlobalCols() << endl;
		else cout << "Kmn: null" << endl;
		if (Kns!=null) cout << "Kns: " << Kns->NumGlobalRows() << " x " << Kns->NumGlobalCols() << endl;
		else cout << "Kns: null" << endl;
		if (Knm!=null) cout << "Knm: " << Knm->NumGlobalRows() << " x " << Knm->NumGlobalCols() << endl;
		else cout << "Knm: null" << endl;
		if (Knn!=null) cout << "Knn: " << Knn->NumGlobalRows() << " x " << Knn->NumGlobalCols() << endl;
		else cout << "Knn: null" << endl;
		cout << "*****************" << endl;
	}
	
	/**********************************************************************/
	/* Split feff into 3 subvectors                                       */
	/**********************************************************************/
	// we want to split f into 3 groups s.m,n
	Epetra_Vector* f1 = NULL;
	Epetra_Vector* f2 = NULL;
	Epetra_Vector* f3 = NULL;
	
	// temporarily we need the group sm
	Epetra_Vector* f12 = NULL;
	
	// do the vector splitting smn -> sm+n -> s+m+n
	LINALG::SplitVector(*feff,*gsmdofs,f12,*gndofs,f3);
	LINALG::SplitVector(*f12,*gsdofrowmap_,f1,*gmdofrowmap_,f2);
	
	// wrap subvectors into RCPs
	RCP<Epetra_Vector> fs = rcp(f1);
	RCP<Epetra_Vector> fm = rcp(f2);
	RCP<Epetra_Vector> fn = rcp(f3);
	
	// output for checking everything
	if(feff->Comm().MyPID()==0)
	{
		cout << endl << "**********" << endl;
		cout << "f:  " << feff->GlobalLength() << endl;
		cout << "fs: " << fs->GlobalLength() << endl;
		cout << "fm: " << fm->GlobalLength() << endl;
		cout << "fn: " << fn->GlobalLength() << endl;
		cout << "**********" << endl;
	}
	
	/**********************************************************************/
	/* Apply basis transformation to K                                    */
	/**********************************************************************/
	// define temporaray RCP mod
	RCP<Epetra_CrsMatrix> mod;
	
	// Kss: nothing to do
	RCP<Epetra_CrsMatrix> Kss_mod = Kss;
	
	// Ksm: add Kss*T(Mbar)
	RCP<Epetra_CrsMatrix> Ksm_mod = LINALG::CreateMatrix(Ksm->RowMap(),100);
	LINALG::Add(*Ksm,false,1.0,*Ksm_mod,1.0);
	mod = LINALG::MatMatMult(*Kss,false,*Mbar,false);
	LINALG::Add(*mod,false,1.0,*Ksm_mod,1.0);
	LINALG::Complete(*Ksm_mod,Ksm->DomainMap(),Ksm->RowMap());

	// Ksn: nothing to do
	RCP<Epetra_CrsMatrix> Ksn_mod = Ksn;
	
	// Kms: add T(Mbar)*Kss
	RCP<Epetra_CrsMatrix> Kms_mod = LINALG::CreateMatrix(Kms->RowMap(),100);
	LINALG::Add(*Kms,false,1.0,*Kms_mod,1.0);
	mod = LINALG::MatMatMult(*Mbar,true,*Kss,false);
	LINALG::Add(*mod,false,1.0,*Kms_mod,1.0);
	LINALG::Complete(*Kms_mod,Kms->DomainMap(),Kms->RowMap());
	
	// Kmm: add Kms*T(Mbar) + T(Mbar)*Ksm + T(Mbar)*Kss*Mbar
	RCP<Epetra_CrsMatrix> Kmm_mod = LINALG::CreateMatrix(Kmm->RowMap(),100);
	LINALG::Add(*Kmm,false,1.0,*Kmm_mod,1.0);
	mod = LINALG::MatMatMult(*Kms,false,*Mbar,false);
	LINALG::Add(*mod,false,1.0,*Kmm_mod,1.0);
	mod = LINALG::MatMatMult(*Mbar,true,*Ksm,false);
	LINALG::Add(*mod,false,1.0,*Kmm_mod,1.0);
	mod = LINALG::MatMatMult(*Mbar,true,*Kss,false);
	mod = LINALG::MatMatMult(*mod,false,*Mbar,false);
	LINALG::Add(*mod,false,1.0,*Kmm_mod,1.0);
	LINALG::Complete(*Kmm_mod,Kmm->DomainMap(),Kmm->RowMap());
	
	// Kmn: add T(Mbar)*Ksn
	RCP<Epetra_CrsMatrix> Kmn_mod = LINALG::CreateMatrix(Kmn->RowMap(),100);
	LINALG::Add(*Kmn,false,1.0,*Kmn_mod,1.0);	
	mod = LINALG::MatMatMult(*Mbar,true,*Ksn,false);
	LINALG::Add(*mod,false,1.0,*Kmn_mod,1.0);
	LINALG::Complete(*Kmn_mod,Kmn->DomainMap(),Kmn->RowMap());
	
	// Kns: nothing to do
	RCP<Epetra_CrsMatrix> Kns_mod = Kns;
	
	// Knm: add Kns*Mbar
	RCP<Epetra_CrsMatrix> Knm_mod = LINALG::CreateMatrix(Knm->RowMap(),100);
	LINALG::Add(*Knm,false,1.0,*Knm_mod,1.0);
	mod = LINALG::MatMatMult(*Kns,false,*Mbar,false);
	LINALG::Add(*mod,false,1.0,*Knm_mod,1.0);
	LINALG::Complete(*Knm_mod,Knm->DomainMap(),Knm->RowMap());
	
	// Knn: nothing to do
	RCP<Epetra_CrsMatrix> Knn_mod = Knn;
	
	// output for checking everything
	if(Kteff->Comm().MyPID()==0)
	{
		cout << endl << "*****************" << endl;
		cout << "K_mod:   " << Kteff->NumGlobalRows() << " x " << Kteff->NumGlobalCols() << endl;
		if (Kss_mod!=null) cout << "Kss_mod: " << Kss_mod->NumGlobalRows() << " x " << Kss_mod->NumGlobalCols() << endl;
		else cout << "Kss_mod: null" << endl;
		if (Ksm_mod!=null) cout << "Ksm_mod: " << Ksm_mod->NumGlobalRows() << " x " << Ksm_mod->NumGlobalCols() << endl;
		else cout << "Ksm_mod: null" << endl;
		if (Ksn_mod!=null) cout << "Ksn_mod: " << Ksn_mod->NumGlobalRows() << " x " << Ksn_mod->NumGlobalCols() << endl;
		else cout << "Ksn_mod: null" << endl;
		if (Kms_mod!=null) cout << "Kms_mod: " << Kms_mod->NumGlobalRows() << " x " << Kms_mod->NumGlobalCols() << endl;
		else cout << "Kms_mod: null" << endl;
		if (Kmm_mod!=null) cout << "Kmm_mod: " << Kmm_mod->NumGlobalRows() << " x " << Kmm_mod->NumGlobalCols() << endl;
		else cout << "Kmm_mod: null" << endl;
		if (Kmn_mod!=null) cout << "Kmn_mod: " << Kmn_mod->NumGlobalRows() << " x " << Kmn_mod->NumGlobalCols() << endl;
		else cout << "Kmn_mod: null" << endl;
		if (Kns_mod!=null) cout << "Kns_mod: " << Kns_mod->NumGlobalRows() << " x " << Kns_mod->NumGlobalCols() << endl;
		else cout << "Kns_mod: null" << endl;
		if (Knm_mod!=null) cout << "Knm_mod: " << Knm_mod->NumGlobalRows() << " x " << Knm_mod->NumGlobalCols() << endl;
		else cout << "Knm_mod: null" << endl;
		if (Knn_mod!=null) cout << "Knn_mod: " << Knn_mod->NumGlobalRows() << " x " << Knn_mod->NumGlobalCols() << endl;
		else cout << "Knn_mod: null" << endl;
		cout << "*****************" << endl;
	}
	
	/**********************************************************************/
	/* Apply basis transformation to f                                    */
	/**********************************************************************/
	// fs: nothing to be done
	RCP<Epetra_Vector> fs_mod = fs;
	
	// fm: add T(Mbar)*fs
	RCP<Epetra_Vector> fm_mod = fm;
	Mbar->Multiply(true,*fs,*fm_mod);

	// fn: nothing to be done
	RCP<Epetra_Vector> fn_mod = fn;
	
	// output for checking everything
	if(feff->Comm().MyPID()==0)
	{
		cout << endl << "**********" << endl;
		cout << "f_mod:  " << feff->GlobalLength() << endl;
		cout << "fs_mod: " << fs_mod->GlobalLength() << endl;
		cout << "fm_mod: " << fm_mod->GlobalLength() << endl;
		cout << "fn_mod: " << fn_mod->GlobalLength() << endl;
		cout << "**********" << endl;
	}
	
	//exit(0);
	/*
#ifdef DEBUG

	// Testing of global matrices and inverse
	 * 
	cout << *D_;
	cout << *M_;
	cout << *g_;
	cout << *N_;
	cout << *T_;
	cout << *invD;
	cout << *Mbar;

	// Testing of LINALG::Transpose
	 * 
	cout << *M_;
	cout << M_->RangeMap();
	cout << M_->DomainMap();
	RCP<Epetra_CrsMatrix> transM = LINALG::Transpose(*M_);
	cout << *transM;
	cout << transM->RangeMap();
	cout << transM->DomainMap();

	// Testing of RemoteIDList function
	cout << *gsnoderowmap_;
	int PID = 0;
	int LID = 0;
	const int GID = 312;
	gsnoderowmap_->RemoteIDList(1,&GID,&PID,&LID);
	cout << "Global ID: " << GID << endl;
	cout << "Local ID : " << LID << endl;
	cout << "Proc ID  : " << PID << endl;

	// Testing of LINALG::SplitMatrix2x2
	RCP<Epetra_CrsMatrix> A11;
	RCP<Epetra_CrsMatrix> A12;
	RCP<Epetra_CrsMatrix> A21;
	RCP<Epetra_CrsMatrix> A22;
	RCP<Epetra_Map> testmap = null;
	RCP<Epetra_Map> testmap2 = null;
	RCP<Epetra_Map> full = rcp(new Epetra_Map(Kteff->DomainMap()));
	LINALG::SplitMatrix2x2(Kteff,full,testmap,full,testmap2,A11,A12,A21,A22);
	
	if(Kteff->Comm().MyPID()==0)
	{
		cout << endl << "*****************" << endl;
		cout << "K:   " << Kteff->NumGlobalRows() << " x " << Kteff->NumGlobalCols() << endl;
		if (A11!=null) cout << "A11: " << A11->NumGlobalRows() << " x " << A11->NumGlobalCols() << endl;
		else cout << "A11: null" << endl;
		if (A12!=null) cout << "A12: " << A12->NumGlobalRows() << " x " << A12->NumGlobalCols() << endl;
		else cout << "A11: null" << endl;
		if (A21!=null) cout << "A21: " << A21->NumGlobalRows() << " x " << A21->NumGlobalCols() << endl;
		else cout << "A21: null" << endl;
		if (A22!=null) cout << "A22: " << A22->NumGlobalRows() << " x " << A22->NumGlobalCols() << endl;
		else cout << "A22: null" << endl;
		cout << "*****************" << endl;
	}
#endif // #ifdef DEBUG
	*/
	
  return;
}

#endif  // #ifdef CCADISCRET
