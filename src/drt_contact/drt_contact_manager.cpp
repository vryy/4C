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

//#define DIMOUTPUT

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::Manager::Manager(DRT::Discretization& discret) :
discret_(discret),
activesetconv_(false)
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
	g_    = LINALG::CreateVector(*gsnoderowmap_,true);
		
	// (re)setup global normal and tangent matrices
	nmatrix_ = rcp(new LINALG::SparseMatrix(*gactiven_,3));
	tmatrix_ = rcp(new LINALG::SparseMatrix(*gactivet_,3));
	
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
void CONTACT::Manager::Evaluate(RCP<LINALG::SparseMatrix>& Kteff,
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
	/* build global matrix N with normal vectors of active nodes          */
	/* and global matrix T with tangent vectors of active nodes           */
	/**********************************************************************/
	for (int i=0; i<(int)interface_.size(); ++i)
	  interface_[i]->AssembleNT(*nmatrix_,*tmatrix_);
	  
	// FillComplete() global matrices N and T
	nmatrix_->Complete(*gactivedofs_,*gactiven_);
	tmatrix_->Complete(*gactivedofs_,*gactivet_);
	
	/**********************************************************************/
	/* Multiply Mortar matrices: M^ = inv(D) * M                          */
	/**********************************************************************/
	RCP<LINALG::SparseMatrix> invD = rcp(new LINALG::SparseMatrix(*dmatrix_));
	RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_,true);
	int err = 0;
	
	// extract diagonal of invD into diag
	invD->ExtractDiagonalCopy(*diag);
	
	// scalar inversion of diagonal values
	err = diag->Reciprocal(*diag);
	if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");
	
	// re-insert inverted diagonal into invD
	err = invD->ReplaceDiagonalValues(*diag);
	if (err>0) dserror("ERROR: ReplaceDiagonalValues: Missing diagonal entry!");
	
	// do the multiplication M^ = inv(D) * M
	mhatmatrix_ = LINALG::Multiply(*invD,false,*mmatrix_,false);
	
	/**********************************************************************/
	/* Split Kteff into 3x3 block matrix                                  */
	/**********************************************************************/
	// we want to split K into 3 groups s,m,n = 9 blocks
	RCP<LINALG::SparseMatrix> Kss, Ksm, Ksn, Kms, Kmm, Kmn, Kns, Knm, Knn;
	
	// temporarily we need the blocks K_sm_sm, K_sm_n, K_n_sm
	// (FIXME: because a direct SplitMatrix3x3 is still missing!) 
	RCP<LINALG::SparseMatrix> Ksmsm, Ksmn, Knsm;
	
	// we also need the combined sm rowmap
	RCP<Epetra_Map> gsmdofs = LINALG::MergeMap(gsdofrowmap_,gmdofrowmap_);
	
	// some temporary RCPs
	RCP<Epetra_Map> tempmap;
	RCP<LINALG::SparseMatrix> tempmtx;
	
	// split into slave/master part + structure part
	LINALG::SplitMatrix2x2(Kteff,gsmdofs,gndofrowmap_,gsmdofs,gndofrowmap_,Ksmsm,Ksmn,Knsm,Knn);
	
	// further splits into slave part + master part
	LINALG::SplitMatrix2x2(Ksmsm,gsdofrowmap_,gmdofrowmap_,gsdofrowmap_,gmdofrowmap_,Kss,Ksm,Kms,Kmm);
	LINALG::SplitMatrix2x2(Ksmn,gsdofrowmap_,gmdofrowmap_,gndofrowmap_,tempmap,Ksn,tempmtx,Kmn,tempmtx);
	LINALG::SplitMatrix2x2(Knsm,gndofrowmap_,tempmap,gsdofrowmap_,gmdofrowmap_,Kns,Knm,tempmtx,tempmtx);
	
	// output for checking everything
#ifdef DIMOUTPUT
	if(Comm().MyPID()==0)
	{
		cout << endl << "*****************" << endl;
		cout << "K:   " << Kteff->EpetraMatrix()->NumGlobalRows() << " x " << Kteff->EpetraMatrix()->NumGlobalCols() << endl;
		if (Kss!=null) cout << "Kss: " << Kss->EpetraMatrix()->NumGlobalRows() << " x " << Kss->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Kss: null" << endl;
		if (Ksm!=null) cout << "Ksm: " << Ksm->EpetraMatrix()->NumGlobalRows() << " x " << Ksm->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Ksm: null" << endl;
		if (Ksn!=null) cout << "Ksn: " << Ksn->EpetraMatrix()->NumGlobalRows() << " x " << Ksn->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Ksn: null" << endl;
		if (Kms!=null) cout << "Kms: " << Kms->EpetraMatrix()->NumGlobalRows() << " x " << Kms->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Kms: null" << endl;
		if (Kmm!=null) cout << "Kmm: " << Kmm->EpetraMatrix()->NumGlobalRows() << " x " << Kmm->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Kmm: null" << endl;
		if (Kmn!=null) cout << "Kmn: " << Kmn->EpetraMatrix()->NumGlobalRows() << " x " << Kmn->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Kmn: null" << endl;
		if (Kns!=null) cout << "Kns: " << Kns->EpetraMatrix()->NumGlobalRows() << " x " << Kns->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Kns: null" << endl;
		if (Knm!=null) cout << "Knm: " << Knm->EpetraMatrix()->NumGlobalRows() << " x " << Knm->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Knm: null" << endl;
		if (Knn!=null) cout << "Knn: " << Knn->EpetraMatrix()->NumGlobalRows() << " x " << Knn->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Knn: null" << endl;
		cout << "*****************" << endl;
	}
#endif // #ifdef DIMOUTPUT
	
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
#ifdef DIMOUTPUT
	if(Comm().MyPID()==0)
	{
		cout << endl << "**********" << endl;
		cout << "f:  " << feff->GlobalLength() << endl;
		cout << "fs: " << fs->GlobalLength() << endl;
		cout << "fm: " << fm->GlobalLength() << endl;
		cout << "fn: " << fn->GlobalLength() << endl;
		cout << "**********" << endl;
	}
#endif // #ifdef DIMOUTPUT
	
	/**********************************************************************/
	/* Apply basis transformation to K                                    */
	/**********************************************************************/
	// define temporary RCP mod
	RCP<LINALG::SparseMatrix> mod;
	
	// Kss: nothing to do
	RCP<LINALG::SparseMatrix> Kssmod = Kss;
	
	// Ksm: add Kss*T(Mbar)
	RCP<LINALG::SparseMatrix> Ksmmod = rcp(new LINALG::SparseMatrix(Ksm->RowMap(),100));
	Ksmmod->Add(*Ksm,false,1.0,1.0);
	mod = LINALG::Multiply(*Kss,false,*mhatmatrix_,false);
	Ksmmod->Add(*mod,false,1.0,1.0);
	Ksmmod->Complete(Ksm->DomainMap(),Ksm->RowMap());

	// Ksn: nothing to do
	RCP<LINALG::SparseMatrix> Ksnmod = Ksn;
	
	// Kms: add T(Mbar)*Kss
	RCP<LINALG::SparseMatrix> Kmsmod = rcp(new LINALG::SparseMatrix(Kms->RowMap(),100));
	Kmsmod->Add(*Kms,false,1.0,1.0);
	mod = LINALG::Multiply(*mhatmatrix_,true,*Kss,false);
	Kmsmod->Add(*mod,false,1.0,1.0);
	Kmsmod->Complete(Kms->DomainMap(),Kms->RowMap());
	
	// Kmm: add Kms*T(Mbar) + T(Mbar)*Ksm + T(Mbar)*Kss*Mbar
	RCP<LINALG::SparseMatrix> Kmmmod = rcp(new LINALG::SparseMatrix(Kmm->RowMap(),100));
	Kmmmod->Add(*Kmm,false,1.0,1.0);
	mod = LINALG::Multiply(*Kms,false,*mhatmatrix_,false);
	Kmmmod->Add(*mod,false,1.0,1.0);
	mod = LINALG::Multiply(*mhatmatrix_,true,*Ksm,false);
	Kmmmod->Add(*mod,false,1.0,1.0);
	mod = LINALG::Multiply(*mhatmatrix_,true,*Kss,false);
	mod = LINALG::Multiply(*mod,false,*mhatmatrix_,false);
	Kmmmod->Add(*mod,false,1.0,1.0);
	Kmmmod->Complete(Kmm->DomainMap(),Kmm->RowMap());
	
	// Kmn: add T(Mbar)*Ksn
	RCP<LINALG::SparseMatrix> Kmnmod = rcp(new LINALG::SparseMatrix(Kmn->RowMap(),100));
	Kmnmod->Add(*Kmn,false,1.0,1.0);	
	mod = LINALG::Multiply(*mhatmatrix_,true,*Ksn,false);
	Kmnmod->Add(*mod,false,1.0,1.0);
	Kmnmod->Complete(Kmn->DomainMap(),Kmn->RowMap());
	
	// Kns: nothing to do
	RCP<LINALG::SparseMatrix> Knsmod = Kns;
	
	// Knm: add Kns*Mbar
	RCP<LINALG::SparseMatrix> Knmmod = rcp(new LINALG::SparseMatrix(Knm->RowMap(),100));
	Knmmod->Add(*Knm,false,1.0,1.0);
	mod = LINALG::Multiply(*Kns,false,*mhatmatrix_,false);
	Knmmod->Add(*mod,false,1.0,1.0);
	Knmmod->Complete(Knm->DomainMap(),Knm->RowMap());
	
	// Knn: nothing to do
	RCP<LINALG::SparseMatrix> Knnmod = Knn;
	
	// output for checking everything
#ifdef DIMOUTPUT
	if(Comm().MyPID()==0)
	{
		cout << endl << "*****************" << endl;
		cout << "Kmod:   " << Kteff->EpetraMatrix()->NumGlobalRows() << " x " << Kteff->EpetraMatrix()->NumGlobalCols() << endl;
		if (Kssmod!=null) cout << "Kssmod: " << Kssmod->EpetraMatrix()->NumGlobalRows() << " x " << Kssmod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Kssmod: null" << endl;
		if (Ksmmod!=null) cout << "Ksmmod: " << Ksmmod->EpetraMatrix()->NumGlobalRows() << " x " << Ksmmod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Ksmmod: null" << endl;
		if (Ksnmod!=null) cout << "Ksnmod: " << Ksnmod->EpetraMatrix()->NumGlobalRows() << " x " << Ksnmod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Ksnmod: null" << endl;
		if (Kmsmod!=null) cout << "Kmsmod: " << Kmsmod->EpetraMatrix()->NumGlobalRows() << " x " << Kmsmod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Kmsmod: null" << endl;
		if (Kmmmod!=null) cout << "Kmmmod: " << Kmmmod->EpetraMatrix()->NumGlobalRows() << " x " << Kmmmod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Kmmmod: null" << endl;
		if (Kmnmod!=null) cout << "Kmnmod: " << Kmnmod->EpetraMatrix()->NumGlobalRows() << " x " << Kmnmod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Kmnmod: null" << endl;
		if (Knsmod!=null) cout << "Knsmod: " << Knsmod->EpetraMatrix()->NumGlobalRows() << " x " << Knsmod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Knsmod: null" << endl;
		if (Knmmod!=null) cout << "Knmmod: " << Knmmod->EpetraMatrix()->NumGlobalRows() << " x " << Knmmod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Knmmod: null" << endl;
		if (Knnmod!=null) cout << "Knnmod: " << Knnmod->EpetraMatrix()->NumGlobalRows() << " x " << Knnmod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Knnmod: null" << endl;
		cout << "*****************" << endl;
	}
#endif // #ifdef DIMOUTPUT
	
	/**********************************************************************/
	/* Apply basis transformation to f                                    */
	/**********************************************************************/
	// fs: nothing to be done
	RCP<Epetra_Vector> fsmod = fs;
	
	// fm: add T(Mbar)*fs
	RCP<Epetra_Vector> fmmod = rcp(new Epetra_Vector(*gmdofrowmap_));
	mhatmatrix_->Multiply(true,*fs,*fmmod);
	fmmod->Update(1.0,*fm,1.0);
	
	// fn: nothing to be done
	RCP<Epetra_Vector> fnmod = fn;
	
	// output for checking everything
#ifdef DIMOUTPUT
	if(Comm().MyPID()==0)
	{
		cout << endl << "**********" << endl;
		cout << "fmod:  " << feff->GlobalLength() << endl;
		cout << "fsmod: " << fsmod->GlobalLength() << endl;
		cout << "fmmod: " << fmmod->GlobalLength() << endl;
		cout << "fnmod: " << fnmod->GlobalLength() << endl;
		cout << "**********" << endl;
	}
#endif // #ifdef DIMOUTPUT
	
	/**********************************************************************/
	/* Split slave quantities into active / inactive                      */
	/**********************************************************************/
	// we want to split Kssmod into 2 groups a,i = 4 blocks
	RCP<LINALG::SparseMatrix> Kaamod, Kaimod, Kiamod, Kiimod;
	
	// we want to split Ksnmod / Ksmmod into 2 groups a,i = 2 blocks
	RCP<LINALG::SparseMatrix> Kanmod, Kinmod, Kammod, Kimmod;
		
	// we will get the i rowmap as a by-product
	RCP<Epetra_Map> gidofs;
		
	// do the splitting
	LINALG::SplitMatrix2x2(Kssmod,gactivedofs_,gidofs,gactivedofs_,gidofs,Kaamod,Kaimod,Kiamod,Kiimod);
	LINALG::SplitMatrix2x2(Ksnmod,gactivedofs_,gidofs,gndofrowmap_,tempmap,Kanmod,tempmtx,Kinmod,tempmtx);
	LINALG::SplitMatrix2x2(Ksmmod,gactivedofs_,gidofs,gmdofrowmap_,tempmap,Kammod,tempmtx,Kimmod,tempmtx);
	
	// output for checking everything
#ifdef DIMOUTPUT
	if(Comm().MyPID()==0)
	{
		cout << endl << "*****************" << endl;
		cout << "Kssmod: " << Kssmod->EpetraMatrix()->NumGlobalRows() << " x " << Kssmod->EpetraMatrix()->NumGlobalCols() << endl;
		if (Kaamod!=null) cout << "Kaamod: " << Kaamod->EpetraMatrix()->NumGlobalRows() << " x " << Kaamod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Kaamod: null" << endl;
		if (Kaimod!=null) cout << "Kaimod: " << Kaimod->EpetraMatrix()->NumGlobalRows() << " x " << Kaimod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Kaimod: null" << endl;
		if (Kiamod!=null) cout << "Kiamod: " << Kiamod->EpetraMatrix()->NumGlobalRows() << " x " << Kiamod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Kiamod: null" << endl;
		if (Kiimod!=null) cout << "Kiimod: " << Kiimod->EpetraMatrix()->NumGlobalRows() << " x " << Kiimod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Kiimod: null" << endl;
		cout << "*****************" << endl;
		
		cout << endl << "*****************" << endl;
		cout << "Ksnmod: " << Ksnmod->EpetraMatrix()->NumGlobalRows() << " x " << Ksnmod->EpetraMatrix()->NumGlobalCols() << endl;
		if (Kanmod!=null) cout << "Kanmod: " << Kanmod->EpetraMatrix()->NumGlobalRows() << " x " << Kanmod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Kanmod: null" << endl;
		if (Kinmod!=null) cout << "Kinmod: " << Kinmod->EpetraMatrix()->NumGlobalRows() << " x " << Kinmod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Kinmod: null" << endl;
		cout << "*****************" << endl;
		
		cout << endl << "*****************" << endl;
		cout << "Ksmmod: " << Ksmmod->EpetraMatrix()->NumGlobalRows() << " x " << Ksmmod->EpetraMatrix()->NumGlobalCols() << endl;
		if (Kammod!=null) cout << "Kammod: " << Kammod->EpetraMatrix()->NumGlobalRows() << " x " << Kammod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Kammod: null" << endl;
		if (Kimmod!=null) cout << "Kimmod: " << Kimmod->EpetraMatrix()->NumGlobalRows() << " x " << Kimmod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "Kimmod: null" << endl;
		cout << "*****************" << endl;
	}
#endif // #ifdef DIMOUTPUT
	
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
#ifdef DIMOUTPUT
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
#endif // #ifdef DIMOUTPUT
	
	// do the multiplications with T-matrix
	RCP<LINALG::SparseMatrix> TKanmod, TKammod, TKaimod, TKaamod;
	RCP<Epetra_Vector> Tfamod;
	
	if(gactivedofs_->NumGlobalElements())
	{
		TKanmod = LINALG::Multiply(*tmatrix_,false,*Kanmod,false);
		TKammod = LINALG::Multiply(*tmatrix_,false,*Kammod,false);
		TKaamod = LINALG::Multiply(*tmatrix_,false,*Kaamod,false);
		
		if (gidofs->NumGlobalElements())
			TKaimod = LINALG::Multiply(*tmatrix_,false,*Kaimod,false);
		
		Tfamod = rcp(new Epetra_Vector(tmatrix_->RowMap()));
		tmatrix_->Multiply(false,*famod,*Tfamod);
	}
	
	// output for checking everything
#ifdef DIMOUTPUT
	if(Comm().MyPID()==0)
	{
		cout << endl << "*****************" << endl;
		if (TKanmod!=null) cout << "TKanmod: " << TKanmod->EpetraMatrix()->NumGlobalRows() << " x " << TKanmod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "TKanmod: null" << endl;
		if (TKammod!=null) cout << "TKammod: " << TKammod->EpetraMatrix()->NumGlobalRows() << " x " << TKammod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "TKammod: null" << endl;
		if (TKaimod!=null) cout << "TKaimod: " << TKaimod->EpetraMatrix()->NumGlobalRows() << " x " << TKaimod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "TKaimod: null" << endl;
		if (TKaamod!=null) cout << "TKaamod: " << TKaamod->EpetraMatrix()->NumGlobalRows() << " x " << TKaamod->EpetraMatrix()->NumGlobalCols() << endl;
		else cout << "TKaamod: null" << endl;
		cout << "*****************" << endl;
		
		cout << endl << "**********" << endl;
		if (Tfamod!=null) cout << "Tfamod: " << Tfamod->GlobalLength() << endl;
		else cout << "Tfamod: null" << endl;
		cout << "**********" << endl;
	}
#endif // #ifdef DIMOUTPUT
	
	/**********************************************************************/
	/* Global setup of Kteffnew, feffnew (including contact)            */
	/**********************************************************************/
	RCP<LINALG::SparseMatrix> Kteffnew = rcp(new LINALG::SparseMatrix(*(discret_.DofRowMap()),81));
	RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*(discret_.DofRowMap()));
	
	// add n / m submatrices to Kteffnew
	Kteffnew->Add(*Knnmod,false,1.0,1.0);
	Kteffnew->Add(*Knmmod,false,1.0,1.0);
	Kteffnew->Add(*Kmnmod,false,1.0,1.0);
	Kteffnew->Add(*Kmmmod,false,1.0,1.0);
	
	// add a / i submatrices to Kteffnew, if existing
	if (Knsmod!=null) Kteffnew->Add(*Knsmod,false,1.0,1.0);
	if (Kmsmod!=null) Kteffnew->Add(*Kmsmod,false,1.0,1.0);
	if (Kinmod!=null) Kteffnew->Add(*Kinmod,false,1.0,1.0);
	if (Kimmod!=null) Kteffnew->Add(*Kimmod,false,1.0,1.0);
	if (Kiimod!=null) Kteffnew->Add(*Kiimod,false,1.0,1.0);
	if (Kiamod!=null) Kteffnew->Add(*Kiamod,false,1.0,1.0);
	
	// add matrix of normals to Kteffnew
	Kteffnew->Add(*nmatrix_,false,1.0,1.0);

	// add submatrices with tangents to Kteffnew, if existing
	if (TKanmod!=null) Kteffnew->Add(*TKanmod,false,1.0,1.0);
	if (TKammod!=null) Kteffnew->Add(*TKammod,false,1.0,1.0);
	if (TKaimod!=null) Kteffnew->Add(*TKaimod,false,1.0,1.0);
	if (TKaamod!=null) Kteffnew->Add(*TKaamod,false,1.0,1.0);
	
	// FillComplete Kteffnew (square)
	Kteffnew->Complete();
	
	// add n / m subvectors to feffnew
	RCP<Epetra_Vector> fnmodexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
	RCP<Epetra_Vector> fmmodexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
	LINALG::Export(*fnmod,*fnmodexp);
	LINALG::Export(*fmmod,*fmmodexp);
	feffnew->Update(1.0,*fnmodexp,1.0,*fmmodexp,1.0);
	
	// add i / Ta subvectors to feffnew, if existing
	RCP<Epetra_Vector> fimodexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
	RCP<Epetra_Vector> Tfamodexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
	if (fimod!=null) LINALG::Export(*fimod,*fimodexp);
	if (Tfamod!=null) LINALG::Export(*Tfamod,*Tfamodexp);
	feffnew->Update(1.0,*fimodexp,1.0,*Tfamodexp,1.0);
	
	// add weighted gap vector to feffnew, if existing
	RCP<Epetra_Vector> gexp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
	if (gact->GlobalLength()) LINALG::Export(*gact,*gexp);
	feffnew->Update(1.0,*gexp,1.0);
	
	/**********************************************************************/
	/* Replace Kteff and feff by Kteffnew and feffnew                   */
	/**********************************************************************/
	Kteff = Kteffnew;
	feff = feffnew;
	LINALG::PrintSparsityToPostscript(*(Kteff->EpetraMatrix()));
	//exit(0);
	
	/**********************************************************************/
	/* Update Lagrange multipliers                                        */
	/**********************************************************************/
	invD->Multiply(false,*fsmod,*z_);

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
	
	return;
}

/*----------------------------------------------------------------------*
 |  Update active set and check for convergence (public)      popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::UpdateActiveSet()
{
	// assume that active set has converged and check for opposite
	activesetconv_=true;
	
	for (int i=0; i<(int)interface_.size(); ++i)
	{
		if (i>0) dserror("ERROR: UpdateActiveSet: Not implemented yet for multiple interfaces!");
		for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
		{
			int gid = interface_[i]->SlaveRowNodes()->GID(j);
			DRT::Node* node = interface_[i]->Discret().gNode(gid);
			if (!node) dserror("ERROR: Cannot find node with gid %",gid);
			CNode* cnode = static_cast<CNode*>(node);
			
			// check nodes of inactive set
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
				
				//cout << "INACTIVE: " << i << " " << j << " " << gid << " " << nincr << " " << (*g_)[g_->Map().LID(gid)] << " " << cnode->HasProj() << endl;
			}
			
			// check nodes of active set
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
				if ((0.5*nz+0.5*nzold)<0)
				{
					cnode->Active() = false;
					activesetconv_ = false;
				}
				
				//cout << "ACTIVE: " << i << " " << j << " " << gid << " " << nz << " " << nzold << " " << 0.5*nz+0.5*nzold << " " << cnode->Getg() << endl;	
			}
			
		}
  }

	// broadcast convergence flag among processors
	int convcheck = 0;
	int localcheck = activesetconv_;
	Comm().SumAll(&localcheck,&convcheck,1);
	if (convcheck!=Comm().NumProc())
		activesetconv_=false;
	
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
	
	// store zold_ Lagrange multipliers if active set converged
	if (activesetconv_) zold_=rcp(new Epetra_Vector(*z_));
	
	if (activesetconv_)
		for (int i=0;i<(int)interface_.size();++i)
			interface_[i]->VisualizeGmsh(interface_[i]->CSegs());
	
	return;
}

/*----------------------------------------------------------------------*
 |  Compute contact forces (public)                           popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::ContactForces(RCP<Epetra_Vector>& fc)
{
	// residuum check
	RCP<Epetra_Vector> testvec1 = rcp(new Epetra_Vector(dmatrix_->RowMap()));
	RCP<Epetra_Vector> testvec2 = rcp(new Epetra_Vector(mmatrix_->DomainMap()));
	dmatrix_->Multiply(false,*z_,*testvec1);
	mmatrix_->Multiply(true,*z_,*testvec2);
	
	RCP<Epetra_Vector> vec1 = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
	RCP<Epetra_Vector> vec2 = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
	LINALG::Export(*testvec1,*vec1);
	LINALG::Export(*testvec2,*vec2);
	
	fc->Update(1.0,*vec1,0.0);
	fc->Update(-1.0,*vec2,1.0);
	
	double testres1 = 0.0;
	double testres2 = 0.0;
	double fcnorm = 0.0;
	testvec1->Norm2(&testres1);
	testvec2->Norm2(&testres2);
	fc->Norm2(&fcnorm);
	//cout << testres1 << " " << testres2 << " " << fcnorm << endl;
		
	return;
}

#endif  // #ifdef CCADISCRET
