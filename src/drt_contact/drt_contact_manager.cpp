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
		
	// (re)setup active global Epetra_Maps
	gactivenodes_ = null;
	gactivedofs_ = null;
	gactiveN_ = null;
	gactiveT_ = null;
	
	// update active global Epetra_Maps
	for (int i=0;i<(int)interface_.size();++i)
	{
		gactivenodes_ = LINALG::MergeMap(gactivenodes_,interface_[i]->ActiveNodes());
		gactivedofs_ = LINALG::MergeMap(gactivedofs_,interface_[i]->ActiveDofs());
		gactiveN_ = LINALG::MergeMap(gactiveN_,interface_[i]->ActiveNDofs());;
		gactiveT_ = LINALG::MergeMap(gactiveT_,interface_[i]->ActiveTDofs());;
	}
	
	// create empty maps, if active set = null
	if (gactivenodes_==null)
	{
		gactivenodes_ = rcp(new Epetra_Map(0,0,Comm()));
		gactivedofs_ = rcp(new Epetra_Map(0,0,Comm()));
		gactiveN_ = rcp(new Epetra_Map(0,0,Comm()));
		gactiveT_ = rcp(new Epetra_Map(0,0,Comm()));
	}
	
	// (re)setup global Mortar Epetra_CrsMatrices and Epetra_Vectors
	D_    = LINALG::CreateMatrix(*gsdofrowmap_,10);
	M_    = LINALG::CreateMatrix(*gsdofrowmap_,100);
	Mbar_ = LINALG::CreateMatrix(*gsdofrowmap_,100);
	g_    = LINALG::CreateVector(*gsnoderowmap_,true);
		
	// (re)setup global normal and tangent matrices
	N_ = LINALG::CreateMatrix(*gactiveN_,3);
	T_ = LINALG::CreateMatrix(*gactiveT_,3);
	
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
	LINALG::Complete(D_);
	LINALG::Complete(M_,*gmdofrowmap_,*gsdofrowmap_);
	
	// export weighted gap vector to gactiveN-map
	RCP<Epetra_Vector> g_act = LINALG::CreateVector(*gactivenodes_,true);
	if (g_act->GlobalLength())
	{
		LINALG::Export(*g_,*g_act);
		g_act->ReplaceMap(*gactiveN_);
	}
	
	/**********************************************************************/
	/* build global matrix N with normal vectors of active nodes          */
	/* and global matrix T with tangent vectors of active nodes           */
	/**********************************************************************/
	for (int i=0; i<(int)interface_.size(); ++i)
	  interface_[i]->Assemble_NT(*N_,*T_);
	  
	// FillComplete() global matrices N and T
	LINALG::Complete(N_,*gactivedofs_,*gactiveN_);
	LINALG::Complete(T_,*gactivedofs_,*gactiveT_);
	
	/**********************************************************************/
	/* Multiply Mortar matrices: M^ = inv(D) * M                          */
	/**********************************************************************/
	RCP<Epetra_CrsMatrix> invD = rcp(new Epetra_CrsMatrix(*D_));
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
	Mbar_ = LINALG::Multiply(invD,false,M_,false);
	
	/**********************************************************************/
	/* Split Kteff into 3x3 block matrix                                  */
	/**********************************************************************/
	// we want to split K into 3 groups s,m,n = 9 blocks
	RCP<Epetra_CrsMatrix> Kss, Ksm, Ksn, Kms, Kmm, Kmn, Kns, Knm, Knn;
	
	// temporarily we need the blocks K_sm_sm, K_sm_n, K_n_sm
	// (FIXME: because a direct SplitMatrix3x3 is still missing!) 
	RCP<Epetra_CrsMatrix> Ksmsm, Ksmn, Knsm;
	
	// we also need the combined sm rowmap
	RCP<Epetra_Map> gsmdofs = LINALG::MergeMap(gsdofrowmap_,gmdofrowmap_);
	
	// we will get the n rowmap as a by-product
	RCP<Epetra_Map> gndofs;
	
	// some temporary RCPs
	RCP<Epetra_Map> tempmap;
	RCP<Epetra_CrsMatrix> tempmtx;
	
	// split into slave/master part + structure part
	LINALG::SplitMatrix2x2(Kteff,gsmdofs,gndofs,Ksmsm,Ksmn,Knsm,Knn);
	
	// further splits into slave part + master part
	LINALG::SplitMatrix2x2(Ksmsm,gsdofrowmap_,gmdofrowmap_,Kss,Ksm,Kms,Kmm);
	LINALG::SplitMatrix2x2(Ksmn,gsdofrowmap_,gmdofrowmap_,gndofs,tempmap,Ksn,tempmtx,Kmn,tempmtx);
	LINALG::SplitMatrix2x2(Knsm,gndofs,tempmap,gsdofrowmap_,gmdofrowmap_,Kns,Knm,tempmtx,tempmtx);
	
	// output for checking everything
	if(Comm().MyPID()==0)
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
	if(Comm().MyPID()==0)
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
	// define temporary RCP mod
	RCP<Epetra_CrsMatrix> mod;
	
	// Kss: nothing to do
	RCP<Epetra_CrsMatrix> Kss_mod = Kss;
	
	// Ksm: add Kss*T(Mbar)
	RCP<Epetra_CrsMatrix> Ksm_mod = LINALG::CreateMatrix(Ksm->RowMap(),100);
	LINALG::Add(Ksm,false,1.0,Ksm_mod,1.0);
	mod = LINALG::Multiply(Kss,false,Mbar_,false);
	LINALG::Add(mod,false,1.0,Ksm_mod,1.0);
	LINALG::Complete(Ksm_mod,Ksm->DomainMap(),Ksm->RowMap());

	// Ksn: nothing to do
	RCP<Epetra_CrsMatrix> Ksn_mod = Ksn;
	
	// Kms: add T(Mbar)*Kss
	RCP<Epetra_CrsMatrix> Kms_mod = LINALG::CreateMatrix(Kms->RowMap(),100);
	LINALG::Add(Kms,false,1.0,Kms_mod,1.0);
	mod = LINALG::Multiply(Mbar_,true,Kss,false);
	LINALG::Add(mod,false,1.0,Kms_mod,1.0);
	LINALG::Complete(Kms_mod,Kms->DomainMap(),Kms->RowMap());
	
	// Kmm: add Kms*T(Mbar) + T(Mbar)*Ksm + T(Mbar)*Kss*Mbar
	RCP<Epetra_CrsMatrix> Kmm_mod = LINALG::CreateMatrix(Kmm->RowMap(),100);
	LINALG::Add(Kmm,false,1.0,Kmm_mod,1.0);
	mod = LINALG::Multiply(Kms,false,Mbar_,false);
	LINALG::Add(mod,false,1.0,Kmm_mod,1.0);
	mod = LINALG::Multiply(Mbar_,true,Ksm,false);
	LINALG::Add(mod,false,1.0,Kmm_mod,1.0);
	mod = LINALG::Multiply(Mbar_,true,Kss,false,Mbar_,false);
	LINALG::Add(mod,false,1.0,Kmm_mod,1.0);
	LINALG::Complete(*Kmm_mod,Kmm->DomainMap(),Kmm->RowMap());
	
	// Kmn: add T(Mbar)*Ksn
	RCP<Epetra_CrsMatrix> Kmn_mod = LINALG::CreateMatrix(Kmn->RowMap(),100);
	LINALG::Add(Kmn,false,1.0,Kmn_mod,1.0);	
	mod = LINALG::Multiply(Mbar_,true,Ksn,false);
	LINALG::Add(mod,false,1.0,Kmn_mod,1.0);
	LINALG::Complete(Kmn_mod,Kmn->DomainMap(),Kmn->RowMap());
	
	// Kns: nothing to do
	RCP<Epetra_CrsMatrix> Kns_mod = Kns;
	
	// Knm: add Kns*Mbar
	RCP<Epetra_CrsMatrix> Knm_mod = LINALG::CreateMatrix(Knm->RowMap(),100);
	LINALG::Add(Knm,false,1.0,Knm_mod,1.0);
	mod = LINALG::Multiply(Kns,false,Mbar_,false);
	LINALG::Add(mod,false,1.0,Knm_mod,1.0);
	LINALG::Complete(Knm_mod,Knm->DomainMap(),Knm->RowMap());
	
	// Knn: nothing to do
	RCP<Epetra_CrsMatrix> Knn_mod = Knn;
	
	// output for checking everything
	if(Comm().MyPID()==0)
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
	RCP<Epetra_Vector> fm_mod = rcp(new Epetra_Vector(*gmdofrowmap_));
	Mbar_->Multiply(true,*fs,*fm_mod);
	fm_mod->Update(1.0,*fm,1.0);
	
	// fn: nothing to be done
	RCP<Epetra_Vector> fn_mod = fn;
	
	// output for checking everything
	if(Comm().MyPID()==0)
	{
		cout << endl << "**********" << endl;
		cout << "f_mod:  " << feff->GlobalLength() << endl;
		cout << "fs_mod: " << fs_mod->GlobalLength() << endl;
		cout << "fm_mod: " << fm_mod->GlobalLength() << endl;
		cout << "fn_mod: " << fn_mod->GlobalLength() << endl;
		cout << "**********" << endl;
	}
	
	/**********************************************************************/
	/* Split slave quantities into active / inactive                      */
	/**********************************************************************/
	// we want to split Kss_mod into 2 groups a,i = 4 blocks
	RCP<Epetra_CrsMatrix> Kaa_mod, Kai_mod, Kia_mod, Kii_mod;
	
	// we want to split Ksn_mod / Ksm_mod into 2 groups a,i = 2 blocks
	RCP<Epetra_CrsMatrix> Kan_mod, Kin_mod, Kam_mod, Kim_mod;
		
	// we will get the i rowmap as a by-product
	RCP<Epetra_Map> gidofs;
		
	// do the splitting
	LINALG::SplitMatrix2x2(Kss_mod,gactivedofs_,gidofs,gactivedofs_,gidofs,Kaa_mod,Kai_mod,Kia_mod,Kii_mod);
	LINALG::SplitMatrix2x2(Ksn_mod,gactivedofs_,gidofs,gndofs,tempmap,Kan_mod,tempmtx,Kin_mod,tempmtx);
	LINALG::SplitMatrix2x2(Ksm_mod,gactivedofs_,gidofs,gmdofrowmap_,tempmap,Kam_mod,tempmtx,Kim_mod,tempmtx);
	
	// output for checking everything
	if(Comm().MyPID()==0)
	{
		cout << endl << "*****************" << endl;
		cout << "Kss_mod: " << Kss_mod->NumGlobalRows() << " x " << Kss_mod->NumGlobalCols() << endl;
		if (Kaa_mod!=null) cout << "Kaa_mod: " << Kaa_mod->NumGlobalRows() << " x " << Kaa_mod->NumGlobalCols() << endl;
		else cout << "Kaa_mod: null" << endl;
		if (Kai_mod!=null) cout << "Kai_mod: " << Kai_mod->NumGlobalRows() << " x " << Kai_mod->NumGlobalCols() << endl;
		else cout << "Kai_mod: null" << endl;
		if (Kia_mod!=null) cout << "Kia_mod: " << Kia_mod->NumGlobalRows() << " x " << Kia_mod->NumGlobalCols() << endl;
		else cout << "Kia_mod: null" << endl;
		if (Kii_mod!=null) cout << "Kii_mod: " << Kii_mod->NumGlobalRows() << " x " << Kii_mod->NumGlobalCols() << endl;
		else cout << "Kii_mod: null" << endl;
		cout << "*****************" << endl;
		
		cout << endl << "*****************" << endl;
		cout << "Ksn_mod: " << Ksn_mod->NumGlobalRows() << " x " << Ksn_mod->NumGlobalCols() << endl;
		if (Kan_mod!=null) cout << "Kan_mod: " << Kan_mod->NumGlobalRows() << " x " << Kan_mod->NumGlobalCols() << endl;
		else cout << "Kan_mod: null" << endl;
		if (Kin_mod!=null) cout << "Kin_mod: " << Kin_mod->NumGlobalRows() << " x " << Kin_mod->NumGlobalCols() << endl;
		else cout << "Kin_mod: null" << endl;
		cout << "*****************" << endl;
		
		cout << endl << "*****************" << endl;
		cout << "Ksm_mod: " << Ksm_mod->NumGlobalRows() << " x " << Ksm_mod->NumGlobalCols() << endl;
		if (Kam_mod!=null) cout << "Kam_mod: " << Kam_mod->NumGlobalRows() << " x " << Kam_mod->NumGlobalCols() << endl;
		else cout << "Kam_mod: null" << endl;
		if (Kim_mod!=null) cout << "Kim_mod: " << Kim_mod->NumGlobalRows() << " x " << Kim_mod->NumGlobalCols() << endl;
		else cout << "Kim_mod: null" << endl;
		cout << "*****************" << endl;
	}
	
	// we want to split fs_mod into 2 groups a,i
	RCP<Epetra_Vector> fa_mod, fi_mod;
	
	// do the vector splitting s -> a+i
	if (!gidofs->NumGlobalElements())
		fa_mod = rcp(new Epetra_Vector(*fs_mod));
	else if (!gactivedofs_->NumGlobalElements())
		fi_mod = rcp(new Epetra_Vector(*fs_mod));
	else
	{
		Epetra_Vector* f4 = NULL;
		Epetra_Vector* f5 = NULL;
		LINALG::SplitVector(*fs_mod,*gactivedofs_,f4,*gidofs,f5);
		RCP<Epetra_Vector> fa_mod = rcp(f4);
		RCP<Epetra_Vector> fi_mod = rcp(f5);
	}
	
	// output for checking everything
	if(Comm().MyPID()==0)
	{
		cout << endl << "**********" << endl;
		if (fs_mod!=null) cout << "fs_mod: " << fs_mod->GlobalLength() << endl;
		else cout << "fs_mod: null" << endl;
		if (fa_mod!=null) cout << "fa_mod: " << fa_mod->GlobalLength() << endl;
		else cout << "fa_mod: null" << endl;
		if (fi_mod!=null) cout << "fi_mod: " << fi_mod->GlobalLength() << endl;
		else cout << "fi_mod: null" << endl;
		cout << "**********" << endl;
	}
	
	// do the multiplications with T-matrix
	RCP<Epetra_CrsMatrix> TKan_mod, TKam_mod, TKai_mod, TKaa_mod;
	RCP<Epetra_Vector> Tfa_mod;
	
	if(gactivedofs_->NumGlobalElements())
	{
		TKan_mod = LINALG::Multiply(T_,false,Kan_mod,false);
		TKam_mod = LINALG::Multiply(T_,false,Kam_mod,false);
		TKaa_mod = LINALG::Multiply(T_,false,Kaa_mod,false);
		
		if (gidofs->NumGlobalElements())
			TKai_mod = LINALG::Multiply(T_,false,Kai_mod,false);
		
		Tfa_mod = rcp(new Epetra_Vector(T_->RowMap()));
		T_->Multiply(false,*fa_mod,*Tfa_mod);
	}
	
	// output for checking everything
	if(Comm().MyPID()==0)
	{
		cout << endl << "*****************" << endl;
		if (TKan_mod!=null) cout << "TKan_mod: " << TKan_mod->NumGlobalRows() << " x " << TKan_mod->NumGlobalCols() << endl;
		else cout << "TKan_mod: null" << endl;
		if (TKam_mod!=null) cout << "TKam_mod: " << TKam_mod->NumGlobalRows() << " x " << TKam_mod->NumGlobalCols() << endl;
		else cout << "TKam_mod: null" << endl;
		if (TKai_mod!=null) cout << "TKai_mod: " << TKai_mod->NumGlobalRows() << " x " << TKai_mod->NumGlobalCols() << endl;
		else cout << "TKai_mod: null" << endl;
		if (TKaa_mod!=null) cout << "TKaa_mod: " << TKaa_mod->NumGlobalRows() << " x " << TKaa_mod->NumGlobalCols() << endl;
		else cout << "TKaa_mod: null" << endl;
		cout << "*****************" << endl;
		
		cout << endl << "**********" << endl;
		if (Tfa_mod!=null) cout << "Tfa_mod: " << Tfa_mod->GlobalLength() << endl;
		else cout << "Tfa_mod: null" << endl;
		cout << "**********" << endl;
	}
	
	/**********************************************************************/
	/* Global setup of Kteff_new, feff_new (including contact)            */
	/**********************************************************************/
	RCP<Epetra_CrsMatrix> Kteff_new = LINALG::CreateMatrix(*(discret_.DofRowMap()),81);
	RCP<Epetra_Vector> feff_new = LINALG::CreateVector(*(discret_.DofRowMap()));
	
	// add n / m submatrices to Kteff_new
	LINALG::Add(Knn_mod,false,1.0,Kteff_new,1.0);
	LINALG::Add(Knm_mod,false,1.0,Kteff_new,1.0);
	LINALG::Add(Kmn_mod,false,1.0,Kteff_new,1.0);
	LINALG::Add(Kmm_mod,false,1.0,Kteff_new,1.0);
	
	// add a / i submatrices to Kteff_new, if existing
	if (Kns_mod!=null) LINALG::Add(Kns_mod,false,1.0,Kteff_new,1.0);
	if (Kms_mod!=null) LINALG::Add(Kms_mod,false,1.0,Kteff_new,1.0);
	if (Kin_mod!=null) LINALG::Add(Kin_mod,false,1.0,Kteff_new,1.0);
	if (Kim_mod!=null) LINALG::Add(Kim_mod,false,1.0,Kteff_new,1.0);
	if (Kii_mod!=null) LINALG::Add(Kii_mod,false,1.0,Kteff_new,1.0);
	if (Kia_mod!=null) LINALG::Add(Kia_mod,false,1.0,Kteff_new,1.0);
	
	// add matrix of normals to Kteff_new
	LINALG::Add(N_,false,1.0,Kteff_new,1.0);

	// add submatrices with tangents to Kteff_new, if existing
	if (TKan_mod!=null) LINALG::Add(TKan_mod,false,1.0,Kteff_new,1.0);
	if (TKam_mod!=null) LINALG::Add(TKam_mod,false,1.0,Kteff_new,1.0);
	if (TKai_mod!=null) LINALG::Add(TKai_mod,false,1.0,Kteff_new,1.0);
	if (TKaa_mod!=null) LINALG::Add(TKaa_mod,false,1.0,Kteff_new,1.0);
	
	// FillComplete Kteff_new (square)
	LINALG::Complete(Kteff_new);
	
	// add n / m subvectors to feff_new
	RCP<Epetra_Vector> fn_mod_exp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
	RCP<Epetra_Vector> fm_mod_exp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
	LINALG::Export(*fn_mod,*fn_mod_exp);
	LINALG::Export(*fm_mod,*fm_mod_exp);
	feff_new->Update(1.0,*fn_mod_exp,1.0,*fm_mod_exp,1.0);
	
	// add i / Ta subvectors to feff_new, if existing
	RCP<Epetra_Vector> fi_mod_exp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
	RCP<Epetra_Vector> Tfa_mod_exp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
	if (fi_mod!=null) LINALG::Export(*fi_mod,*fi_mod_exp);
	if (Tfa_mod!=null) LINALG::Export(*Tfa_mod,*Tfa_mod_exp);
	feff_new->Update(1.0,*fi_mod_exp,1.0,*Tfa_mod_exp,1.0);
	
	// add weighted gap vector to feff_new, if existing
	RCP<Epetra_Vector> g_exp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
	if (g_act->GlobalLength()) LINALG::Export(*g_act,*g_exp);
	feff_new->Update(1.0,*g_exp,1.0);
	
	/**********************************************************************/
	/* Replace Kteff and feff by Kteff_new and feff_new                   */
	/**********************************************************************/
	Kteff = Kteff_new;
	feff = feff_new;
	LINALG::PrintSparsityToPostscript(*Kteff);
	//exit(0);
	
  return;
}

/*----------------------------------------------------------------------*
 |  Transform displacement increment vector (public)          popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::Manager::RecoverDisp(RCP<Epetra_Vector>& disi)
{	
	// extract master displacements from disi
	RCP<Epetra_Vector> disi_m = rcp(new Epetra_Vector(*gmdofrowmap_));
	LINALG::Export(*disi,*disi_m);
	
	// recover slave displacement increments
	RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(Mbar_->RowMap()));
	Mbar_->Multiply(false,*disi_m,*mod);
	
	RCP<Epetra_Vector> mod_exp = rcp(new Epetra_Vector(*(discret_.DofRowMap())));
	LINALG::Export(*mod,*mod_exp);
	disi->Update(1.0,*mod_exp,1.0);
	
	return;
}

#endif  // #ifdef CCADISCRET
