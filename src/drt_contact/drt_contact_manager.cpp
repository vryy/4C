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
void CONTACT::Manager::Evaluate(Epetra_CrsMatrix& Kteff,
																Epetra_Vector& feff)
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

	/*
#ifdef DEBUG
	cout << *D_;
	cout << *M_;
	cout << *g_;
	cout << *N_;
	cout << *T_;
	cout << *invD;
	cout << *Mbar;

	cout << *M_;
	cout << M_->RangeMap();
	cout << M_->DomainMap();
	RCP<Epetra_CrsMatrix> transM = LINALG::Transpose(*M_);
	cout << *transM;
	cout << transM->RangeMap();
	cout << transM->DomainMap();
#endif // #ifdef DEBUG
	*/

  return;
}

#endif  // #ifdef CCADISCRET
