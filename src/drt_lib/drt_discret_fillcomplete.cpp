/*!----------------------------------------------------------------------
\file drt_discret_fillcomplete.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_discret.H"
#include "drt_exporter.H"
#include "drt_dserror.H"




/*----------------------------------------------------------------------*
 |  Finalize construction (public)                           mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Reset(bool killdofs)
{
  filled_ = false;
  if (killdofs)
  {
    havedof_= false;
    for (unsigned i=0; i<dofsets_.size(); ++i)
      dofsets_[i]->Reset();
  }
  dofrowmap_ = null;
  dofcolmap_ = null;
  elerowmap_ = null;
  elecolmap_ = null;
  elerowptr_.clear();
  elecolptr_.clear();
  noderowmap_ = null;
  nodecolmap_ = null;
  noderowptr_.clear();
  nodecolptr_.clear();
  return;
}


/*----------------------------------------------------------------------*
 |  Finalize construction (public)                           mwgee 11/06|
 *----------------------------------------------------------------------*/
int DRT::Discretization::FillComplete(bool assigndegreesoffreedom,
                                      bool initelements,
                                      bool doboundaryconditions)
{
  // set all maps to null
  Reset(assigndegreesoffreedom);

  // (re)build map of nodes noderowmap_, nodecolmap_, noderowptr and nodecolptr
  BuildNodeRowMap();
  BuildNodeColMap();

  // (re)build map of elements elemap_
  BuildElementRowMap();
  BuildElementColMap();

  // (re)construct element -> node pointers
  BuildElementToNodePointers();

  // (re)construct node -> element pointers
  BuildNodeToElementPointers();

#ifdef ELEMINTERFACES
  // bos 12/07
  // (re)construct element -> element pointers for interface-elements 
  BuildElementToElementPointers();
#endif // ELEMINTERFACES
  
  // set the flag indicating Filled()==true
  // as the following methods make use of maps
  // which we just built
  filled_ = true;

  // Assign degrees of freedom to elements and nodes
  if (assigndegreesoffreedom) AssignDegreesOfFreedom(0);

  if (initelements)
  {
    // build the register of elements
    BuildElementRegister();
    // call element routines to initialize
    InitializeElements();
  }
  
  // (Re)build the geometry of the boundary conditions
  if (doboundaryconditions) BoundaryConditionsGeometry();

  return 0;
}


/*----------------------------------------------------------------------*
 |  init elements (public)                                   mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::InitializeElements()
{
  if (!Filled()) dserror("FillComplete was not called");

  map<int,RefCountPtr<ElementRegister> >::iterator fool;
  for (fool=elementregister_.begin(); fool!=elementregister_.end(); ++fool)
  {
    int err = fool->second->Initialize(*this);
    if (err) dserror("Element Initialize returned err=%d",err);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Build elementregister_ (private)                         mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BuildElementRegister()
{
  const int myrank = Comm().MyPID();
  const int numproc = Comm().NumProc();

  // clear any existing data in register
  elementregister_.clear();

  // create a temporary map for communication
  // this maps Element::Type() to ElementRegister
  map<int,RefCountPtr<DRT::ParObject> > tmpmap;

  // loop my row elements and build local register
  vector<int> mygid(500);
  map<int,RefCountPtr<DRT::Element> >::iterator fool;
  map<int,RefCountPtr<DRT::ParObject> >::iterator rcurr;
  int nummygid=0;
  for (fool=element_.begin(); fool!=element_.end(); ++fool)
  {
    if (fool->second->Owner()!=myrank) continue;
    DRT::Element* actele = fool->second.get();
    rcurr = tmpmap.find(actele->Type());
    if (rcurr != tmpmap.end()) continue;
    RefCountPtr<DRT::ElementRegister> tmp = actele->ElementRegister();
    tmpmap[actele->Type()] = tmp;
    if ((int)mygid.size()<=nummygid) mygid.resize(nummygid+100);
    mygid[nummygid] = actele->Type();
    ++nummygid;
  }

  // build the source map
  Epetra_Map sourcemap(-1,nummygid,&mygid[0],0,Comm());

  // build redundant target map
  int numglobalgid = sourcemap.NumGlobalElements();
  vector<int> redgid(numglobalgid);
  for (int i=0; i<numglobalgid; ++i) redgid[i] = 0;
  int position=0;
  for (int proc=0; proc<numproc; ++proc)
  {
    if (proc==myrank)
    {
      for (int i=0; i<nummygid; ++i)
        redgid[position+i] = mygid[i];
      position += nummygid;
    }
    Comm().Broadcast(&position,1,proc);
  }
  vector<int> recvredgid(numglobalgid);
  Comm().SumAll(&redgid[0],&recvredgid[0],numglobalgid);
  Epetra_Map targetmap(-1,numglobalgid,&recvredgid[0],0,Comm());

  // create an exporter and export the tmpmap
  DRT::Exporter exporter(sourcemap,targetmap,Comm());
  exporter.Export(tmpmap);


#if 0
  for (int proc=0; proc<numproc; ++proc)
  {
    if (myrank==proc)
    {
      cout << "Proc " << myrank << endl;
      cout << "tmpmap.size() " << tmpmap.size() << endl;
      for (rcurr=tmpmap.begin(); rcurr!=tmpmap.end(); ++rcurr)
        cout << "rcurr->first " << rcurr->first
             << " rcurr->second->UniqueParObjectId() " << rcurr->second->UniqueParObjectId() << endl;
    }
    fflush(stdout);
    Comm().Barrier();
  }
#endif

  // go through the tmpmap and fill elementregister_
  map<int,RefCountPtr<ElementRegister> >::iterator regcurr;
  for (rcurr=tmpmap.begin(); rcurr!=tmpmap.end(); ++rcurr)
  {
    // cast the parobject to DRT::ElementRegister
    DRT::ElementRegister* elereg = dynamic_cast<DRT::ElementRegister*>(rcurr->second.get());
    if (!elereg) dserror("cast from ParObject to ElementRegister failed");
    // check whether its already in elementregister_
    regcurr = elementregister_.find(elereg->Type());
    if (regcurr!=elementregister_.end()) continue;
    // don't want tmpmap to destroy my ElementRegister class as we
    // will keep it in elementregister_
    // so we release the tmpmap refcountptr and create a new one
    // (this is actually to avoid type conversion)
    rcurr->second.release();
    elementregister_[elereg->Type()] = rcp(elereg);
  }


#if 0
  for (int proc=0; proc<numproc; ++proc)
  {
    if (myrank==proc)
    {
      cout << "Proc " << myrank << endl;
      cout << "elementregister_.size() " << elementregister_.size() << endl;
      for (regcurr=elementregister_.begin(); regcurr!=elementregister_.end(); ++regcurr)
        cout << "regcurr->first " << regcurr->first
             << " regcurr->second->UniqueParObjectId() " << regcurr->second->UniqueParObjectId() << endl;
    }
    fflush(stdout);
    Comm().Barrier();
  }
#endif

  return;
}

/*----------------------------------------------------------------------*
 |  Build noderowmap_ (private)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BuildNodeRowMap()
{
  const int myrank = Comm().MyPID();
  int nummynodes     = 0;
  map<int,RefCountPtr<DRT::Node> >::iterator curr;
  for (curr=node_.begin(); curr != node_.end(); ++curr)
    if (curr->second->Owner() == myrank)
      ++nummynodes;
  vector<int> nodeids(nummynodes);
  noderowptr_.resize(nummynodes);

  int count=0;
  for (curr=node_.begin(); curr != node_.end(); ++curr)
    if (curr->second->Owner() == myrank)
    {
      nodeids[count] = curr->second->Id();
      noderowptr_[count] = curr->second.get();
      ++count;
    }
  if (count != nummynodes) dserror("Mismatch in no. of nodes");
  noderowmap_ = rcp(new Epetra_Map(-1,nummynodes,&nodeids[0],0,Comm()));
  return;
}

/*----------------------------------------------------------------------*
 |  Build nodecolmap_ (private)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BuildNodeColMap()
{
  int nummynodes = (int)node_.size();
  vector<int> nodeids(nummynodes);
  nodecolptr_.resize(nummynodes);

  int count=0;
  map<int,RefCountPtr<DRT::Node> >::iterator curr;
  for (curr=node_.begin(); curr != node_.end(); ++curr)
  {
    nodeids[count] = curr->second->Id();
    nodecolptr_[count] = curr->second.get();
    curr->second->SetLID(count);
    ++count;
  }
  if (count != nummynodes) dserror("Mismatch in no. of nodes");
  nodecolmap_ = rcp(new Epetra_Map(-1,nummynodes,&nodeids[0],0,Comm()));
  return;
}


/*----------------------------------------------------------------------*
 |  Build elerowmap_ (private)                                mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BuildElementRowMap()
{
  const int myrank = Comm().MyPID();
  int nummyeles = 0;
  map<int,RefCountPtr<DRT::Element> >::iterator curr;
  for (curr=element_.begin(); curr != element_.end(); ++curr)
    if (curr->second->Owner()==myrank)
      nummyeles++;
  vector<int> eleids(nummyeles);
  elerowptr_.resize(nummyeles);
  int count=0;
  for (curr=element_.begin(); curr != element_.end(); ++curr)
    if (curr->second->Owner()==myrank)
    {
      eleids[count] = curr->second->Id();
      elerowptr_[count] = curr->second.get();
      ++count;
    }
  if (count != nummyeles) dserror("Mismatch in no. of elements");
  elerowmap_ = rcp(new Epetra_Map(-1,nummyeles,&eleids[0],0,Comm()));
  return;
}

/*----------------------------------------------------------------------*
 |  Build elecolmap_ (private)                                mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BuildElementColMap()
{
  int nummyeles = (int)element_.size();
  vector<int> eleids(nummyeles);
  elecolptr_.resize(nummyeles);
  map<int,RefCountPtr<DRT::Element> >::iterator curr;
  int count=0;
  for (curr=element_.begin(); curr != element_.end(); ++curr)
  {
    eleids[count] = curr->second->Id();
    elecolptr_[count] = curr->second.get();
    curr->second->SetLID(count);
    ++count;
  }
  if (count != nummyeles) dserror("Mismatch in no. of elements");
  elecolmap_ = rcp(new Epetra_Map(-1,nummyeles,&eleids[0],0,Comm()));
  return;
}

/*----------------------------------------------------------------------*
 |  Build ptrs element -> node (private)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BuildElementToNodePointers()
{
  map<int,RefCountPtr<DRT::Element> >::iterator elecurr;
  for (elecurr=element_.begin(); elecurr != element_.end(); ++elecurr)
  {
    bool success = elecurr->second->BuildNodalPointers(node_);
    if (!success)
      dserror("Building element <-> node topology failed");
  }
  return;
}

#ifdef ELEMINTERFACES
/*----------------------------------------------------------------------*
 |  Build ptrs element -> element (private)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BuildElementToElementPointers()
{
  map<int,RefCountPtr<DRT::Element> >::iterator elecurr;
  for (elecurr=element_.begin(); elecurr != element_.end(); ++elecurr)
  {
    bool success = elecurr->second->BuildElementPointers(element_);
    if (!success)
      dserror("Building element <-> element topology failed");
  }
  return;
}
#endif //ELEMINTERFACES

/*----------------------------------------------------------------------*
 |  Build ptrs node -> element (private)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BuildNodeToElementPointers()
{
  map<int,RefCountPtr<DRT::Node> >::iterator nodecurr;
  for (nodecurr=node_.begin(); nodecurr != node_.end(); ++nodecurr)
    nodecurr->second->ClearMyElementTopology();

  map<int,RefCountPtr<DRT::Element> >::iterator elecurr;
  for (elecurr=element_.begin(); elecurr != element_.end(); ++elecurr)
  {
    const int  nnode = elecurr->second->NumNode();
    const int* nodes = elecurr->second->NodeIds();
    for (int j=0; j<nnode; ++j)
    {
      DRT::Node* node = gNode(nodes[j]);
      if (!node) dserror("Node %d is not on this proc %d",j,Comm().MyPID());
      else node->AddElementPtr(elecurr->second.get());
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  set degrees of freedom (public)                          mwgee 03/07|
 *----------------------------------------------------------------------*/
int DRT::Discretization::AssignDegreesOfFreedom(const int start)
{
  if (!Filled()) dserror("Filled()==false");
  if (!NodeRowMap()->UniqueGIDs()) dserror("Nodal row map is not unique");
  if (!ElementRowMap()->UniqueGIDs()) dserror("Element row map is not unique");

  havedof_ = true;
  return dofsets_[currentdofset_]->AssignDegreesOfFreedom(*this,start);
}


#endif  // #ifdef CCADISCRET
