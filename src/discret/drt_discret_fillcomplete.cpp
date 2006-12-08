/*!----------------------------------------------------------------------
\file discret_fillcomplete.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "drt_discret.H"
#include "drt_exporter.H"
#include "drt_dserror.H"




/*----------------------------------------------------------------------*
 |  Finalize construction (public)                           mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Reset()
{
  filled_ = false;
  havedof_=false;
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
int DRT::Discretization::FillComplete()
{
  // set all maps to null
  Reset();  

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
  
  // build the register of elements
  BuildElementRegister();
  
  // set the flag indicating Filled()==true
  filled_ = true;  
  
  return 0;
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
  
  // loop my row elements and build local register
  map<int,RefCountPtr<DRT::Element> >::iterator fool;
  for (fool=element_.begin(); fool!=element_.end(); ++fool)
  {
    if (fool->second->Owner()!=myrank) continue;
    DRT::Element* actele = fool->second.get();
    map<Element::ElementType,RefCountPtr<ElementRegister> >::iterator rcurr;
    rcurr = elementregister_.find(actele->Type());
    if (rcurr != elementregister_.end()) continue;
    elementregister_[actele->Type()] = actele->ElementRegister();
  }


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





#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
