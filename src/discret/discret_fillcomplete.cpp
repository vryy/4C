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

#include "discret.H"
#include "exporter.H"
#include "dserror.H"




/*----------------------------------------------------------------------*
 |  Finalize construction (public)                           mwgee 11/06|
 *----------------------------------------------------------------------*/
int CCADISCRETIZATION::Discretization::FillComplete()
{
  filled_ = false;  

  // (re)build map of nodes noderowmap_ and nodecolmap_
  BuildNodeRowMap();
  BuildNodeColMap();
  
  // (re)build map of elements elemap_
  BuildElementRowMap();
  BuildElementColMap();
  
  // (re)construct element -> node pointers
  BuildElementToNodePointers();

  // (re)construct node -> element pointers
  BuildNodeToElementPointers();
  
  filled_ = true;  
  return 0;
}


/*----------------------------------------------------------------------*
 |  Build noderowmap_ (public)                               mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::BuildNodeRowMap()
{
  const int myrank = Comm().MyPID();
  int nummynodes     = 0;
  int numglobalnodes = 0;
  map<int,RefCountPtr<CCADISCRETIZATION::Node> >::iterator curr;
  for (curr=node_.begin(); curr != node_.end(); ++curr)
    if (curr->second->Owner() == myrank)
      ++nummynodes;
  Comm().SumAll(&nummynodes,&numglobalnodes,1);
  vector<int> nodeids(nummynodes);
  
  int count=0;
  for (curr=node_.begin(); curr != node_.end(); ++curr)
    if (curr->second->Owner() == myrank)
    {
      nodeids[count] = curr->second->Id();
      ++count;
    }
  if (count != nummynodes) dserror("Mismatch in no. of nodes");
  noderowmap_ = rcp(new Epetra_Map(numglobalnodes,nummynodes,&nodeids[0],0,Comm()));
  return;
}

/*----------------------------------------------------------------------*
 |  Build noderowmap_ (public)                               mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::BuildNodeColMap()
{
  int nummynodes     = (int)node_.size();
  int numglobalnodes = 0;
  Comm().SumAll(&nummynodes,&numglobalnodes,1);
  vector<int> nodeids(nummynodes);
  
  int count=0;
  map<int,RefCountPtr<CCADISCRETIZATION::Node> >::iterator curr;
  for (curr=node_.begin(); curr != node_.end(); ++curr)
  {
    nodeids[count] = curr->second->Id();
    ++count;
  }
  if (count != nummynodes) dserror("Mismatch in no. of nodes");
  nodecolmap_ = rcp(new Epetra_Map(numglobalnodes,nummynodes,&nodeids[0],0,Comm()));
  return;
}


/*----------------------------------------------------------------------*
 |  Build elecolmap_ (public)                                mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::BuildElementColMap()
{
  int nummyeles     = (int)element_.size();
  int numglobaleles = 0;
  Comm().SumAll(&nummyeles,&numglobaleles,1);
  
  vector<int> eleids(nummyeles);
  
  map<int,RefCountPtr<CCADISCRETIZATION::Element> >::iterator curr;
  int count=0;
  for (curr=element_.begin(); curr != element_.end(); ++curr)
  {
    eleids[count] = curr->second->Id();
    ++count;
  }
  if (count != nummyeles) dserror("Mismatch in no. of elements");
  elecolmap_ = rcp(new Epetra_Map(numglobaleles,nummyeles,&eleids[0],0,Comm()));
  return;
}

/*----------------------------------------------------------------------*
 |  Build elecolmap_ (public)                                mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::BuildElementRowMap()
{
  const int myrank = Comm().MyPID();
  int nummyeles = 0;
  int numglobaleles = 0;
  map<int,RefCountPtr<CCADISCRETIZATION::Element> >::iterator curr;
  for (curr=element_.begin(); curr != element_.end(); ++curr)
    if (curr->second->Owner()==myrank)
      nummyeles++;

  Comm().SumAll(&nummyeles,&numglobaleles,1);
  
  vector<int> eleids(nummyeles);
  
  int count=0;
  for (curr=element_.begin(); curr != element_.end(); ++curr)
    if (curr->second->Owner()==myrank)
    {
      eleids[count] = curr->second->Id();
      ++count;
    }
  if (count != nummyeles) dserror("Mismatch in no. of elements");
  elerowmap_ = rcp(new Epetra_Map(numglobaleles,nummyeles,&eleids[0],0,Comm()));
  return;
}

/*----------------------------------------------------------------------*
 |  Build ptrs element -> node (public)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::BuildElementToNodePointers()
{
  map<int,RefCountPtr<CCADISCRETIZATION::Element> >::iterator elecurr;
  for (elecurr=element_.begin(); elecurr != element_.end(); ++elecurr)
  {
    bool success = elecurr->second->BuildNodalPointers(node_);
    if (!success)
      dserror("Building element <-> node topology failed");
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Build ptrs node -> element (public)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::BuildNodeToElementPointers()
{
  map<int,RefCountPtr<CCADISCRETIZATION::Node> >::iterator nodecurr;
  for (nodecurr=node_.begin(); nodecurr != node_.end(); ++nodecurr)
    nodecurr->second->element_.resize(0);
  
  map<int,RefCountPtr<CCADISCRETIZATION::Element> >::iterator elecurr;
  for (elecurr=element_.begin(); elecurr != element_.end(); ++elecurr)
  {
    const int  nnode = elecurr->second->NumNode();
    const int* nodes = elecurr->second->NodeIds();
    for (int j=0; j<nnode; ++j)
    {
      CCADISCRETIZATION::Node* node = gNode(nodes[j]);
      if (!node) dserror("Node %d is not on this proc %d",j,Comm().MyPID());
      else
      {
        int size = node->element_.size();
        node->element_.resize(size+1);
        node->element_[size] = elecurr->second.get();
      }
    }
  }
  return;
}





#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
