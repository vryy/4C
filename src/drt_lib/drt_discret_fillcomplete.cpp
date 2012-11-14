/*!----------------------------------------------------------------------
\file drt_discret_fillcomplete.cpp
\brief

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

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

#include "drt_discret.H"
#include "drt_exporter.H"
#include "drt_dserror.H"
#include "drt_parobjectfactory.H"
#include "../linalg/linalg_utils.H"




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

  elerowmap_ = null;
  elecolmap_ = null;
  elerowptr_.clear();
  elecolptr_.clear();
  noderowmap_ = null;
  nodecolmap_ = null;
  noderowptr_.clear();
  nodecolptr_.clear();

  // delete all old geometries that are attached to any conditions
  // as early as possible
  multimap<string,RCP<DRT::Condition> >::iterator fool;
  for (fool=condition_.begin(); fool != condition_.end(); ++fool)
  {
    fool->second->ClearGeometry();
  }

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

  // bos 12/07
  // (re)construct element -> element pointers for interface-elements
  BuildElementToElementPointers();

  // set the flag indicating Filled()==true
  // as the following methods make use of maps
  // which we just built
  filled_ = true;

  // Assign degrees of freedom to elements and nodes
  if (assigndegreesoffreedom) AssignDegreesOfFreedom(0);

  if (initelements)
  {
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

  ParObjectFactory::Instance().InitializeElements( *this );

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
  noderowmap_ = Teuchos::rcp(new Epetra_Map(-1,nummynodes,&nodeids[0],0,Comm()));
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
  nodecolmap_ = Teuchos::rcp(new Epetra_Map(-1,nummynodes,&nodeids[0],0,Comm()));
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
  elerowmap_ = Teuchos::rcp(new Epetra_Map(-1,nummyeles,&eleids[0],0,Comm()));
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
  elecolmap_ = Teuchos::rcp(new Epetra_Map(-1,nummyeles,&eleids[0],0,Comm()));
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
 |  set degrees of freedom (protected)                       mwgee 03/07|
 *----------------------------------------------------------------------*/
int DRT::Discretization::AssignDegreesOfFreedom(int start)
{
  if (!Filled()) dserror("Filled()==false");
  if (!NodeRowMap()->UniqueGIDs()) dserror("Nodal row map is not unique");
  if (!ElementRowMap()->UniqueGIDs()) dserror("Element row map is not unique");

  // Set the havedof flag before dofs are assigned. Some dof set
  // implementations do query the discretization after the assignment has been
  // done and this query demands the havedof flag to be set. An unexpected
  // implicit dependency here.
  havedof_ = true;

  for (unsigned i=0; i<dofsets_.size(); ++i)
    start = dofsets_[i]->AssignDegreesOfFreedom(*this,i,start);
  return start;
}


