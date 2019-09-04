/*---------------------------------------------------------------------*/
/*! \file

\brief Setup of discretization including assignment of degrees of freedom

\level 0

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#include "drt_discret.H"
#include "drt_exporter.H"
#include "drt_dserror.H"
#include "drt_parobjectfactory.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io_pstream.H"



/*----------------------------------------------------------------------*
 |  Finalize construction (public)                           mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Reset(bool killdofs, bool killcond)
{
  filled_ = false;
  if (killdofs)
  {
    havedof_ = false;
    for (unsigned i = 0; i < dofsets_.size(); ++i) dofsets_[i]->Reset();
  }

  elerowmap_ = Teuchos::null;
  elecolmap_ = Teuchos::null;
  elerowptr_.clear();
  elecolptr_.clear();
  noderowmap_ = Teuchos::null;
  nodecolmap_ = Teuchos::null;
  noderowptr_.clear();
  nodecolptr_.clear();

  // delete all old geometries that are attached to any conditions
  // as early as possible
  if (killcond)
  {
    std::multimap<std::string, Teuchos::RCP<DRT::Condition>>::iterator fool;
    for (fool = condition_.begin(); fool != condition_.end(); ++fool)
    {
      fool->second->ClearGeometry();
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Finalize construction (public)                           mwgee 11/06|
 *----------------------------------------------------------------------*/
int DRT::Discretization::FillComplete(
    bool assigndegreesoffreedom, bool initelements, bool doboundaryconditions)
{
  // my processor id
  int myrank = Comm().MyPID();

  // print information to screen
  if (myrank == 0)
  {
    IO::cout(IO::verbose)
        << "\n+--------------------------------------------------------------------+" << IO::endl;
    IO::cout(IO::verbose) << "| FillComplete() on discretization " << std::setw(34) << std::left
                          << Name() << std::setw(1) << std::right << "|" << IO::endl;
  }

  // set all maps to Teuchos::null
  Reset(assigndegreesoffreedom, doboundaryconditions);

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
  if (assigndegreesoffreedom)
  {
    if (myrank == 0)
    {
      IO::cout(IO::verbose)
          << "| AssignDegreesOfFreedom() ...                                       |" << IO::endl;
    }
    AssignDegreesOfFreedom(0);
  }

  // call element routines to initialize
  if (initelements)
  {
    if (myrank == 0)
    {
      IO::cout(IO::verbose)
          << "| InitializeElements() ...                                           |" << IO::endl;
    }
    InitializeElements();
  }

  // (Re)build the geometry of the boundary conditions
  if (doboundaryconditions)
  {
    if (myrank == 0)
    {
      IO::cout(IO::verbose)
          << "| BoundaryConditionsGeometry() ...                                   |" << IO::endl;
    }

    BoundaryConditionsGeometry();
  }

  if (myrank == 0)
  {
    IO::cout(IO::verbose)
        << "+--------------------------------------------------------------------+" << IO::endl;
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  init elements (public)                                   mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::InitializeElements()
{
  if (!Filled()) dserror("FillComplete was not called");

  ParObjectFactory::Instance().InitializeElements(*this);

  return;
}


/*----------------------------------------------------------------------*
 |  Build noderowmap_ (private)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BuildNodeRowMap()
{
  const int myrank = Comm().MyPID();
  int nummynodes = 0;
  std::map<int, Teuchos::RCP<DRT::Node>>::iterator curr;
  for (curr = node_.begin(); curr != node_.end(); ++curr)
    if (curr->second->Owner() == myrank) ++nummynodes;
  std::vector<int> nodeids(nummynodes);
  noderowptr_.resize(nummynodes);

  int count = 0;
  for (curr = node_.begin(); curr != node_.end(); ++curr)
    if (curr->second->Owner() == myrank)
    {
      nodeids[count] = curr->second->Id();
      noderowptr_[count] = curr->second.get();
      ++count;
    }
  if (count != nummynodes) dserror("Mismatch in no. of nodes");
  noderowmap_ = Teuchos::rcp(new Epetra_Map(-1, nummynodes, &nodeids[0], 0, Comm()));
  return;
}

/*----------------------------------------------------------------------*
 |  Build nodecolmap_ (private)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BuildNodeColMap()
{
  int nummynodes = (int)node_.size();
  std::vector<int> nodeids(nummynodes);
  nodecolptr_.resize(nummynodes);

  int count = 0;
  std::map<int, Teuchos::RCP<DRT::Node>>::iterator curr;
  for (curr = node_.begin(); curr != node_.end(); ++curr)
  {
    nodeids[count] = curr->second->Id();
    nodecolptr_[count] = curr->second.get();
    curr->second->SetLID(count);
    ++count;
  }
  if (count != nummynodes) dserror("Mismatch in no. of nodes");
  nodecolmap_ = Teuchos::rcp(new Epetra_Map(-1, nummynodes, &nodeids[0], 0, Comm()));
  return;
}


/*----------------------------------------------------------------------*
 |  Build elerowmap_ (private)                                mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BuildElementRowMap()
{
  const int myrank = Comm().MyPID();
  int nummyeles = 0;
  std::map<int, Teuchos::RCP<DRT::Element>>::iterator curr;
  for (curr = element_.begin(); curr != element_.end(); ++curr)
    if (curr->second->Owner() == myrank) nummyeles++;
  std::vector<int> eleids(nummyeles);
  elerowptr_.resize(nummyeles);
  int count = 0;
  for (curr = element_.begin(); curr != element_.end(); ++curr)
    if (curr->second->Owner() == myrank)
    {
      eleids[count] = curr->second->Id();
      elerowptr_[count] = curr->second.get();
      ++count;
    }
  if (count != nummyeles) dserror("Mismatch in no. of elements");
  elerowmap_ = Teuchos::rcp(new Epetra_Map(-1, nummyeles, &eleids[0], 0, Comm()));
  return;
}

/*----------------------------------------------------------------------*
 |  Build elecolmap_ (private)                                mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BuildElementColMap()
{
  int nummyeles = (int)element_.size();
  std::vector<int> eleids(nummyeles);
  elecolptr_.resize(nummyeles);
  std::map<int, Teuchos::RCP<DRT::Element>>::iterator curr;
  int count = 0;
  for (curr = element_.begin(); curr != element_.end(); ++curr)
  {
    eleids[count] = curr->second->Id();
    elecolptr_[count] = curr->second.get();
    curr->second->SetLID(count);
    ++count;
  }
  if (count != nummyeles) dserror("Mismatch in no. of elements");
  elecolmap_ = Teuchos::rcp(new Epetra_Map(-1, nummyeles, &eleids[0], 0, Comm()));
  return;
}

/*----------------------------------------------------------------------*
 |  Build ptrs element -> node (private)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BuildElementToNodePointers()
{
  std::map<int, Teuchos::RCP<DRT::Element>>::iterator elecurr;
  for (elecurr = element_.begin(); elecurr != element_.end(); ++elecurr)
  {
    bool success = elecurr->second->BuildNodalPointers(node_);
    if (!success) dserror("Building element <-> node topology failed");
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Build ptrs element -> element (private)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BuildElementToElementPointers()
{
  std::map<int, Teuchos::RCP<DRT::Element>>::iterator elecurr;
  for (elecurr = element_.begin(); elecurr != element_.end(); ++elecurr)
  {
    bool success = elecurr->second->BuildElementPointers(element_);
    if (!success) dserror("Building element <-> element topology failed");
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Build ptrs node -> element (private)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::BuildNodeToElementPointers()
{
  std::map<int, Teuchos::RCP<DRT::Node>>::iterator nodecurr;
  for (nodecurr = node_.begin(); nodecurr != node_.end(); ++nodecurr)
    nodecurr->second->ClearMyElementTopology();

  std::map<int, Teuchos::RCP<DRT::Element>>::iterator elecurr;
  for (elecurr = element_.begin(); elecurr != element_.end(); ++elecurr)
  {
    const int nnode = elecurr->second->NumNode();
    const int* nodes = elecurr->second->NodeIds();
    for (int j = 0; j < nnode; ++j)
    {
      DRT::Node* node = gNode(nodes[j]);
      if (!node)
        dserror("Node %d is not on this proc %d", j, Comm().MyPID());
      else
        node->AddElementPtr(elecurr->second.get());
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

  for (unsigned i = 0; i < dofsets_.size(); ++i)
    start = dofsets_[i]->AssignDegreesOfFreedom(*this, i, start);
  return start;
}
