/*!----------------------------------------------------------------------
\file discret_partition.cpp
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
 |  Export nodes owned by a proc (public)                    mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::ExportNodes(const Epetra_Map& newmap)
{
  // destroy all ghosted nodes
  const int myrank = Comm().MyPID();
  map<int,RefCountPtr<CCADISCRETIZATION::Node> >::iterator curr;
  for (curr=node_.begin(); curr!=node_.end(); ++curr)
    if (curr->second->Owner() != myrank)
      node_.erase(curr->first);
  
  // build rowmap of nodes noderowmap_ if it does not exist
  if (noderowmap_==null) BuildNodeRowMap();
  const Epetra_Map& oldmap = *noderowmap_;
  
  // test whether newmap is non-overlapping
  int oldmy = oldmap.NumMyElements();
  int newmy = newmap.NumMyElements();
  int oldglobal=0;
  int newglobal=0;
  Comm().SumAll(&oldmy,&oldglobal,1);
  Comm().SumAll(&newmy,&newglobal,1);
  if (oldglobal != newglobal) dserror("New map is likely not non-overlapping");

  // create an exporter object that will figure out the communication pattern
  CCADISCRETIZATION::Exporter exporter(oldmap,newmap,Comm());
  // Do the communication
  exporter.Export(node_);
  
  // update all ownership flags
  for (curr=node_.begin(); curr!=node_.end(); ++curr)
    curr->second->SetOwner(myrank);
  
  // maps and pointers are no longer correct and need rebuilding
  filled_ = false;  
  return;
}

/*----------------------------------------------------------------------*
 |  Export nodes owned by a proc (public)                    mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::ExportGhostNodes(const Epetra_Map& newmap)
{
  // destroy all ghosted nodes
  const int myrank = Comm().MyPID();
  map<int,RefCountPtr<CCADISCRETIZATION::Node> >::iterator curr;
  for (curr=node_.begin(); curr!=node_.end(); ++curr)
    if (curr->second->Owner() != myrank)
      node_.erase(curr->first);
  
  // build rowmap of nodes noderowmap_ if it does not exist
  if (noderowmap_==null) BuildNodeRowMap();
  const Epetra_Map& oldmap = *noderowmap_;
  
  // test whether all nodes in oldmap are also in newmap, otherwise
  // this would be a change of owner which is not allowed here
  for (int i=0; i<oldmap.NumMyElements(); ++i)
  {
    int gid = oldmap.GID(i);
    if (!(newmap.MyGID(gid))) dserror("Node gid=%d from oldmap is not in newmap",gid);
  }
  
  // create an exporter object that will figure out the communication pattern
  CCADISCRETIZATION::Exporter exporter(oldmap,newmap,Comm());
  // Do the communication
  exporter.Export(node_);
  
  // maps and pointers are no longer correct and need rebuilding
  filled_ = false;  
  return;
}


/*----------------------------------------------------------------------*
 |  Export elements (public)                                 mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::ExportElements(const Epetra_Map& newmap)
{
  // destroy all ghosted elements
  const int myrank = Comm().MyPID();
  map<int,RefCountPtr<CCADISCRETIZATION::Element> >::iterator curr;
  for (curr=element_.begin(); curr!=element_.end(); ++curr)
    if (curr->second->Owner() != myrank)
      element_.erase(curr->first);
  
  // build map of elements elerowmap_ if it does not exist
  if (elerowmap_==null) BuildElementRowMap();
  const Epetra_Map& oldmap = *elerowmap_;
  
  // test whether newmap is non-overlapping
  int oldmy = oldmap.NumMyElements();
  int newmy = newmap.NumMyElements();
  int oldglobal=0;
  int newglobal=0;
  Comm().SumAll(&oldmy,&oldglobal,1);
  Comm().SumAll(&newmy,&newglobal,1);
  if (oldglobal != newglobal) dserror("New map is likely not non-overlapping");
  
  // create an exporter object that will figure out the communication pattern
  CCADISCRETIZATION::Exporter exporter(oldmap,newmap,Comm());
  exporter.Export(element_);
  
  // update ownerships
  for (curr=element_.begin(); curr!=element_.end(); ++curr)
    curr->second->SetOwner(myrank);
  
  // maps and pointers are no longer correct and need rebuilding
  filled_ = false;  
  return;
}


/*----------------------------------------------------------------------*
 |  Export elements (public)                                 mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::ExportGhostElements(const Epetra_Map& newmap)
{
  // destroy all ghosted elements
  const int myrank = Comm().MyPID();
  map<int,RefCountPtr<CCADISCRETIZATION::Element> >::iterator curr;
  for (curr=element_.begin(); curr!=element_.end(); ++curr)
    if (curr->second->Owner() != myrank)
      element_.erase(curr->first);
  
  // build map of elements elerowmap_ if it does not exist
  if (elerowmap_==null) BuildElementRowMap();
  const Epetra_Map& oldmap = *elerowmap_;
  
  // test whether all elements in oldmap are also in newmap
  // Otherwise, this would be a change of owner which is not allowed here
  for (int i=0; i<oldmap.NumMyElements(); ++i)
  {
    int gid = oldmap.GID(i);
    if (!(newmap.MyGID(gid))) dserror("Element gid=%d from oldmap is not in newmap",gid);
  }
  
  // create an exporter object that will figure out the communication pattern
  CCADISCRETIZATION::Exporter exporter(oldmap,newmap,Comm());
  exporter.Export(element_);
  
  // maps and pointers are no longer correct and need rebuilding
  filled_ = false;  
  return;
}

/*----------------------------------------------------------------------*
 |  build nodal graph from discretization (public)           mwgee 11/06|
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_CrsGraph> CCADISCRETIZATION::Discretization::BuildNodeGraph() const
{
  if (!Filled()) dserror("FillComplete() was not called on this discretization");

  // get nodal row map
  const Epetra_Map* noderowmap = NodeRowMap();
  
  // allocate graph
  RefCountPtr<Epetra_CrsGraph> graph = 
                     rcp( new Epetra_CrsGraph(Copy,*noderowmap,108,false));
  
  // iterate all elements on this proc including ghosted ones
  // Note:
  // if a proc stores the appropiate ghosted elements, the resulting
  // graph will be the correct and complete graph of the distributed
  // discretization even if nodes are not ghosted.
  map<int,RefCountPtr<CCADISCRETIZATION::Element> >::const_iterator curr;
  for (curr=element_.begin(); curr!=element_.end(); ++curr)
  {
    const int  nnode   = curr->second->NumNode();
    const int* nodeids = curr->second->NodeIds();
    for (int row=0; row<nnode; ++row)
    {
      const int rownode = nodeids[row];
      if (!noderowmap->MyGID(rownode)) continue;
      for (int col=0; col<nnode; ++col)
      {
        int colnode = nodeids[col];
        int err = graph->InsertGlobalIndices(rownode,1,&colnode);
        if (err<0) dserror("graph->InsertGlobalIndices returned err=%d",err);
      }
    }
  }
  int err = graph->FillComplete();
  if (err) dserror("graph->FillComplete() returned err=%d",err);
  err = graph->OptimizeStorage();
  if (err) dserror("graph->OptimizeStorage() returned err=%d",err);
  return graph;
}


/*----------------------------------------------------------------------*
 |  build element map from discretization (public)           mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::BuildElementGhosting(
                                    const Epetra_Map& noderowmap,
                                    const Epetra_Map& nodecolmap,
                                    RefCountPtr<Epetra_Map>& elerowmap,
                                    RefCountPtr<Epetra_Map>& elecolmap) const
{
  if (!Filled()) dserror("FillComplete() was not called on this discretization");

  // note: 
  // - noderowmap need not match distribution of nodes in this
  //   discretization at all.
  // - noderowmap is a non-overlapping map, that's tested
  int nummy = noderowmap.NumMyElements();
  int numall= 0;
  Comm().SumAll(&nummy,&numall,1);
  if (numall != NumGlobalNodes()) dserror("noderowmap is not a nodal row map or does not match no. of nodes in this discretization");
  
  // build connectivity of elements
  // storing :  element gid
  //            no. of nodes
  //            nodeids
  int stoposize = 2000;
  int count     = 0;
  vector<int> stopo(stoposize);
  map<int,RefCountPtr<CCADISCRETIZATION::Element> >::const_iterator ecurr;
  for (ecurr=element_.begin(); ecurr!=element_.end(); ++ecurr)
  {
    const CCADISCRETIZATION::Element& actele = *(ecurr->second);
    int        gid     = actele.Id();
    int        nnode   = actele.NumNode();
    const int* nodeids = actele.NodeIds();
    if (count+nnode+2>=stoposize)
    {
      stoposize += (count+nnode+2)*100;
      stopo.resize(stoposize);
    }
    stopo[count++] = gid;
    stopo[count++] = nnode;
    for (int j=0; j<nnode; ++j) stopo[count++] = nodeids[j];
  }
  stoposize = count;
  stopo.resize(stoposize);

  const int myrank = Comm().MyPID();
  const int numproc = Comm().NumProc();
  vector<int> rtopo(stoposize);
  vector<bool> mine(50);
  for (int proc=0; proc < numproc; ++proc)
  {
    int size = stoposize;
    Comm().Broadcast(&size,1,proc);
    if (size>(int)rtopo.size()) rtopo.resize(size);
    if (proc==myrank) 
      for (int i=0; i<size; ++i) rtopo[i] = stopo[i];
    Comm().Broadcast(&rtopo[0],size,proc);
    for (int i=0; i<size;)
    {
      const int elegid   = rtopo[i++];
      const int numnod   = rtopo[i++];
      const int* nodeids = &rtopo[i];
      i+= numnod;
      if (numnod>(int)mine.size()) mine.resize(numnod);
      for (int j=0; j<numnod; ++j)
        if (noderowmap.MyGID(nodeids[j]))
          mine[j] = true;
        else
          mine[j] = false;
      int nummine=0;
      for (int j=0; j<numnod; ++j)
        if (mine[j]) ++nummine;
      // - elegid is mine if connected to only my nodes
      // - elegid is mine if connected to some of my nodes
      //   and I own the majority of nodes
      // - elegid is mine if connected to some of my nodes
      //   and there is no majority of nodes but I have the lower rank
      // - elegid is ghosted by me of i own a minority of nodes
      if (nummine==numnod) // I am owner
      {
      }
      else if (nummine=0) // I am not owner and I do not ghost
      {
      }
      else if 
      
      
    }
  }



  return;
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
