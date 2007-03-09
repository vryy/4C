/*!----------------------------------------------------------------------
\file drt_discret_partition.cpp
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
 |  Export nodes owned by a proc (public)                    mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::ExportRowNodes(const Epetra_Map& newmap)
{
  // test whether newmap is non-overlapping
  if (!newmap.UniqueGIDs()) dserror("new map not unique");

  // destroy all ghosted nodes
  const int myrank = Comm().MyPID();
  map<int,RefCountPtr<DRT::Node> >::iterator curr;
  for (curr=node_.begin(); curr!=node_.end(); ++curr)
    if (curr->second->Owner() != myrank)
      node_.erase(curr->first);
  
  // build rowmap of nodes noderowmap_ if it does not exist
  if (noderowmap_==null) BuildNodeRowMap();
  const Epetra_Map& oldmap = *noderowmap_;
  
  // not testing this anymore to allow piecewise exports in jumbo mode input
  //if (newmap.NumGlobalElements() != oldmap.NumGlobalElements())
  //  dserror("Global number of nodes in new and old map does not match");

  // create an exporter object that will figure out the communication pattern
  DRT::Exporter exporter(oldmap,newmap,Comm());
  // Do the communication
  exporter.Export(node_);
  
  // update all ownership flags
  for (curr=node_.begin(); curr!=node_.end(); ++curr)
    curr->second->SetOwner(myrank);
  
  // maps and pointers are no longer correct and need rebuilding
  Reset();

  return;
}

/*----------------------------------------------------------------------*
 |  Export nodes owned by a proc (public)                    mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::ExportColumnNodes(const Epetra_Map& newmap)
{
  // destroy all ghosted nodes
  const int myrank = Comm().MyPID();
  map<int,RefCountPtr<DRT::Node> >::iterator curr;
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
    if (!(newmap.MyGID(gid))) dserror("Proc %d: Node gid=%d from oldmap is not in newmap",myrank,gid);
  }
  
  // create an exporter object that will figure out the communication pattern
  DRT::Exporter exporter(oldmap,newmap,Comm());
  // Do the communication
  exporter.Export(node_);
  
  // maps and pointers are no longer correct and need rebuilding
  Reset();

  return;
}


/*----------------------------------------------------------------------*
 |  Export elements (public)                                 mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::ExportRowElements(const Epetra_Map& newmap)
{
  // destroy all ghosted elements
  const int myrank = Comm().MyPID();
  map<int,RefCountPtr<DRT::Element> >::iterator curr;
  for (curr=element_.begin(); curr!=element_.end(); ++curr)
    if (curr->second->Owner() != myrank)
      element_.erase(curr);

  // build map of elements elerowmap_ if it does not exist
  if (elerowmap_==null) BuildElementRowMap();
  const Epetra_Map& oldmap = *elerowmap_;


  // don't do this test anymore to allow to use this function while
  // the map of elements is still incomplete in the construction phase of
  // the discretization. This is used when creating the discretization in
  // jumbo mode
#if 0
  // test whether newmap is non-overlapping
  int oldmy = oldmap.NumMyElements();
  int newmy = newmap.NumMyElements();
  int oldglobal=0;
  int newglobal=0;
  Comm().SumAll(&oldmy,&oldglobal,1);
  Comm().SumAll(&newmy,&newglobal,1);

  if (oldglobal != newglobal) dserror("New map is likely not non-overlapping");
#endif  
  
  // create an exporter object that will figure out the communication pattern
  DRT::Exporter exporter(oldmap,newmap,Comm());
  exporter.Export(element_);
  
  // update ownerships and kick out everything that's not in newmap
  for (curr=element_.begin(); curr!=element_.end(); ++curr)
    curr->second->SetOwner(myrank);
  
  // maps and pointers are no longer correct and need rebuilding
  Reset();

  return;
}


/*----------------------------------------------------------------------*
 |  Export elements (public)                                 mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::ExportColumnElements(const Epetra_Map& newmap)
{
  // destroy all ghosted elements
  const int myrank = Comm().MyPID();
  map<int,RefCountPtr<DRT::Element> >::iterator curr;
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
    if (!(newmap.MyGID(gid))) dserror("Proc %d: Element gid=%d from oldmap is not in newmap",myrank,gid);
  }
  
  // create an exporter object that will figure out the communication pattern
  DRT::Exporter exporter(oldmap,newmap,Comm());
  exporter.Export(element_);
  
  // maps and pointers are no longer correct and need rebuilding
  Reset();

  return;
}

/*----------------------------------------------------------------------*
 |  build nodal graph from discretization (public)           mwgee 11/06|
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_CrsGraph> DRT::Discretization::BuildNodeGraph() const
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
  map<int,RefCountPtr<DRT::Element> >::const_iterator curr;
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
void DRT::Discretization::BuildElementRowColumn(
                                    const Epetra_Map& noderowmap,
                                    const Epetra_Map& nodecolmap,
                                    RefCountPtr<Epetra_Map>& elerowmap,
                                    RefCountPtr<Epetra_Map>& elecolmap) const
{
  const int myrank = Comm().MyPID();
  const int numproc = Comm().NumProc();

  // note: 
  // - noderowmap need not match distribution of nodes in this
  //   discretization at all.
  // - noderowmap is a non-overlapping map, that's tested
  if (!noderowmap.UniqueGIDs()) dserror("noderowmap is not a unique map");
  
  // find all owners for the overlapping node map
  const int ncnode = nodecolmap.NumMyElements();
  vector<int> cnodeowner(ncnode);
  {
    vector<int> lids(ncnode);
    noderowmap.RemoteIDList(ncnode,nodecolmap.MyGlobalElements(),&cnodeowner[0],&lids[0]);
    lids.clear();
  }
  
  // build connectivity of elements
  // storing :  element gid
  //            no. of nodes
  //            nodeids
  int stoposize = 2000;
  int count     = 0;
  vector<int> stopo(stoposize);
  map<int,RefCountPtr<DRT::Element> >::const_iterator ecurr;
  for (ecurr=element_.begin(); ecurr!=element_.end(); ++ecurr)
  {
    const DRT::Element& actele = *(ecurr->second);
    int        gid     = actele.Id();
    int        nnode   = actele.NumNode();
    const int* nodeids = actele.NodeIds();
    if (count+nnode+2>=stoposize)
    {
      stoposize += (nnode+2)*300;
      stopo.resize(stoposize);
    }
    stopo[count++] = gid;
    stopo[count++] = nnode;
    for (int j=0; j<nnode; ++j) stopo[count++] = nodeids[j];
  }
  stoposize = count;
  stopo.resize(stoposize);

  vector<int> rtopo(stoposize);
  
  // estimate no. of elements equal to no. of nodes
  vector<int> myele(noderowmap.NumMyElements());
  int nummyele=0;
  // estimate no. of ghosted elements much lower
  vector<int> myghostele(noderowmap.NumMyElements()/4);
  int nummyghostele=0;
  
  // loop processors and sort elements into 
  // elements owned by a proc
  // elements ghosted by a proc
  for (int proc=0; proc<numproc; ++proc)
  {
    int size = stoposize;
    Comm().Broadcast(&size,1,proc);
    if (size>(int)rtopo.size()) rtopo.resize(size);
    if (proc==myrank) 
      for (int i=0; i<size; ++i) rtopo[i] = stopo[i];
    Comm().Broadcast(&rtopo[0],size,proc);
    for (int i=0; i<size;)
    {
      const int  elegid  = rtopo[i++];
      const int  numnode = rtopo[i++];
      const int* nodeids = &rtopo[i];
      i += numnode;
          
      // resize arrays
      if (nummyele>=(int)myele.size()) myele.resize(myele.size()+500);
      if (nummyghostele>=(int)myghostele.size()) myghostele.resize(myghostele.size()+500);

      // count nodes I own of this element
      int nummine=0;
      for (int j=0; j<numnode; ++j)
        if (noderowmap.MyGID(nodeids[j]))
          ++nummine;
      
      // if I do not own any of the nodes, it is definitely not my element
      // and I do not ghost it
      if (!nummine) 
        continue;
      
      // check whether I ghost all nodes of this element
      // this is neccessary to be able to own or ghost the element
      for (int j=0; j<numnode; ++j)
        if (!nodecolmap.MyGID(nodeids[j]))
          dserror("I do not have own/ghosted node gid=%d",nodeids[j]);
      
      // find out who owns how many of the nodes
      vector<int> nodeowner(numnode);
      vector<int> numperproc(numproc);
      for (int j=0; j<numproc; ++j) numperproc[j] = 0;
      for (int j=0; j<numnode; ++j)
      {
        const int lid   = nodecolmap.LID(nodeids[j]);
        const int owner = cnodeowner[lid];
        nodeowner[j] = owner;
        numperproc[owner]++;
      }
      
      // the proc with the largest number of nodes owns the element, 
      // all others ghost it
      // if no. of nodes is equal among some procs, 
      // the higher rank owns the element
      int owner   = -1;
      int maxnode = 0;
      for (int j=0; j<numnode; ++j)
      {
        int currentproc = nodeowner[j];
        int ownhowmany  = numperproc[currentproc];
        if (ownhowmany>=maxnode)
        {
          owner   = currentproc;
          maxnode = ownhowmany;
        }
      }
      if (myrank==owner)
      {
        myele[nummyele++] = elegid;
        continue;
      }
      else
      {
        myghostele[nummyghostele++] = elegid;
        continue;
      }
      dserror("Error in logic of element ownerships");
      
    } // for (int i=0; i<size;)
  } // for (int proc=0; proc<numproc; ++proc)

  // at this point we have
  // myele, length nummyele
  // myghostele, length nummyghostele
  myele.resize(nummyele);
  myghostele.resize(nummyghostele);
  
  // allreduced nummyele must match the total no. of elements in this
  // discretization, otherwise we lost some
  // build the rowmap of elements
  elerowmap = rcp(new Epetra_Map(-1,nummyele,&myele[0],0,Comm()));
  if (!elerowmap->UniqueGIDs())
    dserror("Element row map is not unique");

  // build elecolmap
  vector<int> elecol(nummyele+nummyghostele);
  for (int i=0; i<nummyele; ++i) elecol[i] = myele[i];
  for (int i=0; i<nummyghostele; ++i) elecol[nummyele+i] = myghostele[i];
  elecolmap = rcp(new Epetra_Map(-1,nummyghostele+nummyele,
                                 &elecol[0],0,Comm()));

  return;
}




/*----------------------------------------------------------------------*
 |  redistribute discretization (public)                     mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Redistribute(const Epetra_Map& noderowmap,
                                       const Epetra_Map& nodecolmap)
{
  if (!Filled()) dserror("FillComplete() was not called on this discretization");
  
  // build the overlapping and non-overlapping element maps
  RefCountPtr<Epetra_Map> elerowmap;
  RefCountPtr<Epetra_Map> elecolmap;
  BuildElementRowColumn(noderowmap,nodecolmap,elerowmap,elecolmap);
  
  // export nodes and elements to the new maps
  ExportRowNodes(noderowmap);
  ExportColumnNodes(nodecolmap);
  ExportRowElements(*elerowmap);
  ExportColumnElements(*elecolmap);
  
  // these exports have set Filled()=false as all maps are invalid now
  int err = FillComplete();
  if (err) dserror("FillComplete() returned err=%d",err);

  return;
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
