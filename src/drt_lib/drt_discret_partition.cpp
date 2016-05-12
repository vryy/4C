/*!----------------------------------------------------------------------
\file drt_discret_partition.cpp

\brief Partition part of the setup of a discretization

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include <Epetra_FECrsGraph.h>
#include <Epetra_Time.h>

#include "drt_discret.H"
#include "drt_exporter.H"
#include "drt_dserror.H"
#include "drt_utils_metis.H"
#include "drt_dofset_pbc.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*
 |  Export nodes owned by a proc (public)                    mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::ExportRowNodes(const Epetra_Map& newmap, bool killdofs, bool killcond)
{
  // test whether newmap is non-overlapping
  if (!newmap.UniqueGIDs()) dserror("new map not unique");

  // destroy all ghosted nodes
  const int myrank = Comm().MyPID();
  std::map<int,Teuchos::RCP<DRT::Node> >::iterator curr;
  for (curr=node_.begin(); curr!=node_.end();)
  {
    if (curr->second->Owner() != myrank)
      node_.erase(curr++);
    else
      ++curr;
  }

  // build rowmap of nodes noderowmap_ if it does not exist
  if (noderowmap_==Teuchos::null) BuildNodeRowMap();
  const Epetra_Map& oldmap = *noderowmap_;

  // create an exporter object that will figure out the communication pattern
  DRT::Exporter exporter(oldmap,newmap,Comm());

  // Do the communication
  exporter.Export(node_);

  // update all ownership flags
  for (curr=node_.begin(); curr!=node_.end(); ++curr)
    curr->second->SetOwner(myrank);

  // maps and pointers are no longer correct and need rebuilding
  Reset(killdofs,killcond);

  return;
}

/*----------------------------------------------------------------------*
 |  Export nodes owned by a proc (public)                    mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::ExportColumnNodes(const Epetra_Map& newmap, bool killdofs, bool killcond)
{
  // destroy all ghosted nodes
  const int myrank = Comm().MyPID();
  std::map<int,Teuchos::RCP<DRT::Node> >::iterator curr;
  for (curr=node_.begin(); curr!=node_.end();)
  {
    if (curr->second->Owner() != myrank)
    {
      node_.erase(curr++);
    }
    else
    {
      ++curr;
    }
  }
  // build rowmap of nodes noderowmap_ if it does not exist
  if (noderowmap_==Teuchos::null) BuildNodeRowMap();
  const Epetra_Map& oldmap = *noderowmap_;

  // test whether all nodes in oldmap are also in newmap, otherwise
  // this would be a change of owner which is not allowed here
  for (int i=0; i<oldmap.NumMyElements(); ++i)
  {
    int gid = oldmap.GID(i);
    if (!(newmap.MyGID(gid)))
      dserror("Proc %d: Node gid=%d from oldmap is not in newmap",myrank,gid);
  }

  // create an exporter object that will figure out the communication pattern
  DRT::Exporter exporter(oldmap,newmap,Comm());
  // Do the communication
  exporter.Export(node_);

  // maps and pointers are no longer correct and need rebuilding
  Reset(killdofs,killcond);

  return;
}


/*----------------------------------------------------------------------*
 |  Export elements (public)                                 mwgee 02/11|
 *----------------------------------------------------------------------*/
void DRT::Discretization::ProcZeroDistributeElementsToAll(Epetra_Map& target,
                                                          std::vector<int>& gidlist)
{
  const int myrank = Comm().MyPID();

  // proc 0 looks for elements that are to be send to other procs
  int size = (int)gidlist.size();
  std::vector<int> pidlist(size); // gids on proc 0
  int err = target.RemoteIDList(size,&gidlist[0],&pidlist[0],NULL);
  if (err < 0) dserror("Epetra_BlockMap::RemoteIDList returned err=%d",err);

  std::map<int,std::vector<char> > sendmap; // proc to send a set of elements to
  if (!myrank)
  {
    std::map<int,DRT::PackBuffer > sendpb; // proc to send a set of elements to
    for (int i=0; i<size; ++i)
    {
      if (pidlist[i]==myrank or pidlist[i]<0) continue; // do not send to myself
      Element* actele = gElement(gidlist[i]);
      if (!actele) dserror("Cannot find global element %d",gidlist[i]);
      actele->Pack(sendpb[pidlist[i]]);
    }
    for (std::map<int,DRT::PackBuffer >::iterator fool = sendpb.begin(); fool != sendpb.end(); ++fool)
      fool->second.StartPacking();
    for (int i=0; i<size; ++i)
    {
      if (pidlist[i]==myrank or pidlist[i]<0) continue; // do not send to myself
      Element* actele = gElement(gidlist[i]);
      actele->Pack(sendpb[pidlist[i]]);
      element_.erase(actele->Id());
    }
    for (std::map<int,DRT::PackBuffer >::iterator fool = sendpb.begin(); fool != sendpb.end(); ++fool)
      swap(sendmap[fool->first],fool->second());
  }


#ifdef PARALLEL
  // tell everybody who is to receive something
  std::vector<int> receivers;

  for (std::map<int,std::vector<char> >::iterator fool = sendmap.begin(); fool !=sendmap.end(); ++fool)
    receivers.push_back(fool->first);
  size = (int)receivers.size();
  Comm().Broadcast(&size,1,0);
  if (myrank != 0) receivers.resize(size);
  Comm().Broadcast(&receivers[0],size,0);
  int foundme = -1;
  if (myrank != 0)
    for (int i=0; i<size; ++i)
      if (receivers[i]==myrank)
      {
        foundme = i;
        break;
      }


  // proc 0 sends out messages
  int tag = 0;
  DRT::Exporter exporter(Comm());
  std::vector<MPI_Request> request(size);
  if (!myrank)
  {
    for (std::map<int,std::vector<char> >::iterator fool = sendmap.begin(); fool !=sendmap.end(); ++fool)
    {
      exporter.ISend(0,fool->first,&fool->second[0],(int)fool->second.size(),tag,request[tag]);
      tag++;
    }
    if (tag != size) dserror("Number of messages is mixed up");
    // do not delete sendmap until Wait has returned!
  }


  // all other procs listen to message and put element into dis
  if (foundme != -1)
  {
    std::vector<char> recvdata;
    int length = 0;
    int source = -1;
    int tag = -1;
    exporter.ReceiveAny(source,tag,recvdata,length);
    if (source != 0 || tag != foundme)
      dserror("Messages got mixed up");
    // Put received elements into discretization
    std::vector<char>::size_type index = 0;
    while (index < recvdata.size())
    {
      std::vector<char> data;
      ParObject::ExtractfromPack(index,recvdata,data);
      DRT::ParObject* object = DRT::UTILS::Factory(data);
      DRT::Element* ele = dynamic_cast<DRT::Element*>(object);
      if (!ele) dserror("Received object is not an element");
      ele->SetOwner(myrank);
      Teuchos::RCP<DRT::Element> rcpele = Teuchos::rcp(ele);
      AddElement(rcpele);
      //printf("proc %d index %d\n",myrank,index); fflush(stdout);
    }
  }

  // wait for all communication to finish
  if (!myrank)
  {
    for (int i=0; i<size; ++i)
      exporter.Wait(request[i]);
  }
#endif

  Comm().Barrier(); // I feel better this way ;-)
  Reset();
  return;
}

/*----------------------------------------------------------------------*
 |  Export elements (public)                                 mwgee 03/11|
 *----------------------------------------------------------------------*/
void DRT::Discretization::ProcZeroDistributeNodesToAll(Epetra_Map& target)
{
  const int myrank = Comm().MyPID();

#if 0
  Epetra_Time timer(Comm());
  double t1 = timer.ElapsedTime();
#endif

  // proc 0 looks for nodes that are to be distributed
  Reset();
  BuildNodeRowMap();
  const Epetra_Map& oldmap = *noderowmap_;
  int size = oldmap.NumMyElements();
  if (myrank) size = 0;
  std::vector<int> pidlist(size,-1);
  {
    int err = target.RemoteIDList(size,oldmap.MyGlobalElements(),&pidlist[0],NULL);
    if (err) dserror("Epetra_BlockMap::RemoteIDLis returned err=%d",err);
  }

#if 0
  for (int proc=0; proc<Comm().NumProc(); ++proc)
  {
    if (proc==myrank)
    {
      printf("\nProc %d numnode %d\n",myrank,size);
      for (int i=0; i<size; ++i)
        printf("Proc %d gid %d pid %d\n",myrank,oldmap.MyGlobalElements()[i],pidlist[i]);
    }
    fflush(stdout);
    Comm().Barrier();
  }
#endif

#if 0
  double t2 = timer.ElapsedTime();
  if (!myrank) printf("\nTime 1 %10.5e\n",t2-t1); fflush(stdout);
#endif

  std::map<int,std::vector<char> > sendmap;
  if (!myrank)
  {
    std::map<int,DRT::PackBuffer > sendpb;
    for (int i=0; i<size; ++i)
    {
      // proc 0 does not send to itself
      if (pidlist[i]==myrank || pidlist[i]==-1) continue;
      Node* node = gNode(oldmap.MyGlobalElements()[i]);
      if (!node) dserror("Proc 0 cannot find global node %d",oldmap.MyGlobalElements()[i]);
      node->Pack(sendpb[pidlist[i]]);
    }
    for (std::map<int,DRT::PackBuffer >::iterator fool = sendpb.begin(); fool != sendpb.end(); ++fool)
      fool->second.StartPacking();
    for (int i=0; i<size; ++i)
    {
      // proc 0 does not send to itself
      if (pidlist[i]==myrank || pidlist[i]==-1) continue;
      Node* node = gNode(oldmap.MyGlobalElements()[i]);
      node->Pack(sendpb[pidlist[i]]);
      node_.erase(node->Id());
    }
    for (std::map<int,DRT::PackBuffer >::iterator fool = sendpb.begin(); fool != sendpb.end(); ++fool)
      swap(sendmap[fool->first],fool->second());
  }

#if 0
  double t3 = timer.ElapsedTime();
  if (!myrank) printf("Time 2 %10.5e\n",t3-t2); fflush(stdout);
#endif

#ifdef PARALLEL
  // tell everybody who is to receive something
  std::vector<int> receivers;
  for (std::map<int,std::vector<char> >::iterator fool = sendmap.begin(); fool !=sendmap.end(); ++fool)
    receivers.push_back(fool->first);
  size = (int)receivers.size();
  Comm().Broadcast(&size,1,0);
  if (myrank != 0) receivers.resize(size);
  Comm().Broadcast(&receivers[0],size,0);
  int foundme = -1;
  if (myrank != 0)
    for (int i=0; i<size; ++i)
      if (receivers[i]==myrank)
      {
        foundme = i;
        break;
      }

  // proc 0 sends out messages
  int tag = 0;
  DRT::Exporter exporter(Comm());
  std::vector<MPI_Request> request(size);
  if (!myrank)
  {
    for (std::map<int,std::vector<char> >::iterator fool = sendmap.begin(); fool !=sendmap.end(); ++fool)
    {
      exporter.ISend(0,fool->first,&fool->second[0],(int)fool->second.size(),tag,request[tag]);
      tag++;
    }
    if (tag != size) dserror("Number of messages is mixed up");
    // do not delete sendmap until Wait has returned!
  }

  // all other procs listen to message and put node into dis
  if (foundme != -1)
  {
    std::vector<char> recvdata;
    int length = 0;
    int source = -1;
    int tag = -1;
    exporter.ReceiveAny(source,tag,recvdata,length);
    //printf("Proc %d received tag %d length %d\n",myrank,tag,length); fflush(stdout);
    if (source != 0 || tag != foundme)
      dserror("Messages got mixed up");
    // Put received nodes into discretization
    std::vector<char>::size_type index = 0;
    while (index < recvdata.size())
    {
      std::vector<char> data;
      ParObject::ExtractfromPack(index,recvdata,data);
      DRT::ParObject* object = DRT::UTILS::Factory(data);
      DRT::Node* node = dynamic_cast<DRT::Node*>(object);
      if (!node) dserror("Received object is not a node");
      node->SetOwner(myrank);
      Teuchos::RCP<DRT::Node> rcpnode = Teuchos::rcp(node);
      AddNode(rcpnode);
    }
  }


  // wait for all communication to finish
  if (!myrank)
  {
    for (int i=0; i<size; ++i)
      exporter.Wait(request[i]);
  }

#if 0
  Comm().Barrier(); // feel better this way ;-)
  double t4 = timer.ElapsedTime();
  if (!myrank) printf("Time 3 %10.5e\n",t4-t3); fflush(stdout);
#endif


#endif

  Comm().Barrier(); // feel better this way ;-)
  Reset();
  return;
}

/*----------------------------------------------------------------------*
 |  Export elements (public)                                 mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::ExportRowElements(const Epetra_Map& newmap, bool killdofs, bool killcond)
{
  // destroy all ghosted elements
  const int myrank = Comm().MyPID();
  std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
  for (curr=element_.begin(); curr!=element_.end();)
  {
    if (curr->second->Owner() != myrank)
    {
      element_.erase(curr++);
    }
    else
    {
      ++curr;
    }
  }

  // build map of elements elerowmap_ if it does not exist
  if (elerowmap_==Teuchos::null) BuildElementRowMap();
  const Epetra_Map& oldmap = *elerowmap_;

  // create an exporter object that will figure out the communication pattern
  DRT::Exporter exporter(oldmap,newmap,Comm());

  exporter.Export(element_);

  // update ownerships and kick out everything that's not in newmap
  for (curr=element_.begin(); curr!=element_.end(); ++curr)
    curr->second->SetOwner(myrank);

  // maps and pointers are no longer correct and need rebuilding
  Reset(killdofs,killcond);

  return;
}

/*----------------------------------------------------------------------*
 |  Export elements (public)                                 mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::ExportColumnElements(const Epetra_Map& newmap, bool killdofs, bool killcond)
{
  // destroy all ghosted elements
  const int myrank = Comm().MyPID();
  std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
  for (curr=element_.begin(); curr!=element_.end();)
  {
    if (curr->second->Owner() != myrank)
    {
      element_.erase(curr++);
    }
    else
    {
      ++curr;
    }
  }

  // build map of elements elerowmap_ if it does not exist
  if (elerowmap_==Teuchos::null) BuildElementRowMap();
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
  Reset(killdofs,killcond);

  return;
}

/*----------------------------------------------------------------------*
 |  build nodal graph from discretization (public)           mwgee 11/06|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsGraph> DRT::Discretization::BuildNodeGraph() const
{
  if (!Filled()) dserror("FillComplete() was not called on this discretization");

  // get nodal row map
  const Epetra_Map* noderowmap = NodeRowMap();

  // allocate graph
  Teuchos::RCP<Epetra_CrsGraph> graph =
                     Teuchos::rcp( new Epetra_CrsGraph(Copy,*noderowmap,108,false));

  // iterate all elements on this proc including ghosted ones
  // Note:
  // if a proc stores the appropiate ghosted elements, the resulting
  // graph will be the correct and complete graph of the distributed
  // discretization even if nodes are not ghosted.
  std::map<int,Teuchos::RCP<DRT::Element> >::const_iterator curr;
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
                                    Teuchos::RCP<Epetra_Map>& elerowmap,
                                    Teuchos::RCP<Epetra_Map>& elecolmap
                                    ) const
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
  std::vector<int> cnodeowner(ncnode);
  int err = noderowmap.RemoteIDList(ncnode,nodecolmap.MyGlobalElements(),&cnodeowner[0],NULL);
  if (err) dserror("Epetra_BlockMap::RemoteIDLis returned err=%d",err);

  // build connectivity of elements
  // storing :  element gid
  //            no. of nodes
  //            nodeids
  int stoposize = 2000;
  int count     = 0;
  std::vector<int> stopo(stoposize);
  std::map<int,Teuchos::RCP<DRT::Element> >::const_iterator ecurr;
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

  std::vector<int> rtopo(stoposize);

  // communicate number of nodes per proc
  std::vector<int> nodesperproc(numproc);
  int nummynodes = noderowmap.NumMyElements();
  Comm().GatherAll(&nummynodes, &nodesperproc[0], 1);

  // estimate no. of elements equal to no. of nodes
  std::vector<int> myele(nummynodes);
  int nummyele=0;
  // estimate no. of ghosted elements much lower
  std::vector<int> myghostele(nummynodes/4);
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
      // this is necessary to be able to own or ghost the element
      for (int j=0; j<numnode; ++j)
        if (!nodecolmap.MyGID(nodeids[j]))
          dserror("I do not have own/ghosted node gid=%d",nodeids[j]);

      // find out who owns how many of the nodes
      std::vector<int> nodeowner(numnode);
      std::vector<int> numperproc(numproc);
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
      //
      // tie-breaking if number of nodes is equal among some procs:
      // the processor with the smaller number of row nodes owns the element;
      // if still tied, the last node owner with equal number of nodes owns
      // the element
      int owner   = -1;
      int maxnode = 0;
      int minrownodes = noderowmap.NumGlobalElements();
      for (int j=0; j<numnode; ++j)
      {
        int currentproc = nodeowner[j];
        int ownhowmany  = numperproc[currentproc];
        if (ownhowmany > maxnode
            ||
            (ownhowmany == maxnode && nodesperproc[currentproc] <= minrownodes))
        {
          owner   = currentproc;
          maxnode = ownhowmany;
          minrownodes = nodesperproc[currentproc];
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
  elerowmap = Teuchos::rcp(new Epetra_Map(-1,nummyele,&myele[0],0,Comm()));
  if (!elerowmap->UniqueGIDs())
    dserror("Element row map is not unique");

  // build elecolmap
  std::vector<int> elecol(nummyele+nummyghostele);
  for (int i=0; i<nummyele; ++i) elecol[i] = myele[i];
  for (int i=0; i<nummyghostele; ++i) elecol[nummyele+i] = myghostele[i];
  elecolmap = Teuchos::rcp(new Epetra_Map(-1,nummyghostele+nummyele,
                                 &elecol[0],0,Comm()));

  return;
}

/*----------------------------------------------------------------------*
 |  redistribute discretization (public)                     mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Redistribute(const Epetra_Map& noderowmap,
                                       const Epetra_Map& nodecolmap,
                                       bool assigndegreesoffreedom ,
                                       bool initelements           ,
                                       bool doboundaryconditions,
                                       bool killdofs,
                                       bool killcond)
{
  // build the overlapping and non-overlapping element maps
  Teuchos::RCP<Epetra_Map> elerowmap;
  Teuchos::RCP<Epetra_Map> elecolmap;
  BuildElementRowColumn(noderowmap,nodecolmap,elerowmap,elecolmap);

  // export nodes and elements to the new maps
  ExportRowNodes(noderowmap,killdofs,killcond);
  ExportColumnNodes(nodecolmap,killdofs,killcond);
  ExportRowElements(*elerowmap,killdofs,killcond);
  ExportColumnElements(*elecolmap,killdofs,killcond);

  // these exports have set Filled()=false as all maps are invalid now
  int err = FillComplete(assigndegreesoffreedom,initelements,doboundaryconditions);

  if (err) dserror("FillComplete() returned err=%d",err);

  return;
}

/*----------------------------------------------------------------------*
 |  ghost elements according to element column map (public)  rauch 10/13|
 *----------------------------------------------------------------------*/
void DRT::Discretization::ExtendedGhosting(const Epetra_Map& elecolmap,
                                           bool assigndegreesoffreedom,
                                           bool initelements,
                                           bool doboundaryconditions,
                                           bool checkghosting)
{
#ifdef DEBUG
  if(Filled())
  {
    const Epetra_Map* oldelecolmap = ElementColMap();
    // check whether standard ghosting is included in extended ghosting
    for(int i=0; i<oldelecolmap->NumMyElements(); ++i)
    {
      bool hasgid = elecolmap.MyGID(oldelecolmap->GID(i));
      if(!hasgid)
        dserror("standard ghosting of ele %d is not included in extended ghosting", oldelecolmap->GID(i));
    }

    if(checkghosting)
    {
      int diff = elecolmap.NumGlobalElements() - oldelecolmap->NumGlobalElements();
      if (diff==0 and Comm().MyPID()==0)
        dserror("no additional elements have been ghosted");
    }
  }
#endif

  // first export the elements according to the processor local element column maps
  ExportColumnElements(elecolmap);

  // periodic boundary conditions require ghosting of all master and slave nodes,
  // if node of pbc set is contained in list of owned and ghosted elements
  // in case of pbcs, this has to be restored
  bool have_pbc = false;
  Teuchos::RCP<PBCDofSet> pbcdofset = Teuchos::null;
  // map of master nodes and corresponding slave nodes
  std::map<int,std::vector<int> > pbcmap;
  pbcmap.clear();
  // create the inverse map --- slavenode -> masternode
  std::map<int,std::vector<int> > inversenodecoupling;
  inversenodecoupling.clear();
  // map to be filled with new (extended) master nodes and corresponding slave nodes (in col layout)
  Teuchos::RCP<std::map<int,std::vector<int> > > pbcmapnew = Teuchos::rcp(new std::map<int,std::vector<int> >);
  (*pbcmapnew).clear();

  // check for pbcs
  for (int nds = 0; nds<NumDofSets(); nds++)
  {
    pbcdofset = Teuchos::rcp_dynamic_cast<PBCDofSet> (dofsets_[nds]);

    if (pbcdofset!=Teuchos::null)
    {
      have_pbc = true;
      pbcmap = *(pbcdofset->GetCoupledNodes());
      // it is assumed that, if one pbc set is available, all other potential dofsets hold the same layout
      break;
    }
  }

  // if pbcs are available, get master and slave information
  if (have_pbc)
  {
    // communicate all master and slave pairs
    // caution: we build redundant maps here, containing all master nodes
    LINALG::GatherAll(pbcmap,*comm_);

    // and build slave master pairs
    for(std::map<int,std::vector<int> >::iterator curr = pbcmap.begin();
        curr != pbcmap.end();
        ++curr )
    {
      for(unsigned rr=0;rr<curr->second.size();++rr)
        inversenodecoupling[curr->second[rr]].push_back(curr->first);
    }
  }

  // get the node ids of the elements that have to be ghosted and create a proper node column map for their export
  std::set<int> nodes;
  for (int lid=0;lid<elecolmap.NumMyElements();++lid)
  {
    DRT::Element* ele = this->gElement(elecolmap.GID(lid));
    const int* nodeids = ele->NodeIds();
    for(int inode=0; inode<ele->NumNode(); ++inode)
    {
      nodes.insert(nodeids[inode]);

      // for pbcs, take into account all master and slave pairs
      if (have_pbc)
      {
        // is present node a master node?
        std::map<int,std::vector<int> >::iterator foundmaster = pbcmap.find(nodeids[inode]);

        if (foundmaster!=pbcmap.end())
        {
          // also store all corresponding slave nodes in set of col nodes
          for (std::size_t rr=0; rr<foundmaster->second.size(); rr++)
            nodes.insert(foundmaster->second[rr]);

          // add master and corresponding slaves to new col list of master and slave pairs
          std::pair<int,std::vector<int> > mpair(foundmaster->first,foundmaster->second);
          (*pbcmapnew).insert(mpair);
        }
        else
        {
          // is present node a slave node?
          std::map<int,std::vector<int> >::iterator foundslave = inversenodecoupling.find(nodeids[inode]);

          if (foundslave!=inversenodecoupling.end())
          {
            // add corresponding master to set of col nodes
            nodes.insert(foundslave->second[0]);

            // store also all further slave nodes of this master (if multiple pbcs are used)
            for (std::size_t rr=0; rr<pbcmap[foundslave->second[0]].size(); rr++)
              nodes.insert((pbcmap[foundslave->second[0]])[rr]);

            // add master and corresponding slaves to new col list of master and slave pairs
            std::pair<int,std::vector<int> > mpair(foundslave->second[0],pbcmap[foundslave->second[0]]);
            (*pbcmapnew).insert(mpair);
          }
        }
      }
    }
  }

  // transfer master and slave information to pbc dofset
  if (have_pbc)
    pbcdofset->SetCoupledNodes(pbcmapnew);

  std::vector<int> colnodes(nodes.begin(),nodes.end());
  Teuchos::RCP<Epetra_Map> nodecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)colnodes.size(),&colnodes[0],0,Comm()));

  // now ghost the nodes
  ExportColumnNodes(*nodecolmap);

  // these exports have set Filled()=false as all maps are invalid now
  int err = FillComplete(assigndegreesoffreedom,initelements,doboundaryconditions);
  if(err)
    dserror("FillComplete() threw error code %d",err);

  return;
}


/*----------------------------------------------------------------------*
// this is to go away!!!!
 *----------------------------------------------------------------------*/
void DRT::Discretization::SetupGhostingWrongNameDoNotUse(
                                        bool assigndegreesoffreedom ,
                                        bool initelements           ,
                                        bool doboundaryconditions   )
{
  if (Filled())
    dserror("there is really no need to setup ghosting if the discretization is already filled");

  // build the graph ourselves
  std::map<int,std::set<int> > localgraph;
  for (std::map<int,Teuchos::RCP<DRT::Element> >::iterator i=element_.begin();
       i!=element_.end();
       ++i)
  {
    int numnodes = i->second->NumNode();
    const int* nodes = i->second->NodeIds();

    // loop nodes and add this topology to the row in the graph of every node
    for (int n=0; n<numnodes; ++n)
    {
      int nodelid = nodes[n];
      copy(nodes,
           nodes+numnodes,
           inserter(localgraph[nodelid],
                    localgraph[nodelid].begin()));
    }
  }

  // Create node row map. Only the row nodes go there.

  std::vector<int> gids;
  std::vector<int> entriesperrow;

  gids.reserve(localgraph.size());
  entriesperrow.reserve(localgraph.size());

  for (std::map<int,Teuchos::RCP<DRT::Node> >::iterator i=node_.begin();
       i!=node_.end();
       ++i)
  {
    gids.push_back(i->first);
    entriesperrow.push_back(localgraph[i->first].size());
  }

  Epetra_Map rownodes(-1,gids.size(),&gids[0],0,*comm_);

  // Construct FE graph. This graph allows processor off-rows to be inserted
  // as well. The communication issue is solved.

  Teuchos::RCP<Epetra_FECrsGraph> graph = Teuchos::rcp(new Epetra_FECrsGraph(Copy,rownodes,&entriesperrow[0],false));

  gids.clear();
  entriesperrow.clear();

  // Insert all rows into the graph, including the off ones.

  for (std::map<int,std::set<int> >::iterator i=localgraph.begin();
       i!=localgraph.end();
       ++i)
  {
    std::set<int>& rowset = i->second;
    std::vector<int> row;
    row.reserve(rowset.size());
    row.assign(rowset.begin(),rowset.end());
    rowset.clear();

    int err = graph->InsertGlobalIndices(1,&i->first,row.size(),&row[0]);
    if (err<0) dserror("graph->InsertGlobalIndices returned %d",err);
  }

  localgraph.clear();

  // Finalize construction of this graph. Here the communication
  // happens. The ghosting problem is solved at this point.

  int err = graph->GlobalAssemble(rownodes,rownodes);
  if (err) dserror("graph->GlobalAssemble returned %d",err);

  // partition graph using metis
  Epetra_Vector weights(graph->RowMap(),false);
  weights.PutScalar(1.0);
  Teuchos::RCP<Epetra_CrsGraph> gr = DRT::UTILS::PartGraphUsingMetis(*graph,weights);
  graph = Teuchos::null;

  // replace rownodes, colnodes with row and column maps from the graph
  // do stupid conversion from Epetra_BlockMap to Epetra_Map
  const Epetra_BlockMap& brow = gr->RowMap();
  const Epetra_BlockMap& bcol = gr->ColMap();
  Teuchos::RCP<Epetra_Map> noderowmap = Teuchos::rcp(new Epetra_Map(brow.NumGlobalElements(),
                                                  brow.NumMyElements(),
                                                  brow.MyGlobalElements(),
                                                  0,
                                                  *comm_));
  Teuchos::RCP<Epetra_Map> nodecolmap = Teuchos::rcp(new Epetra_Map(bcol.NumGlobalElements(),
                                                  bcol.NumMyElements(),
                                                  bcol.MyGlobalElements(),
                                                  0,
                                                  *comm_));

  gr = Teuchos::null;

  // Redistribute discretization to match the new maps.

  Redistribute(*noderowmap,
               *nodecolmap,
               assigndegreesoffreedom,
               initelements,
               doboundaryconditions);

}

void DRT::Discretization::SetupGhosting(
    bool assigndegreesoffreedom,
    bool initelements,
    bool doboundaryconditions)
{
  if (Filled())
    dserror("there is really no need to setup ghosting if the discretization is already filled");

  // build the graph ourselves
  std::map<int,std::set<int> > localgraph;
  for (std::map<int,Teuchos::RCP<DRT::Element> >::iterator i=element_.begin();
       i!=element_.end();
       ++i)
  {
    int numnodes = i->second->NumNode();
    const int* nodes = i->second->NodeIds();

    // loop nodes and add this topology to the row in the graph of every node
    for (int n=0; n<numnodes; ++n)
    {
      int nodelid = nodes[n];
      copy(nodes,
           nodes+numnodes,
           inserter(localgraph[nodelid],
                    localgraph[nodelid].begin()));
    }
  }

  // Create node row map. Only the row nodes go there.

  std::vector<int> gids;
  std::vector<int> entriesperrow;

  gids.reserve(localgraph.size());
  entriesperrow.reserve(localgraph.size());

  for (std::map<int,Teuchos::RCP<DRT::Node> >::iterator i=node_.begin();
       i!=node_.end();
       ++i)
  {
    gids.push_back(i->first);
    entriesperrow.push_back(localgraph[i->first].size());
  }

  Epetra_Map rownodes(-1,gids.size(),&gids[0],0,*comm_);

  // Construct FE graph. This graph allows processor off-rows to be inserted
  // as well. The communication issue is solved.

  Teuchos::RCP<Epetra_FECrsGraph> graph = Teuchos::rcp(new Epetra_FECrsGraph(Copy,rownodes,&entriesperrow[0],false));

  gids.clear();
  entriesperrow.clear();

  // Insert all rows into the graph, including the off ones.

  for (std::map<int,std::set<int> >::iterator i=localgraph.begin();
       i!=localgraph.end();
       ++i)
  {
    std::set<int>& rowset = i->second;
    std::vector<int> row;
    row.reserve(rowset.size());
    row.assign(rowset.begin(),rowset.end());
    rowset.clear();

    int err = graph->InsertGlobalIndices(1,&i->first,row.size(),&row[0]);
    if (err<0) dserror("graph->InsertGlobalIndices returned %d",err);
  }

  localgraph.clear();

  // Finalize construction of this graph. Here the communication
  // happens. The ghosting problem is solved at this point.

  int err = graph->GlobalAssemble(rownodes,rownodes);
  if (err) dserror("graph->GlobalAssemble returned %d",err);

  // replace rownodes, colnodes with row and column maps from the graph
  // do stupid conversion from Epetra_BlockMap to Epetra_Map
  const Epetra_BlockMap& brow = graph->RowMap();
  const Epetra_BlockMap& bcol = graph->ColMap();
  Teuchos::RCP<Epetra_Map> noderowmap = Teuchos::rcp(new Epetra_Map(brow.NumGlobalElements(),
                                                  brow.NumMyElements(),
                                                  brow.MyGlobalElements(),
                                                  0,
                                                  *comm_));
  Teuchos::RCP<Epetra_Map> nodecolmap = Teuchos::rcp(new Epetra_Map(bcol.NumGlobalElements(),
                                                  bcol.NumMyElements(),
                                                  bcol.MyGlobalElements(),
                                                  0,
                                                  *comm_));

  graph = Teuchos::null;

  // Redistribute discretization to match the new maps.

  Redistribute(*noderowmap,
               *nodecolmap,
               assigndegreesoffreedom,
               initelements,
               doboundaryconditions);

}
