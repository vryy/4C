/*!----------------------------------------------------------------------
\file discret.cpp
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
 |  ctor (public)                                            mwgee 11/06|
 |  comm             (in)  a communicator object                        |
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Discretization::Discretization(RefCountPtr<Epetra_Comm> comm) :
comm_(comm),
filled_(false)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Discretization::Discretization(const CCADISCRETIZATION::Discretization& old)
{
  dserror("Discretization does not have a copy constructor");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Discretization::~Discretization()
{
  return;
}

/*----------------------------------------------------------------------*
 |  Add an element (public)                                  mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::AddElement(RefCountPtr<CCADISCRETIZATION::Element> ele)
{
  filled_ = false;
  element_[ele->Id()] = ele;
  return;
}

/*----------------------------------------------------------------------*
 |  Add a node (public)                                      mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::AddNode(RefCountPtr<CCADISCRETIZATION::Node> node)
{
  filled_ = false;
  node_[node->Id()] = node;
  return;
}

/*----------------------------------------------------------------------*
 |  get element with global id (public)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Element* CCADISCRETIZATION::Discretization::gElement(int gid) const
{
  map<int,RefCountPtr<CCADISCRETIZATION::Element> >:: const_iterator curr = 
    element_.find(gid);
  if (curr == element_.end()) dserror("Element with gobal id gid=%d not stored on this proc",gid);
  else                        return curr->second.get();
  return NULL;
}

/*----------------------------------------------------------------------*
 |  get element with local id (public)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Element* CCADISCRETIZATION::Discretization::lRowElement(int lid) const
{
  if (!Filled()) 
    dserror("CCADISCRETIZATION::Discretization::lRowElement: Filled() != true");
  int gid = ElementRowMap()->GID(lid);
  if (gid<0) dserror("Element with local row index lid=%d not on this proc",lid);
  map<int,RefCountPtr<CCADISCRETIZATION::Element> >:: const_iterator curr = 
    element_.find(gid);
  if (curr == element_.end()) 
    dserror("Cannot find element with row index lid %d gid %d",lid,gid);
  else                        
    return curr->second.get();
  return NULL;
}

/*----------------------------------------------------------------------*
 |  get element with local id (public)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Element* CCADISCRETIZATION::Discretization::lColElement(int lid) const
{
  if (!Filled()) 
    dserror("CCADISCRETIZATION::Discretization::lColElement: Filled() != true");
  int gid = ElementColMap()->GID(lid);
  if (gid<0) dserror("Element with local column index lid=%d not on this proc",lid);
  map<int,RefCountPtr<CCADISCRETIZATION::Element> >:: const_iterator curr = 
    element_.find(gid);
  if (curr == element_.end()) 
    dserror("Cannot find element with column index lid %d gid %d",lid,gid);
  else                        
    return curr->second.get();
  return NULL;
}

/*----------------------------------------------------------------------*
 |  get node with global id (public)                         mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Node* CCADISCRETIZATION::Discretization::gNode(int gid) const
{
  map<int,RefCountPtr<CCADISCRETIZATION::Node> >:: const_iterator curr = 
    node_.find(gid);
  if (curr == node_.end()) dserror("Node with global id gid=%d not stored on this proc",gid);
  else                     return curr->second.get();
  return NULL;
}

/*----------------------------------------------------------------------*
 |  get node with local id (public)                         mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Node* CCADISCRETIZATION::Discretization::lRowNode(int lid) const
{
  if (!Filled()) 
    dserror("CCADISCRETIZATION::Discretization::lRowNode: Filled() != true");
  int gid = NodeRowMap()->GID(lid);
  if (gid<0) return NULL;
  map<int,RefCountPtr<CCADISCRETIZATION::Node> >:: const_iterator curr = 
    node_.find(gid);
  if (curr == node_.end()) dserror("Cannot find node with lid=%d gid=%d",lid,gid);
  else                     return curr->second.get();
  return NULL;
}

/*----------------------------------------------------------------------*
 |  get node with local id (public)                         mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Node* CCADISCRETIZATION::Discretization::lColNode(int lid) const
{
  if (!Filled()) 
    dserror("CCADISCRETIZATION::Discretization::lColNode: Filled() != true");
  int gid = NodeColMap()->GID(lid);
  if (gid<0) return NULL;
  map<int,RefCountPtr<CCADISCRETIZATION::Node> >:: const_iterator curr = 
    node_.find(gid);
  if (curr == node_.end()) dserror("Cannot find node with lid=%d gid=%d",lid,gid);
  else                     return curr->second.get();
  return NULL;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CCADISCRETIZATION::Discretization& dis)
{
  dis.Print(os); 
  return os;
}

/*----------------------------------------------------------------------*
 |  Print discretization (public)                            mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::Print(ostream& os) const
{
  int numglobalelements = 0;
  int numglobalnodes    = 0;
  if (Filled())
  {
    numglobalelements = NumGlobalElements();
    numglobalnodes    = NumGlobalNodes();
  }
  else
  {
    int nummynodes = 0;
    map<int,RefCountPtr<CCADISCRETIZATION::Node> >::const_iterator ncurr;
    for (ncurr=node_.begin(); ncurr != node_.end(); ++ncurr)
      if (ncurr->second->Owner() == Comm().MyPID()) nummynodes++;
    
    int nummyele   = 0;
    map<int,RefCountPtr<CCADISCRETIZATION::Element> >::const_iterator ecurr;
    for (ecurr=element_.begin(); ecurr != element_.end(); ++ecurr)
      if (ecurr->second->Owner() == Comm().MyPID()) nummyele++;
    
    Comm().SumAll(&nummynodes,&numglobalnodes,1);
    Comm().SumAll(&nummyele,&numglobalelements,1);
  }

  // print head
  if (Comm().MyPID()==0)
  {
    os << "--------------------------------------------------\n";
    os << numglobalelements << " Elements " << numglobalnodes << " Nodes (global)\n";
    os << "--------------------------------------------------\n";
    if (Filled())
    os << "Filled() = true\n";
    else
    os << "Filled() = false\n";
    os << "--------------------------------------------------\n";
  }
  // print elements
  for (int proc=0; proc < Comm().NumProc(); ++proc)
  {
    if (proc == Comm().MyPID())
    {
      if ((int)element_.size())
        os << "-------------------------- Proc " << proc << " :\n";
      map<int,RefCountPtr<CCADISCRETIZATION::Element> >:: const_iterator curr;
      for (curr = element_.begin(); curr != element_.end(); ++curr)
        os << *(curr->second) << endl;
      os << endl;
    }
    Comm().Barrier();
  }
  // print nodes
  for (int proc=0; proc < Comm().NumProc(); ++proc)
  {
    if (proc == Comm().MyPID())
    {
      if ((int)node_.size())
        os << "-------------------------- Proc " << proc << " :\n";
      map<int,RefCountPtr<CCADISCRETIZATION::Node> >:: const_iterator curr;
      for (curr = node_.begin(); curr != node_.end(); ++curr)
        os << *(curr->second) << endl;
      os << endl;
    }
    Comm().Barrier();
  }
  return;
}


/*----------------------------------------------------------------------*
 |  node <-> design node topology (public)                   mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::SetDesignEntityIds(Node::OnDesignEntity type, 
                                                           const vector<int>& nfenode, 
                                                           const vector<vector<int> >& fenode)
{
  const int ndentity = nfenode.size();
  for (int i=0; i<ndentity; ++i)
  {
    const int dentityid  = i;
    const int numfenodes = nfenode[i];
    for (int j=0; j<numfenodes; ++j)
    {
      Node* node = gNode(fenode[dentityid][j]);
      if (!node) dserror("Node with gid %d is not on this proc",fenode[dentityid][j]);
      node->SetDesignEntity(type,dentityid);
    }
  }

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
 |  redistribute discretization using metis (public)         mwgee 11/06|
 *----------------------------------------------------------------------*/
int CCADISCRETIZATION::Discretization::DistributeUsingMetis()
{
  if (!Filled()) dserror("FillComplete() was not called on this discretization");

  // we get everything on proc 0 here to do serial metis  
  
  
  
  return 0;
}




#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
