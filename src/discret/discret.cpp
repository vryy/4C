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
#include "dserror.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 |  comm             (in)  a communicator object                        |
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Discretization::Discretization(RefCountPtr<Epetra_Comm> comm) :
comm_(comm),
filled_(false),
elemap_(null),
nodemap_(null)
{
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
 |  get global number of elements (public)                   mwgee 11/06|
 *----------------------------------------------------------------------*/
int CCADISCRETIZATION::Discretization::NumGlobalElements() const
{
  if (Filled()) return ElementMap()->NumGlobalElements();
  else dserror("FillComplete() must be called before call to NumGlobalElements() is possible");
  int local = NumMyElements();
  int global = 0;
  Comm().SumAll(&local,&global,1);
  return global;
}

/*----------------------------------------------------------------------*
 |  get global number of nodes (public)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
int CCADISCRETIZATION::Discretization::NumGlobalNodes() const
{
  if (Filled()) return NodeMap()->NumGlobalElements();
  else dserror("FillComplete() must be called before call to NumGlobalNodes() is possible");
  int local = NumMyNodes();
  int global = 0;
  Comm().SumAll(&local,&global,1);
  return global;
}

/*----------------------------------------------------------------------*
 |  get element with global id (public)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Element* CCADISCRETIZATION::Discretization::gElement(int gid) const
{
  map<int,RefCountPtr<CCADISCRETIZATION::Element> >:: const_iterator curr = 
    element_.find(gid);
  if (curr == element_.end()) return NULL;
  else                        return curr->second.get();
  return NULL;
}

/*----------------------------------------------------------------------*
 |  get element with local id (public)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Element* CCADISCRETIZATION::Discretization::lElement(int lid) const
{
  if (!Filled()) 
    dserror("CCADISCRETIZATION::Discretization::lElement: Filled() != true");
  int gid = ElementMap()->GID(lid);
  if (gid<0) return NULL;
  map<int,RefCountPtr<CCADISCRETIZATION::Element> >:: const_iterator curr = 
    element_.find(gid);
  if (curr == element_.end()) return NULL;
  else                        return curr->second.get();
  return NULL;
}

/*----------------------------------------------------------------------*
 |  get node with global id (public)                         mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Node* CCADISCRETIZATION::Discretization::gNode(int gid) const
{
  map<int,RefCountPtr<CCADISCRETIZATION::Node> >:: const_iterator curr = 
    node_.find(gid);
  if (curr == node_.end()) return NULL;
  else                     return curr->second.get();
  return NULL;
}

/*----------------------------------------------------------------------*
 |  get node with local id (public)                         mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Node* CCADISCRETIZATION::Discretization::lNode(int lid) const
{
  if (!Filled()) 
    dserror("CCADISCRETIZATION::Discretization::lElement: Filled() != true");
  int gid = NodeMap()->GID(lid);
  if (gid<0) return NULL;
  map<int,RefCountPtr<CCADISCRETIZATION::Node> >:: const_iterator curr = 
    node_.find(gid);
  if (curr == node_.end()) return NULL;
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
    int nummynodes = NumMyNodes();
    int nummyele   = NumMyElements();
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
      if (NumMyElements())
        os << "Proc " << proc << " :\n";
      map<int,RefCountPtr<CCADISCRETIZATION::Element> >:: const_iterator curr;
      for (curr = element_.begin(); curr != element_.end(); ++curr)
        os << *(curr->second) << endl;
    }
    Comm().Barrier();
  }
  // print nodes
  for (int proc=0; proc < Comm().NumProc(); ++proc)
  {
    if (proc == Comm().MyPID())
    {
      if (NumMyNodes())
        os << "Proc " << proc << " :\n";
      map<int,RefCountPtr<CCADISCRETIZATION::Node> >:: const_iterator curr;
      for (curr = node_.begin(); curr != node_.end(); ++curr)
        os << *(curr->second) << endl;
    }
    Comm().Barrier();
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Finalize construction (public)                           mwgee 11/06|
 *----------------------------------------------------------------------*/
int CCADISCRETIZATION::Discretization::FillComplete()
{
  filled_ = false;  

  // (re)build map of nodes nodemap_
  BuildNodeMap();
  
  // (re)build map of elements elemap_
  BuildElementMap();
  
  // (re)construct element -> node pointers
  BuildElementToNodePointers();

  // (re)construct node -> element pointers
  BuildNodeToElementPointers();
  
  filled_ = true;  
  return 0;
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
      if (!node)
        dserror("Node %d is not on this proc %d",j,Comm().MyPID());
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
 |  Build nodemap_ (public)                                  mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::BuildNodeMap()
{
  int nummynodes = NumMyNodes();
  int numglobalnodes = 0;
  Comm().SumAll(&nummynodes,&numglobalnodes,1);
  
  vector<int> nodeids(nummynodes);
  
  map<int,RefCountPtr<CCADISCRETIZATION::Node> >::iterator curr;
  int count=0;
  for (curr=node_.begin(); curr != node_.end(); ++curr)
  {
    nodeids[count] = curr->second->Id();
    ++count;
  }
  nodemap_ = rcp(new Epetra_Map(numglobalnodes,nummynodes,&nodeids[0],0,Comm()));
  return;
}


/*----------------------------------------------------------------------*
 |  Build elemap_ (public)                                   mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::BuildElementMap()
{
  int nummyeles = NumMyElements();
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
  elemap_ = rcp(new Epetra_Map(numglobaleles,nummyeles,&eleids[0],0,Comm()));
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









#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
