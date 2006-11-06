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
  int local = NumMyNodes();
  int global = 0;
  Comm().SumAll(&local,&global,1);
  return global;
}

/*----------------------------------------------------------------------*
 |  get element with global id (public)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
const CCADISCRETIZATION::Element* CCADISCRETIZATION::Discretization::Element(int gid) const
{
  map<int,RefCountPtr<CCADISCRETIZATION::Element> >:: const_iterator curr = 
    element_.find(gid);
  if (curr == element_.end()) return NULL;
  else                        return curr->second.get();
  return NULL;
}

/*----------------------------------------------------------------------*
 |  get node with global id (public)                         mwgee 11/06|
 *----------------------------------------------------------------------*/
const CCADISCRETIZATION::Node* CCADISCRETIZATION::Discretization::Node(int gid) const
{
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
  dis.Print(); 
  return os;
}

/*----------------------------------------------------------------------*
 |  Print discretization (public)                            mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::Print() const
{
  int numglobalelements = NumGlobalElements();
  int numglobalnodes    = NumGlobalNodes();

  // print head
  if (Comm().MyPID()==0)
  {
    cout << "--------------------------------------------------\n";
    cout << numglobalelements << " Elements " << numglobalnodes << " Nodes (global)\n";
    cout << "--------------------------------------------------\n";
  }
  
  // print elements
  for (int proc=0; proc < Comm().NumProc(); ++proc)
  {
    if (proc == Comm().MyPID())
    {
      cout << "Proc " << proc << " :\n";
      map<int,RefCountPtr<CCADISCRETIZATION::Element> >:: const_iterator curr;
      for (curr = element_.begin(); curr != element_.end(); ++curr)
        cout << *(curr->second) << endl;
    }
    Comm().Barrier();
  }
  
  // print nodes
    
  
  return;
}

/*----------------------------------------------------------------------*
 |  Finalize construction (public)                           mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::FillComplete()
{
  filled_ = false;  

  // (re)build map of nodes nodemap_
  BuildNodeMap();
  
  // (re)build map of elements elemap_
  BuildElementMap();
  
  // (re)construct element -> nodes pointers
  BuildElementToNodePointers();

  // (re)construct node -> element pointers
  BuildNodeToElementPointers();
  
  filled_ = true;  
  return;
}


/*----------------------------------------------------------------------*
 |  Build ptrs node -> element (public)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::BuildNodeToElementPointers()
{
  map<int,RefCountPtr<CCADISCRETIZATION::Element> >::iterator elecurr;
  for (elecurr=element_.begin(); elecurr != element_.end(); ++elecurr)
  {
    const int  nnode = elecurr->second->NumNode();
    const int* nodes = elecurr->second->NodeIds();
    for (int j=0; j<nnode; ++j)
    {
      const CCADISCRETIZATION::Node* node = Node(nodes[j]);
      if (!node)
      {
        cout << "CCADISCRETIZATION::Discretization::BuildElementToNodePointers:\n"
             << "Node is not on this proc\n"
             << __FILE__ << ":" << __LINE__ << endl;
        exit(EXIT_FAILURE);
      }
      else
        node->AddElementPtr(elecurr->second.get());
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
    {
      cout << "CCADISCRETIZATION::Discretization::BuildElementToNodePointers:\n"
           << "Building element <-> node topology failed\n"
           << __FILE__ << ":" << __LINE__ << endl;
      exit(EXIT_FAILURE);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Build nodemap_ (public)                                  mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Discretization::BuildNodeMap()
{
  int nummynodes = NumMyNodes();
  int numglobalnodes = NumGlobalNodes();
  
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
  int numglobaleles = NumGlobalElements();
  
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












#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
