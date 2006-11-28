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

#include "drt_discret.H"
#include "drt_exporter.H"
#include "drt_dserror.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 |  comm             (in)  a communicator object                        |
 *----------------------------------------------------------------------*/
DRT::Discretization::Discretization(RefCountPtr<Epetra_Comm> comm) :
comm_(comm),
filled_(false)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Discretization::Discretization(const DRT::Discretization& old)
{
  dserror("Discretization does not have a copy constructor");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Discretization::~Discretization()
{
  return;
}

/*----------------------------------------------------------------------*
 |  Add an element (public)                                  mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::AddElement(RefCountPtr<DRT::Element> ele)
{
  filled_ = false;
  element_[ele->Id()] = ele;
  return;
}

/*----------------------------------------------------------------------*
 |  Add a node (public)                                      mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::AddNode(RefCountPtr<DRT::Node> node)
{
  filled_ = false;
  node_[node->Id()] = node;
  return;
}

/*----------------------------------------------------------------------*
 |  get element with global id (public)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Discretization::gElement(int gid) const
{
  map<int,RefCountPtr<DRT::Element> >:: const_iterator curr = 
    element_.find(gid);
  if (curr == element_.end()) dserror("Element with gobal id gid=%d not stored on this proc",gid);
  else                        return curr->second.get();
  return NULL;
}

/*----------------------------------------------------------------------*
 |  get element with local id (public)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Discretization::lRowElement(int lid) const
{
  if (!Filled()) 
    dserror("DRT::Discretization::lRowElement: Filled() != true");
  int gid = ElementRowMap()->GID(lid);
  if (gid<0) dserror("Element with local row index lid=%d not on this proc",lid);
  map<int,RefCountPtr<DRT::Element> >:: const_iterator curr = 
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
DRT::Element* DRT::Discretization::lColElement(int lid) const
{
  if (!Filled()) 
    dserror("DRT::Discretization::lColElement: Filled() != true");
  int gid = ElementColMap()->GID(lid);
  if (gid<0) dserror("Element with local column index lid=%d not on this proc",lid);
  map<int,RefCountPtr<DRT::Element> >:: const_iterator curr = 
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
DRT::Node* DRT::Discretization::gNode(int gid) const
{
  map<int,RefCountPtr<DRT::Node> >:: const_iterator curr = 
    node_.find(gid);
  if (curr == node_.end()) dserror("Node with global id gid=%d not stored on this proc",gid);
  else                     return curr->second.get();
  return NULL;
}

/*----------------------------------------------------------------------*
 |  get node with local id (public)                         mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Node* DRT::Discretization::lRowNode(int lid) const
{
  if (!Filled()) 
    dserror("DRT::Discretization::lRowNode: Filled() != true");
  int gid = NodeRowMap()->GID(lid);
  if (gid<0) return NULL;
  map<int,RefCountPtr<DRT::Node> >:: const_iterator curr = 
    node_.find(gid);
  if (curr == node_.end()) dserror("Cannot find node with lid=%d gid=%d",lid,gid);
  else                     return curr->second.get();
  return NULL;
}

/*----------------------------------------------------------------------*
 |  get node with local id (public)                         mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Node* DRT::Discretization::lColNode(int lid) const
{
  if (!Filled()) 
    dserror("DRT::Discretization::lColNode: Filled() != true");
  int gid = NodeColMap()->GID(lid);
  if (gid<0) dserror("local index lid=%d out of range",lid);
  map<int,RefCountPtr<DRT::Node> >:: const_iterator curr = 
    node_.find(gid);
  if (curr == node_.end()) dserror("Cannot find node with lid=%d gid=%d",lid,gid);
  else                     return curr->second.get();
  return NULL;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const DRT::Discretization& dis)
{
  dis.Print(os); 
  return os;
}

/*----------------------------------------------------------------------*
 |  Print discretization (public)                            mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Print(ostream& os) const
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
    map<int,RefCountPtr<DRT::Node> >::const_iterator ncurr;
    for (ncurr=node_.begin(); ncurr != node_.end(); ++ncurr)
      if (ncurr->second->Owner() == Comm().MyPID()) nummynodes++;
    
    int nummyele   = 0;
    map<int,RefCountPtr<DRT::Element> >::const_iterator ecurr;
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
      map<int,RefCountPtr<DRT::Element> >:: const_iterator curr;
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
      map<int,RefCountPtr<DRT::Node> >:: const_iterator curr;
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
void DRT::Discretization::SetDesignEntityIds(Node::OnDesignEntity type, 
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
