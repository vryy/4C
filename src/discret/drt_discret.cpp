/*!----------------------------------------------------------------------
\file drt_discret.cpp
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
#include "linalg_utils.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 |  comm             (in)  a communicator object                        |
 *----------------------------------------------------------------------*/
DRT::Discretization::Discretization(const string name, RefCountPtr<Epetra_Comm> comm) :
name_(name),
comm_(comm),
filled_(false),
havedof_(false)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Discretization::Discretization(const DRT::Discretization& old) :
name_(old.name_),
state_(old.state_)
{
  comm_ = rcp(old.comm_->Clone());
  Reset();
  
  // deep copy elements
  map<int,RefCountPtr<DRT::Element> >::const_iterator ecurr;
  for (ecurr=old.element_.begin(); ecurr!=old.element_.end(); ++ecurr)
    element_[ecurr->first] = rcp(ecurr->second->Clone());
  
  // deep copy nodes  
  map<int,RefCountPtr<DRT::Node> >::const_iterator ncurr;
  for (ncurr=old.node_.begin(); ncurr!=old.node_.end(); ++ncurr)
    node_[ncurr->first] = rcp(ncurr->second->Clone());

  // do fillcomplete if old was fillcomplete
  if (old.Filled()) FillComplete();
  
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
  element_[ele->Id()] = ele;
  Reset();
  return;
}

/*----------------------------------------------------------------------*
 |  Add a node (public)                                      mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::AddNode(RefCountPtr<DRT::Node> node)
{
  node_[node->Id()] = node;
  Reset();
  return;
}

/*----------------------------------------------------------------------*
 |  get nodal row map (public)                               mwgee 11/06|
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::NodeRowMap() const 
{ 
#ifdef DEBUG
  if (Filled()) return noderowmap_.get();
  else dserror("FillComplete() must be called before call to NodeRowMap()"); 
  return NULL; 
#else
  return noderowmap_.get();
#endif
}

/*----------------------------------------------------------------------*
 |  get nodal column map (public)                            mwgee 11/06|
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::NodeColMap() const 
{ 
#ifdef DEBUG
  if (Filled()) return nodecolmap_.get();
  else dserror("FillComplete() must be called before call to NodeColMap()"); 
  return NULL; 
#else
  return nodecolmap_.get();
#endif
}

/*----------------------------------------------------------------------*
 |  get element row map (public)                             mwgee 11/06|
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::ElementRowMap() const 
{ 
#ifdef DEBUG
  if (Filled()) return elerowmap_.get();
  else dserror("FillComplete() must be called before call to ElementRowMap()"); 
  return NULL; 
#else
  return elerowmap_.get();
#endif
}

/*----------------------------------------------------------------------*
 |  get element column map (public)                          mwgee 11/06|
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::ElementColMap() const 
{ 
#ifdef DEBUG
  if (Filled()) return elecolmap_.get();
  else dserror("FillComplete() must be called before call to ElementColMap()"); 
  return NULL; 
#else
  return elecolmap_.get();
#endif
}

/*----------------------------------------------------------------------*
 |  get global no of elements (public)                       mwgee 11/06|
 *----------------------------------------------------------------------*/
int DRT::Discretization::NumGlobalElements() const 
{ 
#ifdef DEBUG
  if (Filled()) return ElementRowMap()->NumGlobalElements();
  else dserror("FillComplete() must be called before call to NumGlobalElements()");
  return -1; 
#else
  return ElementRowMap()->NumGlobalElements();
#endif
}

/*----------------------------------------------------------------------*
 |  get no of my row elements (public)                       mwgee 11/06|
 *----------------------------------------------------------------------*/
int DRT::Discretization::NumMyRowElements() const 
{ 
#ifdef DEBUG
  if (Filled()) return ElementRowMap()->NumMyElements();
  else dserror("FillComplete() must be called before call to NumMyRowElements()"); 
  return -1;
#else
  return ElementRowMap()->NumMyElements();
#endif
}

/*----------------------------------------------------------------------*
 |  get no of my column elements (public)                    mwgee 11/06|
 *----------------------------------------------------------------------*/
int DRT::Discretization::NumMyColElements() const 
{ 
  if (Filled()) return ElementColMap()->NumMyElements();
  else return (int)element_.size(); 
}

/*----------------------------------------------------------------------*
 |  get global no of nodes (public)                          mwgee 11/06|
 *----------------------------------------------------------------------*/
int DRT::Discretization::NumGlobalNodes() const 
{ 
#ifdef DEBUG
  if (Filled()) return NodeRowMap()->NumGlobalElements();
  else dserror("FillComplete() must be called before call to NumGlobalNodes()");
  return -1; 
#else
  return NodeRowMap()->NumGlobalElements();
#endif
}

/*----------------------------------------------------------------------*
 |  get no of my row nodes (public)                          mwgee 11/06|
 *----------------------------------------------------------------------*/
int DRT::Discretization::NumMyRowNodes() const 
{ 
#ifdef DEBUG
  if (Filled()) return NodeRowMap()->NumMyElements();
  else dserror("FillComplete() must be called before call to NumMyRowNodes()");
  return -1; 
#else
  return NodeRowMap()->NumMyElements();
#endif
}

/*----------------------------------------------------------------------*
 |  get no of my column nodes (public)                       mwgee 11/06|
 *----------------------------------------------------------------------*/
int DRT::Discretization::NumMyColNodes() const 
{ 
  if (Filled()) return NodeColMap()->NumMyElements();
  else return (int)node_.size();
}

/*----------------------------------------------------------------------*
 |  query existance of element (public)                      mwgee 12/06|
 *----------------------------------------------------------------------*/
bool DRT::Discretization::HaveGlobalElement(int gid) const
{
  map<int,RefCountPtr<DRT::Element> >:: const_iterator curr = element_.find(gid);
  if (curr == element_.end()) return false; 
  else                        return true;
}

/*----------------------------------------------------------------------*
 |  get element with global id (public)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Discretization::gElement(int gid) const
{
#ifdef DEBUG
  map<int,RefCountPtr<DRT::Element> >:: const_iterator curr = element_.find(gid);
  if (curr == element_.end()) dserror("Element with gobal id gid=%d not stored on this proc",gid);
  else return curr->second.get();
  return NULL;
#else
  return element_.find(gid)->second.get();
#endif  
}

/*----------------------------------------------------------------------*
 |  query existance of node (public)                         mwgee 12/06|
 *----------------------------------------------------------------------*/
bool DRT::Discretization::HaveGlobalNode(int gid) const
{
  map<int,RefCountPtr<DRT::Node> >:: const_iterator curr = node_.find(gid);
  if (curr == node_.end()) return false; 
  else                     return true;
}

/*----------------------------------------------------------------------*
 |  get node with global id (public)                         mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Node* DRT::Discretization::gNode(int gid) const
{
#ifdef DEBUG
  map<int,RefCountPtr<DRT::Node> >:: const_iterator curr = node_.find(gid);
  if (curr == node_.end()) dserror("Node with global id gid=%d not stored on this proc",gid);
  else                     return curr->second.get();
  return NULL;
#else
  return node_.find(gid)->second.get();
#endif
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
    os << "Discretization: " << Name() << endl;
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





/*----------------------------------------------------------------------*
 |  get dof row map (public)                                 mwgee 12/06|
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::DofRowMap()
{
  if (dofrowmap_ != null) return dofrowmap_.get();
  if (!Filled()) dserror("FillComplete was not called on this discretization");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() not called on this discretization");

  // loop my row nodes and count dofs
  int numnodaldof = 0;
  for (int i=0; i<NumMyRowNodes(); ++i)
    numnodaldof += lRowNode(i)->Dof().NumDof();
    
  // loop my row elements and count dofs
  int numeledof = 0;
  for (int i=0; i<NumMyRowElements(); ++i)
    numeledof += lRowElement(i)->Dof().NumDof();
    
  vector<int> mygid(numnodaldof+numeledof);

  // loop my row nodes and record dofs
  int count=0;
  for (int i=0; i<NumMyRowNodes(); ++i)
  {
    DRT::Node* actnode = lRowNode(i);
    for (int j=0; j<actnode->Dof().NumDof(); ++j)
      mygid[count++] = actnode->Dof()[j];
  }
  
  // loop elements and record dofs
  for (int i=0; i<NumMyRowElements(); ++i)
  {
    DRT::Element* actele = lRowElement(i);
    for (int j=0; j<actele->Dof().NumDof(); ++j)
      mygid[count++] = actele->Dof()[j];
  }
  
  if (count !=  numnodaldof+numeledof)
    dserror("Mismatch in no. of dofs %d <-> %d",count,numnodaldof+numeledof);
    
  dofrowmap_ = rcp(new Epetra_Map(-1,numnodaldof+numeledof,&mygid[0],0,Comm()));
  if (!dofrowmap_->UniqueGIDs()) dserror("Dof row map is not unique");
  
  return dofrowmap_.get();
}


/*----------------------------------------------------------------------*
 |  get dof column map (public)                              mwgee 12/06|
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::DofColMap()
{
  if (dofcolmap_ != null) return dofcolmap_.get();
  if (!Filled()) dserror("FillComplete was not called on this discretization");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() not called on this discretization");

  // loop my column nodes and count dofs
  int numnodaldof = 0;
  for (int i=0; i<NumMyColNodes(); ++i)
    numnodaldof += lColNode(i)->Dof().NumDof();
    
  // loop my column elements and count dofs
  int numeledof = 0;
  for (int i=0; i<NumMyColElements(); ++i)
    numeledof += lColElement(i)->Dof().NumDof();
    
  vector<int> mygid(numnodaldof+numeledof);

  // loop my column nodes and record dofs
  int count=0;
  for (int i=0; i<NumMyColNodes(); ++i)
  {
    DRT::Node* actnode = lColNode(i);
    for (int j=0; j<actnode->Dof().NumDof(); ++j)
      mygid[count++] = actnode->Dof()[j];
  }
  
  // loop elements and record dofs
  for (int i=0; i<NumMyColElements(); ++i)
  {
    DRT::Element* actele = lColElement(i);
    for (int j=0; j<actele->Dof().NumDof(); ++j)
      mygid[count++] = actele->Dof()[j];
  }
  
  if (count !=  numnodaldof+numeledof)
    dserror("Mismatch in no. of dofs %d <-> %d",count,numnodaldof+numeledof);
    
  dofcolmap_ = rcp(new Epetra_Map(-1,numnodaldof+numeledof,&mygid[0],0,Comm()));
  
  return dofcolmap_.get();
}


/*----------------------------------------------------------------------*
 |  set a reference to a data vector (public)                mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::SetState(const string& name,RefCountPtr<const Epetra_Vector> state)
{
  if (!Filled()) dserror("FillComplete() was not called");
  const Epetra_Map* colmap = DofColMap();
  const Epetra_BlockMap& vecmap = state->Map();
  
  // if it's already in column map just set a reference
  if (vecmap.PointSameAs(*colmap))
    state_[name] = state;
  // if it's not in column map export and allocate
  else
  {
    RefCountPtr<Epetra_Vector> tmp = LINALG::CreateVector(*colmap,false);
    LINALG::Export(*state,*tmp);
    state_[name] = tmp;
  }
  return;
}

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
