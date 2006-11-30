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
DRT::Discretization::Discretization(const DRT::Discretization& old) :
filled_(old.Filled())
{
  comm_ = rcp(old.comm_->Clone());
  Reset();
  
  map<int,RefCountPtr<DRT::Element> >::const_iterator ecurr;
  for (ecurr=old.element_.begin(); ecurr!=old.element_.end(); ++ecurr)
    element_[ecurr->first] = rcp(ecurr->second->Clone());
    
  map<int,RefCountPtr<DRT::Node> >::const_iterator ncurr;
  for (ncurr=old.node_.begin(); ncurr!=old.node_.end(); ++ncurr)
    node_[ncurr->first] = rcp(ncurr->second->Clone());

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
 |  get node with global id (public)                         mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Node* DRT::Discretization::gNode(int gid) const
{
#ifdef DEBUG
  map<int,RefCountPtr<DRT::Node> >:: const_iterator curr = 
    node_.find(gid);
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
 |  set degrees of freedom (public)                          mwgee 11/06|
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::DofRowMap() const
{
  // do not reconstruct if it exists
  // dofrowmap_ is destroyed in Reset() which is used in ALL
  // discretization altering methods
  if (dofrowmap_ != null) return dofrowmap_.get();
  if (!Filled()) dserror("Filled()==false");
  
  const int myrank = Comm().MyPID();
  const int numproc = Comm().NumProc();
  
  // loop all my nodes and set number of degrees of freedom to them
  for (int i=0; i<NumMyColNodes(); ++i)
  {
    DRT::Node* actnode = lColNode(i);
    const int numele = actnode->NumElement();
    DRT::Element** myele = actnode->Elements();
    int maxnum=0;
    for (int j=0; j<numele; ++j)
      maxnum = max(maxnum,myele[j]->NumDofPerNode(*actnode));
    actnode->Dof().SetNumDof(maxnum);
  }
  
  // loop all my elements and set number of degrees of freedom to them
  for (int i=0; i<NumMyColElements(); ++i)
  {
    DRT::Element* actele = lColElement(i);
    actele->Dof().SetNumDof(actele->NumDofPerElement());
  }
  
  // count dofs in row nodes and elements
  int nodedofcount=0;
  for (int i=0; i<NumMyRowNodes(); ++i)
    nodedofcount += lRowNode(i)->Dof().NumDof();
  int eledofcount=0;
  for (int i=0; i<NumMyRowElements(); ++i)
    eledofcount += lRowElement(i)->Dof().NumDof();
  
  // communicate these sizes
  vector<int> sendbuff(numproc*2);
  vector<int> numdof(numproc*2);
  for (int i=0; i<(int)sendbuff.size(); ++i) sendbuff[i] = 0;
  sendbuff[myrank*2]   = nodedofcount;
  sendbuff[myrank*2+1] = eledofcount;
  Comm().SumAll(&sendbuff[0],&numdof[0],numproc*2);
  
  // find out where to start numbering dofs
  int start=0;
  for (int i=0; i<myrank*2; ++i)
    start += numdof[i];
  vector<int> myglobaldofs(numdof[myrank*2]+numdof[myrank*2+1]);
  for (int i=0; i<(int)myglobaldofs.size(); ++i)
    myglobaldofs[i] = start+i;
    
  // number my dofs
  map<int,vector<int > > dofpernode;
  map<int,vector<int > > dofperele;
  int count=0;
  for (int i=0; i<NumMyRowNodes(); ++i)
  {
    DRT::Node* actnode = lRowNode(i);
    const int numdof = actnode->Dof().NumDof();
    actnode->Dof().SetDof(&myglobaldofs[count],numdof);
    dofpernode[actnode->Id()].resize(numdof);
    for (int j=0; j<numdof; ++j) 
      dofpernode[actnode->Id()][j] = myglobaldofs[count+j];
    count += numdof;
  }
  if (count != numdof[myrank*2]) 
    dserror("Mismatch in local no. of dofs: %d <-> %d",count,numdof[myrank*2]);
  for (int i=0; i<NumMyRowElements(); ++i)
  {
    DRT::Element* actele = lRowElement(i);
    const int numdof = actele->Dof().NumDof();
    actele->Dof().SetDof(&myglobaldofs[count],numdof);
    dofperele[actele->Id()].resize(numdof);
    for (int j=0; j<numdof; ++j)
      dofperele[actele->Id()][j] = myglobaldofs[count+j];
    count += numdof;
  }
  if (count != numdof[myrank*2]+numdof[myrank*2+1]) 
    dserror("Mismatch in local no. of dofs: %d <-> %d",count,numdof[myrank*2]+numdof[myrank*2+1]);

  // communicate the nodal dofs
  DRT::Exporter exporter(*NodeRowMap(),*NodeColMap(),Comm());
  exporter.Export(dofpernode);

  return NULL;
}





#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
