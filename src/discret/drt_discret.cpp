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
filled_(false),
havedof_(false)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Discretization::Discretization(const DRT::Discretization& old) :
filled_(old.Filled()),
havedof_(old.HaveDofs())
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
  if (old.HaveDofs()) AssignDegreesOfFreedom();
  
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
void DRT::Discretization::AssignDegreesOfFreedom()
{
  if (!Filled()) dserror("Filled()==false");
  if (!NodeRowMap()->UniqueGIDs()) dserror("Nodal row map is not unique");
  if (!ElementRowMap()->UniqueGIDs()) dserror("Element row map is not unique");

  havedof_ = false;
  
  // loop my row nodes and set number of degrees of freedom to them
  for (int i=0; i<NumMyRowNodes(); ++i)
  {
    DRT::Node* actnode = lRowNode(i);
    const int numele = actnode->NumElement();
    DRT::Element** myele = actnode->Elements();
    int maxnum=0;
    for (int j=0; j<numele; ++j)
      maxnum = max(maxnum,myele[j]->NumDofPerNode(*actnode));
    actnode->Dof().SetNumDof(maxnum);
  }
  
  // build a redundant map for nodes
  vector<int> sredundantnodes(NumGlobalNodes());
  for (int i=0; i<NumGlobalNodes(); ++i) sredundantnodes[i] = 0;
  for (int i=0; i<NumMyRowNodes(); ++i)
  {
    const int gid = lRowNode(i)->Id();
    sredundantnodes[gid] = gid;
  }
  vector<int> rredundantnodes(NumGlobalNodes());
  Comm().SumAll(&sredundantnodes[0],&rredundantnodes[0],NumGlobalNodes());
  RefCountPtr<Epetra_Map> rednodemap = 
      rcp(new Epetra_Map(-1,NumGlobalNodes(),&rredundantnodes[0],0,Comm()));
  sredundantnodes.clear();
  rredundantnodes.clear();
  
  // build a map that holds all the node's numdof
  map<int,vector<int> > redundantnodedof;
  for (int i=0; i<NumMyRowNodes(); ++i)
  {
    const int gid = lRowNode(i)->Id();
    const int numdof = lRowNode(i)->Dof().NumDof();
    redundantnodedof[gid].resize(1);
    redundantnodedof[gid][0] = numdof;
  }
  
  // export this map to full overlap (I know this is painful, but how else to do it?))
  {
    DRT::Exporter exporter(*NodeRowMap(),*rednodemap,Comm());
    exporter.Export(redundantnodedof);
  }
  
  // go through the redundant map holding the sizes and assign dofs
  // note that in stl map all gids are ordered ascending
  // so we are numbering dofs in ascending order maintaining the bandwith
  // minimizing property of the node numbering
  int count=0;
  map<int,vector<int> >::iterator fool;
  for (fool=redundantnodedof.begin(); fool!=redundantnodedof.end(); ++fool)
  {
    const int numdof = fool->second[0];
    fool->second.resize(numdof);
    for (int i=0; i<numdof; ++i)
    {
      fool->second[i] = count;
      ++count;
    }
  }
  // element dof numbering starts from count
  const int starteledof = count;
  
  // export the redundant map to column map
  {
    DRT::Exporter exporter(*rednodemap,*NodeColMap(),Comm());
    exporter.Export(redundantnodedof);
  }
  
  // we don't need the redundant map anymore, destroy
  rednodemap = null;
  
  // redundantnodedof now is not redundant any more, better rename for clarity
  map<int,vector<int> >& colnodedofs = redundantnodedof;
  
  // loop my col nodes and assign degrees of freedom
  for (int i=0; i<NumMyColNodes(); ++i)
  {
    DRT::Node* actnode = lColNode(i);
    fool = colnodedofs.find(actnode->Id());
    if (fool == colnodedofs.end()) dserror("Cannot find node gid=%d in colnodedofs",actnode->Id());
    int* dofs = &(fool->second[0]);
    int numdof = (int)fool->second.size();
    actnode->Dof().SetDof(dofs,numdof);
  }
  
  // clear the nodal map
  redundantnodedof.clear();
  
  // Now do all this fun again for the elements
  // loop my row elements and set number of degrees of freedom
  for (int i=0; i<NumMyRowElements(); ++i)
  {
    DRT::Element* actele = lRowElement(i);
    actele->Dof().SetNumDof(actele->NumDofPerElement());
  }
  
  // build a redundant map for elements
  vector<int> sredundanteles(NumGlobalElements());
  for (int i=0; i<NumGlobalElements(); ++i) sredundanteles[i] = 0;
  for (int i=0; i<NumMyRowElements(); ++i)
  {
    const int gid = lRowElement(i)->Id();
    sredundanteles[gid]=gid;
  }
  vector<int> rredundanteles(NumGlobalElements());
  Comm().SumAll(&sredundanteles[0],&rredundanteles[0],NumGlobalElements());
  RefCountPtr<Epetra_Map> redelemap = 
    rcp(new Epetra_Map(-1,NumGlobalElements(),&rredundanteles[0],0,Comm()));
  sredundanteles.clear();
  rredundanteles.clear();
  
  // build a map that holds all the ele's numdof
  map<int,vector<int> > redundanteledof;
  for (int i=0; i<NumMyRowElements(); ++i)
  {
    DRT::Element* actele = lRowElement(i);
    const int gid = actele->Id();
    const int numdof = actele->Dof().NumDof();
    redundanteledof[gid].resize(1);
    redundanteledof[gid][0] = numdof;
  }
  
  // export redundanteledof to full redundance ( I know..., you are welcome to rewrite this method without)
  {
    DRT::Exporter exporter(*ElementRowMap(),*redelemap,Comm());
    exporter.Export(redundanteledof);
  }
  
  // go through the redundant map and assign dofs to elements
  // start with starteledof
  count=starteledof;
  for (fool=redundanteledof.begin(); fool!=redundanteledof.end(); ++fool)
  {
    const int numdof = fool->second[0];
    fool->second.resize(numdof);
    for (int i=0; i<numdof; ++i)
    {
      fool->second[i] = count;
      ++count;
    }
  }
  
  // export the full redundant map to column element map
  {
    DRT::Exporter exporter(*redelemap,*ElementColMap(),Comm());
    exporter.Export(redundanteledof);
  }
  
  // don't need the big redundant map any more
  redelemap = null;
  
  // rename element dof map into what it is now
  map<int,vector<int> >& coleledofs = redundanteledof;
  
  // loop my column elements and set degrees of freedom
  for (int i=0; i<NumMyColElements(); ++i)
  {
    DRT::Element* actele = lColElement(i);
    fool = coleledofs.find(actele->Id());
    if (fool==coleledofs.end()) dserror("Proc %d: Cannot find element gid=%d in coleledofs",Comm().MyPID(),actele->Id());
    int* dofs = &(fool->second[0]);
    int numdof = (int)fool->second.size();
    actele->Dof().SetDof(dofs,numdof);
  }
  
  // clear the element map
  coleledofs.clear();
  
  // set flag indicating that dofs now are present
  havedof_ = true;

  return;
}





#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
