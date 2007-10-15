/*!----------------------------------------------------------------------
\file drt_discret.cpp

\brief a class to manage one discretization

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include <algorithm>

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
havedof_(false),
currentdofset_(0)
{
  dofsets_.push_back(rcp(new DofSet()));
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

  currentdofset_ = old.currentdofset_;
  for (unsigned i=0; i<old.dofsets_.size(); ++i)
    dofsets_.push_back(rcp(new DRT::DofSet::DofSet(*(old.dofsets_[i]))));

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
      {
        os << *(curr->second);
        if (Filled())
        {
          vector<int> dof = Dof(&*(curr->second));
          if (dof.size())
          {
            os << " Dofs ";
            for (unsigned i=0; i<dof.size(); ++i) os << setw(6) << dof[i] << " ";
          }
        }
        os << endl;
      }
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
      {
        os << *(curr->second);
        if (Filled())
        {
          vector<int> dof = Dof(&*(curr->second));
          if (dof.size())
          {
            os << " Dofs ";
            for (unsigned i=0; i<dof.size(); ++i) os << setw(6) << dof[i] << " ";
          }
        }
        os << endl;
      }
      os << endl;
    }
    Comm().Barrier();
  }
  // print conditions
  for (int proc=0; proc<Comm().NumProc(); ++proc)
  {
    if (proc==Comm().MyPID())
    {
      int numcond = condition_.size();
      if (numcond) 
        os << "-------------------------- Proc " << proc << " :\n";
      if (numcond)
      {
        os << numcond << " Conditions:\n";
        map<string,RefCountPtr<Condition> >::const_iterator curr;
        for (curr=condition_.begin(); curr != condition_.end(); ++curr)
        {
          os << curr->first << " ";
          os << *(curr->second) << endl;
        }
      }
      os << endl;
    }
    Comm().Barrier();
  }
  return;
}

/*----------------------------------------------------------------------*
 |  replace the dofset of the discretisation (public)        gammi 05/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::ReplaceDofSet(RefCountPtr<DofSet> newdofset)
{
  havedof_ = false;
  dofsets_[currentdofset_] = newdofset;
  return;
}

/*----------------------------------------------------------------------*
 |  get dof row map (public)                                 mwgee 12/06|
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::DofRowMap()
{
  if (!Filled()) dserror("FillComplete was not called on this discretization");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() not called on this discretization");

  return dofsets_[currentdofset_]->DofRowMap();
}


/*----------------------------------------------------------------------*
 |  get dof column map (public)                              mwgee 12/06|
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::DofColMap()
{
  if (!Filled()) dserror("FillComplete was not called on this discretization");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() not called on this discretization");

  return dofsets_[currentdofset_]->DofColMap();
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
  // This is a rought test, but it might be ok at this place. It is an
  // error anyway to hand in a vector that is not related to our dof
  // maps.
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

/*----------------------------------------------------------------------*
 |  Set a condition of a certain name                          (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
void DRT::Discretization::SetCondition(const string& name,RefCountPtr<Condition> cond)
{
  condition_.insert(pair<string,RefCountPtr<Condition> >(name,cond));
  filled_ = false;
  return;
}

/*----------------------------------------------------------------------*
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
void DRT::Discretization::GetCondition(const string& name,vector<DRT::Condition*>& out) const
{
  const int num = condition_.count(name);
  out.resize(num);
  multimap<string,RefCountPtr<Condition> >::const_iterator startit =
                                         condition_.lower_bound(name);
  multimap<string,RefCountPtr<Condition> >::const_iterator endit =
                                         condition_.upper_bound(name);
  int count=0;
  multimap<string,RefCountPtr<Condition> >::const_iterator curr;
  for (curr=startit; curr!=endit; ++curr)
    out[count++] = curr->second.get();
  if (count != num) dserror("Mismatch in number of conditions found");
  return;
}

/*----------------------------------------------------------------------*
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Condition* DRT::Discretization::GetCondition(const string& name)
{
  multimap<string,RefCountPtr<Condition> >::iterator curr =
                                         condition_.find(name);
  if (curr==condition_.end()) return NULL;
  curr = condition_.lower_bound(name);
  return curr->second.get();
}


/*----------------------------------------------------------------------*
 |  Pack local elements (row map) into buffer                  (public) |
 |                                                          m.kue 02/07 |
 *----------------------------------------------------------------------*/
RefCountPtr<vector<char> > DRT::Discretization::PackMyElements()
{
  if (!Filled()) dserror("FillComplete was not called on this discretization");
  RefCountPtr<vector<char> > block = rcp(new vector<char>);
  for (vector<DRT::Element*>::iterator i=elerowptr_.begin();
       i!=elerowptr_.end();
       ++i)
  {
    vector<char> data;
    (*i)->Pack(data);
    ParObject::AddtoPack(*block,data);
  }
  return block;
}


/*----------------------------------------------------------------------*
 |  Pack local nodes (row map) into buffer                     (public) |
 |                                                          m.kue 02/07 |
 *----------------------------------------------------------------------*/
RefCountPtr<vector<char> > DRT::Discretization::PackMyNodes()
{
  if (!Filled()) dserror("FillComplete was not called on this discretization");
  RefCountPtr<vector<char> > block = rcp(new vector<char>);
  for (vector<DRT::Node*>::iterator i=noderowptr_.begin();
       i!=noderowptr_.end();
       ++i)
  {
    vector<char> data;
    (*i)->Pack(data);
    ParObject::AddtoPack(*block,data);
  }
  return block;
}

/*----------------------------------------------------------------------*
 |  Unpack element buffer and create local elements            (public) |
 |                                                          m.kue 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Discretization::UnPackMyElements(RefCountPtr<vector<char> > e)
{
  int index = 0;
  while (index < static_cast<int>(e->size()))
  {
    vector<char> data;
    ParObject::ExtractfromPack(index,*e,data);
    DRT::ParObject* o = DRT::Utils::Factory(data);
    DRT::Element* ele = dynamic_cast<Element*>(o);
    if (ele == NULL)
    {
      dserror("Failed to build an element from the element data");
    }
    ele->SetOwner(comm_->MyPID());
    AddElement(rcp(ele));
  }
  // in case AddElement forgets...
  Reset();
}

/*----------------------------------------------------------------------*
 |  Unpack nodal buffer and create local nodes                 (public) |
 |                                                          m.kue 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Discretization::UnPackMyNodes(RefCountPtr<vector<char> > e)
{
  int index = 0;
  while (index < static_cast<int>(e->size()))
  {
    vector<char> data;
    ParObject::ExtractfromPack(index,*e,data);
    DRT::ParObject* o = DRT::Utils::Factory(data);
    DRT::Node* n = dynamic_cast<Node*>(o);
    if (n == NULL)
    {
      dserror("Failed to build a node from the node data");
    }
    n->SetOwner(comm_->MyPID());
    AddNode(rcp(n));
  }
  // in case AddNode forgets...
  Reset();
}

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
