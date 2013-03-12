/*!----------------------------------------------------------------------
\file drt_discret.cpp

\brief a class to manage one discretization

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#include <algorithm>
#include <Teuchos_TimeMonitor.hpp>

#include "drt_discret.H"
#include "drt_exporter.H"
#include "drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "drt_dofset_proxy.H"
#include "drt_dofset_subproxy.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 |  comm             (in)  a communicator object                        |
 *----------------------------------------------------------------------*/
DRT::Discretization::Discretization(const std::string name, RCP<Epetra_Comm> comm) :
name_(name),
comm_(comm),
filled_(false),
havedof_(false)
{
  dofsets_.push_back(Teuchos::rcp(new DofSet()));
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Discretization::Discretization(const DRT::Discretization& old) :
name_(old.name_),
state_(old.state_)
{
  comm_ = Teuchos::rcp(old.comm_->Clone());
  Reset();

  // deep copy elements
  std::map<int,RCP<DRT::Element> >::const_iterator ecurr;
  for (ecurr=old.element_.begin(); ecurr!=old.element_.end(); ++ecurr)
    element_[ecurr->first] = Teuchos::rcp(ecurr->second->Clone());

  // deep copy nodes
  std::map<int,RCP<DRT::Node> >::const_iterator ncurr;
  for (ncurr=old.node_.begin(); ncurr!=old.node_.end(); ++ncurr)
    node_[ncurr->first] = Teuchos::rcp(ncurr->second->Clone());

  for (unsigned i=0; i<old.dofsets_.size(); ++i)
    dofsets_.push_back(old.dofsets_[i]->Clone());

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
void DRT::Discretization::AddElement(RCP<DRT::Element> ele)
{
  element_[ele->Id()] = ele;
  Reset();
  return;
}

/*----------------------------------------------------------------------*
 |Calls Reset() on each processor if not Filled() == true on each proc   |
 |                                                            cyron 10/09|
 *----------------------------------------------------------------------*/
void DRT::Discretization::CheckFilledGlobally()
{
  //global filled flag (is true / one if and only if filled_ == true on each processor
  int globalfilled = 0;

  //convert filled_ flag on this procesor  into integer (no Epetra communicator for type bool)
  int localfilled = (int)filled_;

  /*the global filled flag is set to the minimal value of any local filled flag
   * i.e. if on any processor filled_ == false, the flag globalfilled is set to
   * zero*/
  Comm().MinAll(&localfilled,&globalfilled,1);

  //if not Filled() == true on all the processors call Reset()
  if(!globalfilled)
    Reset();

  return;
}

/*----------------------------------------------------------------------*
 |  Add a node (public)                                      mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::AddNode(RCP<DRT::Node> node)
{
  node_[node->Id()] = node;
  Reset();
  return;
}

/*----------------------------------------------------------------------*
 |  delete an node (public)                                  mwgee 10/08|
 *----------------------------------------------------------------------*/
bool DRT::Discretization::DeleteNode(RCP<DRT::Node> node)
{
  std::map<int,RCP<DRT::Node> >::iterator fool = node_.find(node->Id());
  if (fool==node_.end()) return false;
  node_.erase(fool);
  Reset();
  return true;
}

/*----------------------------------------------------------------------*
 |  delete an node (public)                                  mwgee 10/08|
 *----------------------------------------------------------------------*/
bool DRT::Discretization::DeleteNode(const int gid)
{
  std::map<int,RCP<DRT::Node> >::iterator fool = node_.find(gid);
  if (fool==node_.end()) return false;
  node_.erase(fool);
  Reset();
  return true;
}

/*----------------------------------------------------------------------*
 |  delete an element (public)                               mwgee 10/08|
 *----------------------------------------------------------------------*/
bool DRT::Discretization::DeleteElement(RCP<DRT::Element> ele)
{
  std::map<int,RCP<DRT::Element> >::iterator fool = element_.find(ele->Id());
  if (fool==element_.end()) return false;
  element_.erase(fool);
  Reset();
  return true;
}

/*----------------------------------------------------------------------*
 |  delete an element (public)                               mwgee 10/08|
 *----------------------------------------------------------------------*/
bool DRT::Discretization::DeleteElement(const int gid)
{
  std::map<int,RCP<DRT::Element> >::iterator fool = element_.find(gid);
  if (fool==element_.end()) return false;
  element_.erase(fool);
  Reset();
  return true;
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
  std::map<int,RCP<DRT::Element> >:: const_iterator curr = element_.find(gid);
  if (curr == element_.end()) return false;
  else                        return true;
}

/*----------------------------------------------------------------------*
 |  get element with global id (public)                      mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Discretization::gElement(int gid) const
{
#ifdef DEBUG
  std::map<int,RCP<DRT::Element> >:: const_iterator curr = element_.find(gid);
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
  std::map<int,RCP<DRT::Node> >:: const_iterator curr = node_.find(gid);
  if (curr == node_.end()) return false;
  else                     return true;
}

/*----------------------------------------------------------------------*
 |  get node with global id (public)                         mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Node* DRT::Discretization::gNode(int gid) const
{
#ifdef DEBUG
  std::map<int,RCP<DRT::Node> >:: const_iterator curr = node_.find(gid);
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
    std::map<int,RCP<DRT::Node> >::const_iterator ncurr;
    for (ncurr=node_.begin(); ncurr != node_.end(); ++ncurr)
      if (ncurr->second->Owner() == Comm().MyPID()) nummynodes++;

    int nummyele   = 0;
    std::map<int,RCP<DRT::Element> >::const_iterator ecurr;
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
      os << "-------------------------- Proc " << proc << " :\n";
      std::map<int,RCP<DRT::Element> >:: const_iterator curr;
      for (curr = element_.begin(); curr != element_.end(); ++curr)
      {
        os << *(curr->second);
        if (Filled() && HaveDofs())
        {
          std::vector<int> dof = Dof(0,&*(curr->second));
          if (dof.size())
          {
            os << " Dofs ";
            for (unsigned i=0; i<dof.size(); ++i) os << std::setw(6) << dof[i] << " ";
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
      os << "-------------------------- Proc " << proc << " :\n";
      std::map<int,RCP<DRT::Node> >:: const_iterator curr;
      for (curr = node_.begin(); curr != node_.end(); ++curr)
      {
        os << *(curr->second);
        if (Filled() && HaveDofs())
        {
          std::vector<int> dof = Dof(0,&*(curr->second));
          if (dof.size())
          {
            os << " Dofs ";
            for (unsigned i=0; i<dof.size(); ++i) os << std::setw(6) << dof[i] << " ";
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
        std::map<std::string,RCP<Condition> >::const_iterator curr;
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
 |  get dof row map (public)                                 mwgee 12/06|
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::DofRowMap(unsigned nds) const
{
  dsassert(nds<dofsets_.size(),"undefined dof set");
  if (!Filled()) dserror("FillComplete was not called on this discretization");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() not called on this discretization");

  return dofsets_[nds]->DofRowMap();
}


/*----------------------------------------------------------------------*
 |  get dof column map (public)                              mwgee 12/06|
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::DofColMap(unsigned nds) const
{
  dsassert(nds<dofsets_.size(),"undefined dof set");
  if (!Filled()) dserror("FillComplete was not called on this discretization");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() not called on this discretization");

  return dofsets_[nds]->DofColMap();
}


/*----------------------------------------------------------------------*
 |  replace the dofset of the discretisation (public)        gammi 05/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::ReplaceDofSet(unsigned nds, Teuchos::RCP<DofSet> newdofset, bool replaceinstatdofsets)
{
  dsassert(nds<dofsets_.size(),"undefined dof set");
  havedof_ = false;
  if (replaceinstatdofsets)
    newdofset->ReplaceInStaticDofsets(dofsets_[nds]);
  dofsets_[nds] = newdofset;
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::Discretization::AddDofSet(Teuchos::RCP<DofSet> newdofset)
{
  // if we already have our dofs here and we add a properly filled (proxy)
  // DofSet, we do not need (and do not want) to refill.
  havedof_ = havedof_ and newdofset->Filled();
  dofsets_.push_back(newdofset);
  return dofsets_.size()-1;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::DofSet> DRT::Discretization::GetDofSetProxy()
{
  return Teuchos::rcp(new DofSetProxy(&*dofsets_[0]));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::DofSet> DRT::Discretization::GetDofSetProxy(const Epetra_Map* subcolnodes,
                                                              const Epetra_Map* subcoleles )
{
  return Teuchos::rcp(new DofSetSubProxy(&*dofsets_[0],subcolnodes,subcoleles));
}


/*----------------------------------------------------------------------*
 |  replace the dofset of the discretisation (public)        gammi 05/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::ReplaceDofSet(RCP<DofSet> newdofset, bool replaceinstatdofsets)
{
  dsassert(dofsets_.size()==1,"expect just one dof set");
  havedof_ = false;
  if (replaceinstatdofsets)
    newdofset->ReplaceInStaticDofsets(dofsets_[0]);
  dofsets_[0] = newdofset;
  return;
}

/*----------------------------------------------------------------------*
 |  set a reference to a data vector (public)                mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::SetState(unsigned nds,const string& name,RCP<const Epetra_Vector> state)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::SetState");

  if (!HaveDofs()) dserror("FillComplete() was not called");
  const Epetra_Map* colmap = DofColMap(nds);
  const Epetra_BlockMap& vecmap = state->Map();

  if (state_.size()<=nds)
    state_.resize(nds+1);

  // if it's already in column map just set a reference
  // This is a rough test, but it might be ok at this place. It is an
  // error anyway to hand in a vector that is not related to our dof
  // maps.
  if (vecmap.PointSameAs(*colmap))
  {
    state_[nds][name] = state;
  }
  else // if it's not in column map export and allocate
  {
#ifdef DEBUG
    if (not DofRowMap(nds)->SameAs(state->Map()))
    {
      dserror("row map of discretization and state vector %s are different. This is a fatal bug!",name.c_str());
    }
#endif
    RCP<Epetra_Vector> tmp = LINALG::CreateVector(*colmap,false);
    LINALG::Export(*state,*tmp);
    state_[nds][name] = tmp;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Set a condition of a certain name                          (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
void DRT::Discretization::SetCondition(const string& name,RCP<Condition> cond)
{
  condition_.insert(std::pair<std::string,RCP<Condition> >(name,cond));
  filled_ = false;
  return;
}

/*----------------------------------------------------------------------*
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
void DRT::Discretization::GetCondition(const string& name,std::vector<DRT::Condition*>& out) const
{
  const int num = condition_.count(name);
  out.resize(num);
  std::multimap<std::string,RCP<Condition> >::const_iterator startit =
                                         condition_.lower_bound(name);
  std::multimap<std::string,RCP<Condition> >::const_iterator endit =
                                         condition_.upper_bound(name);
  int count=0;
  std::multimap<std::string,RCP<Condition> >::const_iterator curr;
  for (curr=startit; curr!=endit; ++curr)
    out[count++] = curr->second.get();
  if (count != num) dserror("Mismatch in number of conditions found");
  return;
}
/*----------------------------------------------------------------------*
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
void DRT::Discretization::GetCondition(const string& name,std::vector<Teuchos::RCP<Condition> >& out) const
{
  const int num = condition_.count(name);
  out.resize(num);
  std::multimap<std::string,RCP<Condition> >::const_iterator startit =
                                         condition_.lower_bound(name);
  std::multimap<std::string,RCP<Condition> >::const_iterator endit =
                                         condition_.upper_bound(name);
  int count=0;
  std::multimap<std::string,RCP<Condition> >::const_iterator curr;
  for (curr=startit; curr!=endit; ++curr)
    out[count++] = curr->second;
  if (count != num) dserror("Mismatch in number of conditions found");
  return;
}

/*----------------------------------------------------------------------*
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Condition* DRT::Discretization::GetCondition(const string& name) const
{
  std::multimap<std::string,RCP<Condition> >::const_iterator curr =
                                         condition_.find(name);
  if (curr==condition_.end()) return NULL;
  curr = condition_.lower_bound(name);
  return curr->second.get();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::GetConditionNames( std::vector<std::string> & names ) const
{
  std::set<std::string> n;
  for ( std::multimap<std::string,RCP<Condition> >::const_iterator curr=condition_.begin();
        curr!=condition_.end();
        ++curr )
  {
    n.insert( curr->first );
  }
  names.reserve( n.size() );
  names.assign( n.begin(), n.end() );
}


/*----------------------------------------------------------------------*
 |  Pack local elements (row map) into buffer                  (public) |
 |                                                          m.kue 02/07 |
 *----------------------------------------------------------------------*/
RCP<std::vector<char> > DRT::Discretization::PackMyElements() const
{
  if (!Filled()) dserror("FillComplete was not called on this discretization");

  DRT::PackBuffer buffer;

  for (std::vector<DRT::Element*>::const_iterator i=elerowptr_.begin();
       i!=elerowptr_.end();
       ++i)
  {
    DRT::Element * e = *i;
    e->Pack(buffer);
  }

  buffer.StartPacking();

  for (std::vector<DRT::Element*>::const_iterator i=elerowptr_.begin();
       i!=elerowptr_.end();
       ++i)
  {
    DRT::Element * e = *i;
    e->Pack(buffer);
  }

  RCP<std::vector<char> > block = Teuchos::rcp(new std::vector<char>);
  std::swap( *block, buffer() );
  return block;
}


/*----------------------------------------------------------------------*
 |  Pack local nodes (row map) into buffer                     (public) |
 |                                                          m.kue 02/07 |
 *----------------------------------------------------------------------*/
RCP<std::vector<char> > DRT::Discretization::PackMyNodes() const
{
  if (!Filled()) dserror("FillComplete was not called on this discretization");

  DRT::PackBuffer buffer;

  for (std::vector<DRT::Node*>::const_iterator i=noderowptr_.begin();
       i!=noderowptr_.end();
       ++i)
  {
    DRT::Node * n = *i;
    n->Pack(buffer);
  }

  buffer.StartPacking();

  for (std::vector<DRT::Node*>::const_iterator i=noderowptr_.begin();
       i!=noderowptr_.end();
       ++i)
  {
    DRT::Node * n = *i;
    n->Pack(buffer);
  }

  RCP<std::vector<char> > block = Teuchos::rcp(new std::vector<char>);
  std::swap( *block, buffer() );
  return block;
}


/*----------------------------------------------------------------------*
 |  Pack condition into buffer                                 (public) |
 |                                                          a.ger 11/07 |
 *----------------------------------------------------------------------*/
RCP<std::vector<char> > DRT::Discretization::PackCondition(const std::string condname) const
{
  if (!Filled()) dserror("FillComplete was not called on this discretization");

  // get boundary conditions
  std::vector<DRT::Condition*> cond;
  GetCondition(condname,cond);

  DRT::PackBuffer buffer;

  for (std::vector<DRT::Condition*>::const_iterator i = cond.begin();
       i!=cond.end();
       ++i)
  {
    DRT::Condition * c = *i;
    c->Pack(buffer);
  }

  buffer.StartPacking();

  for (std::vector<DRT::Condition*>::const_iterator i = cond.begin();
       i!=cond.end();
       ++i)
  {
    DRT::Condition * c = *i;
    c->Pack(buffer);
  }

  RCP<std::vector<char> > block = Teuchos::rcp(new std::vector<char>);
  std::swap( *block, buffer() );
  return block;
}


/*----------------------------------------------------------------------*
 |  Unpack element buffer and create local elements            (public) |
 |                                                          m.kue 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Discretization::UnPackMyElements(RCP<std::vector<char> > e)
{
  std::vector<char>::size_type index = 0;
  while (index < e->size())
  {
    std::vector<char> data;
    ParObject::ExtractfromPack(index,*e,data);
    DRT::ParObject* o = DRT::UTILS::Factory(data);
    DRT::Element* ele = dynamic_cast<Element*>(o);
    if (ele == NULL)
    {
      dserror("Failed to build an element from the element data");
    }
    ele->SetOwner(comm_->MyPID());
    AddElement(Teuchos::rcp(ele));
  }
  // in case AddElement forgets...
  Reset();
}

/*----------------------------------------------------------------------*
 |  Unpack nodal buffer and create local nodes                 (public) |
 |                                                          m.kue 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Discretization::UnPackMyNodes(RCP<std::vector<char> > e)
{
  std::vector<char>::size_type index = 0;
  while (index < e->size())
  {
    std::vector<char> data;
    ParObject::ExtractfromPack(index,*e,data);
    DRT::ParObject* o = DRT::UTILS::Factory(data);
    DRT::Node* n = dynamic_cast<Node*>(o);
    if (n == NULL)
    {
      dserror("Failed to build a node from the node data");
    }
    n->SetOwner(comm_->MyPID());
    AddNode(Teuchos::rcp(n));
  }
  // in case AddNode forgets...
  Reset();
}

/*----------------------------------------------------------------------*
 |  Unpack buffer and create local condition                   (public) |
 |                                                          a.ger 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Discretization::UnPackCondition(
        const RCP<std::vector<char> > e,
        const std::string condname)
{
  std::vector<char>::size_type index = 0;
  while (index < e->size())
  {
    std::vector<char> data;
    DRT::ParObject::ExtractfromPack(index,*e,data);
    DRT::ParObject* o = DRT::UTILS::Factory(data);
    DRT::Condition* cond = dynamic_cast<DRT::Condition*>(o);
    if (cond == NULL)
    {
      dserror("Failed to build boundary condition from the stored data %s", condname.c_str());
    }
    SetCondition(condname,Teuchos::rcp(cond));
  }
  Reset();
}

