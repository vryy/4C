/*---------------------------------------------------------------------*/
/*! \file

\brief a class to manage one discretization

\level 0

*/
/*---------------------------------------------------------------------*/

#include <algorithm>
#include <Teuchos_TimeMonitor.hpp>
#include <utility>

#include "lib_discret.H"
#include "lib_exporter.H"
#include "utils_exceptions.H"
#include "linalg_utils_sparse_algebra_create.H"
#include "linalg_utils_sparse_algebra_manipulation.H"
#include "lib_dofset_proxy.H"
#include "lib_dofset_pbc.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Discretization::Discretization(const std::string& name, Teuchos::RCP<Epetra_Comm> comm)
    : name_(name), comm_(comm), writer_(Teuchos::null), filled_(false), havedof_(false)
{
  dofsets_.emplace_back(Teuchos::rcp(new DofSet()));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::AddElement(Teuchos::RCP<DRT::Element> ele)
{
  element_[ele->Id()] = ele;
  Reset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::CheckFilledGlobally()
{
  // global filled flag (is true / one if and only if filled_ == true on each processor
  int globalfilled = 0;

  // convert filled_ flag on this procesor  into integer (no Epetra communicator for type bool)
  int localfilled = (int)filled_;

  /*the global filled flag is set to the minimal value of any local filled flag
   * i.e. if on any processor filled_ == false, the flag globalfilled is set to
   * zero*/
  Comm().MinAll(&localfilled, &globalfilled, 1);

  // if not Filled() == true on all the processors call Reset()
  if (!globalfilled) Reset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::AddNode(Teuchos::RCP<DRT::Node> node)
{
  node_[node->Id()] = node;
  Reset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::Discretization::DeleteNode(Teuchos::RCP<DRT::Node> node)
{
  auto it_node = node_.find(node->Id());
  if (it_node == node_.end()) return false;
  node_.erase(it_node);
  Reset();
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::Discretization::DeleteNode(const int gid)
{
  auto it_node = node_.find(gid);
  if (it_node == node_.end()) return false;
  node_.erase(it_node);
  Reset();
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::Discretization::DeleteNodes()
{
  node_.clear();
  Reset();
  CheckFilledGlobally();
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::Discretization::DeleteElements()
{
  element_.clear();
  Reset();
  CheckFilledGlobally();
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::Discretization::DeleteElement(Teuchos::RCP<DRT::Element> ele)
{
  auto it_ele = element_.find(ele->Id());
  if (it_ele == element_.end()) return false;
  element_.erase(it_ele);
  Reset();
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::Discretization::DeleteElement(const int gid)
{
  auto it_ele = element_.find(gid);
  if (it_ele == element_.end()) return false;
  element_.erase(it_ele);
  Reset();
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::Discretization::ClearDiscret()
{
  element_.clear();
  node_.clear();
  condition_.clear();
  Reset();
  CheckFilledGlobally();
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::NodeRowMap() const
{
#ifdef DEBUG
  if (Filled())
    return noderowmap_.get();
  else
    dserror("FillComplete() must be called before call to NodeRowMap()");
  return NULL;
#else
  return noderowmap_.get();
#endif
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::NodeColMap() const
{
#ifdef DEBUG
  if (Filled())
    return nodecolmap_.get();
  else
    dserror("FillComplete() must be called before call to NodeColMap()");
  return NULL;
#else
  return nodecolmap_.get();
#endif
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::ElementRowMap() const
{
#ifdef DEBUG
  if (Filled())
    return elerowmap_.get();
  else
    dserror("FillComplete() must be called before call to ElementRowMap()");
  return NULL;
#else
  return elerowmap_.get();
#endif
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::ElementColMap() const
{
#ifdef DEBUG
  if (Filled())
    return elecolmap_.get();
  else
    dserror("FillComplete() must be called before call to ElementColMap()");
  return NULL;
#else
  return elecolmap_.get();
#endif
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::Discretization::NumGlobalElements() const
{
#ifdef DEBUG
  if (Filled())
    return ElementRowMap()->NumGlobalElements();
  else
    dserror("FillComplete() must be called before call to NumGlobalElements()");
  return -1;
#else
  return ElementRowMap()->NumGlobalElements();
#endif
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::Discretization::NumMyRowElements() const
{
#ifdef DEBUG
  if (Filled())
    return ElementRowMap()->NumMyElements();
  else
    dserror("FillComplete() must be called before call to NumMyRowElements()");
  return -1;
#else
  return ElementRowMap()->NumMyElements();
#endif
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::Discretization::NumMyColElements() const
{
  if (Filled())
    return ElementColMap()->NumMyElements();
  else
    return (int)element_.size();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::Discretization::NumGlobalNodes() const
{
#ifdef DEBUG
  if (Filled())
    return NodeRowMap()->NumGlobalElements();
  else
    dserror("FillComplete() must be called before call to NumGlobalNodes()");
  return -1;
#else
  return NodeRowMap()->NumGlobalElements();
#endif
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::Discretization::NumMyRowNodes() const
{
#ifdef DEBUG
  if (Filled())
    return NodeRowMap()->NumMyElements();
  else
    dserror("FillComplete() must be called before call to NumMyRowNodes()");
  return -1;
#else
  return NodeRowMap()->NumMyElements();
#endif
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::Discretization::NumMyColNodes() const
{
  if (Filled())
    return NodeColMap()->NumMyElements();
  else
    return (int)node_.size();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::Discretization::HaveGlobalElement(const int gid) const
{
  return element_.find(gid) != element_.end();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Discretization::gElement(const int gid) const
{
#ifdef DEBUG
  std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator curr = element_.find(gid);
  if (curr == element_.end())
    dserror("Element with gobal id gid=%d not stored on this proc", gid);
  else
    return curr->second.get();
  return NULL;
#else
  return element_.find(gid)->second.get();
#endif
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::Discretization::HaveGlobalNode(const int gid) const
{
  return node_.find(gid) != node_.end();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Node* DRT::Discretization::gNode(int gid) const
{
#ifdef DEBUG
  std::map<int, Teuchos::RCP<DRT::Node>>::const_iterator curr = node_.find(gid);
  if (curr == node_.end())
    dserror("Node with global id gid=%d not stored on this proc", gid);
  else
    return curr->second.get();
  return NULL;
#else
  return node_.find(gid)->second.get();
#endif
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const DRT::Discretization& dis)
{
  dis.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::Print(std::ostream& os) const
{
  int numglobalelements = 0;
  int numglobalnodes = 0;
  if (Filled())
  {
    numglobalelements = NumGlobalElements();
    numglobalnodes = NumGlobalNodes();
  }
  else
  {
    int nummynodes = 0;
    std::map<int, Teuchos::RCP<DRT::Node>>::const_iterator ncurr;
    for (ncurr = node_.begin(); ncurr != node_.end(); ++ncurr)
      if (ncurr->second->Owner() == Comm().MyPID()) nummynodes++;

    int nummyele = 0;
    std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator ecurr;
    for (ecurr = element_.begin(); ecurr != element_.end(); ++ecurr)
      if (ecurr->second->Owner() == Comm().MyPID()) nummyele++;

    Comm().SumAll(&nummynodes, &numglobalnodes, 1);
    Comm().SumAll(&nummyele, &numglobalelements, 1);
  }

  // print head
  if (Comm().MyPID() == 0)
  {
    os << "--------------------------------------------------\n";
    os << "Discretization: " << Name() << std::endl;
    os << "--------------------------------------------------\n";
    os << numglobalelements << " Elements " << numglobalnodes << " Nodes (global)\n";
    os << "--------------------------------------------------\n";
    if (Filled())
      os << "Filled() = true\n";
    else
      os << "Filled() = false\n";
    os << "--------------------------------------------------\n";
  }
  Comm().Barrier();
  for (int proc = 0; proc < Comm().NumProc(); ++proc)
  {
    if (proc == Comm().MyPID())
    {
      // loop over dofsets
      for (int nds = 0; nds < NumDofSets(); ++nds)
      {
        os << "\n------------------------ Dofset " << nds << " :\n\n";
        // print elements
        {
          os << "-------------------------- Proc " << proc << " :\n";
          std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator curr;
          for (curr = element_.begin(); curr != element_.end(); ++curr)
          {
            os << *(curr->second);
            if (Filled() && HaveDofs())
            {
              std::vector<int> dof = Dof(nds, &*(curr->second));
              if (dof.size())
              {
                os << " Dofs ";
                for (int i : dof) os << std::setw(6) << i << " ";
              }
            }
            os << std::endl;
          }
          os << std::endl;
        }
        // print nodes
        {
          os << "-------------------------- Proc " << proc << " :\n";
          std::map<int, Teuchos::RCP<DRT::Node>>::const_iterator curr;
          for (curr = node_.begin(); curr != node_.end(); ++curr)
          {
            os << *(curr->second);
            if (Filled() && HaveDofs())
            {
              std::vector<int> dof = Dof(nds, &*(curr->second));
              if (dof.size())
              {
                os << " Dofs ";
                for (int i : dof) os << std::setw(6) << i << " ";
              }
            }
            os << std::endl;
          }
          os << std::endl;
        }
      }
      // print conditions
      {
        const unsigned numcond = condition_.size();
        if (numcond) os << "-------------------------- Proc " << proc << " :\n";
        if (numcond)
        {
          os << numcond << " Conditions:\n";
          std::map<std::string, Teuchos::RCP<Condition>>::const_iterator curr;
          for (curr = condition_.begin(); curr != condition_.end(); ++curr)
          {
            os << curr->first << " ";
            os << *(curr->second) << std::endl;
          }
        }
        os << std::endl;
      }
    }
    Comm().Barrier();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::DofRowMap(const unsigned nds) const
{
  dsassert(nds < dofsets_.size(), "undefined dof set");
  if (!Filled()) dserror("FillComplete was not called on this discretization");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() not called on this discretization");

  return dofsets_[nds]->DofRowMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::DofColMap(const unsigned nds) const
{
  dsassert(nds < dofsets_.size(), "undefined dof set");
  if (!Filled()) dserror("FillComplete was not called on this discretization");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() not called on this discretization");

  return dofsets_[nds]->DofColMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::ReplaceDofSet(
    const unsigned nds, Teuchos::RCP<DofSetInterface> newdofset, const bool replaceinstatdofsets)
{
  dsassert(nds < dofsets_.size(), "undefined dof set");
  // if we already have our dofs here and we add a properly filled (proxy)
  // DofSet, we do not need (and do not want) to refill.
  havedof_ = havedof_ and newdofset->Filled() and nds != 0;
  if (replaceinstatdofsets) newdofset->ReplaceInStaticDofsets(dofsets_[nds]);
  dofsets_[nds] = newdofset;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::Discretization::AddDofSet(Teuchos::RCP<DofSetInterface> newdofset)
{
  // if we already have our dofs here and we add a properly filled (proxy)
  // DofSet, we do not need (and do not want) to refill.
  havedof_ = havedof_ and newdofset->Filled();
  dofsets_.push_back(newdofset);
  return static_cast<int>(dofsets_.size() - 1);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::DofSetInterface> DRT::Discretization::GetDofSetProxy(const int nds)
{
  dsassert(nds < (int)dofsets_.size(), "undefined dof set");
  return Teuchos::rcp(new DofSetProxy(&*dofsets_[nds]));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::ReplaceDofSet(
    Teuchos::RCP<DofSetInterface> newdofset, const bool replaceinstatdofsets)
{
  dsassert(dofsets_.size() == 1, "expect just one dof set");
  havedof_ = false;
  if (replaceinstatdofsets) newdofset->ReplaceInStaticDofsets(dofsets_[0]);
  dofsets_[0] = newdofset;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::map<int, std::vector<int>>> DRT::Discretization::GetAllPBCCoupledColNodes()
{
  // check for pbcs
  for (int nds = 0; nds < NumDofSets(); nds++)
  {
    Teuchos::RCP<PBCDofSet> pbcdofset = Teuchos::rcp_dynamic_cast<PBCDofSet>(dofsets_[nds]);

    if (pbcdofset != Teuchos::null)
    {
      // it is assumed that, if one pbc set is available, all other potential dofsets hold the same
      // layout
      return pbcdofset->GetCoupledNodes();
    }
  }

  return Teuchos::rcp(new std::map<int, std::vector<int>>);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::map<int, int>> DRT::Discretization::GetPBCSlaveToMasterNodeConnectivity()
{
  // check for pbcs
  for (int nds = 0; nds < NumDofSets(); nds++)
  {
    Teuchos::RCP<PBCDofSet> pbcdofset = Teuchos::rcp_dynamic_cast<PBCDofSet>(dofsets_[nds]);

    if (pbcdofset != Teuchos::null)
    {
      // it is assumed that, if one pbc set is available, all other potential dofsets hold the same
      // layout
      return pbcdofset->GetSlaveToMasterNodeConnectivity();
    }
  }

  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::SetState(
    const unsigned nds, const std::string& name, Teuchos::RCP<const Epetra_Vector> state)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::SetState");

  if (!HaveDofs()) dserror("FillComplete() was not called");
  const Epetra_Map* colmap = DofColMap(nds);
  const Epetra_BlockMap& vecmap = state->Map();

  if (state_.size() <= nds) state_.resize(nds + 1);

  // if it's already in column map just set a reference
  // This is a rough test, but it might be ok at this place. It is an
  // error anyway to hand in a vector that is not related to our dof
  // maps.
  if (vecmap.PointSameAs(*colmap))
  {
#ifdef DEBUG
    if (not colmap->SameAs(vecmap))
    {
      dserror("col map of discretization and state vector %s are different. This is a fatal bug!",
          name.c_str());
    }
#endif
    // make a copy as in parallel such that no additional RCP points to the state vector
    Teuchos::RCP<Epetra_Vector> tmp = CORE::LINALG::CreateVector(*colmap, false);
    tmp->Update(1.0, *state, 0.0);
    state_[nds][name] = tmp;
  }
  else  // if it's not in column map export and allocate
  {
#ifdef DEBUG
    if (not DofRowMap(nds)->SameAs(state->Map()))
    {
      dserror("row map of discretization and state vector %s are different. This is a fatal bug!",
          name.c_str());
    }
#endif
    Teuchos::RCP<Epetra_Vector> tmp = CORE::LINALG::CreateVector(*colmap, false);

    // this is necessary to find out the number of nodesets in the beginning
    if (stateimporter_.size() <= nds)
    {
      stateimporter_.resize(nds + 1);
      for (unsigned i = 0; i <= nds; ++i) stateimporter_[i] = Teuchos::null;
    }
    // (re)build importer if necessary
    if (stateimporter_[nds] == Teuchos::null or
        not stateimporter_[nds]->SourceMap().SameAs(state->Map()) or
        not stateimporter_[nds]->TargetMap().SameAs(*colmap))
    {
      stateimporter_[nds] = Teuchos::rcp(new Epetra_Import(*colmap, state->Map()));
    }

    // transfer data
    int err = tmp->Import(*state, (*stateimporter_[nds]), Insert);
    if (err) dserror("Export using importer failed for Epetra_Vector: return value = %d", err);

    // save state
    state_[nds][name] = tmp;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::SetCondition(const std::string& name, Teuchos::RCP<Condition> cond)
{
  condition_.insert(std::pair<std::string, Teuchos::RCP<Condition>>(name, cond));
  filled_ = false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::ReplaceConditions(
    const std::string& name, const std::vector<Teuchos::RCP<Condition>>& conds)
{
  if (condition_.count(name) > 0) condition_.erase(name);

  std::vector<Teuchos::RCP<Condition>>::const_iterator cit;
  for (cit = conds.begin(); cit != conds.end(); ++cit)
  {
    // skip null pointers (these conditions will be deleted only and
    // therefore may disappear completely from this discretization)
    if (not cit->is_null())
      condition_.insert(std::pair<std::string, Teuchos::RCP<Condition>>(name, *cit));
  }
  filled_ = false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::GetCondition(
    const std::string& name, std::vector<DRT::Condition*>& out) const
{
  const unsigned num = condition_.count(name);
  out.resize(num);
  unsigned count = 0;

  auto range = condition_.equal_range(name);
  for (auto cond = range.first; cond != range.second; ++cond)
  {
    out[count++] = cond->second.get();
  }
  if (count != num) dserror("Mismatch in number of conditions found");
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::GetCondition(
    const std::string& name, std::vector<Teuchos::RCP<Condition>>& out) const
{
  const unsigned num = condition_.count(name);
  out.resize(num);
  unsigned count = 0;

  auto range = condition_.equal_range(name);
  for (auto cond = range.first; cond != range.second; ++cond)
  {
    out[count++] = cond->second;
  }
  if (count != num) dserror("Mismatch in number of conditions found");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Condition* DRT::Discretization::GetCondition(const std::string& name) const
{
  auto it_cond = condition_.find(name);
  if (it_cond == condition_.end()) return nullptr;
  it_cond = condition_.lower_bound(name);
  return it_cond->second.get();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::GetConditionNames(std::vector<std::string>& names) const
{
  std::set<std::string> n;
  for (const auto& [name, cond] : condition_) n.insert(name);

  names.reserve(n.size());
  names.assign(n.begin(), n.end());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<char>> DRT::Discretization::PackMyElements() const
{
  if (!Filled()) dserror("FillComplete was not called on this discretization");

  DRT::PackBuffer buffer;

  for (auto* ele : elerowptr_) ele->Pack(buffer);

  buffer.StartPacking();

  for (auto* ele : elerowptr_) ele->Pack(buffer);

  auto block = Teuchos::rcp(new std::vector<char>);
  std::swap(*block, buffer());
  return block;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<char>> DRT::Discretization::PackMyNodes() const
{
  if (!Filled()) dserror("FillComplete was not called on this discretization");

  DRT::PackBuffer buffer;

  for (auto* node : noderowptr_) node->Pack(buffer);

  buffer.StartPacking();

  for (auto* node : noderowptr_) node->Pack(buffer);

  auto block = Teuchos::rcp(new std::vector<char>);
  std::swap(*block, buffer());
  return block;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<char>> DRT::Discretization::PackCondition(
    const std::string& condname) const
{
  if (!Filled()) dserror("FillComplete was not called on this discretization");

  // get boundary conditions
  std::vector<DRT::Condition*> cond;
  GetCondition(condname, cond);

  DRT::PackBuffer buffer;

  for (auto* c : cond) c->Pack(buffer);

  buffer.StartPacking();

  for (auto* c : cond) c->Pack(buffer);

  Teuchos::RCP<std::vector<char>> block = Teuchos::rcp(new std::vector<char>);
  std::swap(*block, buffer());
  return block;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::UnPackMyElements(Teuchos::RCP<std::vector<char>> e)
{
  std::vector<char>::size_type index = 0;
  while (index < e->size())
  {
    std::vector<char> data;
    ParObject::ExtractfromPack(index, *e, data);
    DRT::ParObject* o = DRT::UTILS::Factory(data);
    auto* ele = dynamic_cast<Element*>(o);
    if (ele == nullptr)
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
 *----------------------------------------------------------------------*/
void DRT::Discretization::UnPackMyNodes(Teuchos::RCP<std::vector<char>> e)
{
  std::vector<char>::size_type index = 0;
  while (index < e->size())
  {
    std::vector<char> data;
    ParObject::ExtractfromPack(index, *e, data);
    DRT::ParObject* o = DRT::UTILS::Factory(data);
    auto* node = dynamic_cast<Node*>(o);
    if (node == nullptr)
    {
      dserror("Failed to build a node from the node data");
    }
    node->SetOwner(comm_->MyPID());
    AddNode(Teuchos::rcp(node));
  }
  // in case AddNode forgets...
  Reset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::UnPackCondition(
    const Teuchos::RCP<std::vector<char>> e, const std::string& condname)
{
  std::vector<char>::size_type index = 0;
  while (index < e->size())
  {
    std::vector<char> data;
    DRT::ParObject::ExtractfromPack(index, *e, data);
    DRT::ParObject* o = DRT::UTILS::Factory(data);
    auto* cond = dynamic_cast<DRT::Condition*>(o);
    if (cond == nullptr)
    {
      dserror("Failed to build boundary condition from the stored data %s", condname.c_str());
    }
    SetCondition(condname, Teuchos::rcp(cond));
  }
  Reset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::RedistributeState(const unsigned nds, const std::string& name)
{
  // only redistribute if state has been set
  if (HasState(nds, name))
  {
    // get the state and export it to the rowmap to be able to reset the state
    auto statevec = GetState(nds, name);
    auto statevecrowmap = CORE::LINALG::CreateVector(*DofRowMap(nds), true);
    CORE::LINALG::Export(*statevec, *statevecrowmap);

    // now set the state again
    SetState(nds, name, statevecrowmap);
  }
}
