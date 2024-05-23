/*---------------------------------------------------------------------*/
/*! \file

\brief a class to manage one discretization

\level 0

*/
/*---------------------------------------------------------------------*/

#include "4C_lib_discret.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_discretization_dofset_pbc.hpp"
#include "4C_discretization_dofset_proxy.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_TimeMonitor.hpp>

#include <algorithm>
#include <utility>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Discretization::Discretization(const std::string& name, Teuchos::RCP<Epetra_Comm> comm)
    : name_(name), comm_(comm), writer_(Teuchos::null), filled_(false), havedof_(false)
{
  dofsets_.emplace_back(Teuchos::rcp(new CORE::Dofsets::DofSet()));
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
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (Filled())
    return noderowmap_.get();
  else
    FOUR_C_THROW("FillComplete() must be called before call to NodeRowMap()");
  return nullptr;
#else
  return noderowmap_.get();
#endif
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::NodeColMap() const
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (Filled())
    return nodecolmap_.get();
  else
    FOUR_C_THROW("FillComplete() must be called before call to NodeColMap()");
  return nullptr;
#else
  return nodecolmap_.get();
#endif
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::ElementRowMap() const
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (Filled())
    return elerowmap_.get();
  else
    FOUR_C_THROW("FillComplete() must be called before call to ElementRowMap()");
  return nullptr;
#else
  return elerowmap_.get();
#endif
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::ElementColMap() const
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (Filled())
    return elecolmap_.get();
  else
    FOUR_C_THROW("FillComplete() must be called before call to ElementColMap()");
  return nullptr;
#else
  return elecolmap_.get();
#endif
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::Discretization::NumGlobalElements() const
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (Filled())
    return ElementRowMap()->NumGlobalElements();
  else
    FOUR_C_THROW("FillComplete() must be called before call to NumGlobalElements()");
  return -1;
#else
  return ElementRowMap()->NumGlobalElements();
#endif
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::Discretization::NumMyRowElements() const
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (Filled())
    return ElementRowMap()->NumMyElements();
  else
    FOUR_C_THROW("FillComplete() must be called before call to NumMyRowElements()");
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
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (Filled())
    return NodeRowMap()->NumGlobalElements();
  else
    FOUR_C_THROW("FillComplete() must be called before call to NumGlobalNodes()");
  return -1;
#else
  return NodeRowMap()->NumGlobalElements();
#endif
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::Discretization::NumMyRowNodes() const
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (Filled())
    return NodeRowMap()->NumMyElements();
  else
    FOUR_C_THROW("FillComplete() must be called before call to NumMyRowNodes()");
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
#ifdef FOUR_C_ENABLE_ASSERTIONS
  std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator curr = element_.find(gid);
  if (curr == element_.end())
    FOUR_C_THROW("Element with gobal id gid=%d not stored on this proc", gid);
  else
    return curr->second.get();
  return nullptr;
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
#ifdef FOUR_C_ENABLE_ASSERTIONS
  std::map<int, Teuchos::RCP<DRT::Node>>::const_iterator curr = node_.find(gid);
  if (curr == node_.end())
    FOUR_C_THROW("Node with global id gid=%d not stored on this proc", gid);
  else
    return curr->second.get();
  return nullptr;
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
          std::map<std::string, Teuchos::RCP<CORE::Conditions::Condition>>::const_iterator curr;
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
  FOUR_C_ASSERT(nds < dofsets_.size(), "undefined dof set");
  if (!Filled()) FOUR_C_THROW("FillComplete was not called on this discretization");
  if (!HaveDofs()) FOUR_C_THROW("assign_degrees_of_freedom() not called on this discretization");

  return dofsets_[nds]->DofRowMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::DofColMap(const unsigned nds) const
{
  FOUR_C_ASSERT(nds < dofsets_.size(), "undefined dof set");
  if (!Filled()) FOUR_C_THROW("FillComplete was not called on this discretization");
  if (!HaveDofs()) FOUR_C_THROW("assign_degrees_of_freedom() not called on this discretization");

  return dofsets_[nds]->DofColMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::ReplaceDofSet(const unsigned nds,
    Teuchos::RCP<CORE::Dofsets::DofSetInterface> newdofset, const bool replaceinstatdofsets)
{
  FOUR_C_ASSERT(nds < dofsets_.size(), "undefined dof set");
  // if we already have our dofs here and we add a properly filled (proxy)
  // DofSet, we do not need (and do not want) to refill.
  havedof_ = havedof_ and newdofset->Filled() and nds != 0;
  if (replaceinstatdofsets) newdofset->replace_in_static_dofsets(dofsets_[nds]);
  dofsets_[nds] = newdofset;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::Discretization::AddDofSet(Teuchos::RCP<CORE::Dofsets::DofSetInterface> newdofset)
{
  // if we already have our dofs here and we add a properly filled (proxy)
  // DofSet, we do not need (and do not want) to refill.
  havedof_ = havedof_ and newdofset->Filled();
  dofsets_.push_back(newdofset);
  return static_cast<int>(dofsets_.size() - 1);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::Dofsets::DofSetInterface> DRT::Discretization::GetDofSetProxy(const int nds)
{
  FOUR_C_ASSERT(nds < (int)dofsets_.size(), "undefined dof set");
  return Teuchos::rcp(new CORE::Dofsets::DofSetProxy(&*dofsets_[nds]));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::ReplaceDofSet(
    Teuchos::RCP<CORE::Dofsets::DofSetInterface> newdofset, const bool replaceinstatdofsets)
{
  FOUR_C_ASSERT(dofsets_.size() == 1, "expect just one dof set");
  havedof_ = false;
  if (replaceinstatdofsets) newdofset->replace_in_static_dofsets(dofsets_[0]);
  dofsets_[0] = newdofset;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<int, std::vector<int>>* DRT::Discretization::get_all_pbc_coupled_col_nodes()
{
  // check for pbcs
  for (int nds = 0; nds < NumDofSets(); nds++)
  {
    Teuchos::RCP<CORE::Dofsets::PBCDofSet> pbcdofset =
        Teuchos::rcp_dynamic_cast<CORE::Dofsets::PBCDofSet>(dofsets_[nds]);

    if (pbcdofset != Teuchos::null)
    {
      // it is assumed that, if one pbc set is available, all other potential dofsets hold the same
      // layout
      return pbcdofset->GetCoupledNodes();
    }
  }

  return nullptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::map<int, int>> DRT::Discretization::get_pbc_slave_to_master_node_connectivity()
{
  // check for pbcs
  for (int nds = 0; nds < NumDofSets(); nds++)
  {
    Teuchos::RCP<CORE::Dofsets::PBCDofSet> pbcdofset =
        Teuchos::rcp_dynamic_cast<CORE::Dofsets::PBCDofSet>(dofsets_[nds]);

    if (pbcdofset != Teuchos::null)
    {
      // it is assumed that, if one pbc set is available, all other potential dofsets hold the same
      // layout
      return pbcdofset->get_slave_to_master_node_connectivity();
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

  if (!HaveDofs()) FOUR_C_THROW("FillComplete() was not called");
  const Epetra_Map* colmap = DofColMap(nds);
  const Epetra_BlockMap& vecmap = state->Map();

  if (state_.size() <= nds) state_.resize(nds + 1);

  // if it's already in column map just set a reference
  // This is a rough test, but it might be ok at this place. It is an
  // error anyway to hand in a vector that is not related to our dof
  // maps.
  if (vecmap.PointSameAs(*colmap))
  {
#ifdef FOUR_C_ENABLE_ASSERTIONS
    if (not colmap->SameAs(vecmap))
    {
      FOUR_C_THROW(
          "col map of discretization and state vector %s are different. This is a fatal bug!",
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
#ifdef FOUR_C_ENABLE_ASSERTIONS
    if (not DofRowMap(nds)->SameAs(state->Map()))
    {
      FOUR_C_THROW(
          "row map of discretization and state vector %s are different. This is a fatal bug!",
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
    if (err) FOUR_C_THROW("Export using importer failed for Epetra_Vector: return value = %d", err);

    // save state
    state_[nds][name] = tmp;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::SetCondition(
    const std::string& name, Teuchos::RCP<CORE::Conditions::Condition> cond)
{
  condition_.insert(std::pair<std::string, Teuchos::RCP<CORE::Conditions::Condition>>(name, cond));
  filled_ = false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::ReplaceConditions(
    const std::string& name, const std::vector<Teuchos::RCP<CORE::Conditions::Condition>>& conds)
{
  if (condition_.count(name) > 0) condition_.erase(name);

  std::vector<Teuchos::RCP<CORE::Conditions::Condition>>::const_iterator cit;
  for (cit = conds.begin(); cit != conds.end(); ++cit)
  {
    // skip null pointers (these conditions will be deleted only and
    // therefore may disappear completely from this discretization)
    if (not cit->is_null())
      condition_.insert(
          std::pair<std::string, Teuchos::RCP<CORE::Conditions::Condition>>(name, *cit));
  }
  filled_ = false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::GetCondition(
    const std::string& name, std::vector<CORE::Conditions::Condition*>& out) const
{
  const unsigned num = condition_.count(name);
  out.resize(num);
  unsigned count = 0;

  auto range = condition_.equal_range(name);
  for (auto cond = range.first; cond != range.second; ++cond)
  {
    out[count++] = cond->second.get();
  }
  if (count != num) FOUR_C_THROW("Mismatch in number of conditions found");
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::GetCondition(
    const std::string& name, std::vector<Teuchos::RCP<CORE::Conditions::Condition>>& out) const
{
  const unsigned num = condition_.count(name);
  out.resize(num);
  unsigned count = 0;

  auto range = condition_.equal_range(name);
  for (auto cond = range.first; cond != range.second; ++cond)
  {
    out[count++] = cond->second;
  }
  if (count != num) FOUR_C_THROW("Mismatch in number of conditions found");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::Conditions::Condition* DRT::Discretization::GetCondition(const std::string& name) const
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
  if (!Filled()) FOUR_C_THROW("FillComplete was not called on this discretization");

  CORE::COMM::PackBuffer buffer;

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
  if (!Filled()) FOUR_C_THROW("FillComplete was not called on this discretization");

  CORE::COMM::PackBuffer buffer;

  for (auto* node : noderowptr_) node->Pack(buffer);

  buffer.StartPacking();

  for (auto* node : noderowptr_) node->Pack(buffer);

  auto block = Teuchos::rcp(new std::vector<char>);
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
    CORE::COMM::ParObject::ExtractfromPack(index, *e, data);
    CORE::COMM::ParObject* o = CORE::COMM::Factory(data);
    auto* ele = dynamic_cast<Element*>(o);
    if (ele == nullptr)
    {
      FOUR_C_THROW("Failed to build an element from the element data");
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
    CORE::COMM::ParObject::ExtractfromPack(index, *e, data);
    CORE::COMM::ParObject* o = CORE::COMM::Factory(data);
    auto* node = dynamic_cast<Node*>(o);
    if (node == nullptr)
    {
      FOUR_C_THROW("Failed to build a node from the node data");
    }
    node->SetOwner(comm_->MyPID());
    AddNode(Teuchos::rcp(node));
  }
  // in case AddNode forgets...
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

FOUR_C_NAMESPACE_CLOSE
