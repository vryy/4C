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
DRT::Discretization::Discretization(
    const std::string& name, Teuchos::RCP<Epetra_Comm> comm, unsigned int n_dim)
    : name_(name),
      comm_(comm),
      writer_(Teuchos::null),
      filled_(false),
      havedof_(false),
      n_dim_(n_dim)
{
  dofsets_.emplace_back(Teuchos::rcp(new CORE::Dofsets::DofSet()));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::add_element(Teuchos::RCP<CORE::Elements::Element> ele)
{
  element_[ele->Id()] = ele;
  reset();
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
  if (!globalfilled) reset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::AddNode(Teuchos::RCP<CORE::Nodes::Node> node)
{
  node_[node->Id()] = node;
  reset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::Discretization::DeleteNode(Teuchos::RCP<CORE::Nodes::Node> node)
{
  auto it_node = node_.find(node->Id());
  if (it_node == node_.end()) return false;
  node_.erase(it_node);
  reset();
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::Discretization::DeleteNode(const int gid)
{
  auto it_node = node_.find(gid);
  if (it_node == node_.end()) return false;
  node_.erase(it_node);
  reset();
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::Discretization::DeleteNodes()
{
  node_.clear();
  reset();
  CheckFilledGlobally();
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::Discretization::DeleteElements()
{
  element_.clear();
  reset();
  CheckFilledGlobally();
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::Discretization::DeleteElement(Teuchos::RCP<CORE::Elements::Element> ele)
{
  auto it_ele = element_.find(ele->Id());
  if (it_ele == element_.end()) return false;
  element_.erase(it_ele);
  reset();
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::Discretization::DeleteElement(const int gid)
{
  auto it_ele = element_.find(gid);
  if (it_ele == element_.end()) return false;
  element_.erase(it_ele);
  reset();
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::Discretization::ClearDiscret()
{
  element_.clear();
  node_.clear();
  condition_.clear();
  reset();
  CheckFilledGlobally();
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::NodeRowMap() const
{
  FOUR_C_ASSERT(
      Filled(), "fill_complete() must be called before for discretization %s!", name_.c_str());
  return noderowmap_.get();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::NodeColMap() const
{
  FOUR_C_ASSERT(
      Filled(), "fill_complete() must be called before for discretization %s!", name_.c_str());
  return nodecolmap_.get();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::ElementRowMap() const
{
  FOUR_C_ASSERT(
      Filled(), "fill_complete() must be called before for discretization %s!", name_.c_str());
  return elerowmap_.get();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::ElementColMap() const
{
  FOUR_C_ASSERT(
      Filled(), "fill_complete() must be called before for discretization %s!", name_.c_str());
  return elecolmap_.get();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::Discretization::NumGlobalElements() const
{
  FOUR_C_ASSERT(
      Filled(), "fill_complete() must be called before for discretization %s!", name_.c_str());
  return ElementRowMap()->NumGlobalElements();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::Discretization::NumMyRowElements() const
{
  FOUR_C_ASSERT(
      Filled(), "fill_complete() must be called before for discretization %s!", name_.c_str());
  return ElementRowMap()->NumMyElements();
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
  FOUR_C_ASSERT(
      Filled(), "fill_complete() must be called before for discretization %s!", name_.c_str());
  return NodeRowMap()->NumGlobalElements();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::Discretization::NumMyRowNodes() const
{
  FOUR_C_ASSERT(
      Filled(), "fill_complete() must be called before for discretization %s!", name_.c_str());
  return NodeRowMap()->NumMyElements();
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
CORE::Elements::Element* DRT::Discretization::gElement(const int gid) const
{
  std::map<int, Teuchos::RCP<CORE::Elements::Element>>::const_iterator curr = element_.find(gid);
  FOUR_C_ASSERT(
      curr != element_.end(), "Element with global id gid=%d not stored on this proc!", gid);
  return curr->second.get();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::Discretization::HaveGlobalNode(const int gid) const
{
  return node_.find(gid) != node_.end();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::Nodes::Node* DRT::Discretization::gNode(int gid) const
{
  std::map<int, Teuchos::RCP<CORE::Nodes::Node>>::const_iterator curr = node_.find(gid);
  FOUR_C_ASSERT(curr != node_.end(), "Node with global id gid=%d not stored on this proc!", gid);
  return curr->second.get();
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
    std::map<int, Teuchos::RCP<CORE::Nodes::Node>>::const_iterator ncurr;
    for (ncurr = node_.begin(); ncurr != node_.end(); ++ncurr)
      if (ncurr->second->Owner() == Comm().MyPID()) nummynodes++;

    int nummyele = 0;
    std::map<int, Teuchos::RCP<CORE::Elements::Element>>::const_iterator ecurr;
    for (ecurr = element_.begin(); ecurr != element_.end(); ++ecurr)
      if (ecurr->second->Owner() == Comm().MyPID()) nummyele++;

    Comm().SumAll(&nummynodes, &numglobalnodes, 1);
    Comm().SumAll(&nummyele, &numglobalelements, 1);
  }

  // print head
  if (Comm().MyPID() == 0)
  {
    os << "--------------------------------------------------\n";
    os << "discretization: " << Name() << std::endl;
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
          std::map<int, Teuchos::RCP<CORE::Elements::Element>>::const_iterator curr;
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
          std::map<int, Teuchos::RCP<CORE::Nodes::Node>>::const_iterator curr;
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
const Epetra_Map* DRT::Discretization::dof_row_map(const unsigned nds) const
{
  FOUR_C_ASSERT(
      nds < dofsets_.size(), "undefined dof set found in discretization %s!", name_.c_str());
  FOUR_C_THROW_UNLESS(
      Filled(), "fill_complete was not called on discretization %s!", name_.c_str());
  FOUR_C_THROW_UNLESS(
      HaveDofs(), "assign_degrees_of_freedom() not called on discretization %s!", name_.c_str());

  return dofsets_[nds]->dof_row_map();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::Discretization::DofColMap(const unsigned nds) const
{
  FOUR_C_ASSERT(
      nds < dofsets_.size(), "undefined dof set found in discretization %s!", name_.c_str());
  FOUR_C_THROW_UNLESS(
      Filled(), "fill_complete was not called on discretization %s!", name_.c_str());
  FOUR_C_THROW_UNLESS(
      HaveDofs(), "assign_degrees_of_freedom() not called on discretization %s!", name_.c_str());

  return dofsets_[nds]->DofColMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::ReplaceDofSet(const unsigned nds,
    Teuchos::RCP<CORE::Dofsets::DofSetInterface> newdofset, const bool replaceinstatdofsets)
{
  FOUR_C_ASSERT(
      nds < dofsets_.size(), "undefined dof set found in discretization %s!", name_.c_str());
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
  FOUR_C_ASSERT(
      nds < (int)dofsets_.size(), "undefined dof set found in discretization %s!", name_.c_str());
  return Teuchos::rcp(new CORE::Dofsets::DofSetProxy(&*dofsets_[nds]));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::ReplaceDofSet(
    Teuchos::RCP<CORE::Dofsets::DofSetInterface> newdofset, const bool replaceinstatdofsets)
{
  FOUR_C_ASSERT(dofsets_.size() == 1, "Discretization %s expects just one dof set!", name_.c_str());
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
void DRT::Discretization::set_state(
    const unsigned nds, const std::string& name, Teuchos::RCP<const Epetra_Vector> state)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::set_state");

  FOUR_C_THROW_UNLESS(
      HaveDofs(), "fill_complete() was not called for discretization %s!", name_.c_str());
  const Epetra_Map* colmap = DofColMap(nds);
  const Epetra_BlockMap& vecmap = state->Map();

  if (state_.size() <= nds) state_.resize(nds + 1);

  // if it's already in column map just set a reference
  // This is a rough test, but it might be ok at this place. It is an
  // error anyway to hand in a vector that is not related to our dof
  // maps.
  if (vecmap.PointSameAs(*colmap))
  {
    FOUR_C_ASSERT(colmap->SameAs(vecmap),
        "col map of discretization %s and state vector %s are different. This is a fatal bug!",
        name_.c_str(), name.c_str());
    // make a copy as in parallel such that no additional RCP points to the state vector
    Teuchos::RCP<Epetra_Vector> tmp = CORE::LINALG::CreateVector(*colmap, false);
    tmp->Update(1.0, *state, 0.0);
    state_[nds][name] = tmp;
  }
  else  // if it's not in column map export and allocate
  {
    FOUR_C_ASSERT(dof_row_map(nds)->SameAs(state->Map()),
        "row map of discretization %s and state vector %s are different. This is a fatal bug!",
        name_.c_str(), name.c_str());
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
    FOUR_C_THROW_UNLESS(
        !err, "Export using importer failed for Epetra_Vector: return value = %d", err);

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
  FOUR_C_THROW_UNLESS(
      count == num, "Mismatch in number of conditions found in discretization %s!", name_.c_str());
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
  FOUR_C_THROW_UNLESS(
      count == num, "Mismatch in number of conditions found in discretization %s!", name_.c_str());
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
  FOUR_C_THROW_UNLESS(
      Filled(), "fill_complete was not called on discretization %s!", name_.c_str());

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
  FOUR_C_THROW_UNLESS(
      Filled(), "fill_complete was not called on discretization %s!", name_.c_str());

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
    auto* ele = dynamic_cast<CORE::Elements::Element*>(o);
    FOUR_C_THROW_UNLESS(ele != nullptr,
        "Failed to build an element from the element data for discretization %s", name_.c_str());
    ele->SetOwner(comm_->MyPID());
    add_element(Teuchos::rcp(ele));
  }
  // in case add_element forgets...
  reset();
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
    auto* node = dynamic_cast<CORE::Nodes::Node*>(o);
    FOUR_C_THROW_UNLESS(node != nullptr,
        "Failed to build a node from the node data for discretization %s", name_.c_str());
    node->SetOwner(comm_->MyPID());
    AddNode(Teuchos::rcp(node));
  }
  // in case AddNode forgets...
  reset();
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
    auto statevecrowmap = CORE::LINALG::CreateVector(*dof_row_map(nds), true);
    CORE::LINALG::Export(*statevec, *statevecrowmap);

    // now set the state again
    set_state(nds, name, statevecrowmap);
  }
}

FOUR_C_NAMESPACE_CLOSE
