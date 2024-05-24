/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of utils on conditions

\level 1


*/
/*----------------------------------------------------------------------*/

#include "4C_discretization_condition_utils.hpp"

#include "4C_discretization_condition_selector.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace
{
  template <typename Range>
  Teuchos::RCP<Epetra_Map> FillConditionMap(
      const DRT::Discretization& dis, const Range& nodeRange, const std::string& condname)
  {
    std::set<int> condnodeset;

    CORE::Conditions::ConditionSelector conds(dis, condname);

    for (const DRT::Node* node : nodeRange)
    {
      if (conds.ContainsNode(node->Id()))
      {
        condnodeset.insert(node->Id());
      }
    }

    Teuchos::RCP<Epetra_Map> condnodemap = CORE::LINALG::CreateMap(condnodeset, dis.Comm());
    return condnodemap;
  }
}  // namespace

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Conditions::FindConditionedNodes(
    const DRT::Discretization& dis, const std::string& condname, std::vector<int>& nodes)
{
  std::vector<Condition*> conds;
  dis.GetCondition(condname, conds);
  FindConditionedNodes(dis, conds, nodes);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Conditions::FindConditionedNodes(
    const DRT::Discretization& dis, const std::string& condname, std::set<int>& nodeset)
{
  std::vector<Condition*> conds;
  dis.GetCondition(condname, conds);
  FindConditionedNodes(dis, conds, nodeset);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Conditions::FindConditionedNodes(
    const DRT::Discretization& dis, const std::string& condname, std::map<int, DRT::Node*>& nodes)
{
  std::vector<Condition*> conds;
  dis.GetCondition(condname, conds);
  FindConditionedNodes(dis, conds, nodes);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Conditions::FindConditionedNodes(
    const DRT::Discretization& dis, const std::vector<Condition*>& conds, std::vector<int>& nodes)
{
  std::set<int> nodeset;
  const int myrank = dis.Comm().MyPID();
  for (const auto& cond : conds)
  {
    for (const auto node : *cond->GetNodes())
    {
      const int gid = node;
      if (dis.HaveGlobalNode(gid) and dis.gNode(gid)->Owner() == myrank)
      {
        nodeset.insert(gid);
      }
    }
  }

  nodes.reserve(nodeset.size());
  nodes.assign(nodeset.begin(), nodeset.end());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Conditions::FindConditionedNodes(
    const DRT::Discretization& dis, const std::vector<Condition*>& conds, std::set<int>& nodeset)
{
  const int myrank = dis.Comm().MyPID();
  for (auto cond : conds)
  {
    for (int gid : *cond->GetNodes())
    {
      if (dis.HaveGlobalNode(gid) and dis.gNode(gid)->Owner() == myrank)
      {
        nodeset.insert(gid);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Conditions::FindConditionedNodes(const DRT::Discretization& dis,
    const std::vector<Condition*>& conds, std::map<int, DRT::Node*>& nodes)
{
  const int myrank = dis.Comm().MyPID();
  for (auto cond : conds)
  {
    for (int gid : *cond->GetNodes())
    {
      if (dis.HaveGlobalNode(gid) and dis.gNode(gid)->Owner() == myrank)
      {
        nodes[gid] = dis.gNode(gid);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Conditions::FindConditionedNodes(const DRT::Discretization& dis,
    const std::vector<Condition*>& conds, std::map<int, Teuchos::RCP<std::vector<int>>>& nodes,
    bool use_coupling_id)
{
  std::map<int, std::set<int>> nodeset;
  const int myrank = dis.Comm().MyPID();
  for (const auto& cond : conds)
  {
    int id = use_coupling_id ? cond->parameters().Get<int>("coupling id") : 0;
    for (int gid : *cond->GetNodes())
    {
      if (dis.HaveGlobalNode(gid) and dis.gNode(gid)->Owner() == myrank)
      {
        nodeset[id].insert(gid);
      }
    }
  }

  for (const auto& [id, gids] : nodeset)
  {
    nodes[id] = Teuchos::rcp(new std::vector<int>(gids.size()));
    nodes[id]->assign(gids.begin(), gids.end());
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Conditions::FindConditionedNodes(const DRT::Discretization& dis,
    const std::vector<Condition*>& conds, std::map<int, std::map<int, DRT::Node*>>& nodes)
{
  const int myrank = dis.Comm().MyPID();
  for (auto* cond : conds)
  {
    int id = cond->parameters().Get<int>("coupling id");
    for (int gid : *cond->GetNodes())
    {
      if (dis.HaveGlobalNode(gid) and dis.gNode(gid)->Owner() == myrank)
      {
        (nodes[id])[gid] = dis.gNode(gid);
      }
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::Conditions::FindConditionObjects(const DRT::Discretization& dis,
    std::map<int, DRT::Node*>& nodes, std::map<int, Teuchos::RCP<DRT::Element>>& elements,
    const std::string& condname)
{
  int myrank = dis.Comm().MyPID();
  std::vector<Condition*> conds;
  dis.GetCondition(condname, conds);

  FindConditionedNodes(dis, conds, nodes);

  for (auto& cond : conds)
  {
    // get this condition's elements
    std::map<int, Teuchos::RCP<DRT::Element>>& geo = cond->Geometry();
    std::map<int, Teuchos::RCP<DRT::Element>>::iterator iter, pos;
    pos = elements.begin();
    for (iter = geo.begin(); iter != geo.end(); ++iter)
    {
      if (iter->second->Owner() == myrank)
      {
        pos = elements.insert(pos, *iter);
      }
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::Conditions::FindConditionObjects(const DRT::Discretization& dis,
    std::map<int, DRT::Node*>& nodes, std::map<int, DRT::Node*>& gnodes,
    std::map<int, Teuchos::RCP<DRT::Element>>& elements, const std::vector<Condition*>& conds)
{
  FindConditionedNodes(dis, conds, nodes);

  for (const auto& cond : conds)
  {
    // get this condition's elements
    std::map<int, Teuchos::RCP<DRT::Element>>& geo = cond->Geometry();
    std::map<int, Teuchos::RCP<DRT::Element>>::iterator iter, pos;
    pos = elements.begin();
    for (iter = geo.begin(); iter != geo.end(); ++iter)
    {
      // get all elements locally known, including ghost elements
      pos = elements.insert(pos, *iter);
      const int* n = iter->second->NodeIds();
      for (int j = 0; j < iter->second->num_node(); ++j)
      {
        const int gid = n[j];
        if (dis.HaveGlobalNode(gid))
        {
          gnodes[gid] = dis.gNode(gid);
        }
        else
          FOUR_C_THROW("All nodes of known elements must be known. Panic.");
      }
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::Conditions::FindConditionObjects(
    std::map<int, Teuchos::RCP<DRT::Element>>& elements, const std::vector<Condition*>& conds)
{
  for (auto cond : conds)
  {
    // get this condition's elements
    std::map<int, Teuchos::RCP<DRT::Element>>& geo = cond->Geometry();
    std::map<int, Teuchos::RCP<DRT::Element>>::iterator iter, pos;
    pos = elements.begin();
    for (iter = geo.begin(); iter != geo.end(); ++iter)
    {
      // get all elements locally known, including ghost elements
      pos = elements.insert(pos, *iter);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::Conditions::FindConditionObjects(const DRT::Discretization& dis,
    std::map<int, DRT::Node*>& nodes, std::map<int, DRT::Node*>& gnodes,
    std::map<int, Teuchos::RCP<DRT::Element>>& elements, const std::string& condname)
{
  std::vector<Condition*> conds;
  dis.GetCondition(condname, conds);

  FindConditionedNodes(dis, conds, nodes);

  for (const auto& cond : conds)
  {
    // get this condition's elements
    std::map<int, Teuchos::RCP<DRT::Element>>& geo = cond->Geometry();
    std::map<int, Teuchos::RCP<DRT::Element>>::iterator iter, pos;
    pos = elements.begin();
    for (iter = geo.begin(); iter != geo.end(); ++iter)
    {
      // get all elements locally known, including ghost elements
      pos = elements.insert(pos, *iter);
      const int* n = iter->second->NodeIds();
      for (int j = 0; j < iter->second->num_node(); ++j)
      {
        const int gid = n[j];
        if (dis.HaveGlobalNode(gid))
        {
          gnodes[gid] = dis.gNode(gid);
        }
        else
          FOUR_C_THROW("All nodes of known elements must be known. Panic.");
      }
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::Conditions::FindConditionObjects(const DRT::Discretization& dis,
    std::map<int, std::map<int, DRT::Node*>>& nodes,
    std::map<int, std::map<int, DRT::Node*>>& gnodes,
    std::map<int, std::map<int, Teuchos::RCP<DRT::Element>>>& elements, const std::string& condname)
{
  std::vector<Condition*> conds;
  dis.GetCondition(condname, conds);

  FindConditionedNodes(dis, conds, nodes);

  for (auto& cond : conds)
  {
    int id = cond->parameters().Get<int>("coupling id");
    // get this condition's elements
    std::map<int, Teuchos::RCP<DRT::Element>>& geo = cond->Geometry();
    std::map<int, Teuchos::RCP<DRT::Element>>::iterator iter, pos;
    pos = elements[id].begin();
    for (iter = geo.begin(); iter != geo.end(); ++iter)
    {
      // get all elements locally known, including ghost elements
      pos = elements[id].insert(pos, *iter);
      const int* n = iter->second->NodeIds();
      for (int j = 0; j < iter->second->num_node(); ++j)
      {
        const int gid = n[j];
        if (dis.HaveGlobalNode(gid))
        {
          gnodes[id][gid] = dis.gNode(gid);
        }
        else
          FOUR_C_THROW("All nodes of known elements must be known. Panic.");
      }
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::Conditions::FindConditionObjects(const DRT::Discretization& dis,
    std::map<int, Teuchos::RCP<DRT::Element>>& elements, const std::string& condname,
    const int label)
{
  std::vector<Condition*> conds;
  dis.GetCondition(condname, conds);

  bool checklabel = (label >= 0);

  for (auto& cond : conds)
  {
    if (checklabel)
    {
      const int condlabel = cond->parameters().Get<int>("label");

      if (condlabel != label) continue;  // do not consider conditions with wrong label
    }

    // get this condition's elements
    std::map<int, Teuchos::RCP<DRT::Element>>& geo = cond->Geometry();
    std::map<int, Teuchos::RCP<DRT::Element>>::iterator iter, pos;
    pos = elements.begin();
    for (iter = geo.begin(); iter != geo.end(); ++iter)
    {
      // get all elements locally known, including ghost elements
      pos = elements.insert(pos, *iter);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Conditions::FindElementConditions(
    const DRT::Element* ele, const std::string& condname, std::vector<Condition*>& condition)
{
  const DRT::Node* const* nodes = ele->Nodes();

  // We assume the conditions have unique ids. The framework has to provide
  // those.

  // the final set of conditions all nodes of this elements have in common
  std::set<Condition*> fcond;

  // we assume to always have at least one node
  // the first vector of conditions
  std::vector<Condition*> neumcond0;
  nodes[0]->GetCondition(condname, neumcond0);

  // the first set of conditions (copy vector to set)
  std::set<Condition*> cond0;
  std::copy(neumcond0.begin(), neumcond0.end(), std::inserter(cond0, cond0.begin()));


  // loop all remaining nodes
  int iel = ele->num_node();
  for (int inode = 1; inode < iel; ++inode)
  {
    std::vector<Condition*> neumcondn;
    nodes[inode]->GetCondition(condname, neumcondn);

    // the current set of conditions (copy vector to set)
    std::set<Condition*> condn;
    std::copy(neumcondn.begin(), neumcondn.end(), std::inserter(condn, condn.begin()));

    // intersect the first and the current conditions
    std::set_intersection(
        cond0.begin(), cond0.end(), condn.begin(), condn.end(), inserter(fcond, fcond.begin()));

    // make intersection to new starting condition
    cond0.clear();  // ensures that fcond is cleared in the next iteration
    std::swap(cond0, fcond);

    if (cond0.size() == 0)
    {
      // No intersections. Done. empty set is copied into condition-vector
      break;
    }
  }

  condition.clear();
  std::copy(cond0.begin(), cond0.end(), back_inserter(condition));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> CORE::Conditions::ConditionNodeRowMap(
    const DRT::Discretization& dis, const std::string& condname)
{
  return FillConditionMap(dis, dis.MyRowNodeRange(), condname);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> CORE::Conditions::ConditionNodeColMap(
    const DRT::Discretization& dis, const std::string& condname)
{
  return FillConditionMap(dis, dis.MyColNodeRange(), condname);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::set<int>> CORE::Conditions::conditioned_element_map(
    const DRT::Discretization& dis, const std::string& condname)
{
  ConditionSelector conds(dis, condname);

  Teuchos::RCP<std::set<int>> condelementmap = Teuchos::rcp(new std::set<int>());
  const int nummyelements = dis.NumMyColElements();
  for (int i = 0; i < nummyelements; ++i)
  {
    const DRT::Element* actele = dis.lColElement(i);

    const size_t numnodes = actele->num_node();
    const DRT::Node* const* nodes = actele->Nodes();
    for (size_t n = 0; n < numnodes; ++n)
    {
      const DRT::Node* actnode = nodes[n];

      // test if node is covered by condition
      if (conds.ContainsNode(actnode->Id()))
      {
        condelementmap->insert(actele->Id());
      }
    }
  }

  return condelementmap;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool CORE::Conditions::HaveSameNodes(
    const Condition* const condition1, const Condition* const condition2, const bool mustmatch)
{
  // indicates, if both conditions match
  bool matching_conditions = true;

  // get nodes of conditions
  const auto* condition1nodes = condition1->GetNodes();
  const auto* condition2nodes = condition2->GetNodes();

  // simple first check just checks the size
  if (condition1nodes->size() != condition2nodes->size())
  {
    matching_conditions = false;
    if (mustmatch)
    {
      FOUR_C_THROW(
          "Number of nodes that are defined for both conditions do not match! Did you define the "
          "conditions for the same nodesets?");
    }
  }

  // loop over all node global IDs belonging to condition1
  for (auto condition1nodegid : *condition1nodes)
  {
    bool found_node = false;
    // loop over all node global IDs belonging to condition2
    for (auto condition2nodegid : *condition2nodes)
    {
      if (condition1nodegid == condition2nodegid)
      {
        // we found the node, so set foundit to true and continue with next condition1node
        found_node = true;
        continue;
      }
    }
    // throw error if node global ID is not found in condition2
    if (!found_node)
    {
      matching_conditions = false;
      if (mustmatch)
      {
        std::cout << "Node with global ID: " << condition1nodegid
                  << "  which is part of condition: ";
        condition1->Print(std::cout);
        std::cout << " is not part of condition: ";
        condition2->Print(std::cout);
        FOUR_C_THROW(
            "Did you assign those conditions to the same nodeset? Please check your input file and "
            "fix this inconsistency!");
      }
    }
  }

  // when we get here everything is fine
  return matching_conditions;
}

FOUR_C_NAMESPACE_CLOSE
