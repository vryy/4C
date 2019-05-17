/*----------------------------------------------------------------------*/
/*!

\brief Implementation of utils on conditions

\level 1

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#include <map>
#include <set>
#include <string>
#include <vector>
#include <algorithm>

#include "drt_condition_utils.H"
#include "drt_condition_selector.H"
#include "drt_discret_iterator.H"
#include "drt_globalproblem.H"

#include "../drt_lib/drt_utils_parallel.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io_control.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::FindConditionedNodes(
    const DRT::Discretization& dis, std::string condname, std::vector<int>& nodes)
{
  std::vector<DRT::Condition*> conds;
  dis.GetCondition(condname, conds);
  FindConditionedNodes(dis, conds, nodes);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::FindConditionedNodes(
    const DRT::Discretization& dis, std::string condname, std::set<int>& nodeset)
{
  std::vector<DRT::Condition*> conds;
  dis.GetCondition(condname, conds);
  FindConditionedNodes(dis, conds, nodeset);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::FindConditionedNodes(
    const DRT::Discretization& dis, std::string condname, std::map<int, DRT::Node*>& nodes)
{
  std::vector<DRT::Condition*> conds;
  dis.GetCondition(condname, conds);
  FindConditionedNodes(dis, conds, nodes);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::FindConditionedNodes(const DRT::Discretization& dis,
    const std::vector<DRT::Condition*>& conds, std::vector<int>& nodes)
{
  std::set<int> nodeset;
  const int myrank = dis.Comm().MyPID();
  for (unsigned i = 0; i < conds.size(); ++i)
  {
    const std::vector<int>* n = conds[i]->Nodes();
    for (unsigned j = 0; j < n->size(); ++j)
    {
      const int gid = (*n)[j];
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
void DRT::UTILS::FindConditionedNodes(const DRT::Discretization& dis,
    const std::vector<DRT::Condition*>& conds, std::set<int>& nodeset)
{
  const int myrank = dis.Comm().MyPID();
  for (unsigned i = 0; i < conds.size(); ++i)
  {
    const std::vector<int>* n = conds[i]->Nodes();
    for (unsigned j = 0; j < n->size(); ++j)
    {
      const int gid = (*n)[j];
      if (dis.HaveGlobalNode(gid) and dis.gNode(gid)->Owner() == myrank)
      {
        nodeset.insert(gid);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::FindConditionedNodes(const DRT::Discretization& dis,
    const std::vector<DRT::Condition*>& conds, std::map<int, DRT::Node*>& nodes)
{
  const int myrank = dis.Comm().MyPID();
  for (unsigned i = 0; i < conds.size(); ++i)
  {
    const std::vector<int>* n = conds[i]->Nodes();
    for (unsigned j = 0; j < n->size(); ++j)
    {
      const int gid = (*n)[j];
      if (dis.HaveGlobalNode(gid) and dis.gNode(gid)->Owner() == myrank)
      {
        nodes[gid] = dis.gNode(gid);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::FindConditionedNodes(const DRT::Discretization& dis,
    const std::vector<DRT::Condition*>& conds, std::map<int, Teuchos::RCP<std::vector<int>>>& nodes)
{
  std::map<int, std::set<int>> nodeset;
  const int myrank = dis.Comm().MyPID();
  for (unsigned i = 0; i < conds.size(); ++i)
  {
    int id = conds[i]->GetInt("coupling id");
    const std::vector<int>* n = conds[i]->Nodes();
    for (unsigned j = 0; j < n->size(); ++j)
    {
      const int gid = (*n)[j];
      if (dis.HaveGlobalNode(gid) and dis.gNode(gid)->Owner() == myrank)
      {
        nodeset[id].insert(gid);
      }
    }
  }

  std::map<int, std::set<int>>::iterator iter;
  for (iter = nodeset.begin(); iter != nodeset.end(); ++iter)
  {
    nodes[iter->first] = Teuchos::rcp(new std::vector<int>((iter->second).size()));
    nodes[iter->first]->assign((iter->second).begin(), (iter->second).end());
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::FindConditionedNodes(const DRT::Discretization& dis,
    const std::vector<DRT::Condition*>& conds, std::map<int, std::map<int, DRT::Node*>>& nodes)
{
  const int myrank = dis.Comm().MyPID();
  for (unsigned i = 0; i < conds.size(); ++i)
  {
    int id = conds[i]->GetInt("coupling id");
    const std::vector<int>* n = conds[i]->Nodes();
    for (unsigned j = 0; j < n->size(); ++j)
    {
      const int gid = (*n)[j];
      if (dis.HaveGlobalNode(gid) and dis.gNode(gid)->Owner() == myrank)
      {
        (nodes[id])[gid] = dis.gNode(gid);
      }
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::FindConditionObjects(const DRT::Discretization& dis,
    std::map<int, DRT::Node*>& nodes, std::map<int, Teuchos::RCP<DRT::Element>>& elements,
    const std::string& condname)
{
  int myrank = dis.Comm().MyPID();
  std::vector<DRT::Condition*> conds;
  dis.GetCondition(condname, conds);

  FindConditionedNodes(dis, conds, nodes);

  for (unsigned i = 0; i < conds.size(); ++i)
  {
    // get this condition's elements
    std::map<int, Teuchos::RCP<DRT::Element>>& geo = conds[i]->Geometry();
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
void DRT::UTILS::FindConditionObjects(const DRT::Discretization& dis,
    std::map<int, DRT::Node*>& nodes, std::map<int, DRT::Node*>& gnodes,
    std::map<int, Teuchos::RCP<DRT::Element>>& elements, const std::vector<DRT::Condition*>& conds)
{
  FindConditionedNodes(dis, conds, nodes);

  for (size_t i = 0; i < conds.size(); ++i)
  {
    // get this condition's elements
    std::map<int, Teuchos::RCP<DRT::Element>>& geo = conds[i]->Geometry();
    std::map<int, Teuchos::RCP<DRT::Element>>::iterator iter, pos;
    pos = elements.begin();
    for (iter = geo.begin(); iter != geo.end(); ++iter)
    {
      // get all elements locally known, including ghost elements
      pos = elements.insert(pos, *iter);
      const int* n = iter->second->NodeIds();
      for (int j = 0; j < iter->second->NumNode(); ++j)
      {
        const int gid = n[j];
        if (dis.HaveGlobalNode(gid))
        {
          gnodes[gid] = dis.gNode(gid);
        }
        else
          dserror("All nodes of known elements must be known. Panic.");
      }
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::FindConditionObjects(
    std::map<int, Teuchos::RCP<DRT::Element>>& elements, const std::vector<DRT::Condition*>& conds)
{
  for (size_t i = 0; i < conds.size(); ++i)
  {
    // get this condition's elements
    std::map<int, Teuchos::RCP<DRT::Element>>& geo = conds[i]->Geometry();
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
void DRT::UTILS::FindConditionObjects(const DRT::Discretization& dis,
    std::map<int, DRT::Node*>& nodes, std::map<int, DRT::Node*>& gnodes,
    std::map<int, Teuchos::RCP<DRT::Element>>& elements, const std::string& condname)
{
  std::vector<DRT::Condition*> conds;
  dis.GetCondition(condname, conds);

  FindConditionedNodes(dis, conds, nodes);

  for (unsigned i = 0; i < conds.size(); ++i)
  {
    // get this condition's elements
    std::map<int, Teuchos::RCP<DRT::Element>>& geo = conds[i]->Geometry();
    std::map<int, Teuchos::RCP<DRT::Element>>::iterator iter, pos;
    pos = elements.begin();
    for (iter = geo.begin(); iter != geo.end(); ++iter)
    {
      // get all elements locally known, including ghost elements
      pos = elements.insert(pos, *iter);
      const int* n = iter->second->NodeIds();
      for (int j = 0; j < iter->second->NumNode(); ++j)
      {
        const int gid = n[j];
        if (dis.HaveGlobalNode(gid))
        {
          gnodes[gid] = dis.gNode(gid);
        }
        else
          dserror("All nodes of known elements must be known. Panic.");
      }
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::FindConditionObjects(const DRT::Discretization& dis,
    std::map<int, std::map<int, DRT::Node*>>& nodes,
    std::map<int, std::map<int, DRT::Node*>>& gnodes,
    std::map<int, std::map<int, Teuchos::RCP<DRT::Element>>>& elements, const std::string& condname)
{
  std::vector<DRT::Condition*> conds;
  dis.GetCondition(condname, conds);

  FindConditionedNodes(dis, conds, nodes);

  for (unsigned i = 0; i < conds.size(); ++i)
  {
    int id = conds[i]->GetInt("coupling id");
    // get this condition's elements
    std::map<int, Teuchos::RCP<DRT::Element>>& geo = conds[i]->Geometry();
    std::map<int, Teuchos::RCP<DRT::Element>>::iterator iter, pos;
    pos = elements[id].begin();
    for (iter = geo.begin(); iter != geo.end(); ++iter)
    {
      // get all elements locally known, including ghost elements
      pos = elements[id].insert(pos, *iter);
      const int* n = iter->second->NodeIds();
      for (int j = 0; j < iter->second->NumNode(); ++j)
      {
        const int gid = n[j];
        if (dis.HaveGlobalNode(gid))
        {
          gnodes[id][gid] = dis.gNode(gid);
        }
        else
          dserror("All nodes of known elements must be known. Panic.");
      }
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::FindConditionObjects(const DRT::Discretization& dis,
    std::map<int, Teuchos::RCP<DRT::Element>>& elements, const std::string& condname,
    const int label)
{
  std::vector<DRT::Condition*> conds;
  dis.GetCondition(condname, conds);

  bool checklabel = (label >= 0);

  for (unsigned i = 0; i < conds.size(); ++i)
  {
    if (checklabel)
    {
      const int condlabel = conds[i]->GetInt("label");

      if (condlabel != label) continue;  // do not consider conditions with wrong label
    }

    // get this condition's elements
    std::map<int, Teuchos::RCP<DRT::Element>>& geo = conds[i]->Geometry();
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
Teuchos::RCP<Epetra_Map> DRT::UTILS::ConditionElementMap(
    const DRT::Discretization& dis, std::string condname, bool colmap)
{
  std::vector<DRT::Condition*> conds;
  dis.GetCondition(condname, conds);
  std::set<int> elementset;

  if (colmap)
  {
    for (unsigned i = 0; i < conds.size(); ++i)
    {
      std::map<int, Teuchos::RCP<DRT::Element>>& geometry = conds[i]->Geometry();
      std::transform(geometry.begin(), geometry.end(),
          std::inserter(elementset, elementset.begin()),
          LINALG::select1st<std::pair<int, Teuchos::RCP<DRT::Element>>>());
    }
  }
  else
  {
    int myrank = dis.Comm().MyPID();
    for (unsigned i = 0; i < conds.size(); ++i)
    {
      std::map<int, Teuchos::RCP<DRT::Element>>& geometry = conds[i]->Geometry();
      for (std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator iter = geometry.begin();
           iter != geometry.end(); ++iter)
      {
        if (iter->second->Owner() == myrank)
        {
          elementset.insert(iter->first);
        }
      }
    }
  }

  return LINALG::CreateMap(elementset, dis.Comm());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::FindElementConditions(
    const DRT::Element* ele, const std::string& condname, std::vector<DRT::Condition*>& condition)
{
  const DRT::Node* const* nodes = ele->Nodes();

  // We assume the conditions have unique ids. The framework has to provide
  // those.

  // the final set of conditions all nodes of this elements have in common
  std::set<DRT::Condition*> fcond;

  // we assume to always have at least one node
  // the first vector of conditions
  std::vector<DRT::Condition*> neumcond0;
  nodes[0]->GetCondition(condname, neumcond0);

  // the first set of conditions (copy vector to set)
  std::set<DRT::Condition*> cond0;
  std::copy(neumcond0.begin(), neumcond0.end(), std::inserter(cond0, cond0.begin()));


  // loop all remaining nodes
  int iel = ele->NumNode();
  for (int inode = 1; inode < iel; ++inode)
  {
    std::vector<DRT::Condition*> neumcondn;
    nodes[inode]->GetCondition(condname, neumcondn);

    // the current set of conditions (copy vector to set)
    std::set<DRT::Condition*> condn;
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
Teuchos::RCP<Epetra_Map> DRT::UTILS::ConditionNodeRowMap(
    const DRT::Discretization& dis, const std::string& condname)
{
  RowNodeIterator iter(dis);
  return ConditionMap(dis, iter, condname);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> DRT::UTILS::ConditionNodeColMap(
    const DRT::Discretization& dis, const std::string& condname)
{
  ColNodeIterator iter(dis);
  return ConditionMap(dis, iter, condname);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> DRT::UTILS::ConditionMap(const DRT::Discretization& dis,
    const DiscretizationNodeIterator& iter, const std::string& condname)
{
  std::set<int> condnodeset;

  ConditionSelector conds(dis, condname);

  const int numnodes = iter.NumEntries();
  for (int i = 0; i < numnodes; ++i)
  {
    const DRT::Node* node = iter.Entry(i);
    if (conds.ContainsNode(node->Id()))
    {
      condnodeset.insert(node->Id());
    }
  }

  Teuchos::RCP<Epetra_Map> condnodemap = LINALG::CreateMap(condnodeset, dis.Comm());
  return condnodemap;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::set<int>> DRT::UTILS::ConditionedElementMap(
    const DRT::Discretization& dis, const std::string& condname)
{
  ConditionSelector conds(dis, condname);

  Teuchos::RCP<std::set<int>> condelementmap = Teuchos::rcp(new std::set<int>());
  const int nummyelements = dis.NumMyColElements();
  for (int i = 0; i < nummyelements; ++i)
  {
    const DRT::Element* actele = dis.lColElement(i);

    const size_t numnodes = actele->NumNode();
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


/*----------------------------------------------------------------------
 * collects elements by labels (have to be implemented in the           *
 * corresponding condition)                                             *
 *----------------------------------------------------------------------*/
void DRT::UTILS::CollectElementsByConditionLabel(const DRT::Discretization& discret,
    std::map<int, std::set<int>>& elementsByLabel, const std::string& name)
{
  // Reset
  elementsByLabel.clear();
  // get condition
  std::vector<DRT::Condition*> conditions;
  discret.GetCondition(name, conditions);

  // collect elements by xfem coupling label
  for (std::vector<DRT::Condition*>::const_iterator conditer = conditions.begin();
       conditer != conditions.end(); ++conditer)
  {
    const DRT::Condition* condition = *conditer;
    const int label = condition->GetInt("label");
    for (int iele = 0; iele < discret.NumMyColElements(); ++iele)
    {
      // for each element, check, whether all nodes belong to same condition label
      const DRT::Element* element = discret.lColElement(iele);
      int nodecounter = 0;
      for (int inode = 0; inode < element->NumNode(); ++inode)
      {
        const DRT::Node* node = element->Nodes()[inode];
        if (condition->ContainsNode(node->Id())) nodecounter++;
      }
      // if all nodes belong to label, then this element gets a label entry
      if (nodecounter == element->NumNode()) elementsByLabel[label].insert(element->Id());
    }
  }
  int numOfCollectedIds = 0;
  for (std::map<int, std::set<int>>::const_iterator entry = elementsByLabel.begin();
       entry != elementsByLabel.end(); ++entry)
  {
    numOfCollectedIds += entry->second.size();
  }

  if (discret.NumMyColElements() != numOfCollectedIds) dserror("not all elements collected.");
}

/*-----------------------------------------------------------------------*
 * Writes boundary surfaces of a volumetrically coupled problem to file  *
 * 'boundarysurfaces.log' storing the condition-Id as surface-Id. For    *
 * visualisation in gmsh and checking for tetrahedra whose four surfaces *
 * are wrongly contained in the boundary surface of the volumetric       *
 * coupling this file can be used.                         (croth 01/15) *
 *-----------------------------------------------------------------------*/
void DRT::UTILS::WriteBoundarySurfacesVolumeCoupling(
    std::map<std::vector<int>, Teuchos::RCP<DRT::Element>> surfmap, int condID, int numproc,
    int mypid)
{
  if (numproc == 1)
  {
    // Get output prefix
    std::string outputprefix = DRT::Problem::Instance()->OutputControlFile()->NewOutputFileName();
    // Create boundary surface file
    std::ostringstream sf;
    sf << outputprefix << "_boundarysurfaces.log";
    std::string boundarysurffilename;
    boundarysurffilename = sf.str();

    std::ofstream myfile;
    myfile.open(boundarysurffilename.c_str(), std::ios_base::app);
    myfile << "Surfaces in Surfmap for Coupling Condition No. " << condID + 1 << ":\n";
    myfile << "Format: [Node1, Node2, Node3, CondID] \n";
    for (std::map<std::vector<int>, Teuchos::RCP<DRT::Element>>::const_iterator iterat =
             surfmap.begin();
         iterat != surfmap.end(); ++iterat)
    {
      myfile << iterat->first[0] << " " << iterat->first[1] << " " << iterat->first[2] << " "
             << condID + 1 << "\n";
    }
    myfile << "End \n";
    myfile.close();
    std::cout << " Condition " << condID + 1 << " checked and written to file." << std::endl;
  }
  else if (mypid == 0)
  {
    std::cout << " No 'boundarysurfaces.log' written as number of procs = " << numproc
              << " is bigger than 1." << std::endl;
  }
}
