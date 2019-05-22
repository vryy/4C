/*----------------------------------------------------------------------*/
/*!

\brief Split conditions into map extractors

\level 1

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/


#include "drt_condition_selector.H"
#include "drt_dserror.H"
#include "drt_condition.H"
#include "drt_discret.H"

#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::ConditionSelector::ConditionSelector(
    const DRT::Discretization& dis, std::string condname)
    : dis_(dis)
{
  dis.GetCondition(condname, conds_);
  std::sort(conds_.begin(), conds_.end(), DRT::ConditionLess());
}


/*----------------------------------------------------------------------*
 | construct a selector from a given vector of conditions    fang 07/16 |
 *----------------------------------------------------------------------*/
DRT::UTILS::ConditionSelector::ConditionSelector(
    const DRT::Discretization& dis,            //!< discretization
    const std::vector<DRT::Condition*>& conds  //!< given vector of conditions
    )
    : dis_(dis), conds_(conds)
{
  std::sort(conds_.begin(), conds_.end(), DRT::ConditionLess());

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::UTILS::ConditionSelector::SelectDofs(DRT::Node* node, std::set<int>& conddofset)
{
  bool found = false;

  // put all conditioned dofs into conddofset
  if (ContainsNode(node->Id()))
  {
    std::vector<int> dof = Discretization().Dof(0, node);
    for (unsigned k = 0; k < dof.size(); ++k)
    {
      // test for dof position
      if (ContainsDof(dof[k], k))
      {
        conddofset.insert(dof[k]);
        found = true;
      }
    }
  }
  return found;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::UTILS::ConditionSelector::ContainsNode(int ngid)
{
  for (unsigned j = 0; j < conds_.size(); ++j)
  {
    if (conds_[j]->ContainsNode(ngid))
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::UTILS::DirichletSelector::SelectDofs(DRT::Node* node, std::set<int>& conddofset)
{
  bool found = false;
  int ngid = node->Id();

  // The condition vector is sorted by condition type. Thus lesser entity
  // ranks are considered first. The first condition that covers a node gets
  // it.

  const std::vector<DRT::Condition*>& conds = Conditions();
  for (std::vector<DRT::Condition*>::const_iterator i = conds.begin(); i != conds.end(); ++i)
  {
    DRT::Condition& c = **i;
    const std::vector<int>* onoff = c.Get<std::vector<int>>("onoff");
    if (onoff == NULL) dserror("not a valid Dirichlet condition");
    if (c.ContainsNode(ngid))
    {
      std::vector<int> dof = Discretization().Dof(node);
      for (unsigned k = 0; k < dof.size(); ++k)
      {
        if (k > onoff->size()) dserror("not a valid Dirichlet condition");
        if ((*onoff)[k] != 0)
        {
          conddofset.insert(dof[k]);
        }
      }

      // if a node has been covered by one Dirichlet condition do not look for
      // further conditions
      found = true;
      break;
    }
  }
  return found;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::MultiConditionSelector::MultiConditionSelector() : overlapping_(false) {}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::MultiConditionSelector::SetupExtractor(
    const DRT::Discretization& dis, const Epetra_Map& fullmap, LINALG::MultiMapExtractor& extractor)
{
  SetupCondDofSets(dis);

  // Find all non-conditioned dofs by subtracting all conditioned ones.

  std::set<int> otherdofset(
      fullmap.MyGlobalElements(), fullmap.MyGlobalElements() + fullmap.NumMyElements());

  for (unsigned j = 0; j < conddofset_.size(); ++j)
  {
    std::set<int>& conddofset = conddofset_[j];

    std::set<int>::const_iterator conditer;
    for (conditer = conddofset.begin(); conditer != conddofset.end(); ++conditer)
    {
      otherdofset.erase(*conditer);
    }
  }

  // Setup all maps. The "other" map goes first so it becomes the zeroth map
  // of the MultiMapExtractor.

  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  maps.reserve(conddofset_.size() + 1);

  maps.push_back(LINALG::CreateMap(otherdofset, dis.Comm()));
  for (unsigned j = 0; j < conddofset_.size(); ++j)
  {
    std::set<int>& conddofset = conddofset_[j];
    maps.push_back(LINALG::CreateMap(conddofset, dis.Comm()));
  }

  // MultiMapExtractor setup

  extractor.Setup(fullmap, maps);
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::MultiConditionSelector::SetupCondDofSets(const DRT::Discretization& dis)
{
  // we get as many sets as we have selectors
  conddofset_.resize(selectors_.size());

  // for each owned node
  int numrownodes = dis.NumMyRowNodes();
  for (int i = 0; i < numrownodes; ++i)
  {
    DRT::Node* node = dis.lRowNode(i);

    // test each selector
    for (unsigned j = 0; j < selectors_.size(); ++j)
    {
      ConditionSelector& conds = *selectors_[j];

      // if the selector applies, we are done
      if (conds.SelectDofs(node, conddofset_[j]))
        if (!overlapping_) break;
    }
  }
}
