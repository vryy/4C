#ifndef DRT_CONDITION_SELECTOR_H
#define DRT_CONDITION_SELECTOR_H

#include <vector>
#include <string>

#include "drt_condition.H"
#include "drt_discret.H"

namespace DRT
{
namespace UTILS
{

  /// help to find a node inside a condition
  /*!
    We have nodal clouds. Each nodal cloud can have an arbitrary number of
    conditions. And in turn one condition (name) can be assigned to any number
    of nodal clouds. Oftentimes, however, we are not so much concerned with
    nodal clouds. We just want to known if a given node is covered by a
    specific condition. This is what this class is for.

    \author u.kue
    \date 08/09
   */
  class ConditionSelector
  {
  public:

    ConditionSelector(const DRT::Discretization& dis, std::string condname) : dis_(dis), condname_(condname)
    {
      dis.GetCondition(condname, conds_);
    }

    virtual ~ConditionSelector() {}

    /// tell if the node gid is known by any condition of the given name
    bool ContainsNode(int ngid) const
    {
      for (unsigned j=0; j<conds_.size(); ++j)
      {
        if (conds_[j]->ContainsNode(ngid))
        {
          return true;
        }
      }
      return false;
    }

    /// tell if the dof of a node from this condition is covered as well
    virtual bool ContainsDof(int dof, int pos) const { return true; }

  private:
    const DRT::Discretization& dis_;
    std::string condname_;
    std::vector<DRT::Condition*> conds_;
  };


  /// help to find a node inside a condition
  /*!
    In addition to the assignment of nodes we might want to distinguish
    between different dofs. Sometimes a condition applies only to some dofs of
    each node. This is what the ContainsDof() method is for.

    \author u.kue
    \date 08/09
   */
  class NDimConditionSelector : public ConditionSelector
  {
  public:

    NDimConditionSelector(const DRT::Discretization& dis,
                          std::string condname,
                          int startdim,
                          int enddim)
      : ConditionSelector(dis, condname), startdim_(startdim), enddim_(enddim) {}

    virtual bool ContainsDof(int dof, int pos) const { return startdim_<=pos and pos<enddim_; }

  private:
    int startdim_;
    int enddim_;
  };


  /// a collection of ConditionSelector objects used to create a MultiMapExtractor
  /*!
    Oftentimes the dofs of a field need to be split into a set of disjoint
    maps. The MultiMapExtractor class takes care of these splits. However,
    MultiMapExtractor does not do the splitting itself. This is where
    MultiConditionSelector comes in. Here, different ConditionSelector objects
    are collected and afterwards the MultiMapExtractor setup according to the
    conditions.

    MultiConditionSelector is a helper class that is needed during the setup
    only.

    \author u.kue
    \date 08/09
   */
  class MultiConditionSelector
  {
  public:

    /// add a new ConditionSelector
    void AddSelector(RCP<ConditionSelector> s) { selectors_.push_back(s); }

    /// Do the setup
    void SetupExtractor(const DRT::Discretization& dis,
                        Teuchos::RCP<Epetra_Map> fullmap,
                        LINALG::MultiMapExtractor& extractor);

  private:

    /// evaluate the ConditionSelector objects
    void SetupCondDofSets(const DRT::Discretization& dis);

    std::vector<RCP<ConditionSelector> > selectors_;
    std::vector<std::set<int> > conddofset_;
  };

}
}

#endif
