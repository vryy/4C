/*----------------------------------------------------------------------*/
/*! \file

\brief Split conditions into map extrators

\level 1

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#ifndef DRT_CONDITION_SELECTOR_H
#define DRT_CONDITION_SELECTOR_H

#include <vector>
#include <set>
#include <string>

#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>

namespace LINALG
{
  class MultiMapExtractor;
}

namespace DRT
{
  class Discretization;
  class Condition;
  class Node;

  namespace UTILS
  {
    /// Select nodes (and their dofs) that are covered by a condition
    /*!
      We have nodal clouds. Each nodal cloud can have an arbitrary number of
      conditions. And in turn one condition (name) can be assigned to any number
      of nodal clouds. Oftentimes, however, we are not so much concerned with
      nodal clouds. We just want to known if a given node is covered by a
      specific condition. This is what this class is for.

      A stack of ConditionSelector objects is used to build a
      MultiConditionSelector which in turn is used to setup a
      LINALG::MultiMapExtractor object.

      To put it the other way: A LINALG::MultiMapExtractor splits a full (row)
      map into non-overlapping parts. If each part corresponds to a Condition, a
      MultiConditionSelector can be used to setup the LINALG::MultiMapExtractor,
      where each Condition is found by a ConditionSelector object.

      \note The Condition objects know the nodes. The maps contain dofs. The
      ConditionSelector translates between nodes and dofs.

      \see DRT::UTILS::MultiConditionSelector

      \author u.kue
      \date 08/09
     */
    class ConditionSelector
    {
     public:
      /// Construct a selector on the given Discretization for the Condition
      /// with the given name
      ConditionSelector(const DRT::Discretization& dis, std::string condname);

      /// construct a selector from a given vector of conditions
      ConditionSelector(const DRT::Discretization& dis,  //!< discretization
          const std::vector<DRT::Condition*>& conds      //!< given vector of conditions
      );

      /// Destructor
      virtual ~ConditionSelector() {}

      /// select all matching dofs of the node and put them into conddofset
      virtual bool SelectDofs(DRT::Node* node, std::set<int>& conddofset);

      /// tell if the node gid is known by any condition of the given name
      virtual bool ContainsNode(int ngid);

      /// tell if the dof of a node from this condition is covered as well
      virtual bool ContainsDof(int dof, int pos) { return true; }

     protected:
      /// Discretization we are looking at
      const DRT::Discretization& Discretization() const { return dis_; }

      /// all conditions that come by the given name
      const std::vector<DRT::Condition*>& Conditions() const { return conds_; }

     private:
      /// Discretization
      const DRT::Discretization& dis_;

      /// Conditions
      std::vector<DRT::Condition*> conds_;
    };


    /// Select some dofs of the conditioned node
    /*!
      In addition to the assignment of nodes we might want to distinguish
      between different dofs. Sometimes a condition applies only to some dofs of
      each node. This is what the ContainsDof() method is for.

      This selector is used most often. It can be applied e.g. to extract the
      velocity dofs from a fluid node (with velocity and pressure dofs)

      \author u.kue
      \date 08/09
     */
    class NDimConditionSelector : public ConditionSelector
    {
     public:
      NDimConditionSelector(
          const DRT::Discretization& dis, std::string condname, int startdim, int enddim)
          : ConditionSelector(dis, condname), startdim_(startdim), enddim_(enddim)
      {
      }

      /// Contain a dof number only if the dof nodal position is within the
      /// allowed range.
      virtual bool ContainsDof(int dof, int pos) { return startdim_ <= pos and pos < enddim_; }

     private:
      int startdim_;
      int enddim_;
    };


    /// select all dirichlet dofs
    class DirichletSelector : public ConditionSelector
    {
     public:
      DirichletSelector(const DRT::Discretization& dis) : ConditionSelector(dis, "Dirichlet") {}

      /// select all matching dofs of the node and put them into conddofset
      /*!
        Dirichlet condition are special since lower rank conditions override
        higher rank conditions. Thus if a node is covered by e.g. a line
        condition, the node is considered to be selected and not processed by
        surface or volume conditions, even if the line condition selects no dofs
        at all. This way nodes can be freed explicitly by Dirichlet conditions.
       */
      virtual bool SelectDofs(DRT::Node* node, std::set<int>& conddofset);
    };


    /// a collection of ConditionSelector objects used to create a LINALG::MultiMapExtractor
    /*!
      Oftentimes the dofs of a field need to be split into a set of disjoint
      maps. The LINALG::MultiMapExtractor class takes care of these splits. However,
      LINALG::MultiMapExtractor does not do the splitting itself. This is where
      MultiConditionSelector comes in. Here, different ConditionSelector objects
      are collected and afterwards the LINALG::MultiMapExtractor is setup according to the
      conditions.

      MultiConditionSelector is a helper class that is needed during the setup
      only.

      \note Each node ends up in one condition only, independent of whether all
      its dofs are selected or not.

      \author u.kue
      \date 08/09
     */
    class MultiConditionSelector
    {
     public:
      MultiConditionSelector();

      /// add a new ConditionSelector
      /*!
        \note The order of the selector additions determines the slots within
        LINALG::MultiMapExtractor.
       */
      void AddSelector(Teuchos::RCP<ConditionSelector> s) { selectors_.push_back(s); }

      /// Do the setup
      void SetupExtractor(const DRT::Discretization& dis, const Epetra_Map& fullmap,
          LINALG::MultiMapExtractor& extractor);

      /// Activate overlapping
      void SetOverlapping(bool overlapping) { overlapping_ = overlapping; }

     private:
      /// evaluate the ConditionSelector objects
      void SetupCondDofSets(const DRT::Discretization& dis);

      /// condition selectors
      std::vector<Teuchos::RCP<ConditionSelector>> selectors_;

      /// sets of selected dof numbers
      std::vector<std::set<int>> conddofset_;

      /// flag defines if maps are overlapping
      bool overlapping_;
    };

  }  // namespace UTILS
}  // namespace DRT

#endif
