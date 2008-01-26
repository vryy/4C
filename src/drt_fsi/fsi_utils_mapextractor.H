#ifdef CCADISCRET
#ifndef FSI_UTILS_MAPEXTRACTOR_H
#define FSI_UTILS_MAPEXTRACTOR_H

#include "../drt_lib/drt_utils_mapextractor.H"


namespace FSI
{
  namespace UTILS
  {

    /// MapExtractor subclass that remembers the selected nodes
    /*!

    In the FSI context we need to know the nodes at the interface to match
    interfaces of different fields. And since the MapExtractor knowns these
    nodes already, we want to take advantage of this knowledge.

    All we do here is to provide a special Setup() method that remembers the
    nodes (in addition to the dofs) if the supplied ExtractorCond matches. To
    do so a special condition NodeSaveExtractorCond is applied, however this
    is a mere implementation detail.

    \author u.kue
    \date 01/08
    */
    class MapExtractor : public DRT::UTILS::MapExtractor
    {
    public:

      MapExtractor(const Teuchos::RCP<DRT::Discretization> dis) : dis_(dis) {}

      /** \name Setup */
      //@{

      /// setup with implicit storage of nodal GIDs
      template<class C>
      void Setup(const C& condition);

      //@}

      /// the discretization we are working on
      const DRT::Discretization& Discret() const { return *dis_; }

      /// the nodes covered by this condition
      /*!
        This is setup data. We make it available here because others might want
        to use it during their setup.
      */
      const std::vector<int>& Nodes() const { return nodes_; }

    private:

      /// our discretization
      const Teuchos::RCP<DRT::Discretization> dis_;

      /// nodes selected by our condition
      std::vector<int> nodes_;
    };


    /// Internal ExtractorCond class that wraps a real ExtractorCond and remembers the nodal GID on each match.
    template<class C>
    class NodeSaveExtractorCond
    {
    public:

      NodeSaveExtractorCond(const C& c) : c_(c) {}

      /// required function call operator
      bool operator()(const Teuchos::RCP<DRT::Discretization> dis,
                      const DRT::Node* node,
                      const std::vector<int>& dof,
                      int j) const
      {
        if (c_(dis,node,dof,j))
        {
          nodes_.insert(node->Id());
          return true;
        }
        return false;
      }

      const std::set<int>& Nodes() { return nodes_; }

    private:
      C c_;
      mutable std::set<int> nodes_;
    };
  }
}


template<class C>
void FSI::UTILS::MapExtractor::Setup(const C& condition)
{
  NodeSaveExtractorCond<C> save(condition);
  DRT::UTILS::MapExtractor::Setup(dis_,save);

  const std::set<int>& n = save.Nodes();
  nodes_.reserve(n.size());
  nodes_.assign(n.begin(),n.end());
}

#endif
#endif
