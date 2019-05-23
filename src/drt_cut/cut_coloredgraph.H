/*---------------------------------------------------------------------*/
/*!
\file cut_coloredgraph.H

\brief colored graph to create volumecells from facets

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249

*----------------------------------------------------------------------*/


#ifndef CUT_COLOREDGRAPH_H
#define CUT_COLOREDGRAPH_H

#include <map>
#include <list>

#include "cut_utils.H"

namespace GEO
{
  namespace CUT
  {
    namespace COLOREDGRAPH
    {
      class Graph;
      class CycleList;
      class CycleListIterator;

      class ForkFinder
      {
       public:
        ForkFinder(Graph& graph, Graph& used, Graph& cycle, const plain_int_set& free)
            : graph_(graph), used_(used), cycle_(cycle), free_(free)
        {
        }

        bool operator()(const std::pair<const int, plain_int_set>& point);

        Graph& graph_;
        Graph& used_;
        Graph& cycle_;
        const plain_int_set& free_;
      };

      class SingeLineFinder
      {
       public:
        SingeLineFinder(int color_split) : color_split_(color_split) {}

        bool operator()(const std::pair<const int, plain_int_set>& i)
        {
          int p = i.first;
          if (p < color_split_)
          {
            return i.second.size() < 3;
          }
          else
          {
            return i.second.size() < 2;
          }
        }

       private:
        int color_split_;
      };

      /// facet-line graph
      class Graph
      {
       public:
        typedef std::map<int, plain_int_set>::const_iterator const_iterator;

        explicit Graph(int color_split) : color_split_(color_split) {}

        void Add(int row, int col);

        void Add(int p, const plain_int_set& row);

        int FindNext(Graph& used, int p, Graph& cycle, const plain_int_set& free);

        void FindFreeFacets(Graph& graph, Graph& used, plain_int_set& free,
            const std::vector<std::pair<Point*, Point*>>& all_lines, std::vector<int>& split_trace);

        void FindSplitTrace(std::vector<int>& split_trace);

        void GetAll(plain_int_set& all);

        void FixSingleLines();

        void TestClosed();

        void TestSplit();

        void TestFacets();

        void Print() const;

        plain_int_set& at(int p) { return graph_[p]; }

        plain_int_set& operator[](int p) { return graph_[p]; }

        unsigned count(int p) { return graph_.count(p); }

        unsigned size() const { return graph_.size(); }

        const_iterator begin() const { return graph_.begin(); }

        const_iterator end() const { return graph_.end(); }

        void Swap(Graph& other)
        {
          std::swap(graph_, other.graph_);
          std::swap(color_split_, other.color_split_);
        }

        int Split() { return color_split_; }

        void SetSplit(int color_split) { color_split_ = color_split; }

        void Split(Graph& used, plain_int_set& free, Graph& connection,
            const std::vector<int>& split_trace, Graph& c1, Graph& c2, Graph& datagraph);

        bool ContainsTrace(const std::vector<int>& split_trace);

        bool operator==(const Graph& other) const
        {
          return color_split_ == other.color_split_ and graph_ == other.graph_;
        }

        bool operator!=(const Graph& other) const { return not(*this == other); }

        void DumpGraph(const std::string& name);

        void Map(std::map<int, const void*>* index_value_map)
        {
          index_value_map_ = index_value_map;
        };

        // get pointer to the underlying object
        void* GetPointer(int Id) const
        {
          std::map<int, const void*>::iterator it = index_value_map_->find(Id);
          if (it != index_value_map_->end())
          {
            return const_cast<void*>(it->second);
          }
          else
            dserror("Invalid element with id %d", Id);
          return NULL;
        }
        // Split split trace into the closed loops
        void SplitSplittrace(const std::vector<int>& split_trace, Graph& datagraph,
            std::vector<std::vector<int>>& isolated_components);

       private:
        void Fill(const std::vector<int>& split_trace, Graph& used, plain_int_set& free,
            Graph& connection, int seed, Graph& c);

        std::map<int, plain_int_set> graph_;

        /// position where lines start
        int color_split_;

        // mapping between Ids and pointers to the underlying objects
        std::map<int, const void*>* index_value_map_;
      };

      /// One graph that represents the outside facets of one or more volumes
      class Cycle
      {
        friend class CycleListIterator;

       public:
        explicit Cycle(int color_split) : cycle_(color_split) {}

        void Assign(Graph& cycle) { cycle_.Swap(cycle); }

        void Print() const;

        void Split(Graph& used, plain_int_set& free, Graph& connection,
            const std::vector<int>& split_trace, Graph& c1, Graph& c2, Graph& datagraph)
        {
          cycle_.Split(used, free, connection, split_trace, c1, c2, datagraph);
        }

        bool ContainsTrace(const std::vector<int>& split_trace)
        {
          return cycle_.ContainsTrace(split_trace);
        }

        Graph& operator()() { return cycle_; }

       private:
        Graph cycle_;
      };

      class CycleListIterator
      {
       public:
        CycleListIterator(std::list<Cycle>& cycles, std::list<Cycle>::iterator i)
            : cycles_(cycles), i_(i)
        {
          NextActive();
        }

        void NextActive() {}

        CycleListIterator& operator++()
        {
          ++i_;
          NextActive();
          return *this;
        }

        Graph& operator*() { return i_->cycle_; }

        bool operator!=(const CycleListIterator& other) { return i_ != other.i_; }

       private:
        std::list<Cycle>& cycles_;
        std::list<Cycle>::iterator i_;
      };

      /// list of cycles
      class CycleList
      {
       public:
        typedef CycleListIterator iterator;

        void AddPoints(Graph& graph, Graph& used, Graph& cycle, plain_int_set& free,
            const std::vector<std::pair<Point*, Point*>>& all_lines);

        void Print() const;

        unsigned size() const { return cycles_.size(); }

        CycleListIterator begin() { return CycleListIterator(cycles_, cycles_.begin()); }

        CycleListIterator end() { return CycleListIterator(cycles_, cycles_.end()); }

       private:
        void PushBack(Graph& g);

        std::list<Cycle> cycles_;
      };

    }  // namespace COLOREDGRAPH
  }    // namespace CUT
}  // namespace GEO

#endif
