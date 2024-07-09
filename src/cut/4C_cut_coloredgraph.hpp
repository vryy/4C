/*---------------------------------------------------------------------*/
/*! \file

\brief colored graph to create volumecells from facets

\level 3


*----------------------------------------------------------------------*/


#ifndef FOUR_C_CUT_COLOREDGRAPH_HPP
#define FOUR_C_CUT_COLOREDGRAPH_HPP

#include "4C_config.hpp"

#include "4C_cut_utils.hpp"

#include <list>
#include <map>

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  namespace Cut
  {
    namespace ColoredGraph
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

        void add(int row, int col);

        void add(int p, const plain_int_set& row);

        int find_next(Graph& used, int p, Graph& cycle, const plain_int_set& free);

        void find_free_facets(Graph& graph, Graph& used, plain_int_set& free,
            const std::vector<std::pair<Point*, Point*>>& all_lines, std::vector<int>& split_trace);

        void find_split_trace(std::vector<int>& split_trace);

        void get_all(plain_int_set& all);

        void fix_single_lines();

        void test_closed();

        void test_facets();

        void print() const;

        plain_int_set& at(int p) { return graph_[p]; }

        plain_int_set& operator[](int p) { return graph_[p]; }

        unsigned count(int p) { return graph_.count(p); }

        unsigned size() const { return graph_.size(); }

        const_iterator begin() const { return graph_.begin(); }

        const_iterator end() const { return graph_.end(); }

        void swap(Graph& other)
        {
          std::swap(graph_, other.graph_);
          std::swap(color_split_, other.color_split_);
        }

        int split() { return color_split_; }

        void set_split(int color_split) { color_split_ = color_split; }

        void split(Graph& used, plain_int_set& free, Graph& connection,
            const std::vector<int>& split_trace, Graph& c1, Graph& c2, Graph& datagraph);

        bool contains_trace(const std::vector<int>& split_trace);

        bool operator==(const Graph& other) const
        {
          return color_split_ == other.color_split_ and graph_ == other.graph_;
        }

        bool operator!=(const Graph& other) const { return not(*this == other); }

        void dump_graph(const std::string& name);

        void map(std::map<int, const void*>* index_value_map)
        {
          index_value_map_ = index_value_map;
        };

        // get pointer to the underlying object
        void* get_pointer(int Id) const
        {
          std::map<int, const void*>::iterator it = index_value_map_->find(Id);
          if (it != index_value_map_->end())
          {
            return const_cast<void*>(it->second);
          }
          else
            FOUR_C_THROW("Invalid element with id %d", Id);
          return nullptr;
        }
        // Split split trace into the closed loops
        void split_splittrace(const std::vector<int>& split_trace, Graph& datagraph,
            std::vector<std::vector<int>>& isolated_components);

       private:
        void fill(const std::vector<int>& split_trace, Graph& used, plain_int_set& free,
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

        void assign(Graph& cycle) { cycle_.swap(cycle); }

        void print() const;

        void split(Graph& used, plain_int_set& free, Graph& connection,
            const std::vector<int>& split_trace, Graph& c1, Graph& c2, Graph& datagraph)
        {
          cycle_.split(used, free, connection, split_trace, c1, c2, datagraph);
        }

        bool contains_trace(const std::vector<int>& split_trace)
        {
          return cycle_.contains_trace(split_trace);
        }

        Graph& operator()() { return cycle_; }

       private:
        Graph cycle_;
      };

      class CycleListIterator
      {
       public:
        CycleListIterator(std::list<Cycle>& cycles, std::list<Cycle>::iterator i) : i_(i)
        {
          next_active();
        }

        void next_active() {}

        CycleListIterator& operator++()
        {
          ++i_;
          next_active();
          return *this;
        }

        Graph& operator*() { return i_->cycle_; }

        bool operator!=(const CycleListIterator& other) { return i_ != other.i_; }

       private:
        std::list<Cycle>::iterator i_;
      };

      /// list of cycles
      class CycleList
      {
       public:
        typedef CycleListIterator iterator;

        void add_points(Graph& graph, Graph& used, Graph& cycle, plain_int_set& free,
            const std::vector<std::pair<Point*, Point*>>& all_lines);

        void print() const;

        unsigned size() const { return cycles_.size(); }

        CycleListIterator begin() { return CycleListIterator(cycles_, cycles_.begin()); }

        CycleListIterator end() { return CycleListIterator(cycles_, cycles_.end()); }

       private:
        void push_back(Graph& g);

        std::list<Cycle> cycles_;
      };

    }  // namespace ColoredGraph
  }    // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
