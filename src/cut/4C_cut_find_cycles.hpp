/*---------------------------------------------------------------------*/
/*! \file

\brief part of the facetgraph to create facets from lines and edges

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_FIND_CYCLES_HPP
#define FOUR_C_CUT_FIND_CYCLES_HPP

#include "4C_config.hpp"

#include "4C_cut_cycle.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/planar_face_traversal.hpp>

#include <fstream>
#include <map>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  namespace Cut
  {
    class Point;

    namespace Impl
    {
      typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
          boost::property<boost::vertex_name_t, Point*,
              boost::property<boost::vertex_color_t, boost::default_color_type,
                  boost::property<boost::vertex_index_t, int>>>,
          boost::property<boost::edge_index_t, int>>
          graph_t;

      typedef boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
      typedef boost::graph_traits<graph_t>::edge_descriptor edge_t;

      typedef boost::graph_traits<graph_t>::vertex_iterator vertex_iterator;
      typedef boost::graph_traits<graph_t>::edge_iterator edge_iterator;
      typedef boost::graph_traits<graph_t>::adjacency_iterator adjacency_iterator;
      typedef boost::graph_traits<graph_t>::out_edge_iterator out_edge_iterator;

      typedef boost::property_map<graph_t, boost::vertex_name_t>::type name_map_t;
      typedef boost::property_map<graph_t, boost::vertex_color_t>::type color_map_t;
      typedef boost::property_map<graph_t, boost::vertex_index_t>::type vertex_index_map_t;
      typedef boost::property_map<graph_t, boost::edge_index_t>::type edge_index_map_t;

      typedef boost::color_traits<boost::property_traits<color_map_t>::value_type> color_t;

      /// boost::graph visitor that is used to create a spanning tree from a graph
      class SpanningTreeCreator : public boost::default_bfs_visitor
      {
       public:
        SpanningTreeCreator(graph_t& st) : st_(st) {}

        void tree_edge(edge_t e, const graph_t& g)
        {
          vertex_t u = boost::source(e, g);
          vertex_t v = boost::target(e, g);

          if (u > v)
          {
            std::swap(u, v);
          }

          // assume same vertex ids
          boost::add_edge(u, v, st_);
        }

       private:
        graph_t& st_;
      };

      /// boost::graph visitor
      class SpanningTreePathFinder : public boost::default_dfs_visitor
      {
       public:
        SpanningTreePathFinder(vertex_t v, std::vector<vertex_t>& path)
            : done_(false), v_(v), path_(path)
        {
        }

        void discover_vertex(vertex_t u, const graph_t& g)
        {
          if (not done_)
          {
            if (u == v_) done_ = true;
            path_.push_back(u);
          }
        }

        bool operator()(vertex_t u, const graph_t& g) { return done_; }

        void finish_vertex(vertex_t u, const graph_t& g)
        {
          if (not done_)
          {
            if (path_.back() != u) FOUR_C_THROW("confused");
            path_.pop_back();
          }
        }

       private:
        bool done_;
        vertex_t v_;
        std::vector<vertex_t>& path_;
      };

      /// visitor that collects cycles in a planar graph
      struct FaceVisitor : public boost::planar_face_traversal_visitor
      {
        FaceVisitor(name_map_t name_map, std::vector<Cycle>& cycles)
            : name_map_(name_map), cycles_(cycles)
        {
        }

        void begin_face() { cycle.clear(); }

        void end_face()
        {
          cycles_.push_back(Cycle());
          std::swap(cycles_.back(), cycle);
        }

        template <typename Vertex>
        void next_vertex(Vertex v)
        {
          cycle.push_back(name_map_[v]);
        }

        template <typename Edge>
        void next_edge(Edge e)
        {
        }

        name_map_t name_map_;
        Cycle cycle;
        std::vector<Cycle>& cycles_;
      };

      /// filter that selects all edges to a given marker
      struct EdgeFilter
      {
        EdgeFilter() : g_(nullptr), component_(nullptr), c_(0) {}

        EdgeFilter(graph_t& g, const std::vector<int>& component, int c)
            : g_(&g), component_(&component), c_(c)
        {
        }

        template <typename Edge>
        bool operator()(const Edge& e) const
        {
          vertex_t u = boost::source(e, *g_);
          vertex_t v = boost::target(e, *g_);

          vertex_index_map_t vertex_index_map = boost::get(boost::vertex_index, *g_);

          return ((*component_)[vertex_index_map[u]] == c_ and
                  (*component_)[vertex_index_map[v]] == c_);
        }

        graph_t* g_;
        const std::vector<int>* component_;
        int c_;
      };

      typedef boost::filtered_graph<graph_t, EdgeFilter> filtered_graph_t;

      typedef std::vector<vertex_t> cycle_t;

#ifdef CUT_USE_SORTED_VECTOR
      typedef SortedVector<cycle_t*> plain_cycle_set;
      typedef SortedVector<vertex_t> plain_vertix_set;
      typedef SortedVector<std::pair<vertex_t, vertex_t>> plain_graph_edge_set;
#else
      typedef std::set<cycle_t*> plain_cycle_set;
      typedef std::set<vertex_t> plain_vertix_set;
      typedef std::set<std::pair<vertex_t, vertex_t>> plain_graph_edge_set;
#endif

      void find_cycles(graph_t& g, plain_cycle_set& base_cycles);

    }  // namespace Impl
  }    // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
