// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CUT_FIND_CYCLES_HPP
#define FOUR_C_CUT_FIND_CYCLES_HPP

#include "4C_config.hpp"

#include "4C_cut_cycle.hpp"

#include <boost/graph/adjacency_list.hpp>
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

namespace Cut
{
  class Point;

  namespace Impl
  {
    using graph_t = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
        boost::property<boost::vertex_name_t, Point*,
            boost::property<boost::vertex_color_t, boost::default_color_type,
                boost::property<boost::vertex_index_t, int>>>,
        boost::property<boost::edge_index_t, int>>;

    using vertex_t = boost::graph_traits<graph_t>::vertex_descriptor;
    using edge_t = boost::graph_traits<graph_t>::edge_descriptor;

    using vertex_iterator = boost::graph_traits<graph_t>::vertex_iterator;
    using edge_iterator = boost::graph_traits<graph_t>::edge_iterator;
    using adjacency_iterator = boost::graph_traits<graph_t>::adjacency_iterator;
    using out_edge_iterator = boost::graph_traits<graph_t>::out_edge_iterator;

    using name_map_t = boost::property_map<graph_t, boost::vertex_name_t>::type;
    using color_map_t = boost::property_map<graph_t, boost::vertex_color_t>::type;
    using vertex_index_map_t = boost::property_map<graph_t, boost::vertex_index_t>::type;
    using edge_index_map_t = boost::property_map<graph_t, boost::edge_index_t>::type;

    using color_t = boost::color_traits<boost::property_traits<color_map_t>::value_type>;

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

        return (
            (*component_)[vertex_index_map[u]] == c_ and (*component_)[vertex_index_map[v]] == c_);
      }

      graph_t* g_;
      const std::vector<int>* component_;
      int c_;
    };

    using filtered_graph_t = boost::filtered_graph<graph_t, EdgeFilter>;

    using cycle_t = std::vector<vertex_t>;

#ifdef CUT_USE_SORTED_VECTOR
    using plain_cycle_set = SortedVector<cycle_t*>;
    using plain_vertex_set = SortedVector<vertex_t>;
    using plain_graph_edge_set = SortedVector<std::pair<vertex_t, vertex_t>>;
#else
    using plain_cycle_set = std::set<cycle_t*>;
    using plain_vertex_set = std::set<vertex_t>;
    using plain_graph_edge_set = std::set<std::pair<vertex_t, vertex_t>>;
#endif

    void find_cycles(graph_t& g, plain_cycle_set& base_cycles);

  }  // namespace Impl
}  // namespace Cut


FOUR_C_NAMESPACE_CLOSE

#endif
