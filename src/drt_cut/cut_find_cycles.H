#ifndef CUT_FIND_CYCLES_H
#define CUT_FIND_CYCLES_H

#include <fstream>
#include <iostream>

#include <map>
#include <set>
#include <string>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/planar_face_traversal.hpp>

namespace GEO
{
  namespace CUT
  {
    class Point;

    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                                  boost::property<boost::vertex_name_t, Point*,
                                                  boost::property<boost::vertex_color_t, boost::default_color_type,
                                                                  boost::property<boost::vertex_index_t, int> > >,
                                  boost::property<boost::edge_index_t, int> > graph_t;

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

    typedef boost::color_traits<typename boost::property_traits<color_map_t>::value_type> color_t;

    struct face_visitor : public boost::planar_face_traversal_visitor
    {
      face_visitor( name_map_t name_map, std::vector<std::vector<Point*> > & cycles )
        : name_map_( name_map ),
          cycles_( cycles )
      {
      }

      void begin_face()
      {
        cycle.clear();
      }

      void end_face()
      {
        cycles_.push_back( std::vector<Point*>() );
        std::swap( cycles_.back(), cycle );
      }

      template <typename Vertex>
      void next_vertex( Vertex v )
      {
        cycle.push_back( name_map_[v] );
      }

      template <typename Edge>
      void next_edge( Edge e )
      {
      }

      name_map_t name_map_;
      std::vector<Point*> cycle;
      std::vector<std::vector<Point*> > & cycles_;
    };

    struct edge_filter
    {
      edge_filter()
        : g_( NULL ),
          component_( NULL ),
          c_( 0 )
      {
      }

      edge_filter( graph_t & g, const std::vector<int> & component, int c )
        : g_( &g ),
          component_( &component ),
          c_( c )
      {
      }

      template <typename Edge>
      bool operator()( const Edge & e ) const
      {
        vertex_t u = boost::source( e, *g_ );
        vertex_t v = boost::target( e, *g_ );

        vertex_index_map_t vertex_index_map = boost::get( boost::vertex_index, *g_ );

        return ( ( *component_ )[vertex_index_map[u]]==c_ and
                 ( *component_ )[vertex_index_map[v]]==c_ );
      }

      graph_t * g_;
      const std::vector<int> * component_;
      int c_;
    };

    typedef boost::filtered_graph<graph_t,edge_filter> filtered_graph_t;

    typedef std::vector<vertex_t> cycle_t;

    void find_cycles( graph_t & g, std::set<cycle_t*> & base_cycles );

  }
}

#endif
