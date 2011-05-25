
#include <iostream>
#include <iterator>

#include "cut_pointgraph.H"
#include "cut_point.H"
#include "cut_line.H"
#include "cut_side.H"
#include "cut_element.H"
#include "cut_mesh.H"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/planar_face_traversal.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>

GEO::CUT::PointGraph::PointGraph( Mesh & mesh, Element * element, Side * side, bool inner )
  : element_( element ),
    side_( side )
{
  std::vector<int> cycle;
  FillGraph( element, side, cycle );

  if ( graph_.HasSinglePoints() )
  {
#if 1
    graph_.FixSinglePoints();
#else
    graph_.TestClosed();
#endif
  }

#if 0
#ifdef DEBUGCUTLIBRARY
  {
    std::ofstream f( "all_points.plot" );
    graph_.PlotAllPoints( f );
  }
  {
    std::ofstream f( "graph.txt" );
    graph_.Print( f );
  }
#endif
#endif

  graph_.FindCycles( cycle, cycles_, inner );
}

void GEO::CUT::PointGraph::FillGraph( Element * element, Side * side, std::vector<int> & cycle )
{
  const std::vector<Node*> & nodes = side->Nodes();
  const std::vector<Edge*> & edges = side->Edges();
  int end_pos = 0;
  for ( std::vector<Edge*>::const_iterator i=edges.begin(); i!=edges.end(); ++i )
  {
    Edge * e = *i;

    int begin_pos = end_pos;
    end_pos = ( end_pos + 1 ) % nodes.size();

    std::vector<Point*> edge_points;
    e->CutPoint( nodes[begin_pos], nodes[end_pos], edge_points );

    for ( unsigned i=1; i<edge_points.size(); ++i )
    {
      Point * p1 = edge_points[i-1];
      Point * p2 = edge_points[i];
      graph_.AddEdge( p1, p2 );
    }

    for ( std::vector<Point*>::iterator i=edge_points.begin()+1; i!=edge_points.end(); ++i )
    {
      Point * p = *i;
      cycle.push_back( p->Id() );
    }
  }

  const std::vector<Line*> & cut_lines = side->CutLines();

  for ( std::vector<Line*>::const_iterator i=cut_lines.begin(); i!=cut_lines.end(); ++i )
  {
    Line * l = *i;
    graph_.AddEdge( l->BeginPoint(), l->EndPoint() );
    if ( l->IsCut( element ) )
    {
      Point * p1 = l->BeginPoint();
      Point * p2 = l->EndPoint();
      if ( not p1->IsCut( element ) or
           not p2->IsCut( element ) )
      {
        throw std::runtime_error( "point-line inconsistency" );
      }
    }
  }
}


void GEO::CUT::PointGraph::Graph::AddEdge( int row, int col )
{
  graph_[row].insert( col );
  graph_[col].insert( row );
}


void GEO::CUT::PointGraph::Graph::AddEdge( Point * p1, Point * p2 )
{
  all_points_[p1->Id()] = p1;
  all_points_[p2->Id()] = p2;

  AddEdge( p1->Id(), p2->Id() );
}

void GEO::CUT::PointGraph::Graph::Print( std::ostream & stream)
{
  for ( std::map<int, std::set<int> >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
  {
    int p = i->first;
    std::set<int> & row = i->second;
    stream << p << ": ";
    for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
    {
      int p = *i;
      stream << p << " ";
    }
    stream << "\n";
  }
  stream << "\n";
}

void GEO::CUT::PointGraph::Graph::PlotAllPoints( std::ostream & stream )
{
  for ( std::map<int, Point*>::iterator i=all_points_.begin(); i!=all_points_.end(); ++i )
  {
    i->second->Plot( stream );
  }
}

void GEO::CUT::PointGraph::Graph::PlotPoints( Element * element )
{
  for ( std::map<int, Point*>::iterator i=all_points_.begin(); i!=all_points_.end(); ++i )
  {
    Point * p = i->second;
    std::cout << p->Id() << "(" << p->IsCut( element ) << ") ";
  }
  std::cout << "\n";
}

namespace GEO
{
  namespace CUT
  {
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

  }
}

void GEO::CUT::PointGraph::Graph::FindCycles( const std::vector<int> & cycle, std::vector<std::vector<Point*> > & cycles, bool inner )
{
  graph_t g;

  // create boost graph

  name_map_t name_map = boost::get( boost::vertex_name, g );
  edge_index_map_t edge_index_map = boost::get( boost::edge_index, g );

  std::map<int, vertex_t> vertex_map;

  for ( std::map<int, std::set<int> >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
  {
    int n = i->first;

    vertex_t u = add_vertex( g );
    name_map[u] = GetPoint( n );
    vertex_map[n] = u;
  }

  int counter = 0;

  for ( std::map<int, std::set<int> >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
  {
    int u = i->first;
    std::set<int> & row = i->second;

    for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
    {
      int v = *i;
      if ( u < v )
      {
        edge_t e;
        bool inserted;
        boost::tie( e, inserted ) = boost::add_edge( vertex_map[u], vertex_map[v], g );
        if ( inserted )
        {
          edge_index_map[e] = counter;
          counter += 1;
        }
      }
    }
  }

  // Test for planarity - we know it is planar, we just want to compute the
  // planar embedding as a side-effect

  typedef std::vector<edge_t> vec_t;
  std::vector<vec_t> embedding( boost::num_vertices( g ) );
  if ( not boost::boyer_myrvold_planarity_test( boost::boyer_myrvold_params::graph = g,
                                                boost::boyer_myrvold_params::embedding = &embedding[0] ) )
  {
    throw std::runtime_error( "input graph is not planar" );
  }

  face_visitor vis( name_map, cycles );
  boost::planar_face_traversal( g, &embedding[0], vis );
}

void GEO::CUT::PointGraph::Graph::FixSinglePoints()
{
  for ( ;; )
  {
    bool found = false;
    for ( std::map<int, std::set<int> >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
    {
      int p = i->first;
      std::set<int> & row = i->second;
      if ( row.size() < 2 )
      {
        found = true;
        for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
        {
          int p2 = *i;
          std::set<int> & row2 = graph_[p2];
          row2.erase( p );
          if ( row2.size()==0 )
            graph_.erase( p2 );
        }
        graph_.erase( p );
        break;
      }
    }
    if ( not found )
    {
      return;
    }
  }
}

bool GEO::CUT::PointGraph::Graph::HasSinglePoints()
{
  for ( std::map<int, std::set<int> >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
  {
    std::set<int> & row = i->second;
    if ( row.size() < 2 )
    {
      return true;
    }
  }
  return false;
}

