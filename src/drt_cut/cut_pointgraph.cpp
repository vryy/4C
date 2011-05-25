
#include <iostream>
#include <iterator>

#include "cut_pointgraph.H"
#include "cut_point.H"
#include "cut_line.H"
#include "cut_side.H"
#include "cut_element.H"
#include "cut_mesh.H"

#include <boost/graph/copy.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/planar_face_traversal.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/connected_components.hpp>
//#include <boost/graph/subgraph.hpp>
#include <boost/graph/filtered_graph.hpp>

GEO::CUT::PointGraph::PointGraph( Mesh & mesh, Element * element, Side * side, bool inner )
  : element_( element ),
    side_( side )
{
  std::vector<Point*> cycle;
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

  graph_.FindCycles( cycle, inner );
}

void GEO::CUT::PointGraph::FillGraph( Element * element, Side * side, std::vector<Point*> & cycle )
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
      cycle.push_back( p );
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
  }
}

void GEO::CUT::PointGraph::Graph::FindCycles( std::vector<Point*> & cycle, bool inner )
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

  // find unconnected components (main facet(s) and holes)

  std::vector<int> component( boost::num_vertices( g ) );

  int num_comp =
    boost::connected_components( g,
                                 boost::make_iterator_property_map( component.begin(),
                                                                    boost::get( boost::vertex_index, g ) ) );

  // find cycles on each component

  std::sort( cycle.begin(), cycle.end() );

  if ( num_comp == 1 )
  {
    bool main_cycle = FindCycles( g, cycle, main_cycles_ );
    if ( inner and not main_cycle )
    {
      GnuplotDumpCycles( "cycles", main_cycles_ );
      boost::print_graph( g, boost::get( boost::vertex_name, g ) );
      throw std::runtime_error( "cycle needs to contain side edges" );
    }
  }
  else if ( num_comp > 1 )
  {
    for ( int i=0; i<num_comp; ++i )
    {
      typedef boost::filtered_graph<graph_t,edge_filter> filtered_graph_t;
      edge_filter filter( g, component, i );
      filtered_graph_t fg( g, filter );

      std::vector<std::vector<Point*> > filtered_cycles;
#if 1
      bool main_cycle = FindCycles<filtered_graph_t>( fg, cycle, filtered_cycles );
#else
      graph_t cg;
      boost::copy_graph( fg, cg );

      bool main_cycle = FindCycles( cg, cycle, filtered_cycles );
#endif

      if ( main_cycle )
      {
        if ( main_cycles_.size()!=0 )
        {
          throw std::runtime_error( "one set of main cycles only" );
        }
        std::swap( main_cycles_, filtered_cycles );
      }
      else
      {
        hole_cycles_.push_back( std::vector<std::vector<Point*> >() );
        std::swap( hole_cycles_.back(), filtered_cycles );
      }
    }

    if ( inner and main_cycles_.size()==0 )
    {
      throw std::runtime_error( "cycle needs to contain side edges" );
    }
  }
  else
  {
    throw std::runtime_error( "empty graph discovered" );
  }
}

template <typename graph_t>
bool GEO::CUT::PointGraph::Graph::FindCycles( graph_t & g, std::vector<Point*> & cycle, std::vector<std::vector<Point*> > & cycles )
{
  // Test for planarity - we know it is planar, we just want to compute the
  // planar embedding as a side-effect

  typedef std::vector<edge_t> vec_t;
  std::vector<vec_t> embedding( boost::num_vertices( g ) );
  if ( not boost::boyer_myrvold_planarity_test( boost::boyer_myrvold_params::graph = g,
                                                boost::boyer_myrvold_params::embedding = &embedding[0] ) )
  {
    throw std::runtime_error( "input graph is not planar" );
  }

  name_map_t name_map = boost::get( boost::vertex_name, g );

  face_visitor vis( name_map, cycles );
  boost::planar_face_traversal( g, &embedding[0], vis );

//   std::vector<std::vector<Point*> > backup = cycles;

  bool save_first = cycles.size()==2;

  int erase_count = 0;
  for ( std::vector<std::vector<Point*> >::iterator i=cycles.begin(); i!=cycles.end(); )
  {
    std::vector<Point*> & c = *i;
    if ( Equals( cycle, c ) )
    {
      if ( save_first and erase_count == 0 )
      {
        ++i;
      }
      else
      {
        cycles.erase( i );
      }
      erase_count += 1;
    }
    else
    {
      ++i;
    }
  }

  if ( erase_count > ( save_first ? 2 : 1 ) )
  {
    throw std::runtime_error( "more than one back facet" );
  }

  return erase_count != 0;
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

bool GEO::CUT::PointGraph::Graph::Equals( const std::vector<Point*> & sorted, const std::vector<Point*> & test )
{
  if ( sorted.size()!=test.size() )
  {
    return false;
  }

  for ( std::vector<Point*>::const_iterator i=test.begin(); i!=test.end(); ++i )
  {
    Point * p = *i;
    if ( not std::binary_search( sorted.begin(), sorted.end(), p ) )
    {
      return false;
    }
  }
  return true;
}

void GEO::CUT::PointGraph::Graph::GnuplotDumpCycles( const std::string & filename, const std::vector<std::vector<Point*> > & cycles )
{
  int counter = 0;
  for ( std::vector<std::vector<Point*> >::const_iterator i=cycles.begin(); i!=cycles.end(); ++i )
  {
    const std::vector<Point*> & points = *i;

    std::stringstream str;
    str << filename << counter << ".plot";
    std::cout << str.str() << "\n";
    std::ofstream file( str.str().c_str() );

    for ( unsigned i=0; i!=points.size(); ++i )
    {
      Point * p1 = points[i];
      Point * p2 = points[( i+1 ) % points.size()];

      p1->Plot( file );
      p2->Plot( file );
      file << "\n\n";
    }

    counter += 1;
  }
}

