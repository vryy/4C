
#include <iostream>
#include <iterator>

#include <cmath>

#include "cut_pointgraph.H"
#include "cut_point.H"
#include "cut_line.H"
#include "cut_side.H"
#include "cut_element.H"
#include "cut_mesh.H"
#include "cut_find_cycles.H"

#include <boost/graph/copy.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/planar_face_traversal.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/connected_components.hpp>
//#include <boost/graph/subgraph.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graphviz.hpp>

GEO::CUT::PointGraph::PointGraph( Mesh & mesh, Element * element, Side * side, Location location, Strategy strategy )
  : element_( element ),
    side_( side )
{
  std::vector<Point*> cycle;
  FillGraph( element, side, cycle, strategy );

  if ( graph_.HasSinglePoints() )
  {
#if 1
    graph_.FixSinglePoints();
#else
    graph_.TestClosed();
#endif
  }

#if 1
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

  graph_.FindCycles( element, side, cycle, location, strategy );
}

void GEO::CUT::PointGraph::FillGraph( Element * element, Side * side, std::vector<Point*> & cycle, Strategy strategy )
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
    bool element_cut = l->IsCut( element );
    if ( strategy==all_lines or element_cut )
      graph_.AddEdge( l->BeginPoint(), l->EndPoint() );
    if ( element_cut )
    {
      Point * p1 = l->BeginPoint();
      Point * p2 = l->EndPoint();
      if ( not p1->IsCut( element ) or
           not p2->IsCut( element ) )
      {
        std::stringstream str;
        str << "line between " << ( *p1 ) << " and " << ( *p2 ) << " is cut by element, but point cuts are: "
            << p1->IsCut( element ) << " and "
            << p2->IsCut( element );
        throw std::runtime_error( str.str() );
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


bool Equals( const std::vector<Point*> & sorted, const std::vector<Point*> & test )
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

bool FindCycles( graph_t & g, std::vector<Point*> & cycle, std::map<vertex_t, LINALG::Matrix<3,1> > & local, std::vector<std::vector<Point*> > & cycles )
{
  name_map_t name_map = boost::get( boost::vertex_name, g );

  // Initialize the interior edge index
  edge_index_map_t e_index = boost::get( boost::edge_index, g );
  boost::graph_traits<graph_t>::edges_size_type edge_count = 0;
  boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
  for ( boost::tie( ei, ei_end ) = boost::edges( g ); ei != ei_end; ++ei )
    boost::put( e_index, *ei, edge_count++ );

  typedef std::vector<edge_t> vec_t;
  std::vector<vec_t> embedding( boost::num_vertices( g ) );

#if 1
  vertex_iterator vi, vi_end;
  for ( boost::tie( vi, vi_end )=boost::vertices( g ); vi!=vi_end; ++vi )
  {
    const LINALG::Matrix<3,1> & pos = local[*vi];

    std::map<double, vertex_t> arcs;

    adjacency_iterator ai, ai_end;
    for ( boost::tie( ai, ai_end )=boost::adjacent_vertices( *vi, g ); ai!=ai_end; ++ai )
    {
      LINALG::Matrix<3,1> d = local[*ai];
      d.Update( -1, pos, 1 );

      double arc = std::atan2( d( 1 ), d( 0 ) );

      std::map<double, vertex_t>::iterator j = arcs.find( arc );
      if ( j!=arcs.end() )
      {
        throw std::runtime_error( "numeric error: double arc" );
      }

      arcs[arc] = *ai;
    }

    vec_t & em = embedding[*vi];

    for ( std::map<double, vertex_t>::iterator i=arcs.begin(); i!=arcs.end(); ++i )
    {
      out_edge_iterator oi, oi_end;
      for ( boost::tie( oi, oi_end )=boost::out_edges( *vi, g ); oi!=oi_end; ++oi )
      {
        edge_t e = *oi;
        if ( boost::target( e, g )==i->second )
        {
          em.push_back( e );
          break;
        }
      }
    }
  }

#else
  if ( not boost::boyer_myrvold_planarity_test( boost::boyer_myrvold_params::graph = g,
                                                boost::boyer_myrvold_params::embedding = &embedding[0] ) )
  {
    throw std::runtime_error( "input graph is not planar" );
  }
#endif

#if 0
#ifdef DEBUGCUTLIBRARY
  std::cout << "embedding:\n";
  for ( std::vector<vec_t>::iterator i=embedding.begin(); i!=embedding.end(); ++i )
  {
    vec_t & em = *i;
    std::copy( em.begin(), em.end(), std::ostream_iterator<edge_t>( std::cout, " " ) );
    std::cout << "\n";
  }
#endif
#endif

  face_visitor vis( name_map, cycles );
  boost::planar_face_traversal( g, &embedding[0], vis );

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

  }
}

void GEO::CUT::PointGraph::Graph::FindCycles( Element * element, Side * side, std::vector<Point*> & cycle, Location location, Strategy strategy )
{
  graph_t g;

  // create boost graph

  name_map_t name_map = boost::get( boost::vertex_name, g );
  edge_index_map_t edge_index_map = boost::get( boost::edge_index, g );

  std::map<int, vertex_t> vertex_map;

  for ( std::map<int, std::set<int> >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
  {
    int n = i->first;

    Point * p = GetPoint( n );
//     if ( location==element_side or p->IsCut( element ) )
    {
      vertex_t u = add_vertex( g );
      name_map[u] = p;
      vertex_map[n] = u;
    }
  }

  int counter = 0;

  for ( std::map<int, std::set<int> >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
  {
    int u = i->first;

    Point * p1 = GetPoint( u );
    if ( location==element_side or p1->IsCut( element ) )
    {
      std::set<int> & row = i->second;

      for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
      {
        int v = *i;
        Point * p2 = GetPoint( v );
        if ( location==element_side or p2->IsCut( element ) )
        {
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
    }
  }

  if ( strategy==own_lines )
  {
    std::set<cycle_t*> base_cycles;
    find_cycles( g, base_cycles );

    main_cycles_.reserve( base_cycles.size() );

    for ( std::set<cycle_t*>::iterator i=base_cycles.begin(); i!=base_cycles.end(); ++i )
    {
      cycle_t * c = *i;

      main_cycles_.push_back( std::vector<Point*>() );
      std::vector<Point*> & pc = main_cycles_.back();
      pc.reserve( c->size() );

      for ( cycle_t::iterator i=c->begin(); i!=c->end(); ++i )
      {
        vertex_t u = *i;
        pc.push_back( name_map[u] );
      }

      delete c;
    }
  }
  else
  {

    // find local coordinates

    std::map<vertex_t, LINALG::Matrix<3,1> > local;

    vertex_iterator vi, vi_end;
    for ( boost::tie( vi, vi_end )=boost::vertices( g ); vi!=vi_end; ++vi )
    {
      Point * p = name_map[*vi];
      LINALG::Matrix<3,1> xyz( p->X() );
      side->LocalCoordinates( xyz, local[*vi] );
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
      bool main_cycle = GEO::CUT::FindCycles( g, cycle, local, main_cycles_ );
      if ( location==element_side and not main_cycle )
      {
        GnuplotDumpCycles( "cycles", main_cycles_ );
        boost::print_graph( g, boost::get( boost::vertex_name, g ) );

#if 0
        // Output the graph in DOT format
        boost::dynamic_properties dp;
        dp.property( "label", boost::get( boost::vertex_index, g ) );
        std::ofstream out( "side-graph.dot" );
        boost::write_graphviz( out, g, dp, std::string(), boost::get( boost::vertex_index, g ) );
#endif

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

        graph_t cg;
        boost::copy_graph( fg, cg );

        bool main_cycle = GEO::CUT::FindCycles( cg, cycle, local, filtered_cycles );

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

      if ( location==element_side and main_cycles_.size()==0 )
      {
        throw std::runtime_error( "cycle needs to contain side edges" );
      }
    }
    else
    {
      if ( location==element_side )
        throw std::runtime_error( "empty graph discovered" );
    }
  }
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

