
#include <iostream>
#include <iterator>

#include "cut_pointgraph.H"
#include "cut_point.H"
#include "cut_line.H"
#include "cut_side.H"
#include "cut_element.H"
#include "cut_mesh.H"
#include "cut_find_cycles.H"

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

  // find cycles

  std::sort( cycle.begin(), cycle.end() );

  hole_cycles_.reserve( num_comp-1 );

  for ( int i=0; i<num_comp; ++i )
  {
    edge_filter filter( g, component, i );
    filtered_graph_t fg( g, filter );

    graph_t cg;
    boost::copy_graph( fg, cg );

    std::set<cycle_t*> base_cycles;
    find_cycles( cg, base_cycles );

    // see if element nodes included in graph

    name_map_t name_map = boost::get( boost::vertex_name, cg );

    bool main_graph = false;
    vertex_iterator vi, vi_end;
    for ( boost::tie( vi, vi_end )=boost::vertices( cg ); vi!=vi_end; ++vi )
    {
      if ( boost::out_degree( *vi, cg ) > 0 )
      {
        Point * p = name_map[*vi];
        if ( std::binary_search( cycle.begin(), cycle.end(), p ) )
        {
          main_graph = true;
          break;
        }
      }
    }

    std::vector<std::vector<Point*> > * current = NULL;
    if ( main_graph )
    {
      if ( main_cycles_.size()!=0 )
        throw std::runtime_error( "more than one main graph" );
      current = &main_cycles_;
    }
    else
    {
      hole_cycles_.push_back( std::vector<std::vector<Point*> >() );
      current = &hole_cycles_.back();
    }

    // cleanup

    current->reserve( base_cycles.size() );

    for ( std::set<cycle_t*>::iterator i=base_cycles.begin(); i!=base_cycles.end(); ++i )
    {
      cycle_t * c = *i;

      current->push_back( std::vector<Point*>() );
      std::vector<Point*> & pc = current->back();
      pc.reserve( c->size() );

      for ( cycle_t::iterator i=c->begin(); i!=c->end(); ++i )
      {
        vertex_t u = *i;
        Point * p = name_map[u];
        pc.push_back( p );
      }

      delete c;
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

