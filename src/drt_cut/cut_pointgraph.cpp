
#include <iostream>

#include "cut_pointgraph.H"
#include "cut_point.H"
#include "cut_line.H"
#include "cut_side.H"

GEO::CUT::PointGraph::PointGraph( Side * side )
  : side_( side )
{
  std::vector<Point*> cycle;

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
      graph_.AddAll( side_, p1, p2 );
      graph_.AddAll( side_, p2, p1 );
    }

    std::copy( edge_points.begin()+1, edge_points.end(), std::back_inserter( cycle ) );
  }

  graph_.IsValid();

  std::set<Point*> free;
  for ( std::map<Point*, std::set<Point*> >::iterator i=graph_.graph_.begin(); i!=graph_.graph_.end(); ++i )
  {
    Point * p = i->first;
    free.insert( p );
  }
  for ( std::vector<Point*>::iterator i=cycle.begin(); i!=cycle.end(); ++i )
  {
    Point * p = *i;
    free.erase( p );
  }

  AddFacetPoints( cycle, free );
}

void GEO::CUT::PointGraph::AddFacetPoints( std::vector<Point*> & cycle, std::set<Point*> & free )
{
  Graph used;
  used.AddCycle( cycle );

  facet_cycles_.AddFacetPoints( graph_, used, cycle, free );
}

void GEO::CUT::PointGraph::FacetCycleList::AddFacetPoints( Graph & graph, Graph & used, std::vector<Point*> & cycle, std::set<Point*> & free )
{
  facet_cycles_.push_back( FacetCycle() );
  FacetCycle & fc = facet_cycles_.back();
  fc.Assign( cycle );

  fc.Split( graph, used, *this, free );
}

void GEO::CUT::PointGraph::FacetCycle::Split( Graph & graph, Graph & used, FacetCycleList & facet_cycles, std::set<Point*> & free )
{
  std::vector<Point*>::iterator start = std::find_if( cycle_.begin(), cycle_.end(), ForkFinder( graph, used, cycle_, free ) );
  if ( start!=cycle_.end() )
  {
    std::vector<Point*> connection;

    Point * p1 = *start;
    while ( p1!=NULL )
    {
      Point * p2 = graph.FindNext( used, p1, cycle_, free );
      if ( p2 != NULL )
      {
        used.Add( p1, p2 );

        std::vector<Point*>::iterator end = std::find( cycle_.begin(), cycle_.end(), p2 );
        if ( end!=cycle_.end() )
        {
          // split here.

          if ( end < start )
            throw std::runtime_error( "iterator confusion" );

          std::vector<Point*> c1;
          std::vector<Point*> c2;

          std::copy( cycle_.begin(), start, std::back_inserter( c1 ) );
          c1.push_back( *start );
          std::copy( connection.begin(), connection.end(), std::back_inserter( c1 ) );
          std::copy( end, cycle_.end(), std::back_inserter( c1 ) );

          std::copy( start, end, std::back_inserter( c2 ) );
          c2.push_back( *end );
          std::copy( connection.rbegin(), connection.rend(), std::back_inserter( c2 ) );

          for ( std::vector<Point*>::iterator i=connection.begin(); i!=connection.end(); ++i )
          {
            Point * p;
            free.erase( p );
          }

          facet_cycles.AddFacetPoints( graph, used, c1, free );
          facet_cycles.AddFacetPoints( graph, used, c2, free );

          active_ = false;

          p1 = NULL;
        }
        else
        {
          connection.push_back( p2 );
          p1 = p2;
        }
      }
      else
      {
        throw std::runtime_error( "failed to find next point" );
      }
    }
  }
  else
  {
    active_ = true;
  }
}

void GEO::CUT::PointGraph::FacetCycle::Print()
{
  std::cout << active_;
  for ( std::vector<Point*>::iterator i=cycle_.begin(); i!=cycle_.end(); ++i )
  {
    Point * p = *i;
    std::cout << " " << p->Id();
  }
  std::cout << "\n";
}

void GEO::CUT::PointGraph::Graph::Add( Point * p1, Point * p2 )
{
  graph_[p1].insert( p2 );
  graph_[p2].insert( p1 );
}

void GEO::CUT::PointGraph::Graph::AddAll( Side * side, Point * p1, Point * p2 )
{
  std::set<Point*> & row = graph_[p1];
  if ( row.count( p2 )==0 )
  {
    row.insert( p2 );
    std::set<Line*> cut_lines;
    p2->CutLines( side, cut_lines );
    for ( std::set<Line*>::iterator i=cut_lines.begin(); i!=cut_lines.end(); ++i )
    {
      Line * l = *i;
      Point * p3 = l->OtherPoint( p2 );
      if ( p3!=p1 )
      {
        AddAll( side, p2, p3 );
        AddAll( side, p3, p2 );
      }
    }
  }
}

void GEO::CUT::PointGraph::Graph::AddCycle( const std::vector<Point*> & cycle )
{
  unsigned size = cycle.size();
  for ( unsigned i=0; i<size; ++i )
  {
    Point * p1 = cycle[i];
    Point * p2 = cycle[( i+1 ) % size];

    Add( p1, p2 );
  }
}

GEO::CUT::Point * GEO::CUT::PointGraph::Graph::FindNext( Graph & used, Point * point, const std::vector<Point*> & cycle, const std::set<Point*> & free )
{
  std::set<Point*> & row      = graph_[point];
  std::set<Point*> & used_row = used[point];
  for ( std::set<Point*>::iterator i=row.begin(); i!=row.end(); ++i )
  {
    Point * p = *i;
    if ( used_row.count( p ) == 0 )
    {
      if ( free.count( p ) > 0 )
        return p;
      if ( std::find( cycle.begin(), cycle.end(), p )!=cycle.end() )
        return p;
    }
  }
  return NULL;
}

void GEO::CUT::PointGraph::Graph::IsValid()
{
  for ( std::map<Point*, std::set<Point*> >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
  {
    //Point * p = i->first;
    std::set<Point*> & row = i->second;
    if ( row.size() < 2 )
    {
      throw std::runtime_error( "open point in graph" );
    }
  }
}

void GEO::CUT::PointGraph::Graph::Print()
{
  for ( std::map<Point*, std::set<Point*> >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
  {
    Point * p = i->first;
    std::set<Point*> & row = i->second;
    std::cout << p->Id() << ": ";
    for ( std::set<Point*>::iterator i=row.begin(); i!=row.end(); ++i )
    {
      Point * p = *i;
      std::cout << p->Id() << " ";
    }
    std::cout << "\n";
  }
  std::cout << "---- << >> ----\n";
}

void GEO::CUT::PointGraph::FacetCycleList::Print()
{
  for ( std::list<FacetCycle>::iterator i=facet_cycles_.begin(); i!=facet_cycles_.end(); ++i )
  {
    FacetCycle & fc = *i;
    fc.Print();
  }
}
