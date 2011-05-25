
#include <iostream>
#include <iterator>

#include "cut_pointgraph.H"
#include "cut_point.H"
#include "cut_line.H"
#include "cut_side.H"
#include "cut_element.H"
#include "cut_mesh.H"

GEO::CUT::PointGraph::PointGraph( Mesh & mesh, Element * element, Side * side, bool inner )
  : element_( element ),
    side_( side )
{
  std::vector<int> cycle;
  FillGraph( element, side, cycle );

#if 0
#ifdef DEBUGCUTLIBRARY
  if ( side->Id()==3 )
  {
    std::cout << "for element " << element->Id() << "\n";
    graph_.Print();
    std::copy( cycle.begin(), cycle.end(), std::ostream_iterator<int>( std::cout, " " ) );
    std::cout << "\n";
    graph_.PlotAllPoints();
    std::cout << "\n";
  }
#endif
#endif

  if ( graph_.HasSinglePoints() )
  {
#if 1
    graph_.FixSinglePoints();
#else
    graph_.TestClosed();
#endif
  }

  std::set<int> free;
  graph_.GetAll( free );

  for ( std::vector<int>::iterator i=cycle.begin(); i!=cycle.end(); ++i )
  {
    int p = *i;
    free.erase( p );
  }

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

  AddFacetPoints( cycle, free, inner );

#if 0
  std::cout << "GEO::CUT::PointGraph::PointGraph()\n";
  graph_.Print();
  std::copy( cycle.begin(), cycle.end(), std::ostream_iterator<int>( std::cout, " " ) );
  std::cout << "\n";
  graph_.PlotAllPoints();
  std::cout << "\n";
  graph_.PlotPoints( element );
  std::cout << "\n";

  for ( GRAPH::CycleListIterator i=facet_cycles_.begin(); i!=facet_cycles_.end(); ++i )
  {
    const std::vector<int> & points = *i;
    std::copy( points.begin(), points.end(), std::ostream_iterator<int>( std::cout, " " ) );
    std::cout << "\n";
  }

  std::cout << "\n";

  for ( iterator i=begin(); i!=end(); ++i )
  {
    const std::vector<Point*> & points = *i;
    std::copy( points.begin(), points.end(), std::ostream_iterator<Point*>( std::cout, " " ) );
    std::cout << "\n";
  }
#endif
}


void GEO::CUT::PointGraph::AddFacetPoints( std::vector<int> & cycle, std::set<int> & free, bool inner )
{
  GRAPH::Graph used;
  used.AddCycle( cycle );

  facet_cycles_.AddPoints( graph_, used, cycle, free );
  if ( not inner and free.size() > 0 )
  {
    facet_cycles_.AddFreePoints( graph_, used, free );
  }
}


void GEO::CUT::PointGraph::Graph::AddAll( Point * p1, Point * p2 )
{
  all_points_[p1->Id()] = p1;
  all_points_[p2->Id()] = p2;

  Add( p1->Id(), p2->Id() );
}

#if 0

// old version
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
      graph_.AddAll( p1, p2 );
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
    if ( l->IsCut( element ) )
    {
      graph_.AddAll( l->BeginPoint(), l->EndPoint() );
      if ( not l->BeginPoint()->IsCut( element ) or not l->EndPoint()->IsCut( element ) )
        throw std::runtime_error( "point-line inconsistency" );
    }
  }
}

#else

// current version
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
      graph_.AddAll( p1, p2 );
    }

    for ( std::vector<Point*>::iterator i=edge_points.begin()+1; i!=edge_points.end(); ++i )
    {
      Point * p = *i;
      cycle.push_back( p->Id() );
    }
  }

  const std::vector<Line*> & cut_lines = side->CutLines();

  std::vector<Line*> otherlines;

#if 0
  if ( side->Shape()==DRT::Element::quad4 )
  {
    const std::vector<Node*> & nodes = side->Nodes();
    int found = 0;

    for ( std::vector<Line*>::const_iterator i=cut_lines.begin(); i!=cut_lines.end(); ++i )
    {
      Line * l = *i;
      Point * p1 = l->BeginPoint();
      Point * p2 = l->EndPoint();
      if ( ( p1==nodes[0]->point() and p2==nodes[2]->point() ) or
           ( p1==nodes[2]->point() and p2==nodes[0]->point() ) or
           ( p1==nodes[1]->point() and p2==nodes[3]->point() ) or
           ( p1==nodes[3]->point() and p2==nodes[1]->point() ) )
      {
        found += 1;
      }
    }

    if ( found > 1 )
      throw std::runtime_error( "cut line crossing" );
  }
#endif
  for ( std::vector<Line*>::const_iterator i=cut_lines.begin(); i!=cut_lines.end(); ++i )
  {
    Line * l = *i;
    graph_.AddAll( l->BeginPoint(), l->EndPoint() );
    if ( not l->IsCut( element ) )
    {
      otherlines.push_back( l );
    }
    else
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

#if 0
  if ( otherlines.size() > 0 )
  {
    element->GnuplotDump();
    std::cout << "for element " << element->Id() << "\n";
    graph_.Print();
    std::copy( cycle.begin(), cycle.end(), std::ostream_iterator<int>( std::cout, " " ) );
    std::cout << "\n";
    graph_.PlotAllPoints();
    std::cout << "\n";
    for ( std::vector<Line*>::const_iterator i=otherlines.begin(); i!=otherlines.end(); ++i )
    {
      Line * l = *i;
      Point * p1 = l->BeginPoint();
      Point * p2 = l->EndPoint();
      //std::cout << "(" << p1->Id() << "," << p2->Id() << ") ";
      std::cout << "(" << p1->Id()
                << "," << p2->Id()
                << "," << p1->IsCut( element )
                << "," << p2->IsCut( element )
                << "," << l->IsCut( element )
                << ") ";
    }
    std::cout << "\n";
  }
#endif
}

#endif

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
