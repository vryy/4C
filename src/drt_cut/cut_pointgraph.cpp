
#include <iostream>
#include <iterator>

#include "cut_pointgraph.H"
#include "cut_point.H"
#include "cut_line.H"
#include "cut_side.H"
#include "cut_element.H"
#include "cut_mesh.H"

GEO::CUT::PointGraph::PointGraph( Mesh & mesh, Element * element, Side * side, bool inner )
  : side_( side )
{
  std::vector<int> cycle;
  FillGraph( element, side, cycle );

  if ( graph_.HasSinglePoints() )
  {
//     graph_.FindSideLines( mesh, element, side );

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

  AddFacetPoints( cycle, free, inner );
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
void GEO::CUT::PointGraph::Graph::FindSideLines( Mesh & mesh, Element * element, Side * side )
{
  std::vector<int> side_position;
  std::vector<Point*> side_points;
  side_position.reserve( 2 );
  side_points.reserve( 2 );

  std::vector<int> open;

  const std::vector<Side*> & sides = element->Sides();
  for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
  {
    Side * s = *i;

    side_position.clear();
    side_points.clear();

    for ( std::map<int, Point*>::iterator i=all_points_.begin(); i!=all_points_.end(); ++i )
    {
      Point * p = i->second;
      if ( p->IsCut( s ) )
      {
        side_position.push_back( i->first );
        side_points.push_back( p );
      }
    }

    std::cout << ">> sid=" << side->Id()
              << "   #nodes=" << side->Nodes().size()
              << "   #points=" << side_points.size()
              << "   ";
    std::copy( side_position.begin(), side_position.end(), std::ostream_iterator<int>( std::cout, " " ) );

    // find if there is a new line needed and create

    if ( side_position.size() > 2 )
    {
      open.clear();

      for ( std::vector<int>::iterator i=side_position.begin(); i!=side_position.end(); ++i )
      {
        int p = *i;
        std::set<int> & row = at( p );
        if ( row.size() < 2 )
        {
          open.push_back( i-side_position.begin() );
        }
      }

      std::cout << " -- ";
      std::copy( open.begin(), open.end(), std::ostream_iterator<int>( std::cout, " " ) );
      std::cout << "\n";

      for ( std::vector<int>::iterator i=side_position.begin(); i!=side_position.end(); ++i )
      {
        int p = *i;
        Point * point = side_points[i-side_position.begin()];

        std::set<int> & row = at( p );
        for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
        {
          int p2 = *i;
          Point * point2 = GetPoint( p2 );

          //if ( point->Id() < point2->Id() )
          {
            std::cout << "   ";
            point->Print();
            std::cout << "   ";
            point2->Print();
            std::cout << "\n";
          }
        }
      }

      if ( open.size()==2 )
      {
        mesh.NewLine( side_points[open[0]], side_points[open[1]], s, side, element );
        Add( side_position[open[0]], side_position[open[1]] );
      }
    }

    std::cout << "\n";
  }
}
#endif

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
    }
  }
}
