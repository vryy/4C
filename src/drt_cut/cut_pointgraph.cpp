
#include <iostream>

#include "cut_pointgraph.H"
#include "cut_point.H"
#include "cut_line.H"
#include "cut_side.H"

GEO::CUT::PointGraph::PointGraph( Side * side )
  : side_( side )
{
  std::vector<int> cycle;

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

    for ( std::vector<Point*>::iterator i=edge_points.begin()+1; i!=edge_points.end(); ++i )
    {
      Point * p = *i;
      cycle.push_back( p->Id() );
    }
  }

  graph_.FixSinglePoints();

//  graph_.TestClosed();

  std::set<int> free;
  graph_.GetAll( free );

  for ( std::vector<int>::iterator i=cycle.begin(); i!=cycle.end(); ++i )
  {
    int p = *i;
    free.erase( p );
  }

  AddFacetPoints( cycle, free );
}

void GEO::CUT::PointGraph::AddFacetPoints( std::vector<int> & cycle, std::set<int> & free )
{
  GRAPH::Graph used;
  used.AddCycle( cycle );

  facet_cycles_.AddPoints( graph_, used, cycle, free );
}

void GEO::CUT::PointGraph::Graph::AddAll( Side * side, Point * p1, Point * p2 )
{
  all_points_[p1->Id()] = p1;
  all_points_[p2->Id()] = p2;

  std::set<int> & row = at( p1->Id() );
  if ( row.count( p2->Id() )==0 )
  {
    row.insert( p2->Id() );
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
