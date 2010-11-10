
#include "../drt_geometry/intersection_templates.H"

#ifdef QHULL
#undef PI
#ifdef TETGENINCLUDED
#include "../tetgen/tetgen.h"
#else
#include <tetgen.h>
#endif
#undef DOT
#endif

#include "cut_position.H"
#include "cut_position2d.H"
#include "cut_intersection.H"
#include "cut_facet.H"
#include "cut_point_impl.H"
#include "cut_pointcycle.H"
#include "cut_linesegment.H"

#include <string>
#include <stack>

#include "cut_side.H"

void GEO::CUT::Side::FillComplete( Mesh & mesh )
{
  for ( std::vector<Edge*>::iterator i=edges_.begin(); i!=edges_.end(); ++i )
  {
    Edge * edge = *i;
    edge->FillComplete( mesh );
  }
}

GEO::CUT::Edge * GEO::CUT::Side::FindEdge( Point * begin, Point * end )
{
  for ( std::vector<Edge*>::iterator i=edges_.begin(); i!=edges_.end(); ++i )
  {
    Edge * e = *i;
    if ( e->Matches( begin, end ) )
    {
      return e;
    }
  }
  return NULL;
}

bool GEO::CUT::LinearSide::Cut( Mesh & mesh, Side & side, Element * element )
{
  // see if the cut is already there
  bool found_cut = false;
  for ( std::vector<Line*>::iterator i=cut_lines_.begin(); i!=cut_lines_.end(); ++i )
  {
    Line* cut_line = *i;
    if ( cut_line->IsCut( this, &side ) )
    {
      cut_line->AddSide( this );
      cut_line->AddSide( &side );
      cut_line->AddElement( element );
      AddLine( cut_line );
      side.AddLine( cut_line );
      found_cut = true;
    }
  }
  if ( found_cut )
  {
    return true;
  }

  // look at edges for cut
  std::set<Point*, PointPidLess> cuts;
  EdgeCuts( mesh, side, cuts );

  // if edges are cut, create cut line
  if ( cuts.size()>0 )
  {
    if ( cuts.size()==1 )
    {
      // a partial cut with the cut side ending inside the cutted element
      // there could (should?) be a reverse cut
      std::set<Point*, PointPidLess> reverse_cuts;
      side.EdgeCuts( mesh, *this, reverse_cuts );
      reverse_cuts.erase( *cuts.begin() );
      if ( reverse_cuts.size()==1 )
      {
        Line* cut_line = mesh.NewLine( *cuts.begin(), *reverse_cuts.begin(), this, &side, element );
        AddLine( cut_line );
        side.AddLine( cut_line );
        return true;
      }
      else if ( reverse_cuts.size()==0 )
      {
        // Touch of two edges. No lines to create?!
        return false;
      }
      else
      {
        throw std::runtime_error( "most peculiar cut" );
      }
    }
    else if ( cuts.size()==2 )
    {
      // The normal case. A straight cut.
      std::vector<Point*> c( cuts.begin(), cuts.end() );
      Line* cut_line = mesh.NewLine( c[0], c[1], this, &side, element );
      AddLine( cut_line );
      side.AddLine( cut_line );
      return true;
    }
#if 0
    else if ( cuts.size()==3 )
    {
      // A touch between two elements. There are already some lines. But with
      // tree points there is no question, we need three lines. So lets just
      // create them.
      std::vector<Point*> c( cuts.begin(), cuts.end() );

      Line* cut_line = mesh.NewLine( c[0], c[1], this, &side, element );
      AddLine( cut_line );
      side.AddLine( cut_line );

      cut_line = mesh.NewLine( c[0], c[2], this, &side, element );
      AddLine( cut_line );
      side.AddLine( cut_line );

      cut_line = mesh.NewLine( c[1], c[2], this, &side, element );
      AddLine( cut_line );
      side.AddLine( cut_line );

      return true;
    }
#endif
    else
    {
      // If all nodes are catched and nothing else, the cut surface has hit
      // this surface exactly. No need to cut anything. However, the surface
      // might be required for integration.

      bool allonnode = true;
      for ( std::set<Point*, PointPidLess>::iterator i=cuts.begin();
            i!=cuts.end();
            ++i )
      {
        Point * p = *i;
        if ( not p->NodalPoint() )
        {
          allonnode = false;
          break;
        }
      }
      if ( allonnode )
      {
        if ( cuts.size()==Nodes().size() )
        {
          for ( unsigned i=0; i<Nodes().size(); ++i )
          {
            unsigned j = ( i+1 ) % Nodes().size();
            Line* cut_line = mesh.NewLine( Nodes()[i]->point(), Nodes()[j]->point(), this, &side, element );
            AddLine( cut_line );
            side.AddLine( cut_line );
          }
          return true;
        }
        else
        {
          //throw std::runtime_error( "Only nodes touched, but not all nodes covered? How very particular." );
          //
          // Not a cut of this side?!
        }
      }
      else
      {
        // Now we have more than two cuts on the surface. This is a touch, but
        // one that should have been tracked already. I do not have a way to
        // construct all lines here. It is not known which points are to be
        // connected.
        //
        // There should be some lines already. Otherwise we would not get so
        // many cuts between two sides. And thus the line segment
        // reconstruction will recover the line path later on.
        //
        // Maybe this case do never happen at all.
      }
    }
  }
  return false;
}

void GEO::CUT::Side::EdgeCuts( Mesh & mesh, Side & side, std::set<Point*, PointPidLess> & cuts )
{
  for ( std::vector<Edge*>::iterator i=edges_.begin(); i!=edges_.end(); ++i )
  {
    Edge & edge = **i;
    edge.Cut( mesh, side, cuts );
  }
}

void GEO::CUT::LinearSide::AddLine( Line* cut_line )
{
  if ( std::find( cut_lines_.begin(), cut_lines_.end(), cut_line )==cut_lines_.end() )
  {
    cut_lines_.push_back( cut_line );
  }
}

GEO::CUT::Facet * GEO::CUT::LinearSide::FindFacet( const std::vector<Point*> & facet_points )
{
  for ( std::vector<Facet*>::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    if ( f->Equals( facet_points ) )
    {
      return f;
    }
  }
  return NULL;
}

void GEO::CUT::LinearSide::CreateLineSegmentList( Mesh & mesh,
                                                  Element * element,
                                                  std::vector<Teuchos::RCP<LineSegment> > & segments,
                                                  bool inner )
{
  std::set<Line*> cut_lines;

  for ( std::vector<Line*>::const_iterator i=cut_lines_.begin(); i!=cut_lines_.end(); ++i )
  {
    Line * l = *i;
    if ( l->IsCut( element ) )
    {
      cut_lines.insert( l );
    }
  }

  while ( cut_lines.size() )
  {
    segments.push_back( Teuchos::rcp( new LineSegment( mesh, element, this, cut_lines, inner ) ) );
  }
}

Teuchos::RCP<GEO::CUT::LineSegment> GEO::CUT::LinearSide::CreateLineSegment( Mesh & mesh,
                                                                             Element * element )
{
  std::vector<Teuchos::RCP<LineSegment> > segments;

  CreateLineSegmentList( mesh, element, segments, true );

  for ( unsigned i=0; i<segments.size(); ++i )
  {
    if ( segments[i] != Teuchos::null )
    {
      LineSegment & s1 = *segments[i];
      for ( unsigned j=i+1; j<segments.size(); ++j )
      {
        if ( segments[j] != Teuchos::null )
        {
          LineSegment & s2 = *segments[j];
          if ( s1.Combine( mesh, element, this, s2 ) )
          {
            segments[j] = Teuchos::null;
          }
        }
      }
    }
  }

  std::vector<Teuchos::RCP<LineSegment> >::iterator i
    = std::remove_if( segments.begin(), segments.end(),
                      std::bind2nd( std::equal_to<Teuchos::RCP<LineSegment> >(), Teuchos::null ) );
  segments.erase( i, segments.end() );

  if ( segments.size()!=1 )
  {
    throw std::runtime_error( "expect one segment" );
  }
  return segments[0];
}

void GEO::CUT::LinearSide::MakeOwnedSideFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets )
{
  if ( facets_.size()==0 )
  {
    std::vector<Point*> facet_points;

    int end_pos = 0;
    for ( std::vector<Edge*>::const_iterator i=Edges().begin(); i!=Edges().end(); ++i )
    {
      Edge * e = *i;

      int begin_pos = end_pos;
      end_pos = ( end_pos + 1 ) % Nodes().size();

      std::vector<Point*> edge_points;
      e->CutPoint( Nodes()[begin_pos], Nodes()[end_pos], edge_points );

      std::copy( edge_points.begin()+1, edge_points.end(), std::back_inserter( facet_points ) );
    }

    std::set<Point*, PointPidLess> cut_points;
    for ( std::vector<Edge*>::const_iterator i=Edges().begin(); i!=Edges().end(); ++i )
    {
      Edge * e = *i;
      e->CutPoints( this, cut_points );
    }

    if ( cut_points.size()>0 )
    {
      PointCycleList pcl( element, this, facet_points, cut_points );
      pcl.CreateFacets( mesh, this, facets_ );
    }
    else
    {
      // Just a normal side. There might be cut points at the edges that do
      // not concern this side. We have to copy them anyway.

      facets_.push_back( mesh.NewFacet( facet_points, facet_points.size(), this ) );
    }
  }

  std::copy( facets_.begin(), facets_.end(), std::inserter( facets, facets.begin() ) );
}

void GEO::CUT::LinearSide::MakeSideCutFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets )
{
  std::vector<Teuchos::RCP<LineSegment> > segments;
//   CreateLineSegmentList( mesh, element, segments, false );

  std::set<Line*> cut_lines;

  for ( std::vector<Line*>::const_iterator i=cut_lines_.begin(); i!=cut_lines_.end(); ++i )
  {
    Line * l = *i;
    if ( l->IsCut( element ) and
         not OnEdge( l->BeginPoint() ) and
         not OnEdge( l->EndPoint() ) )
    {
      cut_lines.insert( l );
    }
  }

  while ( cut_lines.size() )
  {
    segments.push_back( Teuchos::rcp( new LineSegment( mesh, element, this, cut_lines, false ) ) );
  }

  for ( std::vector<Teuchos::RCP<LineSegment> >::iterator i=segments.begin(); i!=segments.end(); ++i )
  {
    LineSegment & ls = **i;
    if ( ls.IsClosed() )
    {
      const std::vector<Point*> & facet_points = ls.Points();
      Facet * f = FindFacet( facet_points );
      if ( f==NULL )
      {
        // If we have a hole and multiple cuts we have to test which facet the
        // hole belongs to. Not supported now.
        if ( facets_.size()!=1 )
        {
          throw std::runtime_error( "expect side with one (uncut) facet" );
        }
        facets_[0]->AddHole( facet_points );
      }
    }
  }
}

void GEO::CUT::LinearSide::MakeInternalFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets )
{
  std::vector<Teuchos::RCP<LineSegment> > segments;
  CreateLineSegmentList( mesh, element, segments, false );

  if ( segments.size()!=1 )
  {
    throw std::runtime_error( "expect one closed cut" );
  }

  LineSegment & ls = *segments[0];

  if ( not ls.IsClosed() )
  {
    //throw std::runtime_error( "expect one closed cut" );

    // Assume this is a cut along one of our edges. So this side is not
    // responsible.
    return;
  }

  Side * s = ls.OnSide( element );
  if ( s!=NULL )
  {
    const std::vector<Point*> & facet_points = ls.Points();
    Facet * f = s->FindFacet( facet_points );
    if ( f!=NULL )
    {
      f->ExchangeSide( this );
    }
    else
    {
      throw std::runtime_error( "must have matching facet on side" );
    }
  }
  else
  {
    // insert new internal facet
    const std::vector<Point*> & facet_points = ls.Points();
    Facet * f = mesh.NewFacet( facet_points, 0, this );
    facets.insert( f );
  }
}

bool GEO::CUT::Side::OnSide( const std::set<Point*, PointPidLess> & points )
{
  if ( nodes_.size()==points.size() )
  {
    for ( std::vector<Node*>::iterator i=nodes_.begin();
          i!=nodes_.end();
          ++i )
    {
      Node * n = *i;
      if ( points.count( n->point() )==0 )
      {
        return false;
      }
    }
    return true;
  }
  return false;
}

bool GEO::CUT::Side::OnEdge( Point * point )
{
  for ( std::vector<Edge*>::const_iterator i=edges_.begin(); i!=edges_.end(); ++i )
  {
    Edge * e = *i;
    if ( point->IsCut( e ) )
    {
      return true;
    }
  }
  return false;
}

bool GEO::CUT::Side::OnEdge( Line * line )
{
  for ( std::vector<Edge*>::const_iterator i=edges_.begin(); i!=edges_.end(); ++i )
  {
    Edge * e = *i;
    if ( line->OnEdge( e ) )
    {
      return true;
    }
  }
  return false;
}

void GEO::CUT::Side::Print()
{
  std::cout << "[ ";
  for ( std::vector<Edge*>::iterator i=edges_.begin(); i!=edges_.end(); ++i )
  {
    Edge * e = *i;
    e->Print();
    std::cout << " ; ";
  }
  std::cout << " ]";
}

GEO::CUT::Node * GEO::CUT::Side::OnNode( const LINALG::Matrix<3,1> & x )
{
  LINALG::Matrix<3,1> nx;
  for ( std::vector<Node*>::iterator i=nodes_.begin(); i!=nodes_.end(); ++i )
  {
    Node * n = *i;
    n->Coordinates( nx.A() );
    nx.Update( -1, x, 1 );
    if ( nx.Norm2() < RELAXEDTOL )
    {
      return n;
    }
  }
  return NULL;
}

bool GEO::CUT::LinearSide::IsCut()
{
  if ( facets_.size()>1 )
    return true;
  if ( facets_[0]->SideId() > -1 )
    return true;
  return false;
}

void GEO::CUT::LinearSide::ExchangeFacetSide( Side * side )
{
  for ( std::vector<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    f->ExchangeSide( side );
  }
}


void GEO::CUT::ConcreteSide<DRT::Element::tri6>::FillComplete( Mesh & mesh )
{
  if ( subsides_.size()==0 )
  {
    subsides_.reserve( 4 );

    Side::FillComplete( mesh );

    const CellTopologyData * top_data = shards::getCellTopologyData< shards::Triangle<3> >();
    const std::vector<Node*> & nodes = Nodes();

    std::vector<int> nids( 3 );

    nids[0] = nodes[0]->Id();
    nids[1] = nodes[3]->Id();
    nids[2] = nodes[5]->Id();
    subsides_.push_back( mesh.GetSide( Id(), nids, top_data ) );

    nids[0] = nodes[3]->Id();
    nids[1] = nodes[1]->Id();
    nids[2] = nodes[4]->Id();
    subsides_.push_back( mesh.GetSide( Id(), nids, top_data ) );

    nids[0] = nodes[3]->Id();
    nids[1] = nodes[4]->Id();
    nids[2] = nodes[5]->Id();
    subsides_.push_back( mesh.GetSide( Id(), nids, top_data ) );

    nids[0] = nodes[5]->Id();
    nids[1] = nodes[4]->Id();
    nids[2] = nodes[2]->Id();
    subsides_.push_back( mesh.GetSide( Id(), nids, top_data ) );
  }
}

void GEO::CUT::ConcreteSide<DRT::Element::quad8>::FillComplete( Mesh & mesh )
{
  if ( subsides_.size()==0 )
  {
    subsides_.reserve( 4 );

    Side::FillComplete( mesh );

    const CellTopologyData * top_data = shards::getCellTopologyData< shards::Quadrilateral<4> >();
    const std::vector<Node*> & nodes = Nodes();

    // create middle node

    LINALG::Matrix<3,8> xyze;
    Coordinates( xyze );

    LINALG::Matrix<8,1> funct;
    DRT::UTILS::shape_function_2D( funct, 0, 0, DRT::Element::quad8 );

    LINALG::Matrix<3,1> xyz;
    xyz.Multiply( xyze, funct );

    std::set<int> node_nids;
    for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
    {
      Node * n = *i;
      node_nids.insert( n->Id() );
    }
    Node* middle = mesh.GetNode( node_nids, xyz.A() );
    int middle_id = middle->Id();

    std::vector<int> nids( 4 );

    nids[0] = nodes[0]->Id();
    nids[1] = nodes[4]->Id();
    nids[2] = middle_id;
    nids[3] = nodes[7]->Id();
    subsides_.push_back( mesh.GetSide( Id(), nids, top_data ) );

    nids[0] = nodes[4]->Id();
    nids[1] = nodes[1]->Id();
    nids[2] = nodes[5]->Id();
    nids[3] = middle_id;
    subsides_.push_back( mesh.GetSide( Id(), nids, top_data ) );

    nids[0] = middle_id;
    nids[1] = nodes[5]->Id();
    nids[2] = nodes[2]->Id();
    nids[3] = nodes[6]->Id();
    subsides_.push_back( mesh.GetSide( Id(), nids, top_data ) );

    nids[0] = nodes[7]->Id();
    nids[1] = middle_id;
    nids[2] = nodes[6]->Id();
    nids[3] = nodes[3]->Id();
    subsides_.push_back( mesh.GetSide( Id(), nids, top_data ) );
  }
}

void GEO::CUT::ConcreteSide<DRT::Element::quad9>::FillComplete( Mesh & mesh )
{
  if ( subsides_.size()==0 )
  {
    subsides_.reserve( 4 );

    Side::FillComplete( mesh );

    const CellTopologyData * top_data = shards::getCellTopologyData< shards::Quadrilateral<4> >();
    const std::vector<Node*> & nodes = Nodes();

    std::vector<int> nids( 4 );

    nids[0] = nodes[0]->Id();
    nids[1] = nodes[4]->Id();
    nids[2] = nodes[8]->Id();
    nids[3] = nodes[7]->Id();
    subsides_.push_back( mesh.GetSide( Id(), nids, top_data ) );

    nids[0] = nodes[4]->Id();
    nids[1] = nodes[1]->Id();
    nids[2] = nodes[5]->Id();
    nids[3] = nodes[8]->Id();
    subsides_.push_back( mesh.GetSide( Id(), nids, top_data ) );

    nids[0] = nodes[8]->Id();
    nids[1] = nodes[5]->Id();
    nids[2] = nodes[2]->Id();
    nids[3] = nodes[6]->Id();
    subsides_.push_back( mesh.GetSide( Id(), nids, top_data ) );

    nids[0] = nodes[7]->Id();
    nids[1] = nodes[8]->Id();
    nids[2] = nodes[6]->Id();
    nids[3] = nodes[3]->Id();
    subsides_.push_back( mesh.GetSide( Id(), nids, top_data ) );
  }
}


bool GEO::CUT::QuadraticSide::Cut( Mesh & mesh, Side & side, Element * element )
{
  bool success = false;
  for ( std::vector<Side*>::iterator i=subsides_.begin(); i!=subsides_.end(); ++i )
  {
    Side * s = *i;
    if ( s->Cut( mesh, side, element ) )
      success = true;
  }
  return success;
}

void GEO::CUT::QuadraticSide::MakeOwnedSideFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets )
{
  throw std::runtime_error( "not supposed to end up here" );
//   for ( std::vector<Side*>::iterator i=subsides_.begin(); i!=subsides_.end(); ++i )
//   {
//     Side * s = *i;
//     s->MakeOwnedSideFacets( mesh, element, facets );
//   }
}

void GEO::CUT::QuadraticSide::MakeSideCutFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets )
{
  throw std::runtime_error( "not supposed to end up here" );
}

void GEO::CUT::QuadraticSide::MakeInternalFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets )
{
  throw std::runtime_error( "not supposed to end up here" );
//   for ( std::vector<Side*>::iterator i=subsides_.begin(); i!=subsides_.end(); ++i )
//   {
//     Side * s = *i;
//     s->MakeInternalFacets( mesh, element, facets );
//   }
}

bool GEO::CUT::QuadraticSide::IsCut()
{
  throw std::runtime_error( "not supposed to end up here" );
}

void GEO::CUT::ConcreteSide<DRT::Element::tri3>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position2d<DRT::Element::tri3> pos( *this, xyz );
//   if ( not pos.Compute() )
//   {
//     throw std::runtime_error( "global point not within element" );
//   }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::ConcreteSide<DRT::Element::tri6>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position2d<DRT::Element::tri6> pos( *this, xyz );
//   if ( not pos.Compute() )
//   {
//     throw std::runtime_error( "global point not within element" );
//   }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::ConcreteSide<DRT::Element::quad4>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position2d<DRT::Element::quad4> pos( *this, xyz );
//   if ( not pos.Compute() )
//   {
//     throw std::runtime_error( "global point not within element" );
//   }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::ConcreteSide<DRT::Element::quad8>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position2d<DRT::Element::quad8> pos( *this, xyz );
//   if ( not pos.Compute() )
//   {
//     throw std::runtime_error( "global point not within element" );
//   }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::ConcreteSide<DRT::Element::quad9>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position2d<DRT::Element::quad9> pos( *this, xyz );
//   if ( not pos.Compute() )
//   {
//     throw std::runtime_error( "global point not within element" );
//   }
  rst = pos.LocalCoordinates();
}
