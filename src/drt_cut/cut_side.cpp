
//#include "../drt_geometry/intersection_templates.H"

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

bool GEO::CUT::Side::FindCutPoints( Mesh & mesh, Element * element, Side & other )
{
  bool cut = false;
  const std::vector<Edge*> & edges = Edges();
  for ( std::vector<Edge*>::const_iterator i=edges.begin(); i!=edges.end(); ++i )
  {
    Edge * e = *i;
    if ( e->FindCutPoints( mesh, element, *this, other ) )
    {
      cut = true;
    }
  }
  return cut;
}

bool GEO::CUT::Side::FindCutLines( Mesh & mesh, Element * element, Side & other )
{
#if 1
  bool cut = false;
  for ( std::vector<Line*>::iterator i=cut_lines_.begin(); i!=cut_lines_.end(); ++i )
  {
    Line * l = *i;
    if ( l->IsCut( this, &other ) )
    {
      l->AddElement( element );
      other.AddLine( l );
      cut = true;
    }
  }
  if ( cut )
  {
    return true;
  }
#endif

  std::set<Point*> cuts;
  GetCutPoints( element, other, cuts );

  switch ( cuts.size() )
  {
  case 0:
    return false;
  case 1:
  {
    std::set<Point*> reverse_cuts;
    other.GetCutPoints( element, *this, reverse_cuts );
    reverse_cuts.erase( *cuts.begin() );
    if ( reverse_cuts.size()==1 )
    {
      //Line * l =
      mesh.NewLine( *cuts.begin(), *reverse_cuts.begin(), this, &other, element );
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
  case 2:
  {
    // The normal case. A straight cut.
    std::vector<Point*> c( cuts.begin(), cuts.end() );
    //Line * l =
    mesh.NewLine( c[0], c[1], this, &other, element );
    return true;
  }
  default:
  {
    // More that two cut points shows a touch.
    //
    // If all nodes are catched and nothing else, the cut surface has hit this
    // surface exactly. No need to cut anything. However, the surface might be
    // required for integration.

    const std::vector<Node*> & nodes = Nodes();
    if ( cuts.size()==nodes.size() and AllOnNodes( cuts ) )
    {
      for ( unsigned i=0; i<nodes.size(); ++i )
      {
        unsigned j = ( i+1 ) % nodes.size();
        //Line* l =
        mesh.NewLine( nodes[i]->point(), nodes[j]->point(), this, &other, element );
      }
      return true;
    }
    return false;
  }
  }
}

bool GEO::CUT::Side::AllOnNodes( const std::set<Point*> & points )
{
  const std::vector<Node*> & nodes = Nodes();
  for ( std::set<Point*>::const_iterator i=points.begin(); i!=points.end(); ++i )
  {
    Point * p = *i;
    if ( not p->NodalPoint( nodes ) )
    {
      return false;
    }
  }
  return true;
}

void GEO::CUT::Side::GetCutPoints( Element * element, Side & other, std::set<Point*> & cuts )
{
  const std::vector<Edge*> & edges = Edges();
  for ( std::vector<Edge*>::const_iterator i=edges.begin(); i!=edges.end(); ++i )
  {
    Edge * e = *i;
    e->GetCutPoints( element, *this, other, cuts );
  }
}

void GEO::CUT::Side::AddLine( Line* cut_line )
{
  if ( std::find( cut_lines_.begin(), cut_lines_.end(), cut_line )==cut_lines_.end() )
  {
    cut_lines_.push_back( cut_line );
  }
}

GEO::CUT::Facet * GEO::CUT::Side::FindFacet( const std::vector<Point*> & facet_points )
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


void GEO::CUT::Side::GetBoundaryCells( std::set<GEO::CUT::BoundaryCell*> & bcells )
{
  for ( std::vector<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    f->GetBoundaryCells( bcells );
  }
}

void GEO::CUT::Side::CreateLineSegmentList( Mesh & mesh,
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
    LineSegment * ls = new LineSegment( mesh, element, this, cut_lines, inner );
    segments.push_back( Teuchos::rcp( ls ) );
  }
}

void GEO::CUT::Side::CreateLineSegment( Mesh & mesh, Element * element )
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

#if 0
  std::vector<Teuchos::RCP<LineSegment> >::iterator i
    = std::remove_if( segments.begin(), segments.end(),
                      std::bind2nd( std::equal_to<Teuchos::RCP<LineSegment> >(), Teuchos::null ) );
  segments.erase( i, segments.end() );
#endif
}

void GEO::CUT::Side::MakeOwnedSideFacets( Mesh & mesh, const PointLineFilter & filter, std::set<Facet*> & facets )
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
      PointCycleList pcl( filter, this, facet_points, cut_points );
      pcl.CreateFacets( mesh, this, facets_ );
    }
    else
    {
      // Just a normal side. There might be cut points at the edges that do
      // not concern this side. We have to copy them anyway.

      facets_.push_back( mesh.NewFacet( facet_points, this, false ) );
    }
  }

  std::copy( facets_.begin(), facets_.end(), std::inserter( facets, facets.begin() ) );
}

void GEO::CUT::Side::MakeSideCutFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets )
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
        Facet * hole = mesh.NewFacet( facet_points, this, false );
        facets_[0]->AddHole( hole );
      }
    }
  }
}

void GEO::CUT::Side::MakeInternalFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets )
{
  std::vector<Teuchos::RCP<LineSegment> > segments;
  CreateLineSegmentList( mesh, element, segments, false );

  for ( unsigned i=0; i<segments.size(); ++i )
  {
    LineSegment & ls = *segments[i];

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
        f->ExchangeSide( this, true );
      }
      else
      {
        //throw std::runtime_error( "must have matching facet on side" );

        // multiple facets on one cut side within one elemenet: this is a
        // levelset case
        Facet * f = mesh.NewFacet( facet_points, this, true );
        facets.insert( f );
      }
    }
    else
    {
      // insert new internal facet
      const std::vector<Point*> & facet_points = ls.Points();
      Facet * f = mesh.NewFacet( facet_points, this, true );
      facets.insert( f );
    }
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

bool GEO::CUT::Side::HaveCommonEdge( Side & side )
{
  const std::vector<Edge*> & other_edges = side.Edges();
  for ( std::vector<Edge*>::const_iterator i=edges_.begin(); i!=edges_.end(); ++i )
  {
    Edge * e = *i;
    if ( std::find( other_edges.begin(), other_edges.end(), e )!=other_edges.end() )
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
    if ( nx.Norm2() < MINIMALTOL )
    {
      return n;
    }
  }
  return NULL;
}

bool GEO::CUT::Side::IsCut()
{
  if ( facets_.size()>1 )
    return true;
  if ( facets_[0]->SideId() > -1 )
    return true;
  return false;
}

bool GEO::CUT::ConcreteSide<DRT::Element::tri3>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position2d<DRT::Element::tri3> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
  return success;
}

bool GEO::CUT::ConcreteSide<DRT::Element::quad4>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position2d<DRT::Element::quad4> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
  return success;
}
