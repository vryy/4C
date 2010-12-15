
//#include "../drt_geometry/intersection_templates.H"

#include "cut_tetgen.H"
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

bool GEO::CUT::LinearSide::FindCutPoints( Mesh & mesh, LinearElement * element, LinearSide & other )
{
  bool cut = false;
  const std::vector<Edge*> & edges = Edges();
  for ( std::vector<Edge*>::const_iterator i=edges.begin(); i!=edges.end(); ++i )
  {
    ConcreteEdge<DRT::Element::line2> * e = dynamic_cast<ConcreteEdge<DRT::Element::line2>*>( *i );
    if ( e->FindCutPoints( mesh, element, *this, other ) )
    {
      cut = true;
    }
  }
  return cut;
}

bool GEO::CUT::LinearSide::FindCutLines( Mesh & mesh, LinearElement * element, LinearSide & other )
{
#if 0
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

bool GEO::CUT::LinearSide::AllOnNodes( const std::set<Point*> & points )
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

void GEO::CUT::LinearSide::GetCutPoints( Element * element, Side & other, std::set<Point*> & cuts )
{
  const std::vector<Edge*> & edges = Edges();
  for ( std::vector<Edge*>::const_iterator i=edges.begin(); i!=edges.end(); ++i )
  {
    ConcreteEdge<DRT::Element::line2> * e = dynamic_cast<ConcreteEdge<DRT::Element::line2>*>( *i );
    e->GetCutPoints( element, *this, other, cuts );
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


void GEO::CUT::LinearSide::GetBoundaryCells( std::set<GEO::CUT::BoundaryCell*> & bcells )
{
  for ( std::vector<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    f->GetBoundaryCells( bcells );
  }
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
    LineSegment * ls = new LineSegment( mesh, element, this, cut_lines, inner );
    segments.push_back( Teuchos::rcp( ls ) );
  }
}

void GEO::CUT::LinearSide::CreateLineSegment( Mesh & mesh, Element * element )
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

void GEO::CUT::LinearSide::MakeOwnedSideFacets( Mesh & mesh, const PointLineFilter & filter, std::set<Facet*> & facets )
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

bool GEO::CUT::LinearSide::IsCut()
{
  if ( facets_.size()>1 )
    return true;
  if ( facets_[0]->SideId() > -1 )
    return true;
  return false;
}

// void GEO::CUT::LinearSide::ExchangeFacetSide( Side * side, bool cutsurface )
// {
//   for ( std::vector<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
//   {
//     Facet * f = *i;
//     f->ExchangeSide( side, cutsurface );
//   }
// }


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


void GEO::CUT::QuadraticSide::MakeOwnedSideFacets( Mesh & mesh, const PointLineFilter & filter, std::set<Facet*> & facets )
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

void GEO::CUT::QuadraticSide::GetBoundaryCells( std::set<GEO::CUT::BoundaryCell*> & bcells )
{
  for ( std::vector<Side*>::iterator i=subsides_.begin(); i!=subsides_.end(); ++i )
  {
    Side * s = *i;
    s->GetBoundaryCells( bcells );
  }
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

bool GEO::CUT::ConcreteSide<DRT::Element::tri6>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position2d<DRT::Element::tri6> pos( *this, xyz );
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

bool GEO::CUT::ConcreteSide<DRT::Element::quad8>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position2d<DRT::Element::quad8> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
  return success;
}

bool GEO::CUT::ConcreteSide<DRT::Element::quad9>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position2d<DRT::Element::quad9> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
  return success;
}
