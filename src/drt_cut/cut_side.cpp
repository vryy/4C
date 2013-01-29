
#include "cut_position.H"
#include "cut_position2d.H"
#include "cut_intersection.H"
#include "cut_facet.H"
#include "cut_point_impl.H"
#include "cut_pointgraph.H"
#include "cut_creator.H"

#include <string>
#include <stack>


/*-----------------------------------------------------------------------------------------*
      Returns the edge of this side with given begin and end points
 *-----------------------------------------------------------------------------------------*/
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
/*-----------------------------------------------------------------------------------------*
    Calculate the points at which the other side intersects with this considered side
 *-----------------------------------------------------------------------------------------*/
bool GEO::CUT::Side::FindCutPoints( Mesh & mesh, Element * element, Side & other, int recursion )
{
  bool cut = false;
  const std::vector<Edge*> & edges = Edges();
  for ( std::vector<Edge*>::const_iterator i=edges.begin(); i!=edges.end(); ++i )
  {
    Edge * e = *i;
    if ( e->FindCutPoints( mesh, element, *this, other, recursion ) )
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
  // check whether cut lines are already created
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

  // creating cut lines between the given sides for the first time
  PointSet cuts;
  GetCutPoints( element, other, cuts );

  switch ( cuts.size() )
  {
  case 0:
    return false;
  case 1:
  {
    PointSet reverse_cuts;
    other.GetCutPoints( element, *this, reverse_cuts );
    reverse_cuts.erase( *cuts.begin() );
    if ( reverse_cuts.size()==1 )
    {
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
      //throw std::runtime_error( "most peculiar cut" );

      // Assume there is one point too much between the two that
      // matter. Assume all three points are on one line. This is a numerical
      // artefact. But who am I to complain?
      //
      // Try to recover. All points need to be used.

      LINALG::Matrix<3,1> p1( ( *cuts.begin() )->X() );

      std::vector<Point*> rp;
      rp.reserve( 2 );
      rp.assign( reverse_cuts.begin(), reverse_cuts.end() );

      LINALG::Matrix<3,1> p2( rp[0]->X() );
      LINALG::Matrix<3,1> p3( rp[1]->X() );

      LINALG::Matrix<3,1> d;

      d.Update( 1., p2, -1., p1, 0. );
      double dist1 = d.Norm2();

      d.Update( 1., p3, -1., p1, 0. );
      double dist2 = d.Norm2();

      d.Update( 1., p3, -1., p2, 0. );
      double dist3 = d.Norm2();

      if ( dist3 > dist1 and dist3 > dist2 )
      {
        // first point in the middle
        // now this is unexpected
        mesh.NewLine( *cuts.begin(), rp[0], this, &other, element );
        mesh.NewLine( *cuts.begin(), rp[1], this, &other, element );
      }
      else
      {
        if ( dist2 > dist1 )
        {
          mesh.NewLine( *cuts.begin(), rp[1], this, &other, element );
        }
        else
        {
          mesh.NewLine( *cuts.begin(), rp[0], this, &other, element );
        }
        mesh.NewLine( rp[0], rp[1], this, &other, element );
      }
      return true;
    }
  }
  case 2:
  {
    // The normal case. A straight cut.
    std::vector<Point*> c;
    c.reserve( 2 );
    c.assign( cuts.begin(), cuts.end() );
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
        mesh.NewLine( nodes[i]->point(), nodes[j]->point(), this, &other, element );
      }
      return true;
    }
    return other.FindAmbiguousCutLines( mesh, element, *this, cuts );
  }
  }
}

void GEO::CUT::Side::CreateMissingLines( Creator & creator, Element * element )
{
  std::map<Point*, PointSet> pg;

  const std::vector<Line*> & cut_lines = CutLines();
  for ( std::vector<Line*>::const_iterator i=cut_lines.begin(); i!=cut_lines.end(); ++i )
  {
    Line * l = *i;
    if ( l->IsCut( element ) )
    {
      Point * p1 = l->BeginPoint();
      Point * p2 = l->EndPoint();
      pg[p1].insert( p2 );
      pg[p2].insert( p1 );
    }
  }

  if ( pg.size() > 2 )
  {
    // Needs to be a proper cycle. No gaps, no forks.

    std::vector<Point*> open;

    for ( std::map<Point*, PointSet>::iterator i=pg.begin();
          i!=pg.end();
          ++i )
    {
      Point * p = i->first;
      PointSet & row = i->second;
      if ( row.size() < 2 )
      {
        open.push_back( p );
      }
#if 0
      else if ( row.size() > 2 )
      {
        std::stringstream str;
        str << "Fork in line cycle. Maybe one intersection point was not found:\n"
            << "row.size()=" << row.size() << "\n"
          ;
        for ( std::map<Point*, PointSet>::iterator i=pg.begin();
              i!=pg.end();
              ++i )
        {
          str << "  " << ( *i->first ) << " : ";
          std::copy( i->second.begin(), i->second.end(), std::ostream_iterator<Point*>( str, " " ) );
          str << "\n";
        }
        throw std::runtime_error( str.str() );
      }
#endif
    }

#if 0
    if ( open.size() == 2 )
    {
      // handle simple triangle

      Point * p1 = open[0];
      Point * p2 = open[1];

      PointSet & ps1 = pg[p1];
      PointSet & ps2 = pg[p2];

      if ( ps1.size()==1 and ps1.count( p2 )==0 and
           ps2.size()==1 and ps2.count( p1 )==0 )
      {
        creator.NewLine( p1, p2, this, NULL, element );
        open.clear();
      }
    }
#endif

    if ( open.size() > 0 )
    {
      PointSet done;

      std::vector<Point*> open_side_points;
      open_side_points.reserve( 2 );

      const std::vector<Side*> & sides = element->Sides();
      for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
      {
        Side * s = *i;
        open_side_points.clear();

        for ( std::vector<Point*>::iterator i=open.begin(); i!=open.end(); ++i )
        {
          Point * p = *i;
          if ( p->IsCut( s ) )
          {
            open_side_points.push_back( p );
          }
        }

        if ( open_side_points.size()==2 )
        {
          creator.NewLine( open_side_points[0], open_side_points[1], s, this, element );
          done.insert( open_side_points[0] );
          done.insert( open_side_points[1] );
        }
#if 0
        else if ( open_side_points.size() > 0 )
        {
          throw std::runtime_error( "illegal number of open points on element side" );
        }
#endif
      }

      if ( done.size() != open.size() )
      {
#ifdef DEBUGCUTLIBRARY
        element->GnuplotDump();
        element->DumpFacets();

#if 1
  {
    // see if there are any unconnected cut points on any element side

    const std::vector<Side*> & sides = element->Sides();
    for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
    {
      Side * s = *i;
      const PointSet & cutpoints = s->CutPoints();
      for ( PointSet::const_iterator i=cutpoints.begin(); i!=cutpoints.end(); ++i )
      {
        Point * p = *i;
        if ( p->IsCut( this ) )
        {
          pg[p];
        }
      }
    }
  }
#endif

        std::vector<Point*> dummy;
        for ( std::map<Point*, PointSet>::iterator i=pg.begin(); i!=pg.end(); ++i )
        {
          dummy.push_back( i->first );
        }
#endif

        std::stringstream str;
        str << "failed to close open points\npg:\n";
        for ( std::map<Point*, PointSet>::iterator i=pg.begin(); i!=pg.end(); ++i )
        {
          Point * p = i->first;
          PointSet & ps = i->second;
          str << ( *p ) << " : ";
          std::copy( ps.begin(), ps.end(), std::ostream_iterator<Point*>( str, " " ) );
          str << "\n";
        }
        str << "done: ";
        std::copy( done.begin(), done.end(), std::ostream_iterator<Point*>( str, " " ) );
        str << "\n";
        str << "open: ";
        std::copy( open.begin(), open.end(), std::ostream_iterator<Point*>( str, " " ) );
        str << "\n";
        str << "open_side_points: ";
        std::copy( open_side_points.begin(), open_side_points.end(), std::ostream_iterator<Point*>( str, " " ) );
        str << "\n";

        throw std::runtime_error( str.str() );
      }
    }
  }
}

bool GEO::CUT::Side::AllOnNodes( const PointSet & points )
{
  const std::vector<Node*> & nodes = Nodes();
  for ( PointSet::const_iterator i=points.begin(); i!=points.end(); ++i )
  {
    Point * p = *i;
    if ( not p->NodalPoint( nodes ) )
    {
      return false;
    }
  }
  return true;
}

void GEO::CUT::Side::GetCutPoints( Element * element, Side & other, PointSet & cuts )
{
  const std::vector<Edge*> & edges = Edges();
  for ( std::vector<Edge*>::const_iterator i=edges.begin(); i!=edges.end(); ++i )
  {
    Edge * e = *i;
    e->GetCutPoints( element, *this, other, cuts );
  }
}

void GEO::CUT::Side::AddPoint( Point * cut_point )
{
  cut_points_.insert( cut_point );
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

bool GEO::CUT::Side::FindAmbiguousCutLines( Mesh & mesh, Element * element, Side & side, const PointSet & cut )
{
  return false;
}

void GEO::CUT::Side::GetBoundaryCells( plain_boundarycell_set & bcells )
{
  for ( std::vector<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    f->GetBoundaryCells( bcells );
  }
}

/*-----------------------------------------------------------------------------------------------*
                create facets on the background sides of the element
                For all these facets, parent side is an element side
 *-----------------------------------------------------------------------------------------------*/
void GEO::CUT::Side::MakeOwnedSideFacets( Mesh & mesh, Element * element, plain_facet_set & facets )
{
  if ( facets_.size()==0 )
  {
    IMPL::PointGraph pg( mesh, element, this, IMPL::PointGraph::element_side, IMPL::PointGraph::all_lines );

    for ( IMPL::PointGraph::facet_iterator i=pg.fbegin(); i!=pg.fend(); ++i )
    {
      const Cycle & points = *i;

      Facet * f = mesh.NewFacet( points(), this, IsCutSide() );
      if ( f==NULL )
        throw std::runtime_error( "failed to create facet" );
      facets_.push_back( f );
    }

    for ( IMPL::PointGraph::hole_iterator i=pg.hbegin(); i!=pg.hend(); ++i )
    {
      const std::vector<Cycle> & hole = *i;

      // If we have a hole and multiple cuts we have to test which facet the
      // hole belongs to. Not supported now.
      if ( facets_.size()!=1 )
      {
        throw std::runtime_error( "expect side with one (uncut) facet" );
      }

      for ( std::vector<Cycle>::const_iterator i=hole.begin(); i!=hole.end(); ++i )
      {
        const Cycle & points = *i;

        Facet * h = mesh.NewFacet( points(), this, false );
        facets_[0]->AddHole( h );
      }
    }
  }

  std::copy( facets_.begin(), facets_.end(), std::inserter( facets, facets.begin() ) );
}

/*-----------------------------------------------------------------------------------------------*
                     create facets on the cut sides of the element
                     For all these facets, parent side is a cut side
 *-----------------------------------------------------------------------------------------------*/
void GEO::CUT::Side::MakeInternalFacets( Mesh & mesh, Element * element, plain_facet_set & facets )
{
  IMPL::PointGraph pg( mesh, element, this, IMPL::PointGraph::cut_side, IMPL::PointGraph::all_lines );
  for ( IMPL::PointGraph::facet_iterator i=pg.fbegin(); i!=pg.fend(); ++i )
  {
    const Cycle & points = *i;
    MakeInternalFacets( mesh, element, points, facets );
  }
  for ( IMPL::PointGraph::hole_iterator i=pg.hbegin(); i!=pg.hend(); ++i )
  {
    const std::vector<Cycle> & hole = *i;
    for ( std::vector<Cycle>::const_iterator i=hole.begin(); i!=hole.end(); ++i )
    {
      const Cycle & points = *i;
      MakeInternalFacets( mesh, element, points, facets );
    }
  }
}

void GEO::CUT::Side::MakeInternalFacets( Mesh & mesh, Element * element, const Cycle & points, plain_facet_set & facets )
{
  // ignore cycles with all points on one and the same edge
  // ignore cycles with points outside the current element
  if ( not points.IsValid() or not points.IsCut( element ) )
    return;

  Side * s = NULL;

  plain_side_set sides( element->Sides().begin(), element->Sides().end() );
  points.Intersection( sides );

  if ( sides.size()>1 )
  {
    std::stringstream str;
    str << "can touch exactly one element side: "
        << points
        << "found sides:\n";
    std::copy( sides.begin(), sides.end(), std::ostream_iterator<Side*>( str, "\n" ) );
    throw std::runtime_error( str.str() );
  }
  else if ( sides.size()==1 )
  {
    s = *sides.begin();
  }

  if ( s!=NULL )
  {
    Facet * f = s->FindFacet( points() );
    if ( f!=NULL )
    {
      f->ExchangeSide( this, true );
      facets.insert( f );
      facets_.push_back( f );
    }
    else
    {
      //throw std::runtime_error( "must have matching facet on side" );

      // multiple facets on one cut side within one element
      Facet * f = mesh.NewFacet( points(), this, true );
      facets.insert( f );
      facets_.push_back( f );
    }
  }
  else
  {
    // insert new internal facet
    Facet * f = mesh.NewFacet( points(), this, true );
    facets.insert( f );
    facets_.push_back( f );
  }
}

bool GEO::CUT::Side::OnSide( const PointSet & points )
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

bool GEO::CUT::Side::HaveCommonNode( Side & side )
{
  const std::vector<Node*> & other_nodes = side.Nodes();
  for ( std::vector<Node*>::const_iterator i=nodes_.begin(); i!=nodes_.end(); ++i )
  {
    Node * e = *i;
    if ( std::find( other_nodes.begin(), other_nodes.end(), e )!=other_nodes.end() )
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

GEO::CUT::Element * GEO::CUT::Side::CommonElement( Side * other )
{
  plain_element_set intersection;
  std::set_intersection( elements_.begin(), elements_.end(),
                         other->elements_.begin(), other->elements_.end(),
                         std::inserter( intersection, intersection.begin() ) );
  switch ( intersection.size() )
  {
  case 0:
    return NULL;
  case 1:
    return *intersection.begin();
  default:
    throw std::runtime_error( "sides with more than one element in common" );
  }
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
  if ( facets_[0]->OnCutSide() )
    return true;
  return false;
}

/*--------------------------------------------------------------------*
 * is this side closer to the startpoint than the other side?
 * check based on ray-tracing technique
 * set is_closer
 * return if check was successful
 *--------------------------------------------------------------------*/
bool GEO::CUT::Side::IsCloserSide( LINALG::Matrix<3,1>& startpoint_xyz, GEO::CUT::Side* other, bool& is_closer)
{
  // shoot a ray starting from the startpoint through the midpoint of this side
  // and find an intersection point with the other side
  LINALG::Matrix<3,1> ray_point_xyz(true);
  this->SideCenter(ray_point_xyz); // as second point on the ray we define the midpoint of this side

  LINALG::Matrix<2,1> rs(true);
  double line_xi = 0.0;

  // shoot the ray and find a cutpoint with the other side's plane or curved surface space
  bool cut_found = other->RayCut( startpoint_xyz, ray_point_xyz, rs, line_xi);

  if(!cut_found) return false;
  else
  {
    // The main decision if the side lies closer to the start-point than the other side

    //cout << "line_xi " << line_xi << endl;

    if(line_xi > 1.0+TOLERANCE)
    {
      // the first side is closer to the start point than the second side
      is_closer = true;
      return true;
    }
    else if(fabs(line_xi - 1.0) <= TOLERANCE)
    {
      // the found intersection point on the other is the midpoint of the fist side
      // in that case both sides lie within one plane
      // this case is catched in SameNormal afterwards

      // std::cout << "line_xi " << line_xi << endl;
      // std::cout << "check if both sides lie in one plane " << endl;

      return false;
    }
    else if(line_xi < 1.0-TOLERANCE and line_xi > -1.0+TOLERANCE)
    {
      // the most safe check (accept the other side as the nearest side)
      is_closer = false;
      return true;
    }
    else if(fabs(line_xi + 1.0) <= TOLERANCE )
    {
      // the intersection point is the same as the start-point of the ray
      // the other side contains the start-point and the cut-point shared with the original side
      // in that case the line between the start-point and the cut-point lies in the second side
      // then the side is orthogonal to the line
      // this side should be removed in
      std::cout << "line_xi " << line_xi << endl;
      std::cout << "start-point: " << startpoint_xyz << endl;
      std::cout << "side orthogonal ? " << endl; other->Print();

      throw std::runtime_error("IsCloserSide along the ray-tracing line failed! ");

      return false;
    }
    else if(line_xi < -1.0-TOLERANCE)
    {
      // both sides lead to the same result
      is_closer = true; // false would be also okay
      return true;
    }
    else
    {
      // undermined range of local coordinates!

      std::cout << "line_xi " << line_xi << endl;
      std::cout << "cut point found, but the local line coordinates along the ray-tracing line lies in undefined region" << endl;

      throw std::runtime_error("IsCloserSide along the ray-tracing line failed! ");
    }

    return false;
  }

  return false; // return not successful
}


/*--------------------------------------------------------------------*
 * get the global coordinates on side at given local coordinates
 *--------------------------------------------------------------------*/
void GEO::CUT::ConcreteSide<DRT::Element::tri3>::PointAt( double r, double s, LINALG::Matrix<3,1> & xyz)
{
  LINALG::Matrix<3,1> funct(true);
  LINALG::Matrix<3,3> xyz_surface(true);
  this->Coordinates(xyz_surface);

  DRT::UTILS::shape_function_2D( funct, r, s, DRT::Element::tri3 );
  xyz.Multiply(xyz_surface,funct);
}


/*--------------------------------------------------------------------*
 * get the global coordinates on side at given local coordinates
 *--------------------------------------------------------------------*/
void GEO::CUT::ConcreteSide<DRT::Element::quad4>::PointAt( double r, double s, LINALG::Matrix<3,1> & xyz)
{
  LINALG::Matrix<4,1> funct(true);
  LINALG::Matrix<3,4> xyz_surface(true);
  this->Coordinates(xyz_surface);

  DRT::UTILS::shape_function_2D( funct, r, s, DRT::Element::quad4 );
  xyz.Multiply(xyz_surface,funct);
}


/*--------------------------------------------------------------------*
 * get global coordinates of the center of the side
 *--------------------------------------------------------------------*/
void GEO::CUT::ConcreteSide<DRT::Element::tri3>::SideCenter( LINALG::Matrix<3,1> & midpoint )
{
  LINALG::Matrix<2,1> center_rs(DRT::UTILS::getLocalCenterPosition<2>(DRT::Element::tri3));
  PointAt(center_rs(0), center_rs(1), midpoint);
}


/*--------------------------------------------------------------------*
 * get global coordinates of the center of the side
 *--------------------------------------------------------------------*/
void GEO::CUT::ConcreteSide<DRT::Element::quad4>::SideCenter( LINALG::Matrix<3,1> & midpoint )
{
  LINALG::Matrix<2,1> center_rs(DRT::UTILS::getLocalCenterPosition<2>(DRT::Element::quad4));
  PointAt(center_rs(0), center_rs(1), midpoint);
}


/*--------------------------------------------------------------------*
 * lies point with given coordinates within this side?
 *--------------------------------------------------------------------*/
bool GEO::CUT::ConcreteSide<DRT::Element::tri3>::WithinSide( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<2,1> & rs, double & dist, const double & Tol)
{
  Position2d<DRT::Element::tri3> pos( *this, xyz );
  bool success = pos.ComputeDistance();
  if ( not success )
  {
    throw std::runtime_error( "ComputeDistance w.r.t tri3 side not successful" );
  }
  LINALG::Matrix<3,1> rst = pos.LocalCoordinates();

  rs(0)= rst(0);
  rs(1)= rst(1);
  dist = rst(2);

  if(pos.WithinLimitsTol(Tol, false, TOLERANCE))
  {
    return true;
  }

  return false;
}


/*--------------------------------------------------------------------*
 * lies point with given coordinates within this side?
 *--------------------------------------------------------------------*/
bool GEO::CUT::ConcreteSide<DRT::Element::quad4>::WithinSide( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<2,1> & rs, double & dist, const double & Tol)
{
  Position2d<DRT::Element::quad4> pos( *this, xyz );
  bool success = pos.ComputeDistance();
  if ( not success )
  {
    throw std::runtime_error( "ComputeDistance w.r.t quad4 side not successful" );
  }
  LINALG::Matrix<3,1> rst = pos.LocalCoordinates();

  rs(0)= rst(0);
  rs(1)= rst(1);
  dist = rst(2);

  if(pos.WithinLimitsTol(Tol, false, TOLERANCE))
  {
    return true;
  }

  return false;
}


/*--------------------------------------------------------------------*
 * compute the cut of a ray through two points with the 2D space defined by the side
 *--------------------------------------------------------------------*/
bool GEO::CUT::ConcreteSide<DRT::Element::tri3>::RayCut( const LINALG::Matrix<3,1> & p1_xyz, const LINALG::Matrix<3,1> & p2_xyz, LINALG::Matrix<2,1> & rs, double & line_xi)
{

  LINALG::Matrix<3,3> xyze_surfaceElement(true);
  this->Coordinates(xyze_surfaceElement);

  LINALG::Matrix<3,2> xyze_lineElement(true);
  xyze_lineElement(0,0) = p1_xyz(0);
  xyze_lineElement(1,0) = p1_xyz(1);
  xyze_lineElement(2,0) = p1_xyz(2);

  xyze_lineElement(0,1) = p2_xyz(0);
  xyze_lineElement(1,1) = p2_xyz(1);
  xyze_lineElement(2,1) = p2_xyz(2);

  LINALG::Matrix<3,1> xsi(true);

  // do not check for within-limits during the Newton-scheme, since the cut-point is allowed to be not within the side and line
  bool checklimits = false;

  GEO::CUT::KERNEL::ComputeIntersection<DRT::Element::line2, DRT::Element::tri3> ci( xsi, checklimits );
  //GEO::CUT::KERNEL::DebugComputeIntersection<DRT::Element::line2, DRT::Element::tri3> ci( xsi, checklimits );

  // successful line-side intersection
  if ( ci( xyze_surfaceElement, xyze_lineElement ) )
  {
    rs(0)   = xsi(0);
    rs(1)   = xsi(1);
    line_xi = xsi(2);

    return true;
  }

  return false;

}


/*--------------------------------------------------------------------*
 * compute the cut of a ray through two points with the 2D space defined by the side
 *--------------------------------------------------------------------*/
bool GEO::CUT::ConcreteSide<DRT::Element::quad4>::RayCut( const LINALG::Matrix<3,1> & p1_xyz, const LINALG::Matrix<3,1> & p2_xyz, LINALG::Matrix<2,1> & rs, double & line_xi)
{

  LINALG::Matrix<3,4> xyze_surfaceElement(true);
  this->Coordinates(xyze_surfaceElement);

  LINALG::Matrix<3,2> xyze_lineElement(true);
  xyze_lineElement(0,0) = p1_xyz(0);
  xyze_lineElement(1,0) = p1_xyz(1);
  xyze_lineElement(2,0) = p1_xyz(2);

  xyze_lineElement(0,1) = p2_xyz(0);
  xyze_lineElement(1,1) = p2_xyz(1);
  xyze_lineElement(2,1) = p2_xyz(2);

  LINALG::Matrix<3,1> xsi(true);

  // do not check for within-limits during the Newton-scheme, since the cut-point is allowed to be not within the side and line
  bool checklimits = false;

  GEO::CUT::KERNEL::ComputeIntersection<DRT::Element::line2, DRT::Element::quad4> ci( xsi, checklimits );
  //GEO::CUT::KERNEL::DebugComputeIntersection<DRT::Element::line2, DRT::Element::quad4> ci( xsi, checklimits );

  // successful line-side intersection
  if ( ci( xyze_surfaceElement, xyze_lineElement ) )
  {
    rs(0)   = xsi(0);
    rs(1)   = xsi(1);
    line_xi = xsi(2);

    return true;
  }

  return false;

}


/*--------------------------------------------------------------------*
 * Calculates the local coordinates (rst) with respect to the element shape from its global coordinates (xyz), return if successful
 *--------------------------------------------------------------------*/
bool GEO::CUT::ConcreteSide<DRT::Element::tri3>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst, bool allow_dist )
{
  Position2d<DRT::Element::tri3> pos( *this, xyz );
  bool success = pos.Compute(allow_dist);
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
  return success;
}


/*--------------------------------------------------------------------*
 * Calculates the local coordinates (rst) with respect to the element shape from its global coordinates (xyz), return if successful
 *--------------------------------------------------------------------*/
bool GEO::CUT::ConcreteSide<DRT::Element::quad4>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst, bool allow_dist)
{
  Position2d<DRT::Element::quad4> pos( *this, xyz );
  bool success = pos.Compute(allow_dist);
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
  return success;
}


/*--------------------------------------------------------------------*
 * Calculates the normal vector with respect to the element shape at local coordinates xsi
 *--------------------------------------------------------------------*/
void GEO::CUT::ConcreteSide<DRT::Element::tri3>::Normal( const LINALG::Matrix<2,1> & xsi, LINALG::Matrix<3,1> & normal )
{
  // get derivatives at pos
  LINALG::Matrix<3,3> side_xyze( true );
  this->Coordinates(side_xyze);

  LINALG::Matrix<2,3> deriv(true);
  LINALG::Matrix<2,3> A(true);

  DRT::UTILS::shape_function_2D_deriv1(deriv, xsi(0), xsi(1), DRT::Element::tri3);
  A.MultiplyNT( deriv, side_xyze );

  // cross product to get the normal at the point
  normal( 0 ) = A( 0, 1 )*A( 1, 2 ) - A( 0, 2 )*A( 1, 1 );
  normal( 1 ) = A( 0, 2 )*A( 1, 0 ) - A( 0, 0 )*A( 1, 2 );
  normal( 2 ) = A( 0, 0 )*A( 1, 1 ) - A( 0, 1 )*A( 1, 0 );

  double norm = normal.Norm2();
  normal.Scale( 1./norm );
}


/*--------------------------------------------------------------------*
 * Calculates the normal vector with respect to the element shape at local coordinates xsi
 *--------------------------------------------------------------------*/
void GEO::CUT::ConcreteSide<DRT::Element::quad4>::Normal( const LINALG::Matrix<2,1> & xsi, LINALG::Matrix<3,1> & normal )
{
  // get derivatives at pos
  LINALG::Matrix<3,4> side_xyze( true );
  this->Coordinates(side_xyze);

  LINALG::Matrix<2,4> deriv(true);
  LINALG::Matrix<2,3> A(true);

  DRT::UTILS::shape_function_2D_deriv1(deriv, xsi(0), xsi(1), DRT::Element::quad4);
  A.MultiplyNT( deriv, side_xyze );

  // cross product to get the normal at the point
  normal( 0 ) = A( 0, 1 )*A( 1, 2 ) - A( 0, 2 )*A( 1, 1 );
  normal( 1 ) = A( 0, 2 )*A( 1, 0 ) - A( 0, 0 )*A( 1, 2 );
  normal( 2 ) = A( 0, 0 )*A( 1, 1 ) - A( 0, 1 )*A( 1, 0 );

  double norm = normal.Norm2();
  normal.Scale( 1./norm );
}
