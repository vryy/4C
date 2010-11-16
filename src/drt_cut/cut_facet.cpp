
#include <iterator>

#include "../linalg/linalg_gauss.H"

#include "cut_facet.H"
#include "cut_element.H"
#include "cut_mesh.H"

GEO::CUT::Facet::Facet( Mesh & mesh, const std::vector<Point*> & points, Side * side, bool cutsurface )
  : points_( points ),
    parentside_( side ),
    planar_( false ),
    planar_known_( false ),
    position_( cutsurface ? Point::oncutsurface : Point::undecided )
{
  if ( cutsurface )
  {
    for ( std::vector<Point*>::const_iterator i=points.begin(); i!=points.end(); ++i )
    {
      Point * p = *i;
      p->Position( Point::oncutsurface );
    }
  }
  else
  {
    // On multiple cuts there are facets on element sides that belong to an
    // old cut surface. Thus if all nodes are on a cut surface, the facet is
    // as well.
    bool allonsurface = true;
    for ( std::vector<Point*>::const_iterator i=points.begin(); i!=points.end(); ++i )
    {
      Point * p = *i;
      if ( p->Position()!=Point::oncutsurface )
      {
        allonsurface = false;
        break;
      }
    }
    if ( allonsurface )
    {
      // If my side has an id this facet is actually on a cut
      // surface. Otherwise it could be an inside or outside facet. The actual
      // decision does not matter much.
      if ( side->Id()>-1 )
      {
        position_ = Point::oncutsurface;
      }
      else
      {
        position_ = Point::inside;
      }
    }
  }
  for ( std::vector<Point*>::const_iterator i=points.begin(); i!=points.end(); ++i )
  {
    Point * p = *i;
    p->Register( this );
  }
}

int GEO::CUT::Facet::SideId()
{
  return parentside_->Id();
}

#ifdef QHULL
void GEO::CUT::Facet::GenerateTetgen( Mesh & mesh, LinearElement * element, tetgenio::facet & f, int num, int & marker, const std::vector<Point*> & pointlist )
{
  f.numberofpolygons = 1 + holes_.size();
  f.polygonlist = new tetgenio::polygon[f.numberofpolygons];

  f.numberofholes = 0;
  f.holelist = NULL;

  LINALG::Matrix<3,1> x;

  if ( IsPlanar( mesh ) )
  {
    GeneratePolygon( f.polygonlist[0], points_, pointlist, x );

    int pos = 0;
    for ( std::set<std::vector<Point*> >::iterator i=holes_.begin();
          i!=holes_.end();
          ++i )
    {
      const std::vector<Point*> & points = *i;
      GeneratePolygon( f.polygonlist[pos+1], points, pointlist, x );
      pos += 1;
    }
  }
  else
  {
    GeneratePolygon( f.polygonlist[0], triangulation_[num], pointlist, x );
  }

  int sid = SideId();
  switch ( position_ )
  {
  case Point::undecided:
    // There is no node (point) without cut side in this facet. There are no
    // dofs anyway, so assume this is on the inside.
    throw std::runtime_error( "undecided facet position" );
  case Point::inside:
    if ( sid > -1 )
      marker = sid;
    else
      marker = Point::inside;
    break;
  case Point::outside:
    if ( sid > -1 )
      marker = sid;
    else
      marker = Point::outside;
    break;
  case Point::oncutsurface:
    if ( sid > -1 )
      marker = sid;
    else
      throw std::runtime_error( "cannot have facet on cut side without cut side" );
    break;
  }
}
#endif

#ifdef QHULL
void GEO::CUT::Facet::GeneratePolygon( tetgenio::polygon & p,
                             const std::vector<Point*> & points,
                             const std::vector<Point*> & pointlist,
                             LINALG::Matrix<3,1> & mid )
{
  p.numberofvertices = points.size();
  p.vertexlist = new int[p.numberofvertices];

  mid = 0;
  LINALG::Matrix<3,1> x;

  int pos = 0;
  for ( std::vector<Point*>::const_iterator i=points.begin(); i!=points.end(); ++i )
  {
    Point * point = *i;
    point->Coordinates( x.A() );
    mid.Update( 1, x, 1 );
    std::vector<Point*>::const_iterator iter = std::find( pointlist.begin(), pointlist.end(), point );
    if ( iter==pointlist.end() )
    {
      throw std::runtime_error( "facet point not in point list" );
    }
    p.vertexlist[pos] = iter - pointlist.begin();
    pos += 1;
  }

  // assume planar holes!
  mid.Scale( 1./points.size() );
}
#endif

void GEO::CUT::Facet::GetPoints( Mesh & mesh, std::set<Point*, PointPidLess> & points )
{
  if ( IsPlanar( mesh ) )
  {
    std::copy( points_.begin(), points_.end(),
               std::inserter( points, points.begin() ) );

    for ( std::set<std::vector<Point*> >::iterator i=holes_.begin();
          i!=holes_.end();
          ++i )
    {
      const std::vector<Point*> & hole_points = *i;
      std::copy( hole_points.begin(), hole_points.end(),
                 std::inserter( points, points.begin() ) );
    }
  }
  else
  {
    for ( std::vector<std::vector<Point*> >::iterator i=triangulation_.begin();
          i!=triangulation_.end();
          ++i )
    {
      std::vector<Point*> & tri = *i;
      std::copy( tri.begin(), tri.end(),
                 std::inserter( points, points.begin() ) );
    }
  }
}

void GEO::CUT::Facet::CreateLinearElements( Mesh & mesh )
{
  if ( holes_.size()>0 )
    throw std::runtime_error( "cannot have holes if sub-sides are created" );

  std::vector<int> nids;

  if ( IsPlanar( mesh ) )
  {
    switch ( points_.size() )
    {
    case 3:
    {
      GetNodalIds( mesh, points_, nids );
      mesh.CreateTri3( SideId(), nids );
      return;
    }
    case 4:
    {
      GetNodalIds( mesh, points_, nids );
      mesh.CreateQuad4( SideId(), nids );
      return;
    }
    default:
      CreateTriangulation( mesh, points_ );
    }
  }

  for ( std::vector<std::vector<Point*> >::iterator i=triangulation_.begin();
        i!=triangulation_.end();
        ++i )
  {
    const std::vector<Point*> & t = *i;
    nids.clear();
    GetNodalIds( mesh, t, nids );
    mesh.CreateTri3( SideId(), nids );
  }
}

bool GEO::CUT::Facet::IsPlanar( Mesh & mesh )
{
  if ( triangulation_.size()>0 )
  {
    // if there is a triangulation use it
    return false;
  }
  if ( not planar_known_ )
  {
    planar_ = IsPlanar( mesh, points_ );
    planar_known_ = true;
  }
  return planar_;
}

bool GEO::CUT::Facet::IsPlanar( Mesh & mesh, const std::vector<Point*> & points )
{
  if ( points.size() <= 3 )
    return true;

  LINALG::Matrix<3,1> x1;
  LINALG::Matrix<3,1> x2;
  LINALG::Matrix<3,1> x3;

  LINALG::Matrix<3,1> b1;
  LINALG::Matrix<3,1> b2;
  LINALG::Matrix<3,1> b3;

  points[0]->Coordinates( x1.A() );
  points[1]->Coordinates( x2.A() );

  b1.Update( 1, x2, -1, x1, 0 );

  if ( b1.Norm2() < std::numeric_limits<double>::min() )
    throw std::runtime_error( "same point in facet not supported" );

  bool found = false;
  for ( unsigned i=2; i<points.size(); ++i )
  {
    Point * p = points[i];
    p->Coordinates( x3.A() );

    b2.Update( 1, x3, -1, x1, 0 );

    // cross product to get the normal at the point
    b3( 0 ) = b1( 1 )*b2( 2 ) - b1( 2 )*b2( 1 );
    b3( 1 ) = b1( 2 )*b2( 0 ) - b1( 0 )*b2( 2 );
    b3( 2 ) = b1( 0 )*b2( 1 ) - b1( 1 )*b2( 0 );

    if ( b3.Norm2() > MINIMALTOL )
    {
      found = true;
      break;
    }
  }
  if ( not found )
  {
    // all on one line is ok
    return true;
  }

  LINALG::Matrix<3,3> A;
  std::copy( b1.A(), b1.A()+3, A.A() );
  std::copy( b2.A(), b2.A()+3, A.A()+3 );
  std::copy( b3.A(), b3.A()+3, A.A()+6 );

  //std::copy( points.begin(), points.end(), std::ostream_iterator<Point*>( std::cout, "; " ) );
  //std::cout << "\n";

  for ( unsigned i=2; i<points.size(); ++i )
  {
    Point * p = points[i];
    p->Coordinates( x3.A() );

    x3.Update( -1, x1, 1 );

    LINALG::Matrix<3,3> B;
    B = A;
    x2 = 0;
    double det = LINALG::gaussElimination<true, 3>( B, x3, x2 );
    if ( fabs( det ) < LINSOLVETOL )
    {
      throw std::runtime_error( "failed to find point position" );
    }

    if ( fabs( x2( 2 ) ) > MINIMALTOL )
    //if ( fabs( x2( 2 )*det ) > TOLERANCE )
    {
      // there is one point that is not within the plain

      CreateTriangulation( mesh, points );

      return false;
    }
  }
  return true;
}

void GEO::CUT::Facet::CreateTriangulation( Mesh & mesh, const std::vector<Point*> & points )
{
  LINALG::Matrix<3,1> x1;
  LINALG::Matrix<3,1> x2;
  LINALG::Matrix<3,1> x3;

  LINALG::Matrix<3,1> b1;
//   LINALG::Matrix<3,1> b2;
//   LINALG::Matrix<3,1> b3;

  std::vector<Point*> pts( points );
  pts.push_back( points[0] );
  std::vector<std::vector<Point*> > lines;

  for ( unsigned pos = 0; pos<pts.size(); )
  {
    lines.push_back( std::vector<Point*>() );
    std::vector<Point*> * current_line = & lines.back();

    current_line->push_back( pts[pos++] );
    if ( pos==pts.size() )
      break;
    current_line->push_back( pts[pos++] );

    ( *current_line )[0]->Coordinates( x1.A() );
    ( *current_line )[1]->Coordinates( x2.A() );

    b1.Update( 1, x2, -1, x1, 0 );

    if ( b1.Norm2() < std::numeric_limits<double>::min() )
      throw std::runtime_error( "same point in facet not supported" );

    for ( ; pos<pts.size(); ++pos )
    {
      Point * p = pts[pos];
      p->Coordinates( x3.A() );

      x3.Update( -1, x1, 1 );

      double f = x3.Norm2() / b1.Norm2();
      x3.Update( -f, b1, 1 );
      if ( x3.Norm2() > TOLERANCE )
      {
        --pos;
        break;
      }
      current_line->push_back( p );
    }
  }

  if ( lines.size()<4 )
    throw std::runtime_error( "too few lines" );

  // Invent point in middle of surface. Surface needs to be convex to
  // allow all this.

  x1 = 0;
  x2 = 0;

  for ( std::vector<Point*>::const_iterator i=points.begin(); i!=points.end(); ++i )
  {
    Point * p = *i;
    p->Coordinates( x1.A() );
    x2.Update( 1, x1, 1 );
  }

  x2.Scale( 1./points.size() );

  Point * p1 = mesh.NewPoint( x2.A(), NULL, NULL );
  triangulation_.clear();

  for ( std::vector<std::vector<Point*> >::iterator i=lines.begin(); i!=lines.end(); ++i )
  {
    std::vector<Point*> & line = *i;
    triangulation_.push_back( std::vector<Point*>() );
    triangulation_.back().push_back( p1 );
    std::copy( line.begin(), line.end(), std::back_inserter( triangulation_.back() ) );
    //std::copy( line.begin(), line.end(), std::ostream_iterator<Point*>( std::cout, "; " ) );
    //std::cout << "\n";
  }
}

void GEO::CUT::Facet::GetNodalIds( Mesh & mesh, const std::vector<Point*> & points, std::vector<int> & nids )
{
  for ( std::vector<Point*>::const_iterator i=points.begin(); i!=points.end(); ++i )
  {
    Point * p = *i;
    Node * n = p->CutNode();
    if ( n==NULL )
    {
      std::set<int> point_id;
      point_id.insert( p->Id() );
      n = mesh.GetNode( point_id, p->X() );
    }
    nids.push_back( n->Id() );
  }
}

bool GEO::CUT::Facet::Equals( const std::vector<Point*> & facet_points )
{
  if ( points_.size()!=facet_points.size() or holes_.size()>0 )
    return false;

  unsigned size = points_.size();
  unsigned shift = std::find( points_.begin(), points_.end(), facet_points.front() ) - points_.begin();

  bool forward_match = true;
  for ( unsigned i=0; i<size; ++i )
  {
    unsigned j = ( i+shift ) % size;
    if ( points_[j] != facet_points[i] )
    {
      forward_match = false;
      break;
    }
  }
  if ( not forward_match )
  {
    for ( unsigned i=0; i<size; ++i )
    {
      unsigned j = ( shift+size-i ) % size;
      if ( points_[j] != facet_points[i] )
      {
        return false;
      }
    }
  }
  return true;
}

bool GEO::CUT::Facet::IsCutSide( Side * side )
{
  if ( parentside_==side )
    return false;

  int count = 0;
  for ( std::vector<Point*>::iterator i=points_.begin(); i!=points_.end(); ++i )
  {
    Point * p = *i;
    if ( p->IsCut( side ) )
    {
      count += 1;
      if ( count > 1 )
      {
        return true;
      }
    }
  }
  return false;
}

void GEO::CUT::Facet::Position( Point::PointPosition pos )
{
  if ( position_ != pos )
  {
    position_ = pos;
    if ( pos==Point::outside or pos==Point::inside )
    {
      for ( std::vector<Point*>::iterator i=points_.begin(); i!=points_.end(); ++i )
      {
        Point * p = *i;
        if ( p->Position()==Point::undecided )
        {
          p->Position( pos );
        }
      }
    }
  }
}
