
#include <iterator>

#include "../linalg/linalg_gauss.H"

#include "cut_facet.H"
#include "cut_element.H"
#include "cut_mesh.H"
#include "cut_boundarycell.H"
#include "cut_volumecell.H"

GEO::CUT::Facet::Facet( Mesh & mesh, const std::vector<Point*> & points, Side * side, bool cutsurface )
  : points_( points ),
    parentside_( side ),
    planar_( false ),
    planar_known_( false ),
    position_( cutsurface ? Point::oncutsurface : Point::undecided )
{
  FindCornerPoints();

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

int GEO::CUT::Facet::PositionSideId()
{
  int sid = SideId();
  switch ( position_ )
  {
  case Point::undecided:
    throw std::runtime_error( "undecided facet position" );
  case Point::inside:
    if ( sid > -1 )
      return sid;
    else
      return Point::inside;
  case Point::outside:
    if ( sid > -1 )
      return sid;
    else
      return Point::outside;
  case Point::oncutsurface:
    if ( sid > -1 )
      return sid;
    else
      throw std::runtime_error( "cannot have facet on cut side without cut side" );
  }
  throw std::runtime_error( "unhandled position" );
}

void GEO::CUT::Facet::Coordinates( double * x )
{
  for ( std::vector<Point*>::const_iterator i=points_.begin(); i!=points_.end(); ++i )
  {
    Point * p = *i;
    p->Coordinates( x );
    x += 3;
  }
}

void GEO::CUT::Facet::CornerCoordinates( double * x )
{
  FindCornerPoints();
  for ( std::vector<Point*>::const_iterator i=corner_points_.begin(); i!=corner_points_.end(); ++i )
  {
    Point * p = *i;
    p->Coordinates( x );
    x += 3;
  }
}

void GEO::CUT::Facet::GetAllPoints( Mesh & mesh, std::set<Point*> & cut_points )
{
  if ( IsPlanar( mesh ) )
  {
    std::copy( points_.begin(), points_.end(), std::inserter( cut_points, cut_points.begin() ) );
    for ( std::set<Facet*>::iterator i=holes_.begin(); i!=holes_.end(); ++i )
    {
      Facet * h = *i;
      h->GetAllPoints( mesh, cut_points );
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
                 std::inserter( cut_points, cut_points.begin() ) );
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

void GEO::CUT::Facet::AddHole( Facet * hole )
{
  for ( std::vector<Point*>::iterator i=hole->points_.begin(); i!=hole->points_.end(); ++i )
  {
    Point * p = *i;
    p->Register( this );
  }
  holes_.insert( hole );
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
  LINALG::Matrix<3,1> x1;
  LINALG::Matrix<3,1> x2;
  LINALG::Matrix<3,1> x3;

  LINALG::Matrix<3,1> b1;
  LINALG::Matrix<3,1> b2;
  LINALG::Matrix<3,1> b3;

  unsigned i = Normal( points, x1, x2, x3, b1, b2, b3 );
  if ( i==0 )
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

  for ( ++i; i<points.size(); ++i )
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

    if ( fabs( x2( 2 ) ) > PLANARTOL )
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

  Point * p1 = mesh.NewPoint( x2.A(), NULL, NULL /*ParentSide()*/ );
  p1->Position( Position() );
  p1->Register( this );
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

void GEO::CUT::Facet::GetLines( std::map<std::pair<Point*, Point*>, std::set<Facet*> > & lines )
{
#if 0
  // We are interested in the surrounding lines of the facet. Thus the
  // internal lines that result from a triangulation are not wanted
  // here. (There are no other facets connected to any of our internal lines.)
  if ( IsTriangulated() )
  {
    for ( std::vector<std::vector<Point*> >::iterator i=triangulation_.begin();
          i!=triangulation_.end();
          ++i )
    {
      std::vector<Point*> & points = *i;
      GetLines( points, lines );
    }
  }
  else
#endif
  {
    GetLines( points_, lines );

    // add hole lines but do not connect with parent facet
    for ( std::set<Facet*>::iterator i=holes_.begin(); i!=holes_.end(); ++i )
    {
      Facet * hole = *i;
      hole->GetLines( lines );
    }
  }
}

void GEO::CUT::Facet::GetLines( const std::vector<Point*> & points,
                                std::map<std::pair<Point*, Point*>, std::set<Facet*> > & lines )
{
  unsigned length = points.size();
  for ( unsigned i=0; i<length; ++i )
  {
    unsigned j = ( i+1 ) % length;

    Point * p1 = points[i];
    Point * p2 = points[j];

    if ( p1->Id() < p2->Id() )
    {
      lines[std::make_pair( p1, p2 )].insert( this );
    }
    else
    {
      lines[std::make_pair( p2, p1 )].insert( this );
    }
  }
}

bool GEO::CUT::Facet::IsLine( Point * p1, Point * p2 )
{
  if ( IsTriangulated() )
  {
    for ( std::vector<std::vector<Point*> >::iterator i=triangulation_.begin();
          i!=triangulation_.end();
          ++i )
    {
      std::vector<Point*> & points = *i;
      if ( IsLine( points, p1, p2 ) )
        return true;
    }
  }
  else
  {
    if ( IsLine( points_, p1, p2 ) )
      return true;
    for ( std::set<Facet*>::iterator i=holes_.begin(); i!=holes_.end(); ++i )
    {
      Facet * hole = *i;
      if ( hole->IsLine( p1, p2 ) )
        return true;
    }
  }
  return false;
}

bool GEO::CUT::Facet::IsLine( const std::vector<Point*> & points, Point * p1, Point * p2 )
{
  std::vector<Point*>::const_iterator i1 = std::find( points.begin(), points.end(), p1 );
  if ( i1!=points.end() )
  {
    std::vector<Point*>::const_iterator i2 = i1 + 1;
    if ( i2!=points.end() )
    {
      if ( *i2 == p2 )
      {
        return true;
      }
    }
    else
    {
      i2 = points.begin();
      if ( *i2 == p2 )
      {
        return true;
      }
    }
    if ( i1!=points.begin() )
    {
      i2 = i1 - 1;
      if ( *i2 == p2 )
      {
        return true;
      }
    }
    else
    {
      i2 = points.end() - 1;
      if ( *i2 == p2 )
      {
        return true;
      }
    }
  }
  return false;
}

bool GEO::CUT::Facet::Contains( Point * p ) const
{
  if ( IsTriangulated() )
  {
    for ( std::vector<std::vector<Point*> >::const_iterator i=triangulation_.begin();
          i!=triangulation_.end();
          ++i )
    {
      const std::vector<Point*> & points = *i;
      if ( std::find( points.begin(), points.end(), p ) != points.end() )
        return true;
    }
  }
  else
  {
    if ( std::find( points_.begin(), points_.end(), p ) != points_.end() )
      return true;
    for ( std::set<Facet*>::const_iterator i=holes_.begin(); i!=holes_.end(); ++i )
    {
      Facet * hole = *i;
      if ( hole->Contains( p ) )
        return true;
    }
  }
  return false;
}

bool GEO::CUT::Facet::Contains( const std::vector<Point*> & side ) const
{
  for ( std::vector<Point*>::const_iterator i=side.begin(); i!=side.end(); ++i )
  {
    Point * p = *i;
    if ( not Contains( p ) )
    {
      return false;
    }
  }
  return true;
}

bool GEO::CUT::Facet::Touches( Facet * f )
{
  for ( std::vector<Point*>::iterator i=points_.begin(); i!=points_.end(); ++i )
  {
    Point * p = *i;
    if ( f->Contains( p ) )
    {
      return true;
    }
  }
  return false;
}

void GEO::CUT::Facet::Neighbors( Point * p,
                                 const std::set<VolumeCell*> & cells,
                                 const std::set<VolumeCell*> & done,
                                 std::set<VolumeCell*> & connected,
                                 std::set<Element*> & elements )
{
  for ( std::set<VolumeCell*>::iterator i=cells_.begin(); i!=cells_.end(); ++i )
  {
    VolumeCell * c = *i;
    if ( cells.count( c )>0 )
    {
      if ( done.count( c )==0 and
           connected.count( c )==0 and
           elements.count( c->ParentElement() )==0 )
      {
        connected.insert( c );
        elements.insert( c->ParentElement() );
        c->Neighbors( p, cells, done, connected, elements );
      }
    }
  }
}

bool GEO::CUT::Facet::Equals( DRT::Element::DiscretizationType distype )
{
  if ( holes_.size()==0 )
  {
    FindCornerPoints();
    switch ( distype )
    {
    case DRT::Element::quad4:
      return corner_points_.size()==4;
    case DRT::Element::tri3:
      return corner_points_.size()==3;
    default:
      throw std::runtime_error( "unsupported distype requested" );
    }
  }
  return false;
}

unsigned GEO::CUT::Facet::Normal( const std::vector<Point*> & points,
                                  LINALG::Matrix<3,1> & x1,
                                  LINALG::Matrix<3,1> & x2,
                                  LINALG::Matrix<3,1> & x3,
                                  LINALG::Matrix<3,1> & b1,
                                  LINALG::Matrix<3,1> & b2,
                                  LINALG::Matrix<3,1> & b3 )
{
  unsigned pointsize = points.size();
  if ( pointsize < 3 )
    return 0;

  points[0]->Coordinates( x1.A() );
  points[1]->Coordinates( x2.A() );

  b1.Update( 1, x2, -1, x1, 0 );
  b1.Scale( 1./b1.Norm2() );

  if ( b1.Norm2() < std::numeric_limits<double>::min() )
    throw std::runtime_error( "same point in facet not supported" );

  bool found = false;
  unsigned i=2;
  for ( ; i<pointsize; ++i )
  {
    Point * p = points[i];
    p->Coordinates( x3.A() );

    b2.Update( 1, x3, -1, x1, 0 );
    b2.Scale( 1./b2.Norm2() );

    // cross product to get the normal at the point
    b3( 0 ) = b1( 1 )*b2( 2 ) - b1( 2 )*b2( 1 );
    b3( 1 ) = b1( 2 )*b2( 0 ) - b1( 0 )*b2( 2 );
    b3( 2 ) = b1( 0 )*b2( 1 ) - b1( 1 )*b2( 0 );

    if ( b3.Norm2() > PLANARTOL )
    {
      found = true;
      break;
    }
  }
  if ( not found )
  {
    // all on one line, no normal
    return 0;
  }

  b3.Scale( 1./b3.Norm2() );
  return i;
}

void GEO::CUT::Facet::NewTri3Cell( Mesh & mesh, VolumeCell * volume, const std::vector<Point*> & points, std::set<BoundaryCell*> & bcells )
{
//   if ( not Equals( DRT::Element::tri3 ) )
//   {
//     throw std::runtime_error( "cannot create tri3 boundary cell on facet" );
//   }

  BoundaryCell * bc = mesh.NewTri3Cell( volume, this, points );
  bcells.insert( bc );
}

// void GEO::CUT::Facet::NewTri3Cells( Mesh & mesh, VolumeCell * volume, const std::vector<Epetra_SerialDenseMatrix> & xyz, std::set<BoundaryCell*> & bcells )
// {
//   bool mine = bcells_.size()==0;

//   //if ( bcells_.size()==0 )
//   {
// #if 0
//     if ( Equals( DRT::Element::tri3 ) )
//     {
//       NewTri3Cell( mesh, volume, bcells );
//     }
//     else if ( Equals( DRT::Element::quad4 ) )
//     {
//       NewQuad4Cell( mesh, volume, bcells );
//     }
//     else
// #endif
//     {
//       for ( std::vector<Epetra_SerialDenseMatrix>::const_iterator i=xyz.begin();
//             i!=xyz.end();
//             ++i )
//       {
//         const Epetra_SerialDenseMatrix & x = *i;
//         BoundaryCell * bc = mesh.NewTri3Cell( x, volume, this );
//         bcells.insert( bc );
//         if ( mine )
//         {
//           bcells_.insert( bc );
//         }
//       }
//     }
//   }
// }

void GEO::CUT::Facet::NewQuad4Cell( Mesh & mesh, VolumeCell * volume, const std::vector<Point*> & points, std::set<BoundaryCell*> & bcells )
{
//   if ( not Equals( DRT::Element::quad4 ) )
//   {
//     throw std::runtime_error( "cannot create quad4 boundary cell on facet" );
//   }

  BoundaryCell * bc = mesh.NewQuad4Cell( volume, this, points );
  bcells.insert( bc );
}

void GEO::CUT::Facet::GetBoundaryCells( std::set<GEO::CUT::BoundaryCell*> & bcells )
{
  if ( cells_.size()==0 )
    throw std::runtime_error( "no volume cells" );

  VolumeCell * vc = *cells_.begin();

  const std::set<BoundaryCell*> & vbcells = vc->BoundaryCells();
  for ( std::set<BoundaryCell*>::const_iterator i=vbcells.begin();
        i!=vbcells.end();
        ++i )
  {
    BoundaryCell * bc = *i;
    if ( bc->GetFacet()==this )
    {
      bcells.insert( bc );
    }
  }
}

void GEO::CUT::Facet::FindCornerPoints()
{
  if ( corner_points_.size()==0 )
  {
    LINALG::Matrix<3,1> x1;
    LINALG::Matrix<3,1> x2;
    LINALG::Matrix<3,1> x3;
    LINALG::Matrix<3,1> b1;
    LINALG::Matrix<3,1> b2;
    LINALG::Matrix<3,1> b3;

    for ( unsigned i = FindNextCornerPoint( points_, x1, x2, x3, b1, b2, b3, 0 );
          true;
          i = FindNextCornerPoint( points_, x1, x2, x3, b1, b2, b3, i ) )
    {
      Point * p = points_[i];
      if ( corner_points_.size()>0 and corner_points_.front()==p )
        break;
      corner_points_.push_back( p );
    }
  }
}

unsigned GEO::CUT::Facet::FindNextCornerPoint( const std::vector<Point*> & points,
                                               LINALG::Matrix<3,1> & x1,
                                               LINALG::Matrix<3,1> & x2,
                                               LINALG::Matrix<3,1> & x3,
                                               LINALG::Matrix<3,1> & b1,
                                               LINALG::Matrix<3,1> & b2,
                                               LINALG::Matrix<3,1> & b3,
                                               unsigned i )
{
  unsigned pointsize = points.size();
  unsigned j = ( i+1 ) % pointsize;
  if ( pointsize < 3 )
  {
    return j;
  }

  points[i]->Coordinates( x1.A() );
  points[j]->Coordinates( x2.A() );

  b1.Update( 1, x2, -1, x1, 0 );

  double norm = b1.Norm2();
  if ( norm < std::numeric_limits<double>::min() )
    throw std::runtime_error( "same point in facet not supported" );

  b1.Scale( 1./norm );

  if ( b1.Norm2() < std::numeric_limits<double>::min() )
    throw std::runtime_error( "same point in facet not supported" );

  i = j;
  for ( unsigned k=2; k<pointsize; ++k )
  {
    i = ( i+1 ) % pointsize;
    Point * p = points[i];
    p->Coordinates( x3.A() );

    b2.Update( 1, x3, -1, x1, 0 );

    norm = b2.Norm2();
    if ( norm < std::numeric_limits<double>::min() )
      throw std::runtime_error( "same point in facet not supported" );

    b2.Scale( 1./norm );

    // cross product to get the normal at the point
    b3( 0 ) = b1( 1 )*b2( 2 ) - b1( 2 )*b2( 1 );
    b3( 1 ) = b1( 2 )*b2( 0 ) - b1( 0 )*b2( 2 );
    b3( 2 ) = b1( 0 )*b2( 1 ) - b1( 1 )*b2( 0 );

    if ( b3.Norm2() > PLANARTOL )
    {
      // Found. Return last node on this line.
      return ( i+pointsize-1 ) % pointsize;
    }
  }

  // All on one line. Return first and last point.
  if ( j==0 )
  {
    return 0;
  }
  else
  {
    return pointsize-1;
  }
}

void GEO::CUT::Facet::Print( std::ostream & stream )
{
  if ( points_.size() > 0 )
  {
    LINALG::Matrix<3,1> middle;
    LINALG::Matrix<3,1> x;

    middle = 0;
    for ( std::vector<Point*>::iterator i=points_.begin(); i!=points_.end(); ++i )
    {
      Point * p = *i;
      p->Coordinates( x.A() );
      middle.Update( 1, x, 1 );
      //p->Plot( stream );
    }
    middle.Scale( 1./points_.size() );
    for ( unsigned i=0; i<=points_.size(); ++i )
    {
      Point * p = points_[i % points_.size()];
      p->Coordinates( x.A() );
      x.Update( -1, middle, 1 );
      x.Scale( 0.8 );
      x.Update( 1, middle, 1 );
      stream << x( 0, 0 ) << " " << x( 1, 0 ) << " " << x( 2, 0 ) << "\n";
    }
    stream << "\n\n";

    for ( std::set<Facet*>::iterator i=holes_.begin(); i!=holes_.end(); ++i )
    {
      Facet * hole = *i;
      hole->Print( stream );
    }
  }
}

bool GEO::CUT::Facet::IsTriangle( const std::vector<Point*> & tri ) const
{
  if ( tri.size()!=3 )
    throw std::runtime_error( "three points expected" );

  return points_.size()==3 and not IsTriangulated() and not HasHoles() and Contains( tri );
}

bool GEO::CUT::Facet::IsTriangulatedSide( const std::vector<Point*> & tri ) const
{
  if ( tri.size()!=3 )
    throw std::runtime_error( "three points expected" );

  for ( std::vector<std::vector<Point*> >::const_iterator i=triangulation_.begin();
        i!=triangulation_.end();
        ++i )
  {
    const std::vector<Point*> & t = *i;
    bool found = true;
    for ( std::vector<Point*>::const_iterator i=tri.begin(); i!=tri.end(); ++i )
    {
      Point * p = *i;
      if ( std::find( t.begin(), t.end(), p )==t.end() )
      {
        found = false;
        break;
      }
    }
    if ( found )
    {
      return true;
    }
  }
  return false;
}

unsigned GEO::CUT::Facet::NumPoints()
{
  unsigned numpoints = points_.size();
  if ( IsTriangulated() )
    return numpoints + 1;
  for ( std::set<Facet*>::iterator i=holes_.begin(); i!=holes_.end(); ++i )
  {
    Facet * hole = *i;
    numpoints += hole->NumPoints();
  }
  return numpoints;
}
