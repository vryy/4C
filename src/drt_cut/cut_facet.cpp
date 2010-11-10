
#include <iterator>

#include "../linalg/linalg_gauss.H"

#include "cut_facet.H"
#include "cut_element.H"
#include "cut_mesh.H"


int GEO::CUT::Facet::SideId()
{
  return parentside_->Id();
}

#ifdef QHULL
void GEO::CUT::Facet::GenerateTetgen( Element * element, tetgenio::facet & f, int num, int & marker, const std::vector<Point*> & pointlist )
{
  f.numberofpolygons = 1 + holes_.size();
  f.polygonlist = new tetgenio::polygon[f.numberofpolygons];

  f.numberofholes = 0;
  f.holelist = NULL;

  LINALG::Matrix<3,1> x;

  if ( planar_ )
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

  if ( element->OwnedSide( parentside_ ) )
  {
    marker = -1;
  }
  else
  {
    marker = parentside_->Id();
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

void GEO::CUT::Facet::GetPoints( std::set<Point*, PointPidLess> & points )
{
  if ( planar_ )
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

  if ( b1.Norm2()==0 )
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

    if ( fabs( b3.Norm2() ) > TOLERANCE )
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

        if ( b1.Norm2()==0 )
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

      return false;
    }
  }
  return true;
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

