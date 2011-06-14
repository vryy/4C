#ifndef CUT_INTERSECTION_H
#define CUT_INTERSECTION_H

#include <iterator>

#include "../linalg/linalg_gauss.H"

#include "../drt_lib/drt_element.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "cut_tolerance.H"
#include "cut_point.H"
#include "cut_edge.H"
#include "cut_side.H"
#include "cut_mesh.H"

namespace GEO
{
namespace CUT
{

/// find a cut point between side and edge
template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine    = DRT::UTILS::DisTypeToNumNodePerEle<edgetype>::numNodePerElement,
          int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<sidetype>::numNodePerElement
>
class Intersection
{
public:

  Intersection( Mesh & mesh,
                Edge & edge,
                ConcreteSide<sidetype> & side )
    : mesh_( mesh ),
      edge_( edge ),
      side_( side ),
      xsi_( true ),
      x_( true )
  {
    edge_.Coordinates( xyze_lineElement_ );
    side_.Coordinates( xyze_surfaceElement_ );

#if 1
    scale_ = 1;
#else
    scale_ = 0;
    LINALG::Matrix<3,1> d;
    for ( int i=0; i<numNodesSurface; ++i )
    {
      LINALG::Matrix<3,1> x1( &xyze_surfaceElement_( 0, i ),  true );
      LINALG::Matrix<3,1> x2( &xyze_surfaceElement_( 0, ( i+1 )%numNodesSurface ),  true );
      d.Update( 1, x2, -1, x1, 0 );
      scale_ += d.Norm2();
    }
    scale_ /= numNodesSurface;

    xyze_lineElement_.Scale( 1./scale_ );
    xyze_surfaceElement_.Scale( 1./scale_ );
#endif
  }

  bool Intersect( PointSet & cuts );

  bool ComputeCurveSurfaceIntersection();

  bool AtEdge();

  bool SurfaceWithinLimits();

  bool LineWithinLimits();

  LINALG::Matrix<3,1> & FinalPoint();

private:

  double ComputeDistance( Point * p );

  bool ComputeCut( Edge * e1, Edge * e2 );

  /// add cut point that is a node to all edges and sides it touches
  void InsertCut( Node * n, PointSet & cuts )
  {
    cuts.insert( Point::InsertCut( &edge_, &side_, n ) );
  }

  /// test if we are on any edge of the cut side and add those edges
  void TestSideEdges( Point * p, const LINALG::Matrix<3,1> & xsi )
  {
    std::vector<Edge*> edges;
    TestSideEdges( p, xsi, edges );
  }

  void TestSideEdges( Point * p, const LINALG::Matrix<3,1> & xsi, std::vector<Edge*> & edges )
  {
    if ( AtEdge() )
    {
      side_.EdgeAt( xsi( 0, 0 ), xsi( 1, 0 ), edges );
      for ( std::vector<Edge*>::iterator i=edges.begin(); i!=edges.end(); ++i )
      {
        Edge * e = *i;
        p->AddEdge( e );
      }
    }
  }

  Mesh & mesh_;
  Edge & edge_;
  ConcreteSide<sidetype> & side_;

  LINALG::Matrix<3, numNodesSurface> xyze_surfaceElement_;
  LINALG::Matrix<3, numNodesLine> xyze_lineElement_;

  LINALG::Matrix<3,1> xsi_;
  LINALG::Matrix<3,1> x_;

  double scale_;
};


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::Intersect( PointSet & cuts )
{
  bool begin_within_limits = false;
  bool end_within_limits = false;

  LINALG::Matrix<3,1> begin_xsi;
  LINALG::Matrix<3,1> end_xsi;

  if ( ComputeDistance( edge_.BeginNode()->point() ) < TOLERANCE and SurfaceWithinLimits() )
  {
    begin_xsi = xsi_;
    begin_within_limits = true;
    if ( ComputeDistance( edge_.EndNode()->point() ) < TOLERANCE and SurfaceWithinLimits() )
    {
      InsertCut( edge_.BeginNode(), cuts );
      InsertCut( edge_.EndNode()  , cuts );
      TestSideEdges( edge_.BeginNode()->point(), begin_xsi );
      TestSideEdges( edge_.EndNode()  ->point(), xsi_ );
      return true;
    }
  }
  else if ( ComputeDistance( edge_.EndNode()->point() ) < TOLERANCE and SurfaceWithinLimits() )
  {
    end_xsi = xsi_;
    end_within_limits = true;
  }

  if ( ComputeCurveSurfaceIntersection() )
  {
    if ( SurfaceWithinLimits() and LineWithinLimits() )
    {
      FinalPoint();

      if ( fabs( xsi_( 2, 0 )+1 ) < TOLERANCE )
      {
        InsertCut( edge_.BeginNode(), cuts );
        TestSideEdges( edge_.BeginNode()->point(), xsi_ );
        return true;
      }
      else if ( fabs( xsi_( 2, 0 )-1 ) < TOLERANCE )
      {
        InsertCut( edge_.EndNode(), cuts );
        TestSideEdges( edge_.EndNode()->point(), xsi_ );
        return true;
      }
      else
      {
        Node * n = side_.OnNode( x_ );
        if ( n!=NULL )
        {
          InsertCut( n, cuts );
          TestSideEdges( n->point(), xsi_ );
          return true;
        }
        else
        {
          Point * p = Point::NewPoint( mesh_, x_.A(), xsi_( 2, 0 ), &edge_, &side_ );
          TestSideEdges( p, xsi_ );
          cuts.insert( p );
          return true;
        }
      }
    }
  }
  else
  {
    bool success = false;
    Point * nodalpoint = NULL;

    if ( begin_within_limits )
    {
      InsertCut( edge_.BeginNode(), cuts );

      std::vector<Edge*> cut_edges;
      nodalpoint = edge_.BeginNode()->point();
      TestSideEdges( nodalpoint, begin_xsi, cut_edges );

      success = true;
    }
    else if ( end_within_limits )
    {
      InsertCut( edge_.EndNode(), cuts );

      std::vector<Edge*> cut_edges;
      nodalpoint = edge_.EndNode()->point();
      TestSideEdges( nodalpoint, end_xsi, cut_edges );

      success = true;
    }

    // Search all side edges for a cut. We could cross the side.

    const std::vector<Edge*> & side_edges = side_.Edges();
    for ( std::vector<Edge*>::const_iterator i=side_edges.begin(); i!=side_edges.end(); ++i )
    {
      Edge * e = *i;

      if ( nodalpoint!=NULL and nodalpoint->IsCut( e ) )
      {
        // No need to do anything, since the side has a closed cycle of
        // edges. Any matching nodes will be found.
        success = true;
        continue;
      }

      PointSet cut_points;
      edge_.GetCutPoints( e, cut_points );
      if ( cut_points.size() > 0 )
      {
        // Nothing to be done. There are cut points between these edges. We
        // cannot find new ones.
      }
      else if ( ComputeCut( e, &edge_ ) )
      {
        Point * p = Point::NewPoint( mesh_, x_.A(), xsi_( 2 ), &edge_, &side_ );
        p->AddEdge( e );
        cuts.insert( p );
        success = true;
      }
      else
      {
      }
    }

    return success;
  }

  return false;
}


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::ComputeCurveSurfaceIntersection()
{
  LINALG::Matrix<numNodesSurface,1> surfaceFunct;
  LINALG::Matrix<numNodesLine,1> lineFunct;
  LINALG::Matrix<2,numNodesSurface> surfaceDeriv1;
  LINALG::Matrix<1,numNodesLine> lineDeriv1;

  LINALG::Matrix<3,3> A;
  LINALG::Matrix<3,3> B;
  LINALG::Matrix<3,1> b;
  LINALG::Matrix<3,1> dx;

  int off_count = 0;

  // calculate rhs

  xsi_ = 0;

  const double tol = LINSOLVETOL; // / scale_;

  double oldres = 0;
  double oldincr = 0;

// #ifdef DEBUG
//   std::vector<double> res;
//   res.push_back( tol );
//   std::vector<LINALG::Matrix<3,1> > xxx;
// #endif

  DRT::UTILS::shape_function_2D(surfaceFunct, xsi_(0), xsi_(1), sidetype);
  DRT::UTILS::shape_function_1D(lineFunct, xsi_(2), edgetype);

  b = 0;

#if 0
  for ( int inode=0; inode<numNodesSurface; ++inode )
    for( int isd=0; isd<3; ++isd )
    {
      b(isd) -= xyze_surfaceElement_(isd,inode) * surfaceFunct(inode);
    }

  for ( int inode=0; inode<numNodesLine; ++inode )
    for ( int isd=0; isd<3; ++isd )
    {
      b(isd) += xyze_lineElement_(isd,inode) * lineFunct(inode);
    }
#else
  b.Multiply( -1, xyze_surfaceElement_, surfaceFunct, 1 );
  b.Multiply( 1, xyze_lineElement_, lineFunct, 1 );
#endif

  for ( int iter = 0; iter < 20; ++iter )
  {
    // calculate system matrix

    A = 0;

    DRT::UTILS::shape_function_2D_deriv1(surfaceDeriv1, xsi_(0), xsi_(1), sidetype);
    DRT::UTILS::shape_function_1D_deriv1(lineDeriv1, xsi_(2), edgetype);

    for ( int inode=0; inode<numNodesSurface; ++inode )
      for( int isd=0; isd<3; ++isd )
      {
        A(isd,0) += xyze_surfaceElement_(isd,inode) * surfaceDeriv1(0,inode);
        A(isd,1) += xyze_surfaceElement_(isd,inode) * surfaceDeriv1(1,inode);
      }

    for ( int inode=0; inode<numNodesLine; ++inode )
      for ( int isd=0; isd<3; ++isd )
      {
        A(isd,2) -= xyze_lineElement_(isd,inode) * lineDeriv1(0,inode);
      }

    B.Update( 1, A, 0 );

    // solve

    dx = 0;
    double det = LINALG::gaussElimination<true, 3>( A, b, dx );
    //if ( fabs( det ) < tol )
    if ( fabs( det ) == 0 )
    {
      if ( iter==0 )
      {
        // try a different start point
        dx = 1;
      }
      else
      {
        // There is no solution. The line is parallel to the surface.
        return false;
      }
    }

    // update

    xsi_ += dx;

    if ( not SurfaceWithinLimits() or not LineWithinLimits() )
    {
      off_count += 1;
      if ( off_count==3 )
        return false;
    }
    else
    {
      off_count = 0;
    }

    // calculate rhs

    DRT::UTILS::shape_function_2D(surfaceFunct, xsi_(0), xsi_(1), sidetype);
    DRT::UTILS::shape_function_1D(lineFunct, xsi_(2), edgetype);

    b = 0;

#if 0
    for ( int inode=0; inode<numNodesSurface; ++inode )
      for( int isd=0; isd<3; ++isd )
      {
        b(isd) -= xyze_surfaceElement_(isd,inode) * surfaceFunct(inode);
      }

    for ( int inode=0; inode<numNodesLine; ++inode )
      for ( int isd=0; isd<3; ++isd )
      {
        b(isd) += xyze_lineElement_(isd,inode) * lineFunct(inode);
      }
#else
    b.Multiply( -1, xyze_surfaceElement_, surfaceFunct, 1 );
    b.Multiply( 1, xyze_lineElement_, lineFunct, 1 );
#endif

    // test

    double residual = b.Norm2();
    double incr = dx.Norm2();

// #ifdef DEBUG
//     res.push_back( residual );
//     res.push_back( incr );
//     xxx.push_back( xsi_ );
// #endif

    if ( residual < tol and incr < tol )
    {
      return true;
    }
    if ( iter > 0 and off_count == 0 )
    {
      // test for stagnation
      if ( fabs( residual-oldres ) < tol and fabs( incr-oldincr ) < tol )
      {
        return true;
      }
    }
    oldres = residual;
    oldincr = incr;
  }

//   std::cout << "============\n";
//   std::cout << xsi_;
//   std::cout << b;
//   std::cout << dx;
//   std::copy( res.begin(), res.end(), std::ostream_iterator<double>( std::cout, " " ) );
//   std::cout << "\n";

//   LINALG::Matrix<3, numNodesLine> xyze_l( xyze_lineElement_ );
//   xyze_l.Scale( scale_ );
//   std::cout << xyze_l;

//   LINALG::Matrix<3, numNodesSurface> xyze_s( xyze_surfaceElement_ );
//   xyze_s.Scale( scale_ );
//   std::cout << xyze_s;

//   FinalPoint();
//   std::cout << x_;

//   DRT::UTILS::shape_function_2D(surfaceFunct, xsi_(0), xsi_(1), sidetype);

//   b = 0;

//   for ( int inode=0; inode<numNodesSurface; ++inode )
//     for( int isd=0; isd<3; ++isd )
//     {
//       b(isd) += xyze_surfaceElement_(isd,inode) * surfaceFunct(inode);
//     }
//   b.Scale( scale_ );
//   std::cout << b;

//   std::copy( xxx.begin(), xxx.end(), std::ostream_iterator<LINALG::Matrix<3,1> >( std::cout, " " ) );

  // Failed to converge Newton loop! This is bad. We have probably a very slow
  // cut here. But lets assume we are as close to the solution as possible and
  // use it further.
  return true;
}


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
double Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::ComputeDistance( Point * p )
{
  LINALG::Matrix<numNodesSurface,1> surfaceFunct;
  LINALG::Matrix<2,numNodesSurface> surfaceDeriv1;
  LINALG::Matrix<3,numNodesSurface> surfaceDeriv2;

  LINALG::Matrix<3,3> A;
  LINALG::Matrix<3,3> B;
  LINALG::Matrix<3,3> C;
  LINALG::Matrix<3,1> b;
  LINALG::Matrix<3,1> dx;
  LINALG::Matrix<3,1> cross;

  LINALG::Matrix<3,1> px;

  p->Coordinates( px.A() );
  px.Scale( 1./scale_ );

  xsi_ = 0;

  const double tol = LINSOLVETOL / scale_;

#ifdef DEBUG
  std::vector<double> res;
  res.push_back( tol );
#endif

  for ( int iter = 0; iter < 20; ++iter )
  {
    // calculate system matrix
    // calculate rhs

    DRT::UTILS::shape_function_2D(surfaceFunct, xsi_(0), xsi_(1), sidetype);
    DRT::UTILS::shape_function_2D_deriv1(surfaceDeriv1, xsi_(0), xsi_(1), sidetype);
    DRT::UTILS::shape_function_2D_deriv2(surfaceDeriv2, xsi_(0), xsi_(1), sidetype);

    b = px;

    for ( int inode=0; inode<numNodesSurface; ++inode )
      for( int isd=0; isd<3; ++isd )
      {
        b(isd) -= xyze_surfaceElement_(isd,inode) * surfaceFunct(inode);
      }

    A = 0;

    for ( int inode=0; inode<numNodesSurface; ++inode )
      for( int isd=0; isd<3; ++isd )
      {
        A(isd,0) += xyze_surfaceElement_(isd,inode) * surfaceDeriv1(0,inode);
        A(isd,1) += xyze_surfaceElement_(isd,inode) * surfaceDeriv1(1,inode);
      }

    // cross product to get the normal at the point
    cross( 0 ) = A( 1, 0 )*A( 2, 1 ) - A( 2, 0 )*A( 1, 1 );
    cross( 1 ) = A( 2, 0 )*A( 0, 1 ) - A( 0, 0 )*A( 2, 1 );
    cross( 2 ) = A( 0, 0 )*A( 1, 1 ) - A( 1, 0 )*A( 0, 1 );

    // norm normal direction
    double cnorm = cross.Norm2();
    double fact1 = 1./cnorm;
    //double fact2 = -0.5/( cnorm*sqrt( cnorm ) );
    cross.Scale( fact1 );

    C = 0;

    for ( int inode=0; inode<numNodesSurface; ++inode )
      for( int isd=0; isd<3; ++isd )
      {
        C(isd,0) += xyze_surfaceElement_(isd,inode) * surfaceDeriv2(0,inode);
        C(isd,1) += xyze_surfaceElement_(isd,inode) * surfaceDeriv2(1,inode);
        C(isd,2) += xyze_surfaceElement_(isd,inode) * surfaceDeriv2(2,inode);
      }

    b( 0 ) -= cross( 0 )*xsi_( 2 );
    b( 1 ) -= cross( 1 )*xsi_( 2 );
    b( 2 ) -= cross( 2 )*xsi_( 2 );

    A( 0, 2 ) = cross( 0 );
    A( 1, 2 ) = cross( 1 );
    A( 2, 2 ) = cross( 2 );

    B = 0;

    B( 0, 0 ) = ( C( 1, 0 )*A( 2, 1 )+A( 1, 0 )*C( 2, 2 ) - ( C( 2, 0 )*A( 1, 1 )+A( 2, 0 )*C( 1, 2 ) ) )*xsi_( 2 );
    B( 0, 1 ) = ( C( 1, 2 )*A( 2, 1 )+A( 1, 0 )*C( 2, 1 ) - ( C( 2, 2 )*A( 1, 1 )+A( 2, 0 )*C( 1, 1 ) ) )*xsi_( 2 );

    B( 1, 0 ) = ( C( 2, 0 )*A( 0, 1 )+A( 2, 0 )*C( 0, 2 ) - ( C( 0, 0 )*A( 2, 1 )+A( 0, 0 )*C( 2, 2 ) ) )*xsi_( 2 );
    B( 1, 1 ) = ( C( 2, 2 )*A( 0, 1 )+A( 2, 0 )*C( 0, 1 ) - ( C( 0, 2 )*A( 2, 1 )+A( 0, 0 )*C( 2, 1 ) ) )*xsi_( 2 );

    B( 2, 0 ) = ( C( 0, 0 )*A( 1, 1 )+A( 0, 0 )*C( 1, 2 ) - ( C( 1, 0 )*A( 0, 1 )+A( 1, 0 )*C( 0, 2 ) ) )*xsi_( 2 );
    B( 2, 1 ) = ( C( 0, 2 )*A( 1, 1 )+A( 0, 0 )*C( 1, 1 ) - ( C( 1, 2 )*A( 0, 1 )+A( 1, 0 )*C( 0, 1 ) ) )*xsi_( 2 );

    // product rule for normed cross product

    B( 0, 0 ) *= fact1; //+ fact2 * A( 0, 2 );
    B( 0, 1 ) *= fact1; //+ fact2 * A( 0, 2 );

    B( 1, 0 ) *= fact1; //+ fact2 * A( 1, 2 );
    B( 1, 1 ) *= fact1; //+ fact2 * A( 1, 2 );

    B( 2, 0 ) *= fact1; //+ fact2 * A( 2, 2 );
    B( 2, 1 ) *= fact1; //+ fact2 * A( 2, 2 );

    A.Update( 1, B, 1 );

    // test

    double solres = xsi_.Norm2();
    double residual = b.Norm2();

#ifdef DEBUG
    res.push_back( residual );
    res.push_back( solres );
#endif

    if ( ( solres >= std::numeric_limits<double>::min() and
           residual / solres < tol ) or
         residual < tol )
    {
      return fabs( xsi_( 2, 0 ) );
    }

    // solve

    dx = 0;
    double det = LINALG::gaussElimination<true, 3>( A, b, dx );
    //if ( fabs( det ) < LINSOLVETOL )
    if ( det==0 )
    {
      throw std::runtime_error( "intersection: no side-point distance?" );
    }

    // update

    xsi_ += dx;

    if ( not SurfaceWithinLimits() )
      return fabs( xsi_( 2, 0 ) );
  }
#ifdef DEBUG
  std::copy( res.begin(), res.end(), std::ostream_iterator<double>( std::cout, " " ) );
  std::cout << "\n";
#endif
  throw std::runtime_error( "Newton scheme did not converge" );
}


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::ComputeCut( Edge * e1, Edge * e2 )
{
  return e2->ComputeCut( e1, xsi_( 2 ), x_ );
}


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::AtEdge()
{
  switch ( sidetype )
  {
  case DRT::Element::quad4:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
  {
    return ( fabs( xsi_( 0 )+1 )<TOLERANCE or fabs( xsi_( 1 )+1 )<TOLERANCE or
             fabs( xsi_( 0 )-1 )<TOLERANCE or fabs( xsi_( 1 )-1 )<TOLERANCE );
  }
  case DRT::Element::tri3:
  case DRT::Element::tri6:
  {
    return ( fabs( xsi_( 0 )+0 )<TOLERANCE or fabs( xsi_( 1 )+0 )<TOLERANCE or
             fabs( xsi_( 1 )+xsi_( 0 )-1 )<TOLERANCE );
  }
  default:
    throw std::runtime_error( "unsupported side type" );
  }
}


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::SurfaceWithinLimits()
{
  bool intersects = false;
  switch ( sidetype )
  {
  case DRT::Element::quad4:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
  {
    intersects = ( xsi_( 0 ) >= -1-TOLERANCE and xsi_( 1 ) >= -1-TOLERANCE and
                   xsi_( 0 ) <=  1+TOLERANCE and xsi_( 1 ) <=  1+TOLERANCE );
    break;
  }
  case DRT::Element::tri3:
  case DRT::Element::tri6:
  {
    intersects = ( xsi_( 0 ) >= 0-TOLERANCE and xsi_( 1 ) >= 0-TOLERANCE and
                   xsi_( 0 ) <= 1+TOLERANCE and xsi_( 1 ) <= 1+TOLERANCE and
                   xsi_( 1 ) <= ( 1 - xsi_( 0 ) )+TOLERANCE );
    break;
  }
  default:
    throw std::runtime_error( "unsupported side type" );
  }
  return intersects;
}


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::LineWithinLimits()
{
  return xsi_( 2 ) >= -1-TOLERANCE and xsi_( 2 ) <=  1+TOLERANCE;
}


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
LINALG::Matrix<3,1> & Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::FinalPoint()
{
  // get final point
  x_ = 0;
  LINALG::Matrix<numNodesLine,1> lineFunct;
  DRT::UTILS::shape_function_1D(lineFunct, xsi_(2), edgetype);
  for ( int inode=0; inode<numNodesLine; ++inode )
    for ( int isd=0; isd<3; ++isd )
    {
      x_(isd) += xyze_lineElement_(isd,inode) * lineFunct(inode);
    }
  x_.Scale( scale_ );
  return x_;
}

}
}

#endif
