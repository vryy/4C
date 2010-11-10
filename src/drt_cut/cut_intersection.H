#ifndef CUT_INTERSECTION_H
#define CUT_INTERSECTION_H

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

template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine    = DRT::UTILS::DisTypeToNumNodePerEle<edgetype>::numNodePerElement,
          int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<sidetype>::numNodePerElement
>
class Intersection
{
public:

  Intersection( Mesh & mesh,
                ConcreteEdge<edgetype> & edge,
                ConcreteSide<sidetype> & side )
    : mesh_( mesh ),
      edge_( edge ),
      side_( side ),
      xsi_( true ),
      x_( true )
  {
    edge_.Coordinates( xyze_lineElement_ );
    side_.Coordinates( xyze_surfaceElement_ );
  }

  bool Intersect( std::set<Point*, PointPidLess> & cuts );

  bool ComputeCurveSurfaceIntersection();

  bool ComputeDistance( Point * p );

  bool AtEdge();

  bool SurfaceWithinLimits();

  bool LineWithinLimits();

  LINALG::Matrix<3,1> & FinalPoint();

private:

  void InsertCut( Node * n, std::set<Point*, PointPidLess> & cuts, const LINALG::Matrix<3,1> & xsi )
  {
    Point * p = n->point();
    for ( std::set<Edge*>::const_iterator i=n->Edges().begin(); i!=n->Edges().end(); ++i )
    {
      Edge * e = *i;
      p->AddEdge( e );
    }
    p->AddEdge( &edge_ );
    p->AddSide( &side_ );
    TestSideEdges( p, xsi );
    cuts.insert( p );
  }

  void TestSideEdges( Point * p, const LINALG::Matrix<3,1> & xsi )
  {
    if ( AtEdge() )
    {
      std::vector<Edge*> edges;
      side_.EdgeAt( xsi( 0, 0 ), xsi( 1, 0 ), edges );
      for ( std::vector<Edge*>::iterator i=edges.begin(); i!=edges.end(); ++i )
      {
        Edge * e = *i;
        p->AddEdge( e );
      }
    }
  }

  Mesh & mesh_;
  ConcreteEdge<edgetype> & edge_;
  ConcreteSide<sidetype> & side_;

  LINALG::Matrix<3, numNodesSurface> xyze_surfaceElement_;
  LINALG::Matrix<3, numNodesLine> xyze_lineElement_;

  LINALG::Matrix<3,1> xsi_;
  LINALG::Matrix<3,1> x_;
};


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::Intersect( std::set<Point*, PointPidLess> & cuts )
{
  ComputeDistance( edge_.BeginNode()->point() );
  if ( fabs( xsi_( 2, 0 ) ) < TOLERANCE and SurfaceWithinLimits() )
  {
    LINALG::Matrix<3,1> begin_xsi( xsi_ );
    ComputeDistance( edge_.EndNode()->point() );
    if ( fabs( xsi_( 2, 0 ) ) < TOLERANCE and SurfaceWithinLimits() )
    {
      InsertCut( edge_.BeginNode(), cuts, begin_xsi );
      InsertCut( edge_.EndNode()  , cuts, xsi_ );
      return true;
    }
  }

  if ( ComputeCurveSurfaceIntersection() )
  {
    if ( SurfaceWithinLimits() and LineWithinLimits() )
    {
      FinalPoint();

      if ( fabs( xsi_( 2, 0 )+1 ) < TOLERANCE )
      {
        InsertCut( edge_.BeginNode(), cuts, xsi_ );
        return true;
      }
      else if ( fabs( xsi_( 2, 0 )-1 ) < TOLERANCE )
      {
        InsertCut( edge_.EndNode(), cuts, xsi_ );
        return true;
      }
      else
      {
        Node * n = side_.OnNode( x_ );
        if ( n!=NULL )
        {
          InsertCut( n, cuts, xsi_ );
          return true;
        }
        else
        {
          Point * p = mesh_.NewPoint( x_.A(), &edge_, &side_ );
          p->t( &edge_, xsi_( 2, 0 ) );
          TestSideEdges( p, xsi_ );
          cuts.insert( p );
          return true;
        }
      }
    }
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
  LINALG::Matrix<2,numNodesLine> lineDeriv1;

  LINALG::Matrix<3,3> A;
  LINALG::Matrix<3,1> b;
  LINALG::Matrix<3,1> dx;

  // calculate rhs

  xsi_ = 0;

  DRT::UTILS::shape_function_2D(surfaceFunct, xsi_(0), xsi_(1), sidetype);
  DRT::UTILS::shape_function_1D(lineFunct, xsi_(2), edgetype);

  b = 0;

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

    // solve

    dx = 0;
    double det = LINALG::gaussElimination<true, 3>( A, b, dx );
    if ( fabs( det ) < LINSOLVETOL )
    {
      return false;
    }

    // update

    xsi_ += dx;

    // calculate rhs

    DRT::UTILS::shape_function_2D(surfaceFunct, xsi_(0), xsi_(1), sidetype);
    DRT::UTILS::shape_function_1D(lineFunct, xsi_(2), edgetype);

    b = 0;

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

    // test

    double residual = b.Norm2();

    if ( fabs( residual ) < LINSOLVETOL )
    {
      return true;
    }
  }
  return false;
}


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::ComputeDistance( Point * p )
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

  xsi_ = 0;

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

    A.Update( 1, B, 1 );

    // test

    double solres = xsi_.Norm2();
    double residual = b.Norm2();

    if ( ( solres != 0 and fabs( residual / solres ) < LINSOLVETOL ) or
         fabs( residual ) < LINSOLVETOL )
    {
      return true;
    }

    // solve

    dx = 0;
    double det = LINALG::gaussElimination<true, 3>( A, b, dx );
    if ( fabs( det ) < LINSOLVETOL )
    {
      throw std::runtime_error( "no side-point distance?" );
    }

    // update

    xsi_ += dx;
  }
  throw std::runtime_error( "Newton scheme did not converge" );
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
  return x_;
}

}
}

#endif
