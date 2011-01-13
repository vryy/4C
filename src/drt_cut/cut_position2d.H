#ifndef CUT_POSITION2D_H
#define CUT_POSITION2D_H

#include <iterator>

#include "../linalg/linalg_gauss.H"
#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "cut_side.H"
#include "cut_tolerance.H"

namespace GEO
{
namespace CUT
{

template <DRT::Element::DiscretizationType elementtype,
          int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<elementtype>::numNodePerElement
>
class Position2d
{
public:

  Position2d( ConcreteSide<elementtype> & side,
              const LINALG::Matrix<3,1> & xyz )
    : px_( xyz )
  {
    side.Coordinates( xyze_ );
    Scale();
  }

  Position2d( const LINALG::Matrix<3,numNodesElement> & side_xyze,
              const LINALG::Matrix<3,1> & xyz )
    : px_( xyz ),
      xyze_( side_xyze )
  {
    Scale();
  }

  bool Compute();

  bool WithinLimits();

  const LINALG::Matrix<3,1> & LocalCoordinates() { return xsi_; }

  double Distance() const { return xsi_( 2, 0 ); }

private:

  void Scale()
  {
    scale_ = 0;
    LINALG::Matrix<3,1> d;
    for ( int i=0; i<numNodesElement; ++i )
    {
      LINALG::Matrix<3,1> x1( &xyze_( 0, i ),  true );
      LINALG::Matrix<3,1> x2( &xyze_( 0, ( i+1 )%numNodesElement ),  true );
      d.Update( 1, x2, -1, x1, 0 );
      scale_ += d.Norm2();
    }
    scale_ /= numNodesElement;

    px_.Scale( 1./scale_ );
    xyze_.Scale( 1./scale_ );

    shift_ = 0;
    for ( int i=0; i<numNodesElement; ++i )
    {
      LINALG::Matrix<3,1> x1( &xyze_( 0, i ),  true );
      shift_.Update( 1, x1, 1 );
    }
    shift_.Scale( 1./numNodesElement );

    for ( int i=0; i<numNodesElement; ++i )
    {
      LINALG::Matrix<3,1> x1( &xyze_( 0, i ),  true );
      x1.Update( -1, shift_, 1 );
    }
    px_.Update( -1, shift_, 1 );
  }

  LINALG::Matrix<3,1> px_;

  LINALG::Matrix<3,1> xsi_;

  LINALG::Matrix<3, numNodesElement> xyze_;

  double scale_;

  LINALG::Matrix<3,1> shift_;
};


template <DRT::Element::DiscretizationType elementtype,
          int numNodesElement
>
bool Position2d<elementtype,numNodesElement>::Compute()
{
  LINALG::Matrix<numNodesElement,1> surfaceFunct;
  LINALG::Matrix<2,numNodesElement> surfaceDeriv1;
  LINALG::Matrix<3,numNodesElement> surfaceDeriv2;

  LINALG::Matrix<3,3> A;
  LINALG::Matrix<3,3> B;
  LINALG::Matrix<3,3> C;
  LINALG::Matrix<3,1> b;
  LINALG::Matrix<3,1> dx;
  LINALG::Matrix<3,1> cross;

  xsi_ = 0;

  const double tol = LINSOLVETOL; // / scale_;

#ifdef DEBUG
  std::vector<double> res;
  res.push_back( tol );
#endif

  for ( int iter = 0; iter < 20; ++iter )
  {
    // calculate system matrix
    // calculate rhs

    DRT::UTILS::shape_function_2D(surfaceFunct, xsi_(0), xsi_(1), elementtype);
    DRT::UTILS::shape_function_2D_deriv1(surfaceDeriv1, xsi_(0), xsi_(1), elementtype);
    DRT::UTILS::shape_function_2D_deriv2(surfaceDeriv2, xsi_(0), xsi_(1), elementtype);

    b = px_;

    for ( int inode=0; inode<numNodesElement; ++inode )
      for( int isd=0; isd<3; ++isd )
      {
        b(isd) -= xyze_(isd,inode) * surfaceFunct(inode);
      }

    A = 0;

    for ( int inode=0; inode<numNodesElement; ++inode )
      for( int isd=0; isd<3; ++isd )
      {
        A(isd,0) += xyze_(isd,inode) * surfaceDeriv1(0,inode);
        A(isd,1) += xyze_(isd,inode) * surfaceDeriv1(1,inode);
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

    for ( int inode=0; inode<numNodesElement; ++inode )
      for( int isd=0; isd<3; ++isd )
      {
        C(isd,0) += xyze_(isd,inode) * surfaceDeriv2(0,inode);
        C(isd,1) += xyze_(isd,inode) * surfaceDeriv2(1,inode);
        C(isd,2) += xyze_(isd,inode) * surfaceDeriv2(2,inode);
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

    double residual = b.Norm2();

#ifdef DEBUG
    res.push_back( residual );
#endif

    if ( fabs( residual ) < tol )
    {
      return WithinLimits();
    }

    // solve

    dx = 0;
    double det = LINALG::gaussElimination<true, 3>( A, b, dx );
    if ( fabs( det ) < tol )
    {
      throw std::runtime_error( "no side-point distance?" );
    }

    // update

    xsi_ += dx;
  }
#ifdef DEBUG
  std::copy( res.begin(), res.end(), std::ostream_iterator<double>( std::cout, " " ) );
  std::cout << "\n";
#endif
  throw std::runtime_error( "Newton scheme did not converge" );
}


template <DRT::Element::DiscretizationType elementtype,
          int numNodesElement
>
bool Position2d<elementtype,numNodesElement>::WithinLimits()
{
  bool intersects = false;
  switch ( elementtype )
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
  return intersects and fabs( xsi_( 2 ) ) < TOLERANCE;
}

}
}

#endif
