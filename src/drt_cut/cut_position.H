#ifndef CUT_POSITION_H
#define CUT_POSITION_H

#include "../linalg/linalg_gauss.H"
#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "cut_element.H"
#include "cut_tolerance.H"

namespace GEO
{
namespace CUT
{

template <DRT::Element::DiscretizationType elementtype,
          int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<elementtype>::numNodePerElement
>
class Position
{
public:

  Position( ConcreteElement<elementtype> & element,
            Point & point )
  {
    element.Coordinates( xyze_ );
    point.Coordinates( px_.A() );
    Scale();
  }

  Position( ConcreteElement<elementtype> & element,
            const LINALG::Matrix<3,1> & xyz )
    : px_( xyz )
  {
    element.Coordinates( xyze_ );
    Scale();
  }

  bool Compute();

  bool WithinLimits();

  const LINALG::Matrix<3,1> & LocalCoordinates() { return xsi_; }

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
  }

  LINALG::Matrix<3,1> px_;

  LINALG::Matrix<3,1> xsi_;

  LINALG::Matrix<3, numNodesElement> xyze_;

  double scale_;
};


template <DRT::Element::DiscretizationType elementtype,
          int numNodesElement
>
bool Position<elementtype,numNodesElement>::Compute()
{
  LINALG::Matrix<numNodesElement,1> funct;
  LINALG::Matrix<3,numNodesElement> deriv1;

  LINALG::Matrix<3,3> A;
  LINALG::Matrix<3,1> b;
  LINALG::Matrix<3,1> dx;

  xsi_ = 0;

  // calculate rhs

  DRT::UTILS::shape_function_3D( funct, xsi_( 0 ), xsi_( 1 ), xsi_( 2 ), elementtype );

  b = px_;

  for ( int inode=0; inode<numNodesElement; ++inode )
    for( int isd=0; isd<3; ++isd )
    {
      b(isd) -= xyze_(isd,inode) * funct(inode);
    }

  for ( int iter = 0; iter < 20; ++iter )
  {
    // calculate system matrix

    A = 0;

    DRT::UTILS::shape_function_3D_deriv1( deriv1, xsi_( 0 ), xsi_( 1 ), xsi_( 2 ), elementtype );

    for ( int inode=0; inode<numNodesElement; ++inode )
      for( int isd=0; isd<3; ++isd )
      {
        A(isd,0) += xyze_(isd,inode) * deriv1(0,inode);
        A(isd,1) += xyze_(isd,inode) * deriv1(1,inode);
        A(isd,2) += xyze_(isd,inode) * deriv1(2,inode);
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

    DRT::UTILS::shape_function_3D( funct, xsi_( 0 ), xsi_( 1 ), xsi_( 2 ), elementtype );

    b = px_;

    for ( int inode=0; inode<numNodesElement; ++inode )
      for( int isd=0; isd<3; ++isd )
      {
        b(isd) -= xyze_(isd,inode) * funct(inode);
      }

    // test

    double residual = b.Norm2();

    if ( fabs( residual ) < TOLERANCE )
    {
      return WithinLimits();
    }
  }
  return false;
}


template <DRT::Element::DiscretizationType elementtype,
          int numNodesElement
>
bool Position<elementtype,numNodesElement>::WithinLimits()
{
  switch ( elementtype )
  {
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
    return ( xsi_( 0 ) >= -1-TOLERANCE and xsi_( 1 ) >= -1-TOLERANCE and xsi_( 2 ) >= -1-TOLERANCE and
             xsi_( 0 ) <=  1+TOLERANCE and xsi_( 1 ) <=  1+TOLERANCE and xsi_( 2 ) <=  1+TOLERANCE );
  case DRT::Element::tet4:
  case DRT::Element::tet10:
    return ( xsi_( 0 ) >=   -TOLERANCE and xsi_( 1 ) >=   -TOLERANCE and xsi_( 2 ) >=   -TOLERANCE and
             xsi_( 0 )+xsi_( 1 )+xsi_( 2 ) <= 1+TOLERANCE );
  case DRT::Element::pyramid5:
    return ( xsi_( 0 ) >= -1-TOLERANCE and xsi_( 1 ) >= -1-TOLERANCE and xsi_( 2 ) >=   -TOLERANCE and
             xsi_( 0 ) <=  1+TOLERANCE and xsi_( 1 ) <=  1+TOLERANCE and
             ( fabs( xsi_( 0 ) ) > fabs( xsi_( 1 ) ) ?
               ( xsi_( 2 )+fabs( xsi_( 0 ) ) <= 1+TOLERANCE ) :
               ( xsi_( 2 )+fabs( xsi_( 1 ) ) <= 1+TOLERANCE ) )
      );
  case DRT::Element::wedge6:
  case DRT::Element::wedge15:
    return ( xsi_( 0 ) >=   -TOLERANCE and xsi_( 1 ) >=   -TOLERANCE and xsi_( 2 ) >= -1-TOLERANCE and
             xsi_( 2 ) <=  1+TOLERANCE and
             xsi_( 0 )+xsi_( 1 ) <= 1+TOLERANCE );
  default:
    throw std::runtime_error( "unsupported element type" );
  }
  return false;
}

}
}

#endif
