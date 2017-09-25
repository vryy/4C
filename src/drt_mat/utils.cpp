/*----------------------------------------------------------------------*/
/*!
\file utils.cpp

\brief Utility file for the material evaluation
\level 3

\maintainer Michael Hiermeier

\date Sep 20, 2017

*/
/*----------------------------------------------------------------------*/
#include "utils.H"

// initialization of const static members has to be done in the cpp file since
// GCC complains otherwise
const unsigned MAT::IMap::second_[3][3] = { {0, 3, 5},
                                            {3, 1, 4},
                                            {5, 4, 2} };
const unsigned MAT::IMap::fourth_[6][6][4] =
    { { {0,0,0,0}, {0,0,1,1}, {0,0,2,2}, {0,0,0,1}, {0,0,1,2}, {0,0,0,2} },
      { {1,1,0,0}, {1,1,1,1}, {1,1,2,2}, {1,1,0,1}, {1,1,1,2}, {1,1,0,2} },
      { {2,2,0,0}, {2,2,1,1}, {2,2,2,2}, {2,2,0,1}, {2,2,1,2}, {2,2,0,2} },
      { {0,1,0,0}, {0,1,1,1}, {0,1,2,2}, {0,1,0,1}, {0,1,1,2}, {0,1,0,2} },
      { {1,2,0,0}, {1,2,1,1}, {1,2,2,2}, {1,2,0,1}, {1,2,1,2}, {1,2,0,2} },
      { {0,2,0,0}, {0,2,1,1}, {0,2,2,2}, {0,2,0,1}, {0,2,1,2}, {0,2,0,2} } };
template< MAT::Notation type >
const double MAT::VoigtUtils<type>::unscale_fac_[6] = { 1.0, 1.0, 1.0, 0.5, 0.5, 0.5 };
template< MAT::Notation type >
const double MAT::VoigtUtils<type>::scale_fac_[6] = { 1.0, 1.0, 1.0, 2.0, 2.0, 2.0 };

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < MAT::Notation type >
void MAT::VoigtUtils<type>::SymmetricOuterProduct(
    const LINALG::Matrix<3,1>& vec_a,
    const LINALG::Matrix<3,1>& vec_b,
    LINALG::Matrix<6,1>& ab_ba )
{
  std::fill( ab_ba.A(), ab_ba.A()+6, 0.0 );

  LINALG::Matrix<3,3> outer_product;
  outer_product.MultiplyNT( vec_a, vec_b );

  for ( unsigned i=0; i<3; ++i )
    for ( unsigned j=i; j<3; ++j )
      ab_ba( IMap::second_[i][j] ) += outer_product(i,j) + outer_product(j,i);

  // scale off-diagonal values
  ScaleOffDiagonalVals( ab_ba );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < MAT::Notation type >
void MAT::VoigtUtils<type>::MultiplyTensorVector(
    const LINALG::Matrix<6,1>& strain,
    const LINALG::Matrix<3,1>& vec,
    LINALG::Matrix<3,1>& res )
{
  for ( unsigned i=0; i<3; ++i )
    for ( unsigned j=0; j<3; ++j )
    {
      const double fac = UnscaleFactor<type>( IMap::second_[i][j] );
      res(i,0) += strain(IMap::second_[i][j]) * fac * vec(j,0);
    }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < MAT::Notation type >
void MAT::VoigtUtils<type>::PowerOfSymmetricTensor(
    const unsigned pow,
    const LINALG::Matrix<6,1>& strain,
    LINALG::Matrix<6,1>& strain_pow )
{
  std::copy( strain.A(), strain.A()+6, strain_pow.A() );

  if ( pow > 1 )
  {
    // unscale the off-diagonal values
    UnscaleOffDiagonalVals( strain_pow );

    LINALG::Matrix<6,1> prod(false);
    for ( unsigned p=1; p<pow; ++p )
    {
      std::fill( prod.A(), prod.A()+6, 0.0 );

      for ( unsigned i=0; i<3; ++i )
        for ( unsigned j=i; j<3; ++j )
          for ( unsigned k=0; k<3; ++k )
            prod( IMap::second_[i][j], 0 ) += strain_pow( IMap::second_[i][k], 0 )
                * unscale_fac_[IMap::second_[k][j]] * strain( IMap::second_[k][j], 0 );

      std::copy( prod.A(), prod.A()+6, strain_pow.A() );
    }

    // scale the off-diagonal values again
    ScaleOffDiagonalVals( strain_pow );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < MAT::Notation type >
void MAT::VoigtUtils<type>::InverseTensor(
    const LINALG::Matrix<6,1>& strain,
    LINALG::Matrix<6,1>& strain_inv )
{
  LINALG::Matrix<3,3> strain_mat_inv(false);

  for( unsigned i=0; i<3; ++i )
    for ( unsigned j=0; j<3; ++j )
    {
      const double fac = UnscaleFactor<type>( IMap::second_[i][j] );
      strain_mat_inv(i,j) = fac * strain( IMap::second_[i][j] );
    }

  strain_mat_inv.Invert();

  for( unsigned i=0; i<3; ++i )
    for ( unsigned j=i; j<3; ++j )
    {
      const double fac = ScaleFactor<type>( IMap::second_[i][j] );
      strain_inv( IMap::second_[i][j] ) = fac * strain_mat_inv( i, j );
    }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < MAT::Notation type >
void MAT::VoigtUtils<type>::ScaleOffDiagonalVals(
    LINALG::Matrix<6,1>& strain )
{
  for ( unsigned i=3; i<6; ++i )
    strain(i,0) *= ScaleFactor<type>(i);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < MAT::Notation type >
void MAT::VoigtUtils<type>::UnscaleOffDiagonalVals(
    LINALG::Matrix<6,1>& strain )
{
  for ( unsigned i=3; i<6; ++i )
    strain(i,0) *= UnscaleFactor<type>(i);
}

/*----------------------------------------------------------------------------*/
template class MAT::VoigtUtils<MAT::Notation::strain>;
template class MAT::VoigtUtils<MAT::Notation::stress>;
