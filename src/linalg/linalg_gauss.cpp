/*----------------------------------------------------------------------*/
/*!
\file linalg_gauss.cpp

\brief Gauss elimination for small nxn systems

\level 1
\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#include "linalg_gauss.H"
#include <limits>
#include <stdexcept>

namespace LINALG
{

/*!
  \brief computes a Gaussian elimination for a linear system of equations

  \tparam do_piv   (in)    : do_piv = true does pivoting, do_piv = false does not do pivoting
  \tparam dim      (in)    : dimension of the matrix
  \return determinant of system matrix
*/
template<bool do_piv, unsigned dim>
double gaussElimination(
  LINALG::Matrix<dim, dim> & A,  ///< (in)    : system matrix
  LINALG::Matrix<dim, 1>   & b,  ///< (in)    : right-hand-side
  LINALG::Matrix<dim, 1>   & x   ///< (out)   : solution vector
  )
{
  if (dim > 1)
  {
#if 0
    LINALG::Matrix<dim, dim> cA( A );
    LINALG::Matrix<dim, 1>   cb( b );
#endif
    bool changesign = false;
    if ( not do_piv )
    {
      for ( unsigned k=0; k<dim; ++k )
      {
        A(k,k) = 1. / A(k,k);

        for ( unsigned i = k+1; i<dim; ++i )
        {
          A(i,k) *= A(k,k);
          x(i) = A(i,k);

          for ( unsigned j=k+1; j<dim; ++j )
          {
            A(i,j) -= A(i,k) * A(k,j);
          }
        }

        for ( unsigned i=k+1; i<dim; ++i )
        {
          b(i) -= x(i)*b(k);
        }
      }
    }
    else
    {
      for ( unsigned k=0; k<dim; ++k )
      {
        unsigned pivot = k;

        // search for pivot element
        for ( unsigned i=k+1; i<dim; ++i )
        {
          pivot = (fabs(A(pivot,k)) < fabs(A(i,k))) ? i : pivot;
        }

        // exchange pivot row and current row
        if (pivot != k)
        {
          for ( unsigned j=0; j<dim; ++j )
          {
            std::swap( A(k,j), A(pivot,j) );
          }
          std::swap( b(k,0), b(pivot,0) );
          changesign = not changesign;
        }

        if ( fabs( A(k,k) ) < std::numeric_limits<double>::min() )
        {
          return 0;
        }

        A(k,k) = 1./A(k,k);

        for ( unsigned i=k+1; i<dim; ++i )
        {
          A(i,k) *= A(k,k);
          x(i,0) = A(i,k);

          for ( unsigned j=k+1; j<dim; ++j )
          {
            A(i,j) -= A(i,k) * A(k,j);
          }
        }

        for ( unsigned i=k+1; i<dim; ++i )
        {
          b(i,0) -= x(i,0)*b(k,0);
        }
      }
    }

    // back substitution
    x(dim-1,0) = b(dim-1,0) * A(dim-1,dim-1);

    for ( int i=dim-2; i>=0 ; --i )
    {
      for ( int j=dim-1; j>i; --j )
      {
        b(i,0) -= A(i,j)*x(j,0);
      }
      x(i,0) = b(i,0)*A(i,i);
    }

    double det = 1.0;
    for ( unsigned i = 0 ; i < dim; ++i )
      det *= 1.0/A(i,i);

    if(changesign)
      det *= -1.0;

#if 0
    double nx = x.Norm2();
    if ( nx != 0 )
    {
      LINALG::Matrix<dim, 1> res;
      res.Multiply( cA, x );
      res.Update( -1, cb, 1 );
      double nres = res.Norm2();
      if ( fabs( nres/nx ) > 1e-15 )
      {
        std::cout << nres/nx << "\n";
        throw std::runtime_error( "failed to solve linear system" );
      }
    }
#endif

    return det;
  }
  else
  {
    x(0,0) = b(0,0)/A(0,0);
    return x(0,0);
  }
}


template
double gaussElimination<true, 2>(
  LINALG::Matrix<2, 2>   & A,  ///< (in)    : system matrix
  LINALG::Matrix<2, 1>   & b,  ///< (in)    : right-hand-side
  LINALG::Matrix<2, 1>   & x   ///< (out)   : solution vector
  );
template
double gaussElimination<false, 2>(
  LINALG::Matrix<2, 2>   & A,  ///< (in)    : system matrix
  LINALG::Matrix<2, 1>   & b,  ///< (in)    : right-hand-side
  LINALG::Matrix<2, 1>   & x   ///< (out)   : solution vector
  );
template
double gaussElimination<true, 3>(
  LINALG::Matrix<3, 3>   & A,  ///< (in)    : system matrix
  LINALG::Matrix<3, 1>   & b,  ///< (in)    : right-hand-side
  LINALG::Matrix<3, 1>   & x   ///< (out)   : solution vector
  );
template
double gaussElimination<false, 3>(
  LINALG::Matrix<3, 3>   & A,  ///< (in)    : system matrix
  LINALG::Matrix<3, 1>   & b,  ///< (in)    : right-hand-side
  LINALG::Matrix<3, 1>   & x   ///< (out)   : solution vector
  );
template
double gaussElimination<true, 4>(
  LINALG::Matrix<4, 4>   & A,  ///< (in)    : system matrix
  LINALG::Matrix<4, 1>   & b,  ///< (in)    : right-hand-side
  LINALG::Matrix<4, 1>   & x   ///< (out)   : solution vector
  );
template
double gaussElimination<false, 4>(
  LINALG::Matrix<4, 4>   & A,  ///< (in)    : system matrix
  LINALG::Matrix<4, 1>   & b,  ///< (in)    : right-hand-side
  LINALG::Matrix<4, 1>   & x   ///< (out)   : solution vector
  );



/*!
  \brief computes a Gaussian elimination for a linear system of equations after infnorm scaling

  \tparam dim      (in)    : dimension of the matrix
  \return determinant of system matrix
*/
template<unsigned dim>
double scaledGaussElimination(
  LINALG::Matrix<dim, dim> & A,  ///< (in)    : system matrix
  LINALG::Matrix<dim, 1>   & b,  ///< (in)    : right-hand-side
  LINALG::Matrix<dim, 1>   & x   ///< (out)   : solution vector
  )
{
  // infnorm scaling
  for(unsigned i=0; i<dim; ++i)
  {
    // find norm of max entry in row
    double max = std::abs(A(i,0));
    for(unsigned j=1; j<dim; ++j)
    {
      const double norm = std::abs(A(i,j));
      if(norm > max)
        max = norm;
    }

    // close to zero row detected -> matrix does probably not have full rank
    if(max < 1.0e-14)
    {
      return LINALG::gaussElimination<true, dim>( A, b, x );
    }

    // scale row with inv of max entry
    const double scale = 1.0/max;
    for(unsigned j=0; j<dim; ++j)
    {
      A(i,j) *= scale;
    }
    b(i) *= scale;
  }

  // solve scaled system using pivoting
  return LINALG::gaussElimination<true, dim>( A, b, x );
}


template
double scaledGaussElimination<2>(
  LINALG::Matrix<2, 2>   & A,  ///< (in)    : system matrix
  LINALG::Matrix<2, 1>   & b,  ///< (in)    : right-hand-side
  LINALG::Matrix<2, 1>   & x   ///< (out)   : solution vector
  );
template
double scaledGaussElimination<3>(
  LINALG::Matrix<3, 3>   & A,  ///< (in)    : system matrix
  LINALG::Matrix<3, 1>   & b,  ///< (in)    : right-hand-side
  LINALG::Matrix<3, 1>   & x   ///< (out)   : solution vector
  );
template
double scaledGaussElimination<4>(
  LINALG::Matrix<4, 4>   & A,  ///< (in)    : system matrix
  LINALG::Matrix<4, 1>   & b,  ///< (in)    : right-hand-side
  LINALG::Matrix<4, 1>   & x   ///< (out)   : solution vector
  );

}
