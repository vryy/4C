#ifndef LINALG_GAUSS_H
#define LINALG_GAUSS_H

#include <limits>
#include <stdexcept>
#include "../linalg/linalg_fixedsizematrix.H"

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
  typename LINALG::Matrix<dim, dim> & A,  ///< (in)    : system matrix
  typename LINALG::Matrix<dim, 1>   & b,  ///< (in)    : right-hand-side
  typename LINALG::Matrix<dim, 1>   & x   ///< (out)   : solution vector
  )
{
  if (dim > 1)
  {
#if 0
    LINALG::Matrix<dim, dim> cA( A );
    LINALG::Matrix<dim, 1>   cb( b );
#endif
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


/*!
  \brief computes a Gaussian elimination for a linear system of equations after row scaling (if necessary)

  \tparam dim      (in)    : dimension of the matrix
  \return determinant of system matrix
*/
template<unsigned dim>
double scaledGaussElemination(
  typename LINALG::Matrix<dim, dim> & A,  ///< (in)    : system matrix
  typename LINALG::Matrix<dim, 1>   & b,  ///< (in)    : right-hand-side
  typename LINALG::Matrix<dim, 1>   & x,  ///< (out)   : solution vector
  const double tol = 1.0e-3               ///< (in)    : tolerance if row scaling is applied
  )
{
  // scale rows of the system with inverse of the diagonal entry of A if entry is small
  double min = A(0,0);
  for(unsigned d=1; d<dim; ++d)
    min = std::min(min, A(d,d));
  if(min < tol)
  {
    for(unsigned j=0; j<dim; ++j)
    {
      const double scale = 1.0/A(j,j);
      for(unsigned i=0; i<dim; ++i)
      {
        A(j,i) *= scale;
      }
      b(j) *= scale;
    }
  }

  // no pivoting applied due to scaling
  return LINALG::gaussElimination<false, dim>( A, b, x );
}

}

#endif
