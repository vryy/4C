#ifndef LINALG_GAUSS_H
#define LINALG_GAUSS_H


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
          pivot = (fabs(A(pivot,pivot)) < fabs(A(i,k))) ? i : pivot;
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

        if ( A(k,k)==0 )
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

    return det;
  }
  else
  {
    x(0,0) = b(0,0)/A(0,0);
    return x(0,0);
  }
}

}

#endif
