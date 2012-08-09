/*!
\file intersection_math.H

\brief collection of math tools for the interface computation of two
       curved meshes

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/


#ifndef INTERSECTION_MATH_H_
#define INTERSECTION_MATH_H_


#include "../linalg/linalg_fixedsizematrix.H"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialDenseVector.h"


namespace GEO
{


    //! tolerance used for residuals in Newton-methods
    const double                 TOL14 = 1e-14;

    //! tolerance used for residuals in Newton-methods
    const double                 TOL13 = 1e-13;
    const double                 TOL12 = 1e-12;

    //! tolerance used for residuals in Newton-methods
    const double                 TOL10 = 1e-10;

    //! named tolerance for easy search/grep
    const double                 TOL7 = 1e-7;
    //! named tolerance for easy search/grep
    const double                 TOL6 = 1e-6;
    //! named tolerance for easy search/grep
    const double                 TOL5 = 1e-5;
    //! named tolerance for easy search/grep
    const double                 TOL4 = 1e-4;
    //! named tolerance for easy search/grep
    const double                 TOL3 = 1e-3;
    //! named tolerance for easy search/grep
    const double                 TOL2 = 1e-2;
    //! named tolerance for easy search/grep
    const double                 TOLPLUS8 = 1e8;
    //! named tolerance for easy search/grep
    const double                 TOLPLUS15 = 1e15;
    //! large number to start computations of nearest distance in tree
    const double                 LARGENUMBER = 1e30;



    inline double XSIGN(const double a, const double b)
    {
        return (b >= 0.0 ? fabs(a) : -fabs(a));
    }
    
 
  

    /*!
    \brief calculate a*a
    \param a
    \return a^2
     */
    template <typename T>
    inline T sqr(const T a)
    {
        return a * a;
    }



    /*!
    \brief  computes the Theorem of Pythagoras
            (a^2 + b^2)^(1/2)
            (modified from NUMERICAL RECIPES)

    \param a    (in)        : a
    \param b    (in)        : b
    \return result of the Theorem of Pythagoras
    */
    inline double pythagoras(
        const double  a,
        const double  b
        )
    {
        //note: sqrt(0.0) is properly defined and returns 0.0 as expected
      //      no check for zero needed
        return sqrt(sqr(a)+sqr(b));
    }


    /*!
    \brief sum up the absolute value of all entries of a vector
    \param v (in) : LINALG::Matrix<3,1>
    \return sum
    */
    inline double SumOfFabsEntries(
        const LINALG::Matrix<3,1>&     v
        )
    {
        return fabs(v(0)) + fabs(v(1)) + fabs(v(2));
    }



    /*!
    \brief computes the cross product of two LINALG::Matrix<3,1> a x b
    \param a (in) : arbitrary vector
    \param b (in) : arbitrary vector
    \return cross product of two LINALG::Matrix<3,1> a x b
    */
    LINALG::Matrix<3,1> computeCrossProduct(
        const LINALG::Matrix<3,1>& a,
        const LINALG::Matrix<3,1>& b
        );


  
/*-------------------------------------------------------------------*
| 	VON MIR GESCHRIEBEN                                             |
*--------------------------------------------------------------------*/

    /*!
        \brief  Given a matrix A[1..m][1..n], this method computes its
                singular value decomposition, A =U ·W ·V T .
                The matrix U replaces a on output. The diagonal matrix of
                singular values W is output as a vector W[1..n].
                The matrix V (not the transpose V T ) is output as V[1..n][1..n].
                (modified from NUMERICAL RECIPES)

        \param A    (in/out)        : system matrix / matrix U
        \param W    (out)           : diagonal matrix stored in a vector
        \param V    (out)           : V matrix (not its transpose)
        */
    
    //Selbe Funktion wie svdcmp nur SerialDenseMatrizen als Eingabeparameter

 void svdcmpSerialDense(
     Epetra_SerialDenseMatrix&  A,
     Epetra_SerialDenseMatrix&  W,
     Epetra_SerialDenseMatrix&  V
         );

  
  /*!
  \brief  Tests the singular value decomposition

  \param A    (in)        : system matrix
  \param U    (in)        : decomposed matrix U
  \param W    (in)        : diagonal matrix W stored in  vector
  \param V    (in)        : decomposed matrix V
  \param dim  (in)        : dimension
  */
  void test_svdcmp(
      LINALG::Matrix<3,3>&   A,
      LINALG::Matrix<3,3>&   U,
      LINALG::Matrix<3,1>&   W,
      LINALG::Matrix<3,3>&   V,
      int dim
      );


  /*!
  \brief computes a Gaussian elimination for a linear system of equations

  \tparam do_piv   (in)    : do_piv = true does pivoting, do_piv = false does not do pivoting
  \tparam dim      (in)    : dimension of the matrix
  \tparam order    (in)    : ???
  \return true if matrix is not singular , false if matrix is singular
  */
  template<bool do_piv, int dim, class M, class V1, class V2>
  inline bool gaussElimination(
          M&                    A,  ///< (in)    : system matrix
          V1&                   b,  ///< (in)    : right-hand-side
          V2&                   x,  ///< (out)   : solution vector
          const double          tol
          )
  {
      bool solution = true;

      //printf("A(0,) = %8e %8e %8e | %8e\n", A(0,0),A(0,1),A(0,2),b(0));
      //printf("A(1,) = %8e %8e %8e | %8e\n", A(1,0),A(1,1),A(1,2),b(1));
      //printf("A(2,) = %8e %8e %8e | %8e\n\n", A(2,0),A(2,1),A(2,2),b(2));
      if(dim > 1)
      {

          if (!do_piv)
          {
              for (int k=0;k<dim;k++)
              {
                  A(k,k)=1./A(k,k);

                  for (int i=k+1;i<dim;i++)
                  {
                      A(i,k) *= A(k,k);
                      x(i) = A(i,k);

                      for (int j=k+1;j<dim;j++)
                      {
                          A(i,j) -= A(i,k) * A(k,j);
                      }
                  }

                  for (int i=k+1;i<dim;i++)
                  {
                      b(i) -= x(i)*b(k);
                  }
              }
          }
          else
          {
              for (int k=0;k<dim;k++)
              {
                  int pivot = k;
                  /* search for pivot element */
                  for (int i=k+1;i<dim;i++)
                  {
                      pivot = (fabs(A(pivot,pivot)) < fabs(A(i,k))) ? i : pivot;
                  }
                  /* copy pivot row to current row */
                  if (pivot != k)
                  {
                      LINALG::Matrix<4,1> tmp;    // check changed
                      for (int j=0;j<dim;j++)
                          tmp(j) = A(pivot,j);
                      tmp(dim) = b(pivot);
                      for (int j=0;j<dim;j++)
                          A(pivot,j) = A(k,j);
                      b(pivot) = b(k);
                      for (int j=0;j<dim;j++)
                          A(k,j) = tmp(j);
                      b(k) = tmp(dim);
                  }

                  A(k,k) = 1./A(k,k);
                  //printf("inf_diag = %8e\n", A(k,k));
                  //fflush(NULL);

                  for (int i=k+1;i<dim;i++)
                  {
                      A(i,k) *= A(k,k);
                      x(i) = A(i,k);

                      for (int j=k+1;j<dim;j++)
                      {
                          A(i,j) -= A(i,k) * A(k,j);
                      }
                  }

                  for (int i=k+1;i<dim;i++)
                  {
                      b(i) -= x(i)*b(k);
                  }
                  //printf("A(0,) = %8e %8e %8e | %8e\n", A(0,0),A(0,1),A(0,2),b(0));
                  //printf("A(1,) = %8e %8e %8e | %8e\n", A(1,0),A(1,1),A(1,2),b(1));
                  //printf("A(2,) = %8e %8e %8e | %8e\n\n", A(2,0),A(2,1),A(2,2),b(2));
              }
          }


          /* backward substitution */
          x(dim-1) = b(dim-1) * A(dim-1,dim-1);

          for (int i=dim-2;i>=0;i--)
          {
              for (int j=dim-1;j>i;j--)
              {
                  b(i) -= A(i,j)*x(j);
              }
              x(i) = b(i)*A(i,i);
          }

          //for (i=0;i<dim;i++)
          //    printf("%8e ",x(i));
          //printf("\n");

          double det = 1.0;
          for(int i = 0 ; i < dim; i++)
              det *= 1.0/A(i,i);
          //printf("matrix is singular A1 = %f, A2 = %f, A3 = %f, det = %f\n ", 1/A(0,0), 1/A(1,1), 1/A(2,2), fabs(det) );
          //printf("det = %f\n ", fabs(det) );

          if(fabs(det) < tol)
          {
              solution = false;
              //printf("matrix is singular A1 = %f, A2 = %f, A3 = %f, det = %f\n ", 1/A(0,0), 1/A(1,1), 1/A(2,2), det );
          }

      }
      else
      {
          if( fabs(A(0,0)) < tol)
          {
              printf("singular \n");
              solution = false;
          }
          x(0) = b(0)/A(0,0);

          printf("x = %f\n", x(0));
      }

      return solution;
  }


  /*----------------------------------------------------------------------*
   |  ML:     solves a linear system of equation               u.may 02/08|
   |          with help of a Gaussian Elimination provide by Epetra       |
   *----------------------------------------------------------------------*/
  template<int dim, class MA, class Vb, class Vx>
    inline bool gaussEliminationEpetraT(
          MA&   A,
          Vb&   b,
          Vx&   x,
          const double tol)
  {
      bool solution = true;

      // view on hopefully column arrays
      /*Epetra_SerialDenseMatrix A_Epetra(Copy, A.data(), A.columns(), A.rows(), A.columns());
      Epetra_SerialDenseVector b_Epetra(Copy, b.data(), b.rows());
      Epetra_SerialDenseVector x_Epetra(dim);
      */
      Epetra_SerialDenseMatrix A_Epetra(dim,dim);
      Epetra_SerialDenseVector b_Epetra(dim);
      Epetra_SerialDenseVector x_Epetra(dim);

      // copy
      for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
          A_Epetra(i,j) = A(i,j);

      for(int i = 0; i < dim; i++)
        b_Epetra(i) = b(i);

      Epetra_SerialDenseSolver ge;
      ge.SetMatrix(A_Epetra);
      ge.SetVectors(x_Epetra, b_Epetra);


      // Lu factorization
      ge.Factor();
      Epetra_SerialDenseMatrix* FactoredMa = ge.FactoredMatrix();

      if(!ge.Factored())
          return false;

      // check if singular by computing the determinate
      double det = 1.0;
      for(int i = 0; i < dim; i++)
      {
          det = det*((*FactoredMa)(i,i) );
          //printf("det =  %f\n", det);
      }

      //printf("matrix is singular A1 = %f, A2 = %f, A3 = %f, det = %f\n ", 1/A(0,0), 1/A(1,1), fabs(det) );
      if(fabs(det) < tol) //TODO: check! was TOL7
      {
          solution = false;
          //printf("matrix is singular A1 = %f, A2 = %f, A3 = %f, det = %f\n ", 1/A(0,0), 1/A(1,1), det );
      }
      else
      {
          ge.Solve();
      }

      //for(int i = 0; i < dim; i++)
      //    printf("X = %f, x_epetra = %f\n", x(i), x_Epetra(i));

      // copy solution
      for(int i = 0; i < dim; i++)
        x(i) = x_Epetra(i);


      return solution;
  }


}


#endif  //INTERSECTION_MATH_H_

