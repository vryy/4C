/*!
\file intersection.H

\brief collection of math tools for the interface computation of two
       curved meshes

<pre>
Maintainer: Ursula Mayer
</pre>
*/

#ifdef CCADISCRET
#ifndef INTERSECTION_MATH_H_
#define INTERSECTION_MATH_H_

#include "intersection_service.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_element.H"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RefCountPtr.hpp"


namespace XFEM
{

    /*!
    \brief  Given a matrix A[1..m][1..n], this method computes its 
            singular value decomposition, A =U ·W ·V T .
            The matrix U replaces a on output. The diagonal matrix of 
            singular values W is output as a vector W[1..n]. 
            The matrix V (not the transpose V T ) is output as V[1..n][1..n].
            (modified from NUMERICAL RECIPES) 
                   
    \param A    (in/out)        : system matrix / matrix U
    \param W    (out)           : diagonal matrix stored in a vector
    \param W    (out)           : V matrix (not its transpose)
    */       
    void svdcmp(
        Epetra_SerialDenseMatrix&  A, 
        Epetra_SerialDenseVector&  W, 
        Epetra_SerialDenseMatrix&  V,
        const int n,
        const int m
        );
        
     
    /*!
    \brief  Solves a linear system of equations with help of
            a singular value decomposition. The system mytrix
            may be singular !!
                   
    \param U    (in/out)        : system matrix / matrix U
    \param b    (in)            : right-hand-side
    \param x    (out)           : solution vector
    \param dim  (in)            : dimension of the square system matrix
    */          
    bool solveLinearSystemWithSVD(
        Epetra_SerialDenseMatrix&   U,
        const Epetra_SerialDenseVector&   b,
        Epetra_SerialDenseVector&   x,
        const int dim
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
        Epetra_SerialDenseMatrix&   A,
        Epetra_SerialDenseMatrix&   U,
        Epetra_SerialDenseVector&   W,
        Epetra_SerialDenseMatrix&   V,
        int dim
        );
        
        
    /*!
    \brief computes a Gaussian elimination for a linear system of equations

    \param A        (in)    : system matrix
    \param b        (in)    : right-hand-side
    \param x        (out)   : solution vector
    \param do_piv   (in)    : do_piv = true does pivoting, do_piv = false does not do pivoting
    \param dim      (in)    : dimension of the matrix
    \return true if matrix is not singular , false if matrix is singular
    */
    bool gaussElimination(  Epetra_SerialDenseMatrix&   A,
                            Epetra_SerialDenseVector&   b,
                            Epetra_SerialDenseVector&   x,
                            const bool                  do_piv,
                            const int                   dim,
                            const int                   order);
    
        
};


#endif  //INTERSECTION_MATH_H_
#endif  // #ifdef CCADISCRET


