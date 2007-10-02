/*!
\file intersection.H

\brief collection of math tools for the interface computation of two
       curved meshes

<pre>
Maintainer: Ursula Mayer
</pre>
*/

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef INTERSECTION_MATH_H_
#define INTERSECTION_MATH_H_


#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_element.H"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RefCountPtr.hpp"


using namespace std;
using namespace DRT;
using namespace DRT::Utils;
using namespace Teuchos;


namespace XFEM
{


    /*!
    \brief Returns the sum of two Epetra_SerialDenseVectors

    \param v1 (in) : arbitrary Epetra_SerialDenseVector
    \param v2 (in) : arbitrary Epetra_SerialDenseVector
    \return sum of two Epetra_SerialDenseVectors
    */
    Epetra_SerialDenseVector addTwoVectors(     
        const Epetra_SerialDenseVector&   v1,
        const Epetra_SerialDenseVector&   v2
        );      
    
    
    /*!
    \brief Returns the sum of two vector<double>

    \param v1 (in) :    arbitrary vector<double>
    \param v2 (in) :    arbitrary vector<double>
    \return sum of two vector<double>
    */
    vector<double> addTwoVectors(   
        const vector<double>&   v1,
        const vector<double>&   v2
        );
    
                                          
    /*!
    \brief Returns the difference of Epetra_SerialDenseVectors v1 - v2
    \param v1 (in) : Epetra_SerialDenseVector      
    \param v2 (in) : Epetra_SerialDenseVector
    \return difference of Epetra_SerialDenseVectors v1 - v2
    */
    Epetra_SerialDenseVector subtractsTwoVectors(   
        const Epetra_SerialDenseVector&   v1,
        const Epetra_SerialDenseVector&   v2
        );       
    
                                        
    /*!
    \brief Returns the difference of two vector<double> v1 - v2
    \param v1 (in) : arbitrary vector<double>       
    \param v2 (in) : arbitrary vector<double>
    \return difference of two vector<double> v1 - v2
    */
    vector<double> subtractsTwoVectors( 
        const vector <double>&   v1,
        const vector <double>&   v2
        );
    
             
    /*!
    \brief computes the cross product of two Epetra_SerialDenseVectors a x b
    \param a (in) : arbitrary Epetra_SerialDenseVector     
    \param b (in) : arbitrary Epetra_SerialDenseVector
    \return cross product of two Epetra_SerialDenseVectors a x b
    */
    Epetra_SerialDenseVector computeCrossProduct(   
        const Epetra_SerialDenseVector& a,
        const Epetra_SerialDenseVector& b
        );   
    
                                                                                           
    /*!
    \brief compute absolute value (L2-norm) of a Epetra_SerialDenseVector
    \param a    (in) : arbitrary Epetra_SerialDenseVector  
    \return absolute value
    */                                            
    double computeAbsoluteValue( 
        const Epetra_SerialDenseVector&   a
        );
    
                   
    /*!
    \brief compute a scalar vector multiplication
    \param s (in) : scalar  
    \param v (in) : Epetra_SerialDenseVector
    \return Epetra_SerialDenseVector resulting from scalar vector multiplication
    */       
    Epetra_SerialDenseVector computeScalarVectorMultiplication( 
        const double                        s, 
        const Epetra_SerialDenseVector&     v
        ); 
    
     
    /*!
    \brief normalizes a Epetra_SerialDenseVector
    \param v (in) : Epetra_SerialDenseVector
    \return normalized Epetra_SerialDenseVector 
    */   
    Epetra_SerialDenseVector normalizeVector(   
        const Epetra_SerialDenseVector&     v
        );


    /*!
    \brief  computes the Theorem of Pythagoras
            (a^2 + b^2)^(1/2)
            (modified from NUMERICAL RECIPES)
                   
    \param a    (in)        : a
    \param b    (in)        : b
    \return result of the Theorem of Pythagoras
    */    
    double pythagoras(
        double  a, 
        double  b
        );
        
       
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
        int n,
        int m
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
        Epetra_SerialDenseVector&   b,
        Epetra_SerialDenseVector&   x,
        int dim
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
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET


