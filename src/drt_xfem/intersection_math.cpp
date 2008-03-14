/*!----------------------------------------------------------------------
\file intersection_math.cpp

\brief collection of math tools for the interface determination of trv1o meshes

    ML      math library for the interface computation
 
    
<pre>
Maintainer: Ursula Mayer
</pre>
*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "intersection_service.H"
#include "intersection_math.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_utils_local_connectivity_matrices.H"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialDenseMatrix.h"
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <math.h>

using namespace XFEM;
using namespace DRT::UTILS;


void XFEM::test_svdcmp(
    BlitzMat&   A,
    BlitzMat&   U,
    BlitzVec&   W,
    BlitzMat&   V,
    int dim)
{

    BlitzMat H1(3,3);
    BlitzMat H2(3,3);
    
    printf("W U\n");
    for(int i = 0 ; i < dim; i++)
    {
        printf("W = %f\t", W(i));
        for(int j = 0 ; j < dim; j++)
        {
            printf("U = %f\t", U(i,j));
        }   
        printf("\n");
    }
    printf("\n");
    
    for(int i = 0 ; i < dim; i++)
    {
        for(int j = 0 ; j < dim; j++)
        {
            H1(i,j) = U(i,j)*W(j);
        }
    }
    
    printf("H1\n");
    for(int i = 0 ; i < dim; i++)
    {
        for(int j = 0 ; j < dim; j++)
        {
            printf("H1 = %f\t", H1(i,j));
        }   
        printf("\n");
    }
    printf("\n");
    
    printf("V\n");
    for(int i = 0 ; i < dim; i++)
    {
        for(int j = 0 ; j < dim; j++)
        {
            printf("V = %f\t", V(i,j));
        }   
        printf("\n");
    }
    printf("\n");
    
    for(int i = 0 ; i < dim; i++)
    {
        for(int j = 0 ; j < dim; j++)
        {
            H2(i,j) = 0.0;
            for(int k = 0 ; k < dim; k++)
            {
                H2(i,j) += H1(i,k)*V(j,k);
            }
        }
    }
    
    printf("system matrix\n");
    for(int i = 0 ; i < dim; i++)
    {
        for(int j = 0 ; j < dim; j++)
        {
            printf("A = %f\t", A(i,j));
        }   
        printf("\n");
    }
    printf("\n");
    
    printf("system matrix SVD\n");
    for(int i = 0 ; i < dim; i++)
    {
        for(int j = 0 ; j < dim; j++)
        {
            printf("H2 = %f\t", H2(i,j));
        }   
        printf("\n");
    }
    printf("\n");

}


/*----------------------------------------------------------------------*
 |  ML:     solves a linear system of equation               u.may 02/08|
 |          with help of a Gaussian Elimination provide by Epetra       |
 *----------------------------------------------------------------------*/
bool XFEM::gaussEliminationEpetra(
        BlitzMat&   A,
        BlitzVec&   b,
        BlitzVec&   x)
{
    bool solution = true;
    int dim = 2;

    // view on hopefullx column major blitz arrays
    Epetra_SerialDenseMatrix A_Epetra(Copy, A.data(), A.columns(), A.rows(), A.columns());
    Epetra_SerialDenseVector b_Epetra(Copy, b.data(), b.rows());
    Epetra_SerialDenseVector x_Epetra(2); 
    
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
        det = det*((*FactoredMa)[i][i]);
        //printf("det =  %f\n", det);
    }
    
    //printf("matrix is singular A1 = %f, A2 = %f, A3 = %f, det = %f\n ", 1/A(0,0), 1/A(1,1), fabs(det) );
    if(fabs(det) < TOL7)
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
    
     x(0) = x_Epetra(0);
     x(1) = x_Epetra(1);

    return solution;
}



#endif  // #ifdef CCADISCRET

