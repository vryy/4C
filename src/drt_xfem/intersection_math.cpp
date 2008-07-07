/*!----------------------------------------------------------------------
\file intersection_math.cpp

\brief collection of math tools for the interface determination of trv1o meshes

    ML      math library for the interface computation
 
    
<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "intersection_service.H"
#include "intersection_math.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "Epetra_SerialDenseMatrix.h"
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <math.h>



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




#endif  // #ifdef CCADISCRET

