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
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_element.H"
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

using namespace XFEM;
using namespace DRT::UTILS;


#define                 SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
//inline double SIGN(double a, double b)
//{
//    return ((b) >= 0.0 ? fabs(a) : -fabs(a));
//}

static double           maxarg1,maxarg2;
#define                 FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
//inline double FMAX(double maxarg1, double maxarg2)
//{
//    return ((maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2));
//}

static int              iminarg1,iminarg2;
#define                 IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))
//inline double IMIN(double iminarg1, double iminarg2)
//{
//    return ((iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2));
//}

/*----------------------------------------------------------------------*
 |    computes the singular value decomposition             u.may 09/07 |
 |    of a matrix A (modified from NUMERICAL RECIPES)                   |
 *----------------------------------------------------------------------*/  
void XFEM::svdcmp(
    Epetra_SerialDenseMatrix&  A, 
    Epetra_SerialDenseVector&  W, 
    Epetra_SerialDenseMatrix&  V,
    const int n,
    const int m)
{
    int flag,i,its,j,jj,k,l,nm;
    double anorm,c,f,g,h,s,scale,x,y,z;
    Epetra_SerialDenseVector rv1(n);


    //Householder reduction to bidiagonal form.
    g=scale=anorm=0.0;

    for (i=0;i<n;i++) {
        l=i+1;
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i < m) {
            for (k=i;k<m;k++) scale += fabs(A[k][i]);
                if (scale) {
                    for (k=i;k<m;k++) {
                        A[k][i] /= scale;
                        s += A[k][i]*A[k][i];
                    }
                    f=A[i][i];
                    g = -SIGN(sqrt(s),f);
                    h=f*g-s;
                    A[i][i]=f-g;
                    for (j=l;j<n;j++) {
                        for (s=0.0,k=i;k<m;k++) s += A[k][i]*A[k][j];
                        f=s/h; 
                        for (k=i;k<m;k++) A[k][j] += f*A[k][i];
                    }
                    for (k=i;k<m;k++) A[k][i] *= scale;
                }
            }
            W[i]=scale*g;
            g=s=scale=0.0;
            if (i < m && i != (n-1)) {
                for (k=l;k<n;k++) scale += fabs(A[i][k]);
                if (scale) {
                    for (k=l;k<n;k++) {
                        A[i][k] /= scale;
                        s += A[i][k]*A[i][k];
                    }
                    f=A[i][l];
                    g = -SIGN(sqrt(s),f);
                    h=f*g-s;
                    A[i][l]=f-g;
                    for (k=l;k<n;k++) rv1[k]=A[i][k]/h;
                    for (j=l;j<m;j++) {
                        for (s=0.0,k=l;k<n;k++) s += A[j][k]*A[i][k];
                        for (k=l;k<n;k++) A[j][k] += s*rv1[k];
                    }
                    for (k=l;k<n;k++) A[i][k] *= scale;
                }
            }
            anorm=FMAX(anorm,(fabs(W[i])+fabs(rv1[i])));
        }
        //Accumulation of right-hand transformations.
        for (i=(n-1);i>=0;i--) {
            if (i < n) {
                if (g) {
                    //Double division to avoid possible underflow.
                    for (j=l;j<n;j++)
                        V[j][i]=(A[i][j]/A[i][l])/g;
                    for (j=l;j<n;j++) {
                        for (s=0.0,k=l;k<n;k++) s += A[i][k]*V[k][j];
                        for (k=l;k<n;k++) V[k][j] += s*V[k][i];
                    }
                }
                for (j=l;j<n;j++) V[i][j]=V[j][i]=0.0;
            }
            V[i][i]=1.0;
            g=rv1[i];
            l=i;
        }
        //Accumulation of left-hand transformations.
        for (i=IMIN((m-1),(n-1));i>=0;i--) {
            l=i+1;
            g=W[i];
            for (j=l;j<n;j++) A[i][j]=0.0;
            if (g) {
                g=1.0/g;
                for (j=l;j<n;j++) {
                    for (s=0.0,k=l;k<m;k++) s += A[k][i]*A[k][j];
                    f=(s/A[i][i])*g;
                    for (k=i;k<m;k++) A[k][j] += f*A[k][i];
                }
                for (j=i;j<m;j++) A[j][i] *= g;
            } else for (j=i;j<m;j++) A[j][i]=0.0;
            ++A[i][i];
        }
        //Diagonalization of the bidiagonal form: Loop over
        for (k=(n-1);k>=0;k--) {
            //singular values, and over allowed iterations.
            for (its=1;its<=30;its++) {
                flag=1;
                //Test for splitting.
                for (l=k;l>=0;l--) {
                    //Note that rv1[1] is always zero.
                    nm=l-1;
                    if ((double)(fabs(rv1[l])+anorm) == anorm) {
                        flag=0;
                        break;
                    }
                    if ((double)(fabs(W[nm])+anorm) == anorm) break;
                }
                if (flag) {
                    //Cancellation of rv1[l], if l > 1.
                    c=0.0;
                    s=1.0;
                    for (i=l;i<=k;i++) {
                        f=s*rv1[i];
                        rv1[i]=c*rv1[i];
                        if ((double)(fabs(f)+anorm) == anorm) break;
                        g=W[i];
                        h=XFEM::pythagoras(f,g);
                        W[i]=h;
                        h=1.0/h;
                        c=g*h;
                        s = -f*h;
                        for (j=0;j<m;j++) {
                            y=A[j][nm];
                            z=A[j][i];
                            A[j][nm]=y*c+z*s;
                            A[j][i]=z*c-y*s;
                        }
                    }
                }
                z=W[k];
                //Convergence.
                if (l == k) {
                    //Singular value is made nonnegative.
                    if (z < 0.0) {
                        W[k] = -z;
                        for (j=0;j<n;j++) V[j][k] = -V[j][k];
                    }
                    break;
                }
                if (its == 30) dserror("no convergence in 30 svdcmp iterations");
                //Shift from bottom 2-by-2 minor.
                x=W[l];
                nm=k-1;
                y=W[nm];
                g=rv1[nm];
                h=rv1[k];
                f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
                g=XFEM::pythagoras(f,1.0);
                f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
                //Next QR transformation:
                c=s=1.0;
                for (j=l;j<=nm;j++) {
                    i=j+1;
                    g=rv1[i];
                    y=W[i];
                    h=s*g;
                    g=c*g;
                    z=XFEM::pythagoras(f,h);
                    rv1[j]=z;
                    c=f/z;
                    s=h/z;
                    f=x*c+g*s;
                    g = g*c-x*s;
                    h=y*s;
                    y *= c;
                    for (jj=0;jj<n;jj++) {
                        x=V[jj][j];
                        z=V[jj][i];
                        V[jj][j]=x*c+z*s;
                        V[jj][i]=z*c-x*s;
                    }
                z=XFEM::pythagoras(f,h);
                //Rotation can be arbitrary if z = 0.
                W[j]=z;
                if (z) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
                for (jj=0;jj<m;jj++) {
                    y=A[jj][j];
                    z=A[jj][i];
                    A[jj][j]=y*c+z*s;
                    A[jj][i]=z*c-y*s;
                }
            }
            rv1[l]=0.0;
            rv1[k]=f;
            W[k]=x;
        }
    }
}





bool XFEM::solveLinearSystemWithSVD(
    Epetra_SerialDenseMatrix&   U,
    const Epetra_SerialDenseVector&   b,
    Epetra_SerialDenseVector&   x,
    const int dim)
{
    bool nonsingular = true;
    double svdtemp = 0.0;;
    Epetra_SerialDenseMatrix V(dim, dim);
    Epetra_SerialDenseMatrix A(dim, dim);
    Epetra_SerialDenseVector W(dim);
    
    // initialize vectors and matrices
    for(int  i = 0; i < dim; i++ )
    {
        x[i] = 0.0;
        W[i] = 0.0;
        for(int j = 0; j < dim; j++)
        {
            V[i][j] = 0.0;
        }
    }
    
    A = U;
    
    
    svdcmp( U, W, V, dim, dim );
    
    //test_svdcmp( A, U, W, V, 3);
       
    
    for(int i = 0; i < dim; i++ )
    {
        for(int  k = 0; k < dim; k++ )
        {
            svdtemp = 0.0;
            for(int  j = 0; j < dim; j++ )
            {
                if( fabs(W[j]) > 1e-7 )
                    svdtemp += V[i][j]*U[k][j] / W[j];
                  
            }
            x[i] += svdtemp * b[k];
        }
    }
    
    for(int  j = 0; j < dim; j++ )
        if( fabs(W[j]) <= 1e-7 )
        {
             nonsingular = false;
             break;
        }
    
   Epetra_SerialDenseVector b1(3);
 /*  for(int i = 0; i < dim; i++ )
    {
        b1[i] = 0.0;
        for(int  k = 0; k < dim; k++ )
        {  
            b1[i] += A[i][k]*x[k];
        }
        printf("b = %f\t  b1 = %f\n", b[i], b1[i]);
    }
    printf("\n");
    
    for(int i = 0; i < dim; i++ )
        printf("x = %f\t", x[i]);
        
    printf("\n");
    */
    return nonsingular;
}


void XFEM::test_svdcmp(
    Epetra_SerialDenseMatrix&   A,
    Epetra_SerialDenseMatrix&   U,
    Epetra_SerialDenseVector&   W,
    Epetra_SerialDenseMatrix&   V,
    int dim)
{

    Epetra_SerialDenseMatrix H1(3,3);
    Epetra_SerialDenseMatrix H2(3,3);
    
    printf("W U\n");
    for(int i = 0 ; i < dim; i++)
    {
        printf("W = %f\t", W[i]);
        for(int j = 0 ; j < dim; j++)
        {
            printf("U = %f\t", U[i][j]);
        }   
        printf("\n");
    }
    printf("\n");
    
    for(int i = 0 ; i < dim; i++)
    {
        for(int j = 0 ; j < dim; j++)
        {
            H1[i][j] = U[i][j]*W[j];
        }
    }
    
    printf("H1\n");
    for(int i = 0 ; i < dim; i++)
    {
        for(int j = 0 ; j < dim; j++)
        {
            printf("H1 = %f\t", H1[i][j]);
        }   
        printf("\n");
    }
    printf("\n");
    
    printf("V\n");
    for(int i = 0 ; i < dim; i++)
    {
        for(int j = 0 ; j < dim; j++)
        {
            printf("V = %f\t", V[i][j]);
        }   
        printf("\n");
    }
    printf("\n");
    
    for(int i = 0 ; i < dim; i++)
    {
        for(int j = 0 ; j < dim; j++)
        {
            H2[i][j] = 0.0;
            for(int k = 0 ; k < dim; k++)
            {
                H2[i][j] += H1[i][k]*V[j][k];
            }
        }
    }
    
    printf("system matrix\n");
    for(int i = 0 ; i < dim; i++)
    {
        for(int j = 0 ; j < dim; j++)
        {
            printf("A = %f\t", A[i][j]);
        }   
        printf("\n");
    }
    printf("\n");
    
    printf("system matrix SVD\n");
    for(int i = 0 ; i < dim; i++)
    {
        for(int j = 0 ; j < dim; j++)
        {
            printf("H2 = %f\t", H2[i][j]);
        }   
        printf("\n");
    }
    printf("\n");

}



/*----------------------------------------------------------------------*
 |  ML:     solves a linear system of equation               u.may 06/07|
 |          with help of a Gaussian Elimination                         |
 *----------------------------------------------------------------------*/
bool XFEM::gaussElimination(    
    Epetra_SerialDenseMatrix&   A,
    Epetra_SerialDenseVector&   b,
    Epetra_SerialDenseVector&   x,
    const bool                  do_piv,
    const int                   dim,    
    const int                   order)
{
    bool solution = true;

    //printf("A[0][] = %8e %8e %8e | %8e\n", A[0][0],A[0][1],A[0][2],b[0]);
    //printf("A[1][] = %8e %8e %8e | %8e\n", A[1][0],A[1][1],A[1][2],b[1]);
    //printf("A[2][] = %8e %8e %8e | %8e\n\n", A[2][0],A[2][1],A[2][2],b[2]);
if(dim > 1)
{

    if (!do_piv) {
        for (int k=0;k<dim;k++)
        {
            A[k][k]=1./A[k][k];

            for (int i=k+1;i<dim;i++)
            {
                A[i][k] = A[i][k] * A[k][k];
                x[i] = A[i][k];

                for (int j=k+1;j<dim;j++)
                {
                    A[i][j] = A[i][j] - A[i][k] * A[k][j];
                }
            }

            for (int i=k+1;i<dim;i++)
            {
                b[i]=b[i]-x[i]*b[k];
            }
        }
    }
    else {
        for (int k=0;k<dim;k++)
        {
            int pivot = k;
            /* search for pivot element */
            for (int i=k+1;i<dim;i++)
            {
                pivot = (fabs(A[pivot][pivot]) < fabs(A[i][k])) ? i : pivot;
            }
            /* copy pivot row to current row */
            if (pivot != k)
            {
                double tmp[4];    // check changed
                for (int j=0;j<dim;j++)
                    tmp[j] = A[pivot][j];
                tmp[dim] = b[pivot];
                for (int j=0;j<dim;j++)
                    A[pivot][j] = A[k][j];
                b[pivot] = b[k];
                for (int j=0;j<dim;j++)
                    A[k][j] = tmp[j];
                b[k] = tmp[dim];
            }

            A[k][k]=1./A[k][k];
            //printf("inf_diag = %8e\n", A[k][k]);
            //fflush(NULL);

            for (int i=k+1;i<dim;i++)
            {
                A[i][k] = A[i][k] * A[k][k];
                x[i] = A[i][k];

                for (int j=k+1;j<dim;j++)
                {
                    A[i][j] = A[i][j] - A[i][k] * A[k][j];
                }
            }

            for (int i=k+1;i<dim;i++)
            {
                b[i]=b[i]-x[i]*b[k];
            }
            //printf("A[0][] = %8e %8e %8e | %8e\n", A[0][0],A[0][1],A[0][2],b[0]);
            //printf("A[1][] = %8e %8e %8e | %8e\n", A[1][0],A[1][1],A[1][2],b[1]);
            //printf("A[2][] = %8e %8e %8e | %8e\n\n", A[2][0],A[2][1],A[2][2],b[2]);
        }
    }


    /* backward substitution */
    x[dim-1]=b[dim-1]*A[dim-1][dim-1];

    for (int i=dim-2;i>=0;i--)
    {
        for (int j=dim-1;j>i;j--)
        {
            b[i]=b[i]-A[i][j]*x[j];
        }
        x[i]=b[i]*A[i][i];
    }

    //for (i=0;i<dim;i++)
    //    printf("%8e ",x[i]);
    //printf("\n");
    
    double det = 1.0;
    for(int i = 0 ; i < dim; i++)
        det = det * (1.0/A(i,i));
    //printf("matrix is singular A1 = %f, A2 = %f, A3 = %f, det = %f\n ", 1/A(0,0), 1/A(1,1), 1/A(2,2), fabs(det) );
    if(fabs(det) < TOL7 && order == 1)
    {
        solution = false;
        //printf("matrix is singular A1 = %f, A2 = %f, A3 = %f, det = %f\n ", 1/A(0,0), 1/A(1,1), 1/A(2,2), det );
    }
    
}
else
{
    if( fabs(A[0][0]) < TOL7)
    {
        printf("singular \n");
        solution = false;
    }
    x[0] = b[0]/A[0][0];
    
    printf("x = %f\n", x[0]);
}
    
 
   
    return solution;
}




#endif  // #ifdef CCADISCRET

