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

#ifdef CCADISCRET
#ifndef INTERSECTION_MATH_H_
#define INTERSECTION_MATH_H_

#include "intersection_service.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_element.H"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_SerialDenseSolver.h"


namespace XFEM
{


inline double XSIGN(const double a, const double b)
{
    return (b >= 0.0 ? fabs(a) : -fabs(a));
}

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
template <int n, int m>
static void svdcmp(
        blitz::TinyMatrix<double,m,n>&  A, 
        blitz::TinyVector<double,n>&  W, 
        blitz::TinyMatrix<double,n,n>&  V)
{
    static blitz::TinyVector<double,n> rv1;

    //Householder reduction to bidiagonal form.
    double g = 0.0;
    double scale = 0.0;
    double anorm = 0.0;
    for (int i=0;i<n;++i)
    {
        const int l=i+1;
        rv1(i)=scale*g;
        g=scale=0.0;
        if (i < m) 
        {
            for (int k=i;k<m;k++)
            {
                scale += fabs(A(k,i));
            }
            if (scale != 0.0)
            {
                double s = 0.0;
                for (int k=i;k<m;++k)
                {
                    A(k,i) /= scale;
                    s += pow(A(k,i), 2);
                }
                const double f=A(i,i);
                g = -XSIGN(sqrt(s),f);
                const double h=f*g-s;
                A(i,i)=f-g;
                for (int j=l;j<n;++j)
                {
                    double s = 0.0;
                    for (int k=i;k<m;++k)
                        s += A(k,i)*A(k,j);
                    const double f=s/h; 
                    for (int k=i;k<m;++k)
                        A(k,j) += f*A(k,i);
                }
                for (int k=i;k<m;++k)
                    A(k,i) *= scale;
            }
        }
        W(i)=scale*g;
        g=scale=0.0;
        if (i < m && i != (n-1))
        {
            for (int k=l;k<n;k++)
            {
                scale += fabs(A(i,k));
            }
            if (scale)
            {
                double s = 0.0;
                for (int k=l;k<n;k++)
                {
                    A(i,k) /= scale;
                    s += pow(A(i,k),2);
                }
                const double f=A(i,l);
                g = -XSIGN(sqrt(s),f);
                const double h=f*g-s;
                const double h_inv=1.0/h;
                A(i,l)=f-g;
                for (int k=l;k<n;k++)
                    rv1(k)=A(i,k)*h_inv;
                for (int j=l;j<m;++j)
                {
                    double s = 0.0;
                    for (int k=l;k<n;++k)
                        s += A(j,k)*A(i,k);
                    for (int k=l;k<n;++k)
                        A(j,k) += s*rv1(k);
                }
                for (int k=l;k<n;k++)
                    A(i,k) *= scale;
            }
        }
        anorm=std::max(anorm,(fabs(W(i))+fabs(rv1(i))));
    }
    //Accumulation of right-hand transformations.
    for (int i=(n-1);i>=0;i--)
    {
        const int l=i+1;
        if (i < n)
        {
            if (g)
            {
                //Double division to avoid possible underflow.
                const double g_inv = 1.0/g;
                for (int j=l;j<n;j++)
                    V(j,i)=(A(i,j)/A(i,l))*g_inv;
                for (int j=l;j<n;j++)
                {
                    double s = 0.0;
                    for (int k=l;k<n;k++) 
                        s += A(i,k)*V(k,j);
                    for (int k=l;k<n;k++) 
                        V(k,j) += s*V(k,i);
                }
            }
            for (int j=l;j<n;j++) 
                V(i,j)=V(j,i)=0.0;
        }
        V(i,i)=1.0;
        g=rv1(i);
    }
    //Accumulation of left-hand transformations.
    for (int i=std::min((m-1),(n-1));i>=0;i--)
    {
        const int l=i+1;
        const double g=W(i);
        for (int j=l;j<n;++j)
            A(i,j)=0.0;
        if (g)
        {
            const double g_inv=1.0/g;
            for (int j=l;j<n;++j)
            {
                double s = 0.0;
                for (int k=l;k<m;++k)
                    s += A(k,i)*A(k,j);
                const double f=(s/A(i,i))*g_inv;
                for (int k=i;k<m;++k)
                    A(k,j) += f*A(k,i);
            }
            for (int j=i;j<m;++j) 
                A(j,i) *= g_inv;
        } 
        else 
            for (int j=i;j<m;++j) 
                A(j,i)=0.0;
        ++A(i,i);
    }
    //Diagonalization of the bidiagonal form: Loop over
    for (int k=(n-1);k>=0;k--)
    {
        //singular values, and over allowed iterations.
        for (int its=1;its<=30;its++)
        {
            //Test for splitting.
            bool flag=true;
            int nm;
            int l = k;
            for (l=k;l>=0;l--)
            {
                //Note that rv1(1) is always zero.
                nm=l-1;
                if ((double)(fabs(rv1(l))+anorm) == anorm)
                {
                    flag=false;
                    break;
                }
                if ((double)(fabs(W(nm))+anorm) == anorm)
                    break;
            }
            if (flag)
            {
                //Cancellation of rv1(l), if l > 1.
                double c=0.0;
                double s=1.0;
                for (int i=l;i<=k;i++)
                {
                    const double f=s*rv1(i);
                    rv1(i) *= c;
                    if ((double)(fabs(f)+anorm) == anorm)
                        break;
                    const double g=W(i);
                    const double h=XFEM::pythagoras(f,g);
                    W(i)=h;
                    const double h_inv=1.0/h;
                    c=g*h_inv;
                    s=-f*h_inv;
                    for (int j=0;j<m;j++)
                    {
                        const double y=A(j,nm);
                        const double z=A(j,i);
                        A(j,nm)=y*c+z*s;
                        A(j,i)=z*c-y*s;
                    }
                }
            }
            const double z=W(k);
            //Convergence.
            if (l == k)
            {
                //Singular value is made nonnegative.
                if (z < 0.0)
                {
                    W(k) = -z;
                    for (int j=0;j<n;j++)
                        V(j,k) = -V(j,k);
                }
                break;
            }
            if (its == 30)
                dserror("no convergence in 30 svdcmp iterations");
            //Shift from bottom 2-by-2 minor.
            double x=W(l);
            const int km=k-1;
            const double y=W(km);
            double g=rv1(km);
            double h=rv1(k);
            double f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=XFEM::pythagoras(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+XSIGN(g,f)))-h))/x;
            //Next QR transformation:
            double c = 1.0;
            double s = 1.0;
            for (int j=l;j<=km;j++)
            {
                const int i=j+1;
                double g=rv1(i);
                double y=W(i);
                double h=s*g;
                g=c*g;
                double z=XFEM::pythagoras(f,h);
                rv1(j)=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g = g*c-x*s;
                h=y*s;
                y *= c;
                for (int jj=0;jj<n;jj++)
                {
                    const double x=V(jj,j);
                    const double z=V(jj,i);
                    V(jj,j)=x*c+z*s;
                    V(jj,i)=z*c-x*s;
                }
                z=XFEM::pythagoras(f,h);
                //Rotation can be arbitrary if z = 0.
                W(j)=z;
                if (z)
                {
                    const double z_inv=1.0/z;
                    c=f*z_inv;
                    s=h*z_inv;
                }
                f=c*g+s*y;
                x=c*y-s*g;
                for (int jj=0;jj<m;jj++)
                {
                    const double y=A(jj,j);
                    const double z=A(jj,i);
                    A(jj,j)=y*c+z*s;
                    A(jj,i)=z*c-y*s;
                }
            }
            rv1(l)=0.0;
            rv1(k)=f;
            W(k)=x;
        }
    }
}
        
     
    /*!
    \brief  Solves a linear system of equations with help of
            a singular value decomposition. The system matrix
            may be singular !!
                   
    \tparam dim  (in)            : dimension of the square system matrix
    */
    template<int dim>
    inline bool solveLinearSystemWithSVD(
        blitz::TinyMatrix<double,dim,dim>&    U,  ///< (in/out)        : system matrix / matrix U
        const blitz::TinyVector<double,dim>&  b,  ///< (in)            : right-hand-side
        blitz::TinyVector<double,dim>&        x,  ///< (out)           : solution vector
        const double                          tol ///< (in)            : tolerance
        )
    {
        static blitz::TinyMatrix<double,dim,dim> V;
        static blitz::TinyMatrix<double,dim,dim> A;
        static blitz::TinyVector<double,dim>     W;
        
        // initialize vectors and matrices
        W = 0.0;
        V = 0.0;
        
        A = U;
        
        svdcmp<dim,dim>( U, W, V);
        
        //test_svdcmp( A, U, W, V, 3);
           
        x = 0.0;
        for(int i = 0; i < dim; i++ )
        {
            for(int  k = 0; k < dim; k++ )
            {
                double svdtemp = 0.0;
                for(int  j = 0; j < dim; j++ )
                {
                    if( fabs(W(j)) > tol )
                        svdtemp += V(i,j)*U(k,j) / W(j);
                }
                x(i) += svdtemp * b(k);
            }
        }
        
        bool nonsingular = true;
        for(int  j = 0; j < dim; j++ )
            if( fabs(W(j)) <= tol )
            {
                 nonsingular = false;
                 break;
            }
        
       /*  
       BlitzVec b1(3);
       for(int i = 0; i < dim; i++ )
        {
            b1(i) = 0.0;
            for(int  k = 0; k < dim; k++ )
            {  
                b1(i) += A(i,k)*x(k);
            }
            printf("b = %f\t  b1 = %f\n", b(i), b1(i));
        }
        printf("\n");
        
        for(int i = 0; i < dim; i++ )
            printf("x = %f\t", x(i));
            
        printf("\n");
        */
        return nonsingular;
    }
     
      
    /*!
    \brief  Tests the singular value decomposition
                   
    \param A    (in)        : system matrix 
    \param U    (in)        : decomposed matrix U
    \param W    (in)        : diagonal matrix W stored in  vector
    \param V    (in)        : decomposed matrix V
    \param dim  (in)        : dimension
    */            
    void test_svdcmp(
            BlitzMat&   A,
            BlitzMat&   U,
            BlitzVec&   W,
            BlitzMat&   V,
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
                        blitz::TinyVector<double,4> tmp;    // check changed
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

        // view on hopefullx column major blitz arrays
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
#endif  // #ifdef CCADISCRET


