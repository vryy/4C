/*!----------------------------------------------------------------------
\file interface_service.cpp

\brief collection of math tools for the interface determination of trv1o meshes

    ML      math library for the interface computation
 
    
<pre>
Maintainer: Ursula Mayer
</pre>
*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "intersection_service.H"
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
using namespace DRT::Utils;


static double           sqrarg;
#define                 SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define                 SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static double           maxarg1,maxarg2;
#define                 FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
                        (maxarg1) : (maxarg2))
static int              iminarg1,iminarg2;
#define                 IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
                        (iminarg1) : (iminarg2))


/*----------------------------------------------------------------------*
 |  ML:     adds two Epetra_SerialDenseVector                u.may 06/07|
 *----------------------------------------------------------------------*/
Epetra_SerialDenseVector XFEM::addTwoVectors(   
    const Epetra_SerialDenseVector&   v1,
    const Epetra_SerialDenseVector&   v2)
{   
    Epetra_SerialDenseVector vResult(v1.Length());
    
    if(v1.Length() != v2.Length())
        dserror("both vectors need to have the same size\n"); 

    for(int i = 0; i < v1.Length(); i++)
        vResult[i] = v1[i] + v2[i];
 
    return vResult;
}
    
    
   
/*----------------------------------------------------------------------*
 |  ML:     adds two vector<double>                          u.may 06/07|
 *----------------------------------------------------------------------*/
vector<double> XFEM::addTwoVectors(
    const vector<double>&   v1,
    const vector<double>&   v2)
{   
    vector<double> vResult(v1.size());
    
    if(v1.size() != v2.size())
        dserror("both vectors need to have the same size\n"); 

    for(unsigned int i = 0; i < v1.size(); i++)
        vResult[i] = v1[i] + v2[i];
 
    return vResult;
}
    
    

/*----------------------------------------------------------------------*
 |  ML:     subtracts one Epetra_SerialDenseVector from   u.may 06/07   |
 |          another Epetra_SerialDenseVector.                           |
 |          The result is stored in v1                                  |
 *----------------------------------------------------------------------*/
Epetra_SerialDenseVector XFEM::subtractsTwoVectors( 
    const Epetra_SerialDenseVector& v1,
    const Epetra_SerialDenseVector& v2)
{   
    
    Epetra_SerialDenseVector vResult(v1.Length());
    
    if(v1.Length() != v2.Length())
        dserror("both vectors need to have the same size\n"); 

    for(int i = 0; i < v1.Length(); i++)
        vResult[i] = v1[i] - v2[i];
 
    return vResult;
}

    

/*----------------------------------------------------------------------*
 |  ML :    subtracts one vector<double> from another        u.may 06/07|
 |          vector<double> . The result is stored in v1.                |
 *----------------------------------------------------------------------*/
vector<double> XFEM::subtractsTwoVectors(   
    const vector <double>& v1,
    const vector <double>& v2)
{   
    vector <double>  vResult(v1.size());
    
    if(v1.size() != v2.size())
        dserror("both vectors need to have the same size\n"); 

    for(unsigned int i = 0; i < v1.size(); i++)
        vResult[i] = v1[i] - v2[i];
 
    return vResult;
}



/*----------------------------------------------------------------------*
 |  ML:     computes the cross product                       u.may 08/07|
 |          of 2 Epetra_SerialDenseVector c = a x b                     |
 *----------------------------------------------------------------------*/  
Epetra_SerialDenseVector XFEM::computeCrossProduct(
    const Epetra_SerialDenseVector& a,
    const Epetra_SerialDenseVector& b)
{
    Epetra_SerialDenseVector c(3);
   
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    
    return c;
}



/*----------------------------------------------------------------------*
 |  ML:     normalizes a Epetra_SerialDenseVector            u.may 08/07|
 *----------------------------------------------------------------------*/  
void XFEM::normalizeVector(   
    Epetra_SerialDenseVector&     v)
{
    const double norm = v.Norm2();
    v.Scale(1.0/norm);
    return;
}



/*----------------------------------------------------------------------*
 |    theorem of pythagoras                                 u.may 09/07 |
 |    computes ( a^2 + b^2 ) ^ (1/2)                                    |
 |    (modified from NUMERICAL RECIPES)                                 |
 *----------------------------------------------------------------------*/  
double XFEM::pythagoras(
    const double  a, 
    const double  b)
{
    const double absa=fabs(a);
    const double absb=fabs(b);
    if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
    else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}


#endif  // #ifdef CCADISCRET
