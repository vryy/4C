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
#include "intersection_math.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_element.H"
//#include <ctype.h>
//#include <errno.h>
//#include <stdio.h>
//#include <string.h>
//#include <limits.h>
//#include <stdlib.h>
//#include <assert.h>
//#include <math.h>

using namespace XFEM;
using namespace DRT::Utils;


static double           sqrarg;
#define                 SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


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


/*----------------------------------------------------------------------*
 |  GM:     checks if a certain element is a                 u.may 06/07|
 |          volume  element with help of the discretization type        |
 *----------------------------------------------------------------------*/  
bool XFEM::checkIfVolumeElement(
        const DRT::Element::DiscretizationType distype)
{
    bool isVolume = false;
    
    if( distype == DRT::Element::hex8  ||
        distype == DRT::Element::hex20 ||
        distype == DRT::Element::hex27 ||
        distype == DRT::Element::tet4  ||
        distype == DRT::Element::tet10  )
    {
        isVolume = true;        
    }
    return isVolume;
}


/*----------------------------------------------------------------------*
 |  GM:     checks if a certain element is a                 u.may 06/07|
 |          with help of the discretization type                        |
 *----------------------------------------------------------------------*/  
bool XFEM::checkIfSurfaceElement(
        const DRT::Element::DiscretizationType distype)
{
    bool isSurface = false;
    
    if( distype == DRT::Element::quad4 ||
        distype == DRT::Element::quad8 ||
        distype == DRT::Element::quad9 ||
        distype == DRT::Element::tri3  ||
        distype == DRT::Element::tri6  )
    {
        isSurface = true;       
    }
    return isSurface;
}


/*----------------------------------------------------------------------*
 |  GM:     checks if a certain element is a                 u.may 06/07|
 |          line element with help of the discretization type           |
 *----------------------------------------------------------------------*/  
bool XFEM::checkIfLineElement(
        const DRT::Element::DiscretizationType distype)
{
    bool isLine = false;
    
    if( distype == DRT::Element::line2 ||
        distype == DRT::Element::line3 )
    {
        isLine = true;      
    }
    return isLine;
}


/*----------------------------------------------------------------------*
 |  ICS:    checks if a node is within an XAABB               u.may 06/07|
 *----------------------------------------------------------------------*/
bool XFEM::isNodeWithinXAABB(    
    const std::vector<double>&         node,
    const Epetra_SerialDenseMatrix&    XAABB)
{
    bool isWithin = true;
    /*
    for (int dim=0; dim<3; dim++)
    {
        double diffMin = node[dim] - XAABB(dim,0);
        double diffMax = XAABB(dim,1) - node[dim];
        
        if((diffMin < -tol)||(diffMax < -tol)) //check again !!!!!      
            isWithin = false;
    }
    */
    for (int dim=0; dim<3; dim++)
    {
        double diffMin = XAABB(dim,0) - TOL7;
        double diffMax = XAABB(dim,1) + TOL7;
        
       // printf("nodal value =  %f, min =  %f, max =  %f\n", node[dim], diffMin, diffMax);
        
        if((node[dim] < diffMin)||(node[dim] > diffMax)) //check again !!!!!   
        {
            isWithin = false;
            break;
        }
    }
   
    return isWithin;
}


/*----------------------------------------------------------------------*
 |  ICS:    checks if a node is within an XAABB               u.may 06/07|
 *----------------------------------------------------------------------*/
bool XFEM::isLineWithinXAABB(    
    const std::vector<double>&         node1,
    const std::vector<double>&         node2,
    const Epetra_SerialDenseMatrix&    XAABB)
{
    bool isWithin = true;
    int dim = -1;
    
    for(dim=0; dim<3; dim++)
        if(fabs(node1[dim]-node2[dim]) > TOL7)
            break;
    
    for(int i = 0; i < 3; i++)
    {
        if(i != dim)
        {
            double min = XAABB(i,0) - TOL7;
            double max = XAABB(i,1) + TOL7;
   
            if((node1[i] < min)||(node1[i] > max))
                isWithin = false;
            
        }
        if(!isWithin)
            break;
    }
        
    if(isWithin && dim > -1)
    {
        isWithin = false;
        double min = XAABB(dim,0) - TOL7;
        double max = XAABB(dim,1) + TOL7;
                            
        if( ((node1[dim] < min) && (node2[dim] > max)) ||  
            ((node2[dim] < min) && (node1[dim] > max)) )
            isWithin = true;
    }
    return isWithin;
}


/*----------------------------------------------------------------------*
 |  CLI:    checks if a node is within a given element       u.may 06/07|   
 *----------------------------------------------------------------------*/
bool XFEM::checkNodeWithinElement(  
    DRT::Element*                       element,
    const Epetra_SerialDenseVector&     x)
{

    bool nodeWithinElement = true;
    int iter = 0;
    const int dim = getDimension(element);
    const int maxiter = 20;
    double residual = 1.0;
    
    Epetra_SerialDenseMatrix A(dim,dim);
    Epetra_SerialDenseVector b(dim);
    Epetra_SerialDenseVector dx(dim);
    Epetra_SerialDenseVector xsi(dim);
            
    updateRHSForNWE( dim, b, xsi, x, element);
   
    while(residual > TOL14)
    {   
        updateAForNWE( dim, A, xsi, element);
   
        if(!gaussElimination(A, b, dx, true, dim, 1))
        {
            nodeWithinElement = false;
            break;
        }   
        
        xsi = addTwoVectors(xsi,dx);
        updateRHSForNWE(dim, b, xsi, x, element);
        residual = b.Norm2();
        iter++; 
        
        if(iter >= maxiter)
        {   
            nodeWithinElement = false;
            break;
        }   
    }
    
    //printf("iter = %d\n", iter);
    //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi[0],xsi[1],xsi[2], residual, TOL14);
    
    for(int i=0; i<dim; i++)
        if( (fabs(xsi[i])-1.0) > TOL7)     
        {    
            nodeWithinElement = false;
            break;
        }
        
    return nodeWithinElement;
}


/*----------------------------------------------------------------------*
 |  CLI:    checks if a node is within a given element       u.may 06/07|   
 *----------------------------------------------------------------------*/
bool XFEM::checkNodeWithinElement(  
    DRT::Element*                       element,
    const Epetra_SerialDenseVector&     x,
    Epetra_SerialDenseVector&           xsi)
{

    bool nodeWithinElement = true;
    int iter = 0;
    const int dim = getDimension(element);
    const int maxiter = 20;
    double residual = 1.0;
    
    Epetra_SerialDenseMatrix A(dim,dim);
    Epetra_SerialDenseVector b(dim);
    Epetra_SerialDenseVector dx(dim);
    
    xsi.Scale(0.0);
            
    updateRHSForNWE( dim, b, xsi, x, element);
   
    while(residual > TOL14)
    {   
        updateAForNWE( dim, A, xsi, element);
   
        if(!gaussElimination(A, b, dx, true, dim, 1))
        {
            nodeWithinElement = false;
            break;
        }   
        
        xsi = addTwoVectors(xsi,dx);
        updateRHSForNWE(dim, b, xsi, x, element);
        residual = b.Norm2();
        iter++; 
        
        if(iter >= maxiter)
        {   
            nodeWithinElement = false;
            break;
        }   
    }
    
    //printf("iter = %d\n", iter);
    //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi[0],xsi[1],xsi[2], residual, TOL14);
    
    for(int i=0; i<dim; i++)
        if( (fabs(xsi[i])-1.0) > TOL7)     
        {    
            nodeWithinElement = false;
            break;
        }
        
    return nodeWithinElement;
}


/*----------------------------------------------------------------------*
 |  CLI:    updates the Jacobi matrix for the computation    u.may 06/07|
 |          if a node is in a given element                             |
 *----------------------------------------------------------------------*/
void XFEM::updateAForNWE(   
    const int                   dim,
    Epetra_SerialDenseMatrix&   A,
    Epetra_SerialDenseVector&   xsi,
    DRT::Element*               element)                                                  
{   
    const int numNodes = element->NumNode();
    Epetra_SerialDenseMatrix deriv1(dim, numNodes);
    
    A.Scale(0.0);
   
    switch(dim)
    {
        case 1:
        {
            shape_function_1D_deriv1(deriv1, xsi[0], element->Shape());
            break;
        }
        case 2:
        {
            shape_function_2D_deriv1(deriv1, xsi[0], xsi[1], element->Shape());
            break;
        }
        case 3:
        {
            shape_function_3D_deriv1(deriv1, xsi[0], xsi[1], xsi[2], element->Shape());
            break;
        }
        default:
            dserror("dimension of the element is not correct");
    }
    
    for(int j=0; j<numNodes; j++) 
    {
        DRT::Node* node = element->Nodes()[j];
        for(int i=0; i<dim; i++)
        {
            double nodalCoord = node->X()[i];
            for(int k=0; k<dim; k++)
                A[i][k] += nodalCoord * deriv1[j][k];
        }
    }      
}


/*----------------------------------------------------------------------*
 |  CLI:    updates the rhs for the computation if a         u.may 06/07|
 |          node is in a given element                                  |
 *----------------------------------------------------------------------*/
void XFEM::updateRHSForNWE( 
    const int                           dim,
    Epetra_SerialDenseVector&           b,
    Epetra_SerialDenseVector&           xsi,
    const Epetra_SerialDenseVector&     x,
    DRT::Element*                       element)                                                  
{
    const int numNodes = element->NumNode();
    Epetra_SerialDenseVector funct(numNodes);
      
    b.Scale(0.0);
     
    switch(dim)
    {
        case 1:
        {
            shape_function_1D(funct, xsi[0], element->Shape());
            break;
        }
        case 2:
        {
            shape_function_2D(funct, xsi[0], xsi[1], element->Shape());
            break;
        }
        case 3:
        {
            shape_function_3D(funct, xsi[0], xsi[1], xsi[2], element->Shape());
            break;
        }
        default:
            dserror("dimension of the element is not correct");
    }
    
    for(int j=0; j<numNodes; j++)
    {
        DRT::Node* node = element->Nodes()[j];
        for(int i=0; i<dim; i++)
            b[i] += (-1.0) * node->X()[i] * funct[j];
    }
      
    for(int i=0; i<dim; i++)
        b[i] += x[i];
}



#endif  // #ifdef CCADISCRET
