/*!----------------------------------------------------------------------
\file intersection_service.cpp

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

using namespace XFEM;
using namespace DRT::UTILS;


/*----------------------------------------------------------------------*
 |  ML:     computes the cross product                       u.may 08/07|
 |          of 2 vectors c = a x b                                      |
 *----------------------------------------------------------------------*/  
BlitzVec XFEM::computeCrossProduct(
    const BlitzVec& a,
    const BlitzVec& b)
{
    BlitzVec c(3);
   
    c(0) = a(1)*b(2) - a(2)*b(1);
    c(1) = a(2)*b(0) - a(0)*b(2);
    c(2) = a(0)*b(1) - a(1)*b(0);
    
    return c;
}


///*----------------------------------------------------------------------*
// |  GM:     checks if a certain element is a                 u.may 06/07|
// |          volume  element with help of the discretization type        |
// *----------------------------------------------------------------------*/  
//bool XFEM::checkIfVolumeElement(
//        const DRT::Element::DiscretizationType distype)
//{
//    bool isVolume = false;
//    
//    if( distype == DRT::Element::hex8  ||
//        distype == DRT::Element::hex20 ||
//        distype == DRT::Element::hex27 ||
//        distype == DRT::Element::tet4  ||
//        distype == DRT::Element::tet10  )
//    {
//        isVolume = true;        
//    }
//    return isVolume;
//}


///*----------------------------------------------------------------------*
// |  GM:     checks if a certain element is a                 u.may 06/07|
// |          with help of the discretization type                        |
// *----------------------------------------------------------------------*/  
//bool XFEM::checkIfSurfaceElement(
//        const DRT::Element::DiscretizationType distype)
//{
//    bool isSurface = false;
//    
//    if( distype == DRT::Element::quad4 ||
//        distype == DRT::Element::quad8 ||
//        distype == DRT::Element::quad9 ||
//        distype == DRT::Element::tri3  ||
//        distype == DRT::Element::tri6  )
//    {
//        isSurface = true;       
//    }
//    return isSurface;
//}


///*----------------------------------------------------------------------*
// |  GM:     checks if a certain element is a                 u.may 06/07|
// |          line element with help of the discretization type           |
// *----------------------------------------------------------------------*/  
//bool XFEM::checkIfLineElement(
//        const DRT::Element::DiscretizationType distype)
//{
//    bool isLine = false;
//    
//    if( distype == DRT::Element::line2 ||
//        distype == DRT::Element::line3 )
//    {
//        isLine = true;      
//    }
//    return isLine;
//}



/*----------------------------------------------------------------------*
 | GM:      transforms a node in element coordinates       u.may 07/07|
 |          into current coordinates                                    |
 *----------------------------------------------------------------------*/  
BlitzVec XFEM::elementToCurrentCoordinates(   
        const DRT::Element* element,
        const BlitzVec&     eleCoord) 
{
    const int numNodes = element->NumNode();
    BlitzVec funct(numNodes);
    BlitzVec physCoord(3);
    
    switch(getDimension(element->Shape()))
    {
        case 1:
        {
            shape_function_1D(funct,eleCoord(0), element->Shape());
            break;
        }
        case 2:
        {
            shape_function_2D(funct, eleCoord(0), eleCoord(1), element->Shape());
            break;
        }
        case 3:
        {
            shape_function_3D(funct, eleCoord(0), eleCoord(1), eleCoord(2), element->Shape());
            break;
        }
        default:
            dserror("dimension of the element is not correct");
    }          
        
    physCoord = 0.0;
    
    for(int i=0; i<numNodes; i++)
    {
        const DRT::Node* node = element->Nodes()[i];
        for(int j=0; j<3; j++)
        	physCoord(j) += node->X()[j] * funct(i);
    }
    
    return physCoord;
}


/*----------------------------------------------------------------------*
 | GM:      transforms a node in element coordinates       u.may 07/07|
 |          into current coordinates                                    |
 *----------------------------------------------------------------------*/  
void XFEM::elementToCurrentCoordinatesInPlace(   
        const DRT::Element* element,
        BlitzVec&           eleCoord) 
{
    const int numNodes = element->NumNode();
    BlitzVec funct(numNodes);
    dsassert(eleCoord.size() == 3, "inplace coordinate transfer only in 3d!");
    
    switch(getDimension(element->Shape()))
    {
        case 1:
        {
            shape_function_1D(funct,eleCoord(0), element->Shape());
            break;
        }
        case 2:
        {
            shape_function_2D(funct, eleCoord(0), eleCoord(1), element->Shape());
            break;
        }
        case 3:
        {
            shape_function_3D(funct, eleCoord(0), eleCoord(1), eleCoord(2), element->Shape());
            break;
        }
        default:
            dserror("dimension of the element is not correct");
    }          
        
    eleCoord = 0.0;
    for(int i=0; i<numNodes; i++)
    {
        const DRT::Node* node = element->Nodes()[i];
        for(int j=0; j<3; j++)
            eleCoord(j) += node->X()[j] * funct(i);
    }
}

/*----------------------------------------------------------------------*
 | GM:  transforms a node in current coordinates            u.may 07/07 |
 |      into element coordinates                                        |
 *----------------------------------------------------------------------*/  
BlitzVec XFEM::currentToElementCoordinatesExact(   
    const DRT::Element*               element, 
    const BlitzVec&                   x) 
{
    dsassert(x.size() == 3, "current coordinates have to be 3!");
    const int dim = DRT::UTILS::getDimension(element->Shape());
    
    BlitzVec xsi(dim);
    currentToElementCoordinates(element, x, xsi);
   
    // rounding 1 and -1 to be exact for the CDT
    for(int j = 0; j < 3; j++)
    {
        if( fabs((fabs(xsi(j))-1.0)) < TOL7 &&  xsi(j) < 0)    xsi(j) = -1.0;
        if( fabs((fabs(xsi(j))-1.0)) < TOL7 &&  xsi(j) > 0)    xsi(j) =  1.0;      
    }
    return xsi;
}     


/*!
\brief updates the system matrix at the corresponding element coordinates for the 
       computation if a node in current coordinates lies within an element 

\param dim               (in)       : dimension of the problem
\param A                 (out)      : system matrix
\param xsi               (in)       : vector of element coordinates
\param element           (in)       : element 
*/  
template <DRT::Element::DiscretizationType DISTYPE>
static inline void updateAForNWE(
    BlitzMat&                           A,
    const BlitzVec&                     xsi,
    const BlitzMat&                     xyze
    )                                                  
{   
    const int numNodes = DRT::UTILS::getNumberOfElementNodes<DISTYPE>();
    const int dim = DRT::UTILS::getDimension<DISTYPE>();
    static BlitzMat deriv1(dim,numNodes);
    shape_function_deriv1<DISTYPE,BlitzVec>(xsi, deriv1);
    
    blitz::firstIndex isd;
    blitz::secondIndex jsd;
    blitz::thirdIndex inode;
    A = blitz::sum(xyze(isd,inode) * deriv1(jsd,inode), inode);
}



/*!
\brief updates the rhs at the corresponding element coordinates for the 
       computation whether a node in current coordinates lies within an element 

\param dim               (in)       : dimension of the problem
\param b                 (out)      : right-hand-side       
\param xsi               (in)       : vector of element coordinates
\param x                 (in)       : node in current coordinates
\param element           (in)       : element
*/
template <DRT::Element::DiscretizationType DISTYPE>
static inline void updateRHSForNWE(
        BlitzVec&           b,
        const BlitzVec&     xsi,
        const BlitzVec&     x,
        const BlitzMat&     xyze)                                                  
{
    const int numNodes = DRT::UTILS::getNumberOfElementNodes<DISTYPE>();
    static BlitzVec funct(numNodes);
    shape_function<DISTYPE,BlitzVec>(xsi, funct);
    
    blitz::firstIndex isd;
    blitz::secondIndex inode;
    b = x(isd) - blitz::sum(xyze(isd,inode) * funct(inode), inode);
}


/*!
\brief transforms a point in current coordinates to a point
       in element coordinates with respect to a given element       
       The nonlinear system of equation is solved with help of the Newton-method.
       Fast templated version

\param element              (in)        : element 
\param x                    (in)        : node in current coordinates (x, y, z)
\param xsi                  (inout)     : node in element coordinates
*/  
template <DRT::Element::DiscretizationType DISTYPE, int dim>
static inline bool currentToElementCoordinatesT(
    const DRT::Element*                 element,
    const BlitzVec&                     x,
    BlitzVec&                           xsi)
{
    dsassert(element->Shape() == DISTYPE, "this is a bug in currentToElementCoordinatesT!");
    bool nodeWithinElement = true;
    //const int dim = DRT::UTILS::getDimension<DISTYPE>();
    const int maxiter = 20;
    double residual = 1.0;
    
    static BlitzMat A(dim,dim);
    static BlitzVec b(dim);
    static BlitzVec dx(dim);
    
    static BlitzMat xyze(3, DRT::UTILS::getNumberOfElementNodes<DISTYPE>());
    fillPositionArrayBlitzT<DISTYPE>(element, xyze);
    
    // initial guess
    xsi = 0.0;
            
    updateRHSForNWE<DISTYPE>(b, xsi, x, xyze);
    
    int iter = 0;
    while(residual > TOL14)
    {   
        updateAForNWE<DISTYPE>( A, xsi, xyze);
   
        if(!gaussElimination<true, dim, 1>(A, b, dx))
        {
            nodeWithinElement = false;
            break;
        }   
        
        xsi += dx;
        updateRHSForNWE<DISTYPE>(b, xsi, x, xyze);
        
        residual = Norm2(b);
        iter++; 
        
        if(iter >= maxiter)
        {   
            nodeWithinElement = false;
            break;
        }   
    }
    return nodeWithinElement;
}

/*----------------------------------------------------------------------*
 | GM:  transforms a node in current coordinates            u.may 12/07 |
 |      into element coordinates                                        | 
 *----------------------------------------------------------------------*/
bool XFEM::currentToElementCoordinates(  
    const DRT::Element*                 element,
    const BlitzVec&                     x,
    BlitzVec&                           xsi)
{
    bool nodeWithinElement = false;
    switch (element->Shape())
    {
    case DRT::Element::hex8:
        nodeWithinElement = currentToElementCoordinatesT<DRT::Element::hex8, 3>(element, x, xsi);
        break;
    case DRT::Element::hex20:
        nodeWithinElement = currentToElementCoordinatesT<DRT::Element::hex20, 3>(element, x, xsi);
        break;
    case DRT::Element::hex27:
        nodeWithinElement = currentToElementCoordinatesT<DRT::Element::hex27, 3>(element, x, xsi);
        break;
    case DRT::Element::line2:
        nodeWithinElement = currentToElementCoordinatesT<DRT::Element::line2, 1>(element, x, xsi);
        break;
    case DRT::Element::line3:
        nodeWithinElement = currentToElementCoordinatesT<DRT::Element::line3, 1>(element, x, xsi);
        break;
    default:
        dserror("add your distype to this switch!");
        nodeWithinElement = false;
    }   
    //printf("iter = %d\n", iter);
    //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi(0),xsi(1),xsi(2), residual, TOL14);
    return nodeWithinElement;
}


/*----------------------------------------------------------------------*
 |  ICS:    checks if a position is within an XAABB          u.may 06/07|
 *----------------------------------------------------------------------*/
bool XFEM::isPositionWithinXAABB(    
    const BlitzVec&                    pos,
    const BlitzMat&                    XAABB)
{
    const int nsd = 3;
    bool isWithin = true;
    for (int isd=0; isd<nsd; isd++)
    {
        const double diffMin = XAABB(isd,0) - TOL7;
        const double diffMax = XAABB(isd,1) + TOL7;
        
       // printf("nodal value =  %f, min =  %f, max =  %f\n", node[dim], diffMin, diffMax);
        
        if((pos(isd) < diffMin)||(pos(isd) > diffMax)) //check again !!!!!   
        {
            isWithin = false;
            break;
        }
    }
   
    return isWithin;
}


/*----------------------------------------------------------------------*
 |  ICS:    checks if a pos is within an XAABB               u.may 06/07|
 *----------------------------------------------------------------------*/
bool XFEM::isLineWithinXAABB(    
    const BlitzVec&                    pos1,
    const BlitzVec&                    pos2,
    const BlitzMat&                    XAABB)
{
    const int nsd = 3;
    bool isWithin = true;
    int isd = -1;
    
    for(isd=0; isd<nsd; isd++)
        if(fabs(pos1(isd)-pos2(isd)) > TOL7)
            break;
    
    for(int ksd = 0; ksd < nsd; ksd++)
    {
        if(ksd != isd)
        {
            const double min = XAABB(ksd,0) - TOL7;
            const double max = XAABB(ksd,1) + TOL7;
   
            if((pos1(ksd) < min)||(pos1(ksd) > max))
                isWithin = false;
            
        }
        if(!isWithin)
            break;
    }
        
    if(isWithin && isd > -1)
    {
        isWithin = false;
        const double min = XAABB(isd,0) - TOL7;
        const double max = XAABB(isd,1) + TOL7;
                            
        if( ((pos1(isd) < min) && (pos2(isd) > max)) ||  
            ((pos2(isd) < min) && (pos1(isd) > max)) )
            isWithin = true;
    }
    return isWithin;
}


/*----------------------------------------------------------------------*
 |  CLI:    checks if a position is within a given element   u.may 06/07|   
 *----------------------------------------------------------------------*/
bool XFEM::checkPositionWithinElement(  
    const DRT::Element*                 element,
    const BlitzVec&                     x)
{
    BlitzVec xsi(getDimension(element->Shape()));
    bool nodeWithinElement = currentToElementCoordinates(element, x, xsi);
    //printf("iter = %d\n", iter);
    //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi(0),xsi(1),xsi(2), residual, TOL14);
    
    nodeWithinElement = checkPositionWithinElementParameterSpace(xsi,element->Shape());
//    for(int i=0; i<dim; i++)
//        if( (fabs(xsi(i))-1.0) > TOL7)     
//        {    
//            nodeWithinElement = false;
//            break;
//        }
        
    return nodeWithinElement;
}


/*----------------------------------------------------------------------*
 |  CLI:    checks if a position is within a given mesh      a.ger 12/07|   
 *----------------------------------------------------------------------*/
bool XFEM::PositionWithinDiscretization(  
    const RCP<DRT::Discretization>      dis,
    const BlitzVec&                     x)
{
    bool nodeWithinMesh = false;
    
    // loop all elements on this processor
    for (int i=0; i<dis->NumMyRowElements(); ++i)
    {
        const DRT::Element* ele = dis->lRowElement(i);
        if (isPositionWithinXAABB(x, computeFastXAABB(ele)))
        {
            nodeWithinMesh = checkPositionWithinElement(ele, x);
            if (nodeWithinMesh)
                break;
        }
    }

    // TODO: in parallel, we have to ask all processors, whether there is any match!!!!
#ifdef PARALLEL
    dserror("not implemented, yet");
#endif
    return nodeWithinMesh;
}

/*----------------------------------------------------------------------*
 |  CLI:    checks if a position is within condition-enclosed region      a.ger 12/07|   
 *----------------------------------------------------------------------*/
bool XFEM::PositionWithinCondition(
        const BlitzVec&                     x_in,
        const int                           xfem_condition_label, 
        const RCP<DRT::Discretization>      cutterdis
    )
{
    const int nsd = 3;
    static BlitzVec x(nsd);
//    for (int isd = 0; isd < nsd; ++isd) {
//        x(isd) = x_in(isd);
//    }
    x = x_in;

    // TODO: use label to identify the surface/xfem condition
    
    bool nodeWithinMesh = false;
    // loop all elements on this processor
    for (int i=0; i<cutterdis->NumMyRowElements(); ++i)
    {
        const DRT::Element* ele = cutterdis->lRowElement(i);
        if (isPositionWithinXAABB(x, computeFastXAABB(ele)))
        {
            nodeWithinMesh = checkPositionWithinElement(ele, x);
            if (nodeWithinMesh)
                break;
        }
    }

    // TODO: in parallel, we have to ask all processors, whether there is any match!!!!
#ifdef PARALLEL
    dserror("not implemented, yet");
#endif
    return nodeWithinMesh;
}


/*----------------------------------------------------------------------*
 |  RQI:    searches the nearest point on a surface          u.may 02/08|
 |          element for a given point in physical coordinates           |
 *----------------------------------------------------------------------*/
bool XFEM::searchForNearestPointOnSurface(
    const DRT::Element*                   	surfaceElement,
    const BlitzVec&     	            	physCoord,
    BlitzVec&                       		eleCoord,
    BlitzVec&           	              	normal,
    double&									distance)
{
	distance = -1.0;
	normal = 0;
	
	eleCoord = CurrentToSurfaceElementCoordinates(surfaceElement, physCoord);
	
	const bool pointWithinElement = checkPositionWithinElementParameterSpace(eleCoord, surfaceElement->Shape());
	
	if(pointWithinElement)
	{
		const BlitzVec x_surface_phys(elementToCurrentCoordinates(surfaceElement, eleCoord));
		normal = x_surface_phys - physCoord;
		distance = sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2));
	}

	return pointWithinElement;
}



/*----------------------------------------------------------------------*
 |  RQI:    compute element coordinates from a given point              |
 |          in the 3-dim physical space lies on a given surface element |
 *----------------------------------------------------------------------*/
BlitzVec XFEM::CurrentToSurfaceElementCoordinates(
        const DRT::Element*        	surfaceElement,
        const BlitzVec&     		physCoord
   	    )
{
	bool nodeWithinElement = true;
    
    BlitzVec eleCoord(2);
    eleCoord = 0.0;
    
    blitz::firstIndex i;    // Placeholder for the first blitz array index
    blitz::secondIndex j;   // Placeholder for the second blitz array index
   								
  	const int maxiter = 20;
  	int iter = 0;
  	while(iter < maxiter)
    {   
  	    iter++;
  	    
        // compute Jacobian, f and b
  	    const BlitzMat Jacobi(updateJacobianForMap3To2(eleCoord, surfaceElement));
        const BlitzVec F(updateFForMap3To2(eleCoord, physCoord, surfaceElement));
        BlitzVec b(blitz::sum(-Jacobi(j,i)*F(j),j));
     
        const double residual = sqrt(b(0)*b(0)+b(1)*b(1));
        if (residual < TOL14)
        {   
            nodeWithinElement = true;
            break;
        }  
  	    
        // compute system matrix A
        BlitzMat A(2,2,blitz::ColumnMajorArray<2>());
  		updateAForMap3To2(A, Jacobi, F, eleCoord, surfaceElement);
  	
  		BlitzVec dx(2);
  		dx = 0;
        if(!gaussEliminationEpetra(A, b, dx))
        {
            nodeWithinElement = false;
            break;
        }   
           
        eleCoord += dx;
    }
  
    //printf("iter = %d\n", iter);
    //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi[0],xsi[1],xsi[2], residual, TOL14);
	return eleCoord;
}


/*
 * checks if a position in element coordinates lies within a certain surfaceElement
 */  
bool XFEM::checkPositionWithinElementParameterSpace(
        const BlitzVec&                        eleCoord,
        const DRT::Element::DiscretizationType distype
        )
{
    if (distype != DRT::Element::line2 and
            distype != DRT::Element::line3 and
            distype != DRT::Element::quad4 and 
            distype != DRT::Element::quad8 and 
            distype != DRT::Element::quad9 and
            distype != DRT::Element::hex8 and
            distype != DRT::Element::hex20 and 
            distype != DRT::Element::hex27)
        dserror("function only defined for rectangular element types at the moment");
    
    bool nodeWithinElement = true;
    
    // loop over r and s (local coordinates)
    for(int i=0; i<DRT::UTILS::getDimension(distype); i++)
        if( (fabs(eleCoord(i))-1.0) > TOL7)     
        {    
            nodeWithinElement = false;
            break;
        }
    
    return nodeWithinElement;
}


/*----------------------------------------------------------------------*
 |  RQI:    updates the Jacobian for the computation         u.may 01/08|
 |          whether a point in the 3-dim physical space lies            |
 | 			on a surface element 										|                 
 *----------------------------------------------------------------------*/
BlitzMat XFEM::updateJacobianForMap3To2(   
    const BlitzVec&	                xsi,
    const DRT::Element*             surfaceElement)                                                  
{   
    BlitzMat Jacobi(3,2,blitz::ColumnMajorArray<2>());
    Jacobi = 0;
    
    const int numNodes = surfaceElement->NumNode();
    const BlitzMat deriv1(shape_function_2D_deriv1(xsi(0), xsi(1), surfaceElement->Shape()));
    for(int inode=0; inode<numNodes; inode++) 
    {
        const double* x = surfaceElement->Nodes()[inode]->X();
        for(int isd=0; isd<3; isd++)
        {
            for(int jsd=0; jsd<2; jsd++)
            {
            	Jacobi(isd, jsd) += x[isd] * deriv1(jsd,inode);
            }
        }
    } 
    return Jacobi;
}



/*----------------------------------------------------------------------*
 |  RQI:    updates the system of nonlinear equations        u.may 01/08|
 |          for the computation, whether a point in the 3-dim physical 	|
 |			space lies on a surface element 							|                 
 *----------------------------------------------------------------------*/
BlitzVec XFEM::updateFForMap3To2(   
        const BlitzVec&	    xsi,
        const BlitzVec&	    x,
        const DRT::Element* surfaceElement)                                                  
{   
    BlitzVec F(3);
    F = 0.0;
    
    const int numNodes = surfaceElement->NumNode();
    const BlitzVec funct(shape_function_2D(xsi(0), xsi(1), surfaceElement->Shape()));
    for(int inode=0; inode<numNodes; inode++) 
    {
        const double* coord = surfaceElement->Nodes()[inode]->X();
        for(int isd=0; isd<3; ++isd)
     		F(isd) += coord[isd] * funct(inode);
    }   
    
    F -= x;
    return F;
}



/*----------------------------------------------------------------------*
 |  RQI:    updates the system matrix			             u.may 01/08|
 |          for the computation, whether a point in the 3-dim physical 	|
 |			space lies on a surface element 							|                 
 *----------------------------------------------------------------------*/
void XFEM::updateAForMap3To2(   
        BlitzMat& 		            A,
        const BlitzMat& 	        Jacobi,
        const BlitzVec&             F,
        const BlitzVec&	            xsi,
        const DRT::Element*         surfaceElement
        )                                                  
{   
    const int numNodes = surfaceElement->NumNode();
    blitz::firstIndex i;    // Placeholder for the first blitz array index
    blitz::secondIndex j;   // Placeholder for the second blitz array index
    blitz::thirdIndex k;    // Placeholder for the third blitz array index
    
    const BlitzMat deriv2(shape_function_2D_deriv2(xsi(0), xsi(1), surfaceElement->Shape()));
	blitz::Array<double, 3> tensor3Ord(3,2,2, blitz::ColumnMajorArray<3>());		
    tensor3Ord = 0.0;
   
	for(int inode=0; inode<numNodes; inode++) 	
	{
		const double* x = surfaceElement->Nodes()[inode]->X();
		for(int isd=0; isd<3; ++isd)
		{
			const double nodalCoord = x[isd];
			for(int jsd=0; jsd<2; ++jsd)
				for(int ksd=0; ksd<2; ++ksd)
					tensor3Ord(isd, jsd, ksd) += nodalCoord * deriv2(jsd,inode);
		}
	}
	
	A = blitz::sum(Jacobi(k,i) * Jacobi(k,j), k) + blitz::sum(F(k)*tensor3Ord(k,i,j),k);
}



/*----------------------------------------------------------------------*
 |  ICS:    computes an extended axis-aligned bounding box   u.may 06/07|
 |          XAABB for a given element                                   |
 *----------------------------------------------------------------------*/
BlitzMat XFEM::computeFastXAABB( 
    const DRT::Element* element)
{
    const int nsd = 3;
    BlitzMat XAABB(nsd, 2);
    
    // first node
    const double* pos = element->Nodes()[0]->X();
    for(int dim=0; dim<nsd; ++dim)
    {
        XAABB(dim, 0) = pos[dim] - TOL7;
        XAABB(dim, 1) = pos[dim] + TOL7;
    }
    // remaining node
    for(int i=1; i<element->NumNode(); ++i)
    {
        const double* posEle = element->Nodes()[i]->X();
        for(int dim=0; dim<nsd; dim++)
        {
            XAABB(dim, 0) = std::min( XAABB(dim, 0), posEle[dim] - TOL7);
            XAABB(dim, 1) = std::max( XAABB(dim, 1), posEle[dim] + TOL7);
        }
    }
    
    double maxDistance = fabs(XAABB(0,1) - XAABB(0,0));
    for(int dim=1; dim<nsd; ++dim)
       maxDistance = std::max(maxDistance, fabs(XAABB(dim,1)-XAABB(dim,0)) );
    
    // subtracts half of the maximal distance to minX, minY, minZ
    // adds half of the maximal distance to maxX, maxY, maxZ 
    const double halfMaxDistance = 0.5*maxDistance;
    for(int dim=0; dim<nsd; ++dim)
    {
        XAABB(dim, 0) -= halfMaxDistance;
        XAABB(dim, 1) += halfMaxDistance;
    }   
    
    /*
    printf("\n");
    printf("axis-aligned bounding box:\n minX = %f\n minY = %f\n minZ = %f\n maxX = %f\n maxY = %f\n maxZ = %f\n", 
              XAABB(0,0), XAABB(1,0), XAABB(2,0), XAABB(0,1), XAABB(1,1), XAABB(2,1));
    printf("\n");
    */
    
    return XAABB;
}


/*----------------------------------------------------------------------*
 |  ICS:    checks if two XAABB's intersect                  u.may 06/07|
 *----------------------------------------------------------------------*/
bool XFEM::intersectionOfXAABB(  
    const BlitzMat&     cutterXAABB, 
    const BlitzMat&     xfemXAABB)
{
    
  /*====================================================================*/
  /* bounding box topology*/
  /*--------------------------------------------------------------------*/
  /* parameter coordinates (x,y,z) of nodes
   * node 0: (minX, minY, minZ)
   * node 1: (maxX, minY, minZ)
   * node 2: (maxX, maxY, minZ)
   * node 3: (minX, maxY, minZ)
   * node 4: (minX, minY, maxZ)
   * node 5: (maxX, minY, maxZ)
   * node 6: (maxX, maxY, maxZ)
   * node 7: (minX, maxY, maxZ)
   * 
   *                      z
   *                      |           
   *             4========|================7
   *           //|        |               /||
   *          // |        |              //||
   *         //  |        |             // ||
   *        //   |        |            //  ||
   *       //    |        |           //   ||
   *      //     |        |          //    ||
   *     //      |        |         //     ||
   *     5=========================6       ||
   *    ||       |        |        ||      ||
   *    ||       |        o--------||---------y
   *    ||       |       /         ||      ||
   *    ||       0------/----------||------3
   *    ||      /      /           ||     //
   *    ||     /      /            ||    //
   *    ||    /      /             ||   //
   *    ||   /      /              ||  //
   *    ||  /      /               || //
   *    || /      x                ||//
   *    ||/                        ||/
   *     1=========================2
   *
   */
  /*====================================================================*/
    
    bool intersection =  false;
    std::vector < BlitzVec > nodes(8, BlitzVec(3));
    
    nodes[0](0) = cutterXAABB(0,0); nodes[0](1) = cutterXAABB(1,0); nodes[0](2) = cutterXAABB(2,0); // node 0   
    nodes[1](0) = cutterXAABB(0,1); nodes[1](1) = cutterXAABB(1,0); nodes[1](2) = cutterXAABB(2,0); // node 1
    nodes[2](0) = cutterXAABB(0,1); nodes[2](1) = cutterXAABB(1,1); nodes[2](2) = cutterXAABB(2,0); // node 2
    nodes[3](0) = cutterXAABB(0,0); nodes[3](1) = cutterXAABB(1,1); nodes[3](2) = cutterXAABB(2,0); // node 3
    nodes[4](0) = cutterXAABB(0,0); nodes[4](1) = cutterXAABB(1,0); nodes[4](2) = cutterXAABB(2,1); // node 4
    nodes[5](0) = cutterXAABB(0,1); nodes[5](1) = cutterXAABB(1,0); nodes[5](2) = cutterXAABB(2,1); // node 5
    nodes[6](0) = cutterXAABB(0,1); nodes[6](1) = cutterXAABB(1,1); nodes[6](2) = cutterXAABB(2,1); // node 6
    nodes[7](0) = cutterXAABB(0,0); nodes[7](1) = cutterXAABB(1,1); nodes[7](2) = cutterXAABB(2,1); // node 7
    
    for (int i = 0; i < 8; i++)
        if(isPositionWithinXAABB(nodes[i], xfemXAABB))
        {
            intersection = true;
            break;
        }
    
    if(!intersection)
    {
        for (int i = 0; i < 12; i++)
        {
            const int index1 = eleNodeNumbering_hex27_lines[i][0];
            const int index2 = eleNodeNumbering_hex27_lines[i][1];
            if(isLineWithinXAABB(nodes[index1], nodes[index2], xfemXAABB))
            {
                intersection = true;
                break;
            }
        }
    }
    
    if(!intersection)
    {
        nodes[0](0) = xfemXAABB(0,0);   nodes[0](1) = xfemXAABB(1,0);   nodes[0](2) = xfemXAABB(2,0);   // node 0   
        nodes[1](0) = xfemXAABB(0,1);   nodes[1](1) = xfemXAABB(1,0);   nodes[1](2) = xfemXAABB(2,0);   // node 1
        nodes[2](0) = xfemXAABB(0,1);   nodes[2](1) = xfemXAABB(1,1);   nodes[2](2) = xfemXAABB(2,0);   // node 2
        nodes[3](0) = xfemXAABB(0,0);   nodes[3](1) = xfemXAABB(1,1);   nodes[3](2) = xfemXAABB(2,0);   // node 3
        nodes[4](0) = xfemXAABB(0,0);   nodes[4](1) = xfemXAABB(1,0);   nodes[4](2) = xfemXAABB(2,1);   // node 4
        nodes[5](0) = xfemXAABB(0,1);   nodes[5](1) = xfemXAABB(1,0);   nodes[5](2) = xfemXAABB(2,1);   // node 5
        nodes[6](0) = xfemXAABB(0,1);   nodes[6](1) = xfemXAABB(1,1);   nodes[6](2) = xfemXAABB(2,1);   // node 6
        nodes[7](0) = xfemXAABB(0,0);   nodes[7](1) = xfemXAABB(1,1);   nodes[7](2) = xfemXAABB(2,1);   // node 7
    
        for (int i = 0; i < 8; i++)
            if(isPositionWithinXAABB(nodes[i], cutterXAABB))
            {
                intersection = true;
                break;
            }
    }   
    
    if(!intersection)
    {
        for (int i = 0; i < 12; i++)
        {
            const int index1 = eleNodeNumbering_hex27_lines[i][0];
            const int index2 = eleNodeNumbering_hex27_lines[i][1];
            if(isLineWithinXAABB(nodes[index1], nodes[index2], cutterXAABB))
            {
                intersection = true;
                break;
            }
        }
    }
    return intersection;
}


#endif  // #ifdef CCADISCRET
