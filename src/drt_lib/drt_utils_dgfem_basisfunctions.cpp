/*!
 \file drt_utils_dgfem_basisfunctions.cpp

 \brief Provide function values of basis-functions used in STDG-FEM

 Provided are 1D, 2D and 3D shape functions with special functions to be used 
 for the FEM discretization in time.
 
 Use is made of Legendre polynomial basis-functions using the Rodrigues representation 
 for which the $p$-th order polynomial is given by $P_p=\frac{1}{2^p p!}\frac{d^p}{dx^p}(x^2-1)^p$
 and the first realizations are given by $1$,$x$,$(3x^2-1)/2$,$\ldots$.
 
 \author Fedderik van der Bos
 bos@lnm.mw.tum.de
 http://www.lnm.mw.tum.de
 089 - 289-15253
 */
#ifdef CCADISCRET

#include "drt_element.H"
#include "drt_utils.H"
#include "drt_dserror.H"

//
// 3D Legendre basis-functions
//
void DRT::UTILS::DGFEM_BasisFunction_3D(
        blitz::Array<double, 1>& funct,
        const double& r,
        const double& s,
        const double& t,
        const DRT::Element::DiscretizationType& shape,
        const int& npol)
{
    switch (shape)
    {
    case DRT::Element::hex8:
    	switch (npol)
    	{
    	case 1:
    		funct(0)=1.0;
    		
    	case 4:
    		funct(0)=1.0;
    		
    		funct(1)=r;
    		funct(2)=s;
    		funct(3)=t;
    		
    	case 10:
    		funct(0)=1.0;
    		
    		funct(1)=r;
    		funct(2)=s;
    		funct(3)=t;
    		
    		funct(4)=0.5*(3.0*r*r-1.0); 	//P2(r)
    		funct(5)=r*s;					//P1(r)*P1(s)
    		funct(6)=r*t;					//P1(r)*P1(t)
    		funct(7)=0.5*(3.0*s*s-1.0);		//P2(s)
    		funct(8)=s*t;					//P1(s)*P1(t)
    		funct(9)=0.5*(3.0*t*t-1.0);		//P2(t)
    		
    	case 20:
    		funct(0)=1.0;
    		
    		funct(1)=r;
    		funct(2)=s;
    		funct(3)=t;
    		
    		funct(4)=0.5*(3.0*r*r-1.0); 	//P2(r)
    		funct(5)=r*s;					//P1(r)*P1(s)
    		funct(6)=r*t;					//P1(r)*P1(t)
    		funct(7)=0.5*(3.0*s*s-1.0);		//P2(s)
    		funct(8)=s*t;					//P1(s)*P1(t)
    		funct(9)=0.5*(3.0*t*t-1.0);		//P2(t)

    		funct(10)=0.5*(5.0*r*r*r-3.0*r);		//P3(r)
    		funct(11)=0.5*(3.0*r*r-1.0)*s;			//P2(r)*P1(s)
    		funct(12)=0.5*(3.0*r*r-1.0)*t;			//P2(r)*P1(t)
    		funct(13)=r*0.5*(3.0*s*s-1.0);			//P1(r)*P2(s)
    		funct(14)=r*s*t;						//P1(r)*P1(s)*P1(t)
    		funct(15)=r*0.5*(3.0*t*t-1.0);			//P1(r)*P2(t)
    		funct(16)=0.5*(5.0*s*s*s-3.0*s);		//P3(s)
    		funct(17)=0.5*(3.0*s*s-1.0)*t;			//P2(s)*P1(t)
    		funct(18)=s*0.5*(3.0*t*t-1.0);			//P1(s)*P2(t)
    		funct(19)=0.5*(5.0*t*t*t-3.0*t);		//P3(t)

    	default:
    		dserror("Unknown number of polynomials\n");
    	}

    default:
        dserror("shape unknown\n");
    } /* end switch(distype) */

    return;
}

//
// first natural derivative of 3D Legendre basis-functions
//
void DRT::UTILS::DGFEM_BasisFunction_3D_deriv1(
        blitz::Array<double, 2>& deriv1,
        const double& r,
        const double& s,
        const double& t,
        const DRT::Element::DiscretizationType& shape,
        const int& npol)
{
	deriv1=0.0;
	
    switch (shape)
    {
    case DRT::Element::hex8:
    	switch (npol)
    	{
    	case 1:
    		//Nothing to do
    		
    	case 4:
    		deriv1(0,1)=1.0;
    		deriv1(1,2)=1.0;
    		deriv1(2,3)=1.0;

    	case 10:
    		deriv1(0,1)=1.0;
    		deriv1(1,2)=1.0;
    		deriv1(2,3)=1.0;

    		deriv1(0,4)=3.0*r;
    		deriv1(0,5)=s;
    		deriv1(1,5)=r;
    		deriv1(0,6)=t;
    		deriv1(2,6)=r;
    		deriv1(1,7)=3.0*s;
    		deriv1(1,8)=s;
    		deriv1(2,8)=r;
    		deriv1(2,9)=3.0*t;
    		
    	case 20:
    		deriv1(0,1)=1.0;
    		deriv1(1,2)=1.0;
    		deriv1(2,3)=1.0;

    		deriv1(0,4)=3.0*r;
    		deriv1(0,5)=s;
    		deriv1(1,5)=r;
    		deriv1(0,6)=t;
    		deriv1(2,6)=r;
    		deriv1(1,7)=3.0*s;
    		deriv1(1,8)=s;
    		deriv1(2,8)=r;
    		deriv1(2,9)=3.0*t;

    		deriv1(0,10)=0.5*(15.0*r*r-3.0);	
    		deriv1(0,11)=3.0*r*s;
    		deriv1(1,11)=0.5*(3.0*r*r-1.0);
    		deriv1(0,12)=3.0*r*t;
    		deriv1(2,12)=0.5*(3.0*r*r-1.0);
    		deriv1(0,13)=0.5*(3.0*s*s-1.0);
    		deriv1(1,13)=r*3.0*s;
    		deriv1(0,14)=s*t;
    		deriv1(1,14)=r*t;
    		deriv1(2,14)=r*s;
    		deriv1(0,15)=0.5*(3.0*t*t-1.0);
    		deriv1(2,15)=r*3.0*t;
    		deriv1(1,16)=0.5*(15.0*s*s-3.0);
    		deriv1(1,17)=3.0*s*t;
    		deriv1(2,17)=0.5*(3.0*s*s-1.0);
    		deriv1(1,18)=0.5*(3.0*t*t-1.0);
    		deriv1(2,18)=s*3.0*t;
    		deriv1(2,19)=0.5*(15.0*t*t-3.0);
    		
    	default:
    		dserror("Unknown number of polynomials\n");
    	}

    default:
        dserror("shape unknown\n");
    } /* end switch(distype) */

    return;
}

//
// second natural derivative of 3D Legendre basis-functions
//
void DRT::UTILS::DGFEM_BasisFunction_3D_deriv2(
        blitz::Array<double, 2>& deriv2,
        const double& r,
        const double& s,
        const double& t,
        const DRT::Element::DiscretizationType& shape,
        const int& npol)
{
    const int drdr = 0;
    const int dsds = 1;
    const int dtdt = 2;
    const int drds = 3;
    const int drdt = 4;
    const int dsdt = 5;
	
	deriv2=0.0;
	
    switch (shape)
    {
    case DRT::Element::hex8:
    	switch (npol)
    	{
    	case 1:
    		//Nothing to do
    		
    	case 4:
    		//Nothing to do

    	case 10:
    		deriv2(drdr,4)=3.0;
    		deriv2(dsds,7)=3.0;
    		deriv2(dtdt,9)=3.0;
    		
    	case 20:
    		deriv2(drdr,4)=3.0;
    		deriv2(drds,5)=1.0;
    		deriv2(drdt,6)=1.0;
    		deriv2(dsds,7)=3.0;
    		deriv2(dsdt,8)=1.0;
    		deriv2(dtdt,9)=3.0;

    		deriv2(drdr,10)=15.0*r;	
    		deriv2(drdr,11)=3.0*s;
    		deriv2(drds,11)=3.0*r;
    		deriv2(drdr,12)=3.0*t;
    		deriv2(drdt,12)=3.0*r;
    		deriv2(dsds,13)=3.0*r;
    		deriv2(drds,13)=3.0*s;
    		deriv2(drds,14)=t;
    		deriv2(drdt,14)=s;
    		deriv2(dsdt,14)=r;
    		deriv2(dtdt,15)=3.0*r;
    		deriv2(drdt,15)=3.0*t;
    		deriv2(dsds,16)=15.0*s;
    		deriv2(dsds,17)=3.0*t;
    		deriv2(dsdt,17)=3.0*s;
    		deriv2(dsdt,18)=3.0*t;
    		deriv2(dtdt,18)=3.0*s;
    		deriv2(dtdt,19)=15.0*t;
    		
    	default:
    		dserror("Unknown number of polynomials\n");
    	}

    default:
        dserror("shape unknown\n");
    } /* end switch(distype) */

    return;
}



//
// 2D Legendre basis-functions
//
void DRT::UTILS::DGFEM_BasisFunction_2D(
        blitz::Array<double, 1>& funct,
        const double& r,
        const double& s,
        const DRT::Element::DiscretizationType& shape,
        const int& npol)
{
    switch (shape)
    {
    case DRT::Element::quad4:
    	switch (npol)
    	{
    	case 1:
    		funct(0)=1.0;
    		
    	case 3:
    		funct(0)=1.0;
    		
    		funct(1)=r;
    		funct(2)=s;
    		
    	case 6:
    		funct(0)=1.0;
    		
    		funct(1)=r;
    		funct(2)=s;
    		
    		funct(3)=0.5*(3.0*r*r-1.0); 	//P2(r)
    		funct(4)=r*s;					//P1(r)*P1(s)
    		funct(5)=0.5*(3.0*s*s-1.0);		//P2(s)
    		
    	case 10:
    		funct(0)=1.0;
    		
    		funct(1)=r;
    		funct(2)=s;
    		
    		funct(3)=0.5*(3.0*r*r-1.0); 	//P2(r)
    		funct(4)=r*s;					//P1(r)*P1(s)
    		funct(5)=0.5*(3.0*s*s-1.0);		//P2(s)

    		funct(6)=0.5*(5.0*r*r*r-3.0*r);		//P3(r)
    		funct(7)=0.5*(3.0*r*r-1.0)*s;		//P2(r)*P1(s)
    		funct(8)=r*0.5*(3.0*s*s-1.0);		//P1(r)*P2(s)
    		funct(9)=0.5*(5.0*s*s*s-3.0*s);		//P3(s)

    	default:
    		dserror("Unknown number of polynomials\n");
    	}

    default:
        dserror("shape unknown\n");
    } /* end switch(distype) */

    return;
}

//
// first natural derivative of 2D Legendre basis-functions
//
void DRT::UTILS::DGFEM_BasisFunction_2D_deriv1(
        blitz::Array<double, 2>& deriv1,
        const double& r,
        const double& s,
        const DRT::Element::DiscretizationType& shape,
        const int& npol)
{
	deriv1=0.0;
	
    switch (shape)
    {
    case DRT::Element::quad4:
    	switch (npol)
    	{
    	case 1:
    		//Nothing to do
    		
    	case 3:
    		deriv1(0,1)=1.0;
    		deriv1(1,2)=1.0;

    	case 6:
    		deriv1(0,1)=1.0;
    		deriv1(1,2)=1.0;

    		deriv1(0,3)=3.0*r;
    		deriv1(0,4)=s;
    		deriv1(1,4)=r;
    		deriv1(1,5)=3.0*s;
    		
    	case 10:
    		deriv1(0,1)=1.0;
    		deriv1(1,2)=1.0;

    		deriv1(0,3)=3.0*r;
    		deriv1(0,4)=s;
    		deriv1(1,4)=r;
    		deriv1(1,5)=3.0*s;

    		deriv1(0,6)=0.5*(15.0*r*r-3.0);	
    		deriv1(0,7)=3.0*r*s;
    		deriv1(1,7)=0.5*(3.0*r*r-1.0);
    		deriv1(0,8)=0.5*(3.0*s*s-1.0);
    		deriv1(1,8)=r*3.0*s;
    		deriv1(1,9)=0.5*(15.0*s*s-3.0);
    		
    	default:
    		dserror("Unknown number of polynomials\n");
    	}

    default:
        dserror("shape unknown\n");
    } /* end switch(distype) */

    return;
}

//
// second natural derivative of 2D Legendre basis-functions
//
void DRT::UTILS::DGFEM_BasisFunction_2D_deriv2(
        blitz::Array<double, 2>& deriv2,
        const double& r,
        const double& s,
        const DRT::Element::DiscretizationType& shape,
        const int& npol)
{
    const int drdr = 0;
    const int dsds = 1;
    const int drds = 2;
	
	deriv2=0.0;
	
    switch (shape)
    {
    case DRT::Element::quad4:
    	switch (npol)
    	{
    	case 1:
    		//Nothing to do
    		
    	case 3:
    		//Nothing to do

    	case 6:
    		deriv2(drdr,3)=3.0;
    		deriv2(dsds,5)=3.0;
    		
    	case 10:
    		deriv2(drdr,3)=3.0;
    		deriv2(dsds,5)=3.0;

    		deriv2(drdr,6)=15.0*r;	
    		deriv2(drdr,7)=3.0*s;
    		deriv2(drds,7)=3.0*r;
    		deriv2(dsds,8)=3.0*r;
    		deriv2(drds,8)=3.0*s;
    		deriv2(dsds,9)=15.0*s;
    		
    	default:
    		dserror("Unknown number of polynomials\n");
    	}

    default:
        dserror("shape unknown\n");
    } /* end switch(distype) */

    return;
}

//
// 1D Legendre basis-functions
//
void DRT::UTILS::DGFEM_BasisFunction_1D(
        blitz::Array<double, 1>& funct,
        const double& r,
        const DRT::Element::DiscretizationType& shape,
        const int& npol)
{
    switch (shape)
    {
    case DRT::Element::line2:
    	switch (npol)
    	{
    	case 1:
    		funct(0)=1.0;
    		
    	case 2:
    		funct(0)=1.0;
    		funct(1)=r;
    		
    	case 3:
    		funct(0)=1.0;
    		funct(1)=r;
    		funct(2)=0.5*(3.0*r*r-1.0); 	//P2(r)
    		
    	case 4:
    		funct(0)=1.0;
    		funct(1)=r;
    		funct(2)=0.5*(3.0*r*r-1.0); 	//P2(r)
    		funct(3)=0.5*(5.0*r*r*r-3.0*r);	//P3(r)

    	default:
    		dserror("Unknown number of polynomials\n");
    	}

    default:
        dserror("shape unknown\n");
    } /* end switch(distype) */

    return;
}

//
// first natural derivative of 1D Legendre basis-functions
//
void DRT::UTILS::DGFEM_BasisFunction_1D_deriv1(
        blitz::Array<double, 1>& deriv1,
        const double& r,
        const DRT::Element::DiscretizationType& shape,
        const int& npol)
{
	deriv1=0.0;
	
    switch (shape)
    {
    case DRT::Element::line2:
    	switch (npol)
    	{
    	case 1:
    		//Nothing to do
    		
    	case 2:
    		deriv1(1)=1.0;

    	case 3:
    		deriv1(1)=1.0;
    		deriv1(2)=3.0*r;
    		
    	case 4:
    		deriv1(1)=1.0;
    		deriv1(2)=3.0*r;
    		deriv1(3)=0.5*(15.0*r*r-3.0);	
    		
    	default:
    		dserror("Unknown number of polynomials\n");
    	}

    default:
        dserror("shape unknown\n");
    } /* end switch(distype) */

    return;
}

//
// second natural derivative of 1D Legendre basis-functions
//
void DRT::UTILS::DGFEM_BasisFunction_1D_deriv2(
        blitz::Array<double, 1>& deriv2,
        const double& r,
        const DRT::Element::DiscretizationType& shape,
        const int& npol)
{
	deriv2=0.0;
	
    switch (shape)
    {
    case DRT::Element::line2:
    	switch (npol)
    	{
    	case 1:
    		//Nothing to do
    		
    	case 2:
    		//Nothing to do

    	case 3:
    		deriv2(2)=3.0;
    		
    	case 4:
    		deriv2(2)=3.0;
    		deriv2(3)=15.0*r;	
    		
    	default:
    		dserror("Unknown number of polynomials\n");
    	}

    default:
        dserror("shape unknown\n");
    } /* end switch(distype) */

    return;
}

#endif  // #ifdef CCADISCRET
