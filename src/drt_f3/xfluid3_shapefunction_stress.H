/*----------------------------------------------------------------------*/
/*!
\file xfluid3_shapefunction_stress.H

\brief rotuines to assemble element stiffness matrix and vector

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET


#ifndef XFLUID_SHAPEFUNCTION_STRESS_H_
#define XFLUID_SHAPEFUNCTION_STRESS_H_


typedef blitz::Array<double,1> BlitzVec; ///< a short-cut for a one-dimensional BLITZ vector
typedef blitz::Array<double,2> BlitzMat; ///< a short-cut for a two-dimensional BLITZ matrix


namespace DRT
{
namespace Elements
{
namespace FLUID
{
/*!
\brief Returns the number of nodes of the discretization type
       for each discretization type
       - template version
     
\return number of nodes of the discretization type
*/  
template <DRT::Element::DiscretizationType DISTYPE>
int getNumberOfStressDofs()
{
   switch(DISTYPE)
   {
       case DRT::Element::hex8:
           return 1;
           break;
       case DRT::Element::hex20:
       case DRT::Element::hex27:
           return 8;
           break;
       case DRT::Element::tet4:
           return 1;
           break;
       case DRT::Element::tet10:
           return 4;
           break;   
       default:
           dserror("discretization type not yet implemented");     
   }
};


/*!
\brief Fill a vector with with 3D stress shape function
*/
template <DRT::Element::DiscretizationType DISTYPE>
BlitzVec stress_shape_function_3D(
        const double r,
        const double s,
        const double t
        )
{
    
    BlitzVec funct(getNumberOfStressDofs<DISTYPE>());
    
    const double Q18 = 1.0/8.0;
    const double Q12 = 1.0/2.0;

    switch (DISTYPE)
    {
    case DRT::Element::hex8:
    {

        funct(0)=1.0;

        break;
    }
    
    case DRT::Element::hex20:
    case DRT::Element::hex27:
    {
        const double rp=1.0+r;
        const double rm=1.0-r;
        const double sp=1.0+s;
        const double sm=1.0-s;
        const double tp=1.0+t;
        const double tm=1.0-t;

        funct(0)=Q18*rm*sm*tm;
        funct(1)=Q18*rp*sm*tm;
        funct(2)=Q18*rp*sp*tm;
        funct(3)=Q18*rm*sp*tm;
        funct(4)=Q18*rm*sm*tp;
        funct(5)=Q18*rp*sm*tp;
        funct(6)=Q18*rp*sp*tp;
        funct(7)=Q18*rm*sp*tp;

        break;
    }
    case DRT::Element::tet4:
    {
        funct(0)= 1.0;
        break;
    }
    case DRT::Element::tet10:
    {
        const double t1=1.0-r-s-t;
        const double t2=r;
        const double t3=s;
        const double t4=t;

        funct(0)= t1;
        funct(1)= t2;
        funct(2)= t3;
        funct(3)= t4;
        break;
    }
    case DRT::Element::wedge6:
    {
        funct(0)=1.0;
        break;
    }
    case DRT::Element::wedge15:
    {
        const double t3=1.0-r-s;

        funct(0)=Q12*r*(1-t);
        funct(1)=Q12*s*(1-t);
        funct(2)=Q12*t3*(1-t);
        funct(3)=Q12*r*(1+t);
        funct(4)=Q12*s*(1+t);
        funct(5)=Q12*t3*(1+t);
        break;
    }
    case DRT::Element::pyramid5:
    {
        funct(0)=1.0;
        break;
    }
    default:
        dserror("distyp unknown\n");
    }

    return funct;
}

/*!
\brief Fill a vector with with first 3D stress shape function derivatives
*/
template <DRT::Element::DiscretizationType DISTYPE>
BlitzMat stress_shape_function_3D_deriv1(
        const double& r,
        const double& s,
        const double& t
        )
{
    
    BlitzMat deriv1(3,getNumberOfStressDofs<DISTYPE>());
    
    const double Q18 = 1.0/8.0;

    switch (DISTYPE)
    {
    case DRT::Element::hex8:
    {
        deriv1(0, 0) = 0.0;
        deriv1(1, 0) = 0.0;
        deriv1(2, 0) = 0.0;
        break;
    }
    case DRT::Element::hex20:
    case DRT::Element::hex27:
    {
        const double rp=1.0+r;
        const double rm=1.0-r;
        const double sp=1.0+s;
        const double sm=1.0-s;
        const double tp=1.0+t;
        const double tm=1.0-t;

        deriv1(0, 0)=-Q18*sm*tm;
        deriv1(1, 0)=-Q18*tm*rm;
        deriv1(2, 0)=-Q18*rm*sm;

        deriv1(0, 1)= Q18*sm*tm;
        deriv1(1, 1)=-Q18*tm*rp;
        deriv1(2, 1)=-Q18*rp*sm;

        deriv1(0, 2)= Q18*sp*tm;
        deriv1(1, 2)= Q18*tm*rp;
        deriv1(2, 2)=-Q18*rp*sp;

        deriv1(0, 3)=-Q18*sp*tm;
        deriv1(1, 3)= Q18*tm*rm;
        deriv1(2, 3)=-Q18*rm*sp;

        deriv1(0, 4)=-Q18*sm*tp;
        deriv1(1, 4)=-Q18*tp*rm;
        deriv1(2, 4)= Q18*rm*sm;

        deriv1(0, 5)= Q18*sm*tp;
        deriv1(1, 5)=-Q18*tp*rp;
        deriv1(2, 5)= Q18*rp*sm;

        deriv1(0, 6)= Q18*sp*tp;
        deriv1(1, 6)= Q18*tp*rp;
        deriv1(2, 6)= Q18*rp*sp;

        deriv1(0, 7)=-Q18*sp*tp;
        deriv1(1, 7)= Q18*tp*rm;
        deriv1(2, 7)= Q18*rm*sp;
        break;
    }
    default:
        dserror("distyp unknown\n");
    } /* end switch(distype) */

    return deriv1;
}
}
}
}



#endif /*XFLUID_SHAPEFUNCTION_STRESS_H_*/
#endif
