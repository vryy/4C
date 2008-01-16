/*----------------------------------------------------------------------*/
/*!
\file xfluid3_interpolation.H

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


#ifndef XFLUID3_INTERPOLATION_H
#define XFLUID3_INTERPOLATION_H


namespace XFLUID
{

enum ElementFormulation
{
    QQQQundefined, ///< dummy value to identify undefined behaviour
    Q2Q1,      ///< Bi-quadratic velocity, linear continuous pressure (Taylor Hood) 
    Q2Pm1,     ///< Bi-quadratic velocity, linear discontinuous pressure
    Q2Q2,      ///< cont. quadratic velocity and pressure
    Q1Q1P0,    ///< cont. linear velocity and pressure, discont. constant stress
    Q2Q2Pm1,
    Q2Pm1Pm1   ///< cont. quadratic velocity, discontinous pressure and stress (P_-1)
};

/*!
\brief Returns the distype for stress interpolation
       for each parent element discretization type
     
\return distype for stress interpolation
*/  
inline DRT::Element::DiscretizationType getStressInterpolationType(ElementFormulation interpol)
{
    //  if interpol names a interpolation type, we can return the appropriate stress interpolation
    switch (interpol)
    {
        case Q1Q1P0:
            return DRT::Element::point1;
            break;
        case Q2Q2Pm1:
            return DRT::Element::tet4;
            break;
        case Q2Pm1Pm1:
            return DRT::Element::tet4;
            break;
        default:
            dserror("discretization type not yet implemented");
            return DRT::Element::dis_none;
    }
};

inline ElementFormulation getElementFormulation(DRT::Element::DiscretizationType distype)
{
    switch (distype)
    {
    case DRT::Element::hex8:
        return XFLUID::Q1Q1P0;
        break;
    case DRT::Element::hex27:
        return XFLUID::Q2Q2Pm1;
        break;
    default:
        dserror("not defined!");
        return XFLUID::QQQQundefined;
    }
};

}


#endif /*XFLUID3_INTERPOLATION_H*/
#endif
