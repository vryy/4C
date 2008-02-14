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

#include "../drt_lib/drt_element.H"


namespace XFLUID
{

//enum InterpolationType
//{
//    P0,        ///< constant
//    P1,        ///< complete linear
//    P2,        ///< complete quadratic
//    Q1,        ///< bi-/tri-linear
//    Q2,        ///< bi-/tri-quadratic
//    Pm1,       ///< discontinuous complete linear
//    Pm2,       ///< discontinuous complete quadratic
//    Qm1,       ///< discontinuous bi-/tri-linear
//    Qm2        ///< discontinuous bi-/tri-quadratic
//};


//enum ElementFormulation
//{
//    // default parameter
//    QQQ_undefined,
//    
//    // inf-sup stable velocity pressure formulations
//    Q2Q1,      ///< bi-/tri-quadratic velocity, bi-/tri-linear continuous pressure (Taylor Hood) 
//    Q2Pm1,     ///< bi-/tri-quadratic velocity, linear discontinuous pressure
//
//    // stabilized velocity-pressure formulation
//    P1P1,      ///< cont. linear velocity and pressure
//    P2P2,      ///< cont. quadratic velocity and pressure
//    Q1Q1,      ///< cont. bi-/tri-linear velocity and pressure
//    Q2Q2,      ///< cont. bi-/tri-quadratic velocity and pressure
//    
//    // inf-sup stable stress formulations (pressure is not inf-sup stable -> stabilized as always)
//    Qm1Q1Q1,   ///< discont. bi-/tri-linear stress, cont. bi-/tri-linear velocity and pressure
//    Qm2Q2Q2,   ///< discont. bi-/tri-quadratic stress, cont. bi-/tri-quadratic velocity and pressure
//    Pm2Q2Q2,   // ??? is the stress interpolation sufficient high?
//    
//    // inf-sup unstable -- need stress stabilization
//    P0Q1Q1,    ///< discont. constant stress, cont. linear velocity and pressure
//    Pm1Q1Q1,   ///< discont. linear stress, cont. bi-linear velocity and pressure
//    P0Q2Q2,
//    Pm1Q2Q2,
//    Qm1Q2Q2,
//
//    // 
//    Pm1Q2Pm1   ///< cont. quadratic velocity, discontinous pressure and stress (P_-1)
//    
//};


/*!
\brief Returns the distype for stress interpolation
       for each parent element discretization type in 3d
     
\return distype for stress interpolation
*/
template <DRT::Element::DiscretizationType veldistype>
inline DRT::Element::DiscretizationType getStressInterpolationType3D()
{
    DRT::Element::DiscretizationType distype = DRT::Element::dis_none;
    switch (veldistype)
    {
        case DRT::Element::tet4:
            distype = DRT::Element::point1;
            break;
        case DRT::Element::hex8:
            distype = DRT::Element::hex8;
            break;
        case DRT::Element::hex20:
            distype = DRT::Element::hex20;
            break;
        case DRT::Element::hex27:
            distype = DRT::Element::hex27;
            break;
        default:
            dserror("no stressdistype defined for velocity distype");
    }
    return distype;
}


/*!
\brief Returns the distype for stress interpolation
       for each parent element discretization type in 3d
     
\return distype for stress interpolation
*/  
inline DRT::Element::DiscretizationType getStressInterpolationType3D(DRT::Element::DiscretizationType veldistype)
{
    DRT::Element::DiscretizationType distype = DRT::Element::dis_none;
    switch (veldistype)
    {
        case DRT::Element::tet4:
            distype = getStressInterpolationType3D<DRT::Element::tet4>();
            break;
        case DRT::Element::hex8:
            distype = getStressInterpolationType3D<DRT::Element::hex8>();
            break;
        case DRT::Element::hex20:
            distype = getStressInterpolationType3D<DRT::Element::hex20>();
            break;
        case DRT::Element::hex27:
            distype = getStressInterpolationType3D<DRT::Element::hex27>();
            break;
        default:
            dserror("no stressdistype defined for velocity distype");
    }
    return distype;
}

///*!
//\brief Returns the distype for stress interpolation
//       for each parent element discretization type in 3d
//     
//\return distype for stress interpolation
//*/  
//inline DRT::Element::DiscretizationType getStressInterpolationType3D(ElementFormulation interpol)
//{
//    //  if interpol names a interpolation type, we can return the appropriate stress interpolation
//    DRT::Element::DiscretizationType distype = DRT::Element::dis_none;
//    switch (interpol)
//    {
//        case P1P1:
//            distype = DRT::Element::dis_none;
//            break;
//        case P2P2:
//            distype = DRT::Element::dis_none;
//            break;
//        case Q1Q1:
//            distype = DRT::Element::dis_none;
//            break;
//        case Q2Q2:
//            distype = DRT::Element::dis_none;
//            break;
//        case P0Q1Q1:
//            distype = DRT::Element::point1;
//            break;
//        case Pm1Q1Q1:
//            distype = DRT::Element::tet4;
//            break;
//        case Qm1Q1Q1:
//            distype = DRT::Element::hex8;
//            break;
//        case P0Q2Q2:
//            distype = DRT::Element::point1;
//            break;
//        case Qm1Q2Q2:
//            distype = DRT::Element::hex8;
//            break;
//        case Qm2Q2Q2:
//            distype = DRT::Element::hex27;
//            break;
//        case Pm1Q2Q2:
//            distype = DRT::Element::tet4;
//            break;
//        case Pm2Q2Q2:
//            distype = DRT::Element::tet10;
//            break;
//        case Pm1Q2Pm1:
//            distype = DRT::Element::tet4;
//            break;
//        default:
//            dserror("rule not defined");
//    }
//    return distype;
//}

//inline ElementFormulation getElementFormulation(DRT::Element::DiscretizationType distype)
//{
//    XFLUID::ElementFormulation form = XFLUID::QQQ_undefined;
//    switch (distype)
//    {
//    case DRT::Element::tet4:
//        //form = XFLUID::Qm1Q1Q1;
//        form = XFLUID::P1P1;
//        break;
//    case DRT::Element::tet10:
//        //form = XFLUID::Qm1Q1Q1;
//        form = XFLUID::P2P2;
//        break;
//    case DRT::Element::hex8:
//        //form = XFLUID::Qm1Q1Q1;
//        form = XFLUID::Q1Q1;
//        break;
//    case DRT::Element::hex27:
//        //form = XFLUID::Qm2Q2Q2;
//        //form = XFLUID::Qm1Q2Q2;
//        //form = XFLUID::Pm2Q2Q2;
//        form = XFLUID::Q2Q2;
//        break;
//    default:
//        dserror("not defined!");
//    }
//    return form;
//}

}


#endif /*XFLUID3_INTERPOLATION_H*/
#endif
