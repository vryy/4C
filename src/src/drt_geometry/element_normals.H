/*!----------------------------------------------------------------------
\file element_normals.H

\brief computes element normals

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ELEMENT_NORMALS_H_
#define ELEMENT_NORMALS_H_


// #include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"


namespace GEO
{


/*!
\brief calculate surface element normal vector in gausspoint (r,s)

\param deriv              (in)        : shape functions derivative in (r,s) -> (2,numnode)
\param xyze               (in)        : element nodal coordinates in (x,y,z) space -> (3,numnode)
\param unrm               (out)       : unit normal vector at (r,s)
*/
template<DRT::Element::DiscretizationType DISTYPE, class M1, class M2>
inline void computeNormalToSurfaceElementViaShape(
        const M1&                 deriv,
        const M2&                 xyze,
        LINALG::Matrix<3,1>&      unrm)
{
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // compute dXYZ / drs
    static LINALG::Matrix<3,2> dxyzdrs;
    dxyzdrs.Clear();
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < numnode; k++)
                dxyzdrs(i,j) += xyze(i,k)*deriv(j,k);

    /* compute covariant metric tensor G for surface element
    **                        | g11   g12 |
    **                    G = |           |
    **                        | g12   g22 |
    ** where (o denotes the inner product, xyz a vector)
    **
    **       dXYZ   dXYZ          dXYZ   dXYZ          dXYZ   dXYZ
    ** g11 = ---- o ----    g12 = ---- o ----    g22 = ---- o ----
    **        dr     dr            dr     ds            ds     ds
    */
    static LINALG::Matrix<2,2> metrictensor;
    metrictensor.MultiplyTN(dxyzdrs,dxyzdrs);

    unrm(0) = dxyzdrs(1,0) * dxyzdrs(2,1) - dxyzdrs(2,0) * dxyzdrs(1,1);
    unrm(1) = dxyzdrs(2,0) * dxyzdrs(0,1) - dxyzdrs(0,0) * dxyzdrs(2,1);
    unrm(2) = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(1,0) * dxyzdrs(0,1);

    // scale
    unrm.Scale(1.0 / unrm.Norm2());

    return;
}



/*!
\brief calculate element normal in gausspoint (r,s)
\param xyze_surfaceElement   (in)   nodal position array (3,numnode)
\param eleCoord              (in)   node in element coordinates (r, s)
\param unitnormalvec         (out)  unit normal vector
*/
template<DRT::Element::DiscretizationType DISTYPE, class M>
inline void computeNormalToSurfaceElementT(
    const M&                      xyze_surfaceElement,
    const LINALG::Matrix<2,1>&    eleCoord,
    LINALG::Matrix<3,1>&          unitnormalvec)
{
  const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
  static LINALG::Matrix<2,numnode> deriv_boundary;
  DRT::UTILS::shape_function_2D_deriv1(deriv_boundary,eleCoord(0),eleCoord(1),DISTYPE);

  computeNormalToSurfaceElementViaShape<DISTYPE>(deriv_boundary, xyze_surfaceElement, unitnormalvec);
  return;
}



/*!
\brief calculate element normal in gausspoint (r,s)
\param surfaceElement        (in) :  surface element
\param xyze_surfaceElement   (in) :  nodal position array (3,numnode)
\param eleCoord              (in) :  node in element coordinates (r, s)
\param unitnormalvec         (out):  unit normal vector
*/
template <class M>
inline void computeNormalToSurfaceElement(
    const DRT::Element::DiscretizationType  distype,
    const M&                                xyze_surfaceElement,
    const LINALG::Matrix<2,1>&              eleCoord,
    LINALG::Matrix<3,1>&                    unitnormalvec)
{
  switch (distype)
  {
  case DRT::Element::quad4:
    computeNormalToSurfaceElementT<DRT::Element::quad4>(xyze_surfaceElement, eleCoord, unitnormalvec);
    break;
  case DRT::Element::quad8:
    computeNormalToSurfaceElementT<DRT::Element::quad8>(xyze_surfaceElement, eleCoord, unitnormalvec);
    break;
  case DRT::Element::quad9:
    computeNormalToSurfaceElementT<DRT::Element::quad9>(xyze_surfaceElement, eleCoord, unitnormalvec);
    break;
  case DRT::Element::tri3:
    computeNormalToSurfaceElementT<DRT::Element::tri3>(xyze_surfaceElement, eleCoord, unitnormalvec);
    break;
  case DRT::Element::tri6:
    computeNormalToSurfaceElementT<DRT::Element::tri6>(xyze_surfaceElement, eleCoord, unitnormalvec);
    break;
  default:
    cout << DRT::DistypeToString(distype) << " not found!" << endl;
    dserror("please add your distype!");
  }
  return;
}


} // namespace GEO


#endif /*ELEMENT_NORMALS_H_NORMALS_H_*/
#endif /*CCADISCRET*/
