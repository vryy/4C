/*!----------------------------------------------------------------------
\file element_normals.H

\brief computes element normals

\level 2

<pre>
\maintainer Ager Christoph
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
</pre>
*----------------------------------------------------------------------*/


#ifndef ELEMENT_NORMALS_H_
#define ELEMENT_NORMALS_H_


#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"


namespace GEO
{
  /*!
  \brief calculate surface element normal vector in gausspoint (r,s)

  \param deriv              (in)        : shape functions derivative in (r,s) -> (2,numnode)
  \param xyze               (in)        : element nodal coordinates in (x,y,z) space -> (3,numnode)
  \param unrm               (out)       : unit normal vector at (r,s)
  */
  template <DRT::Element::DiscretizationType DISTYPE, class M1, class M2>
  inline void computeNormalToSurfaceElementViaShape(
      const M1& deriv, const M2& xyze, LINALG::Matrix<3, 1>& unrm)
  {
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // compute dXYZ / drs
    static LINALG::Matrix<3, 2> dxyzdrs;
    dxyzdrs.Clear();
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 2; j++)
        for (int k = 0; k < numnode; k++) dxyzdrs(i, j) += xyze(i, k) * deriv(j, k);

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
    static LINALG::Matrix<2, 2> metrictensor;
    metrictensor.MultiplyTN(dxyzdrs, dxyzdrs);

    unrm(0) = dxyzdrs(1, 0) * dxyzdrs(2, 1) - dxyzdrs(2, 0) * dxyzdrs(1, 1);
    unrm(1) = dxyzdrs(2, 0) * dxyzdrs(0, 1) - dxyzdrs(0, 0) * dxyzdrs(2, 1);
    unrm(2) = dxyzdrs(0, 0) * dxyzdrs(1, 1) - dxyzdrs(1, 0) * dxyzdrs(0, 1);

    // scale
    unrm.Scale(1.0 / unrm.Norm2());

    return;
  }


  /*!
  \brief calculate surface element normal vector in gausspoint (r,s)

  \param deriv              (in)        : shape functions derivative in (r,s) -> (2,numnode)
  \param xyze               (in)        : element nodal coordinates in (x,y,z) space -> (3,numnode)
  \param unrm               (out)       : unit normal vector at (r,s)
  \param utang1             (out)       : unit tangential vector at (r,s)
  \param utang2             (out)       : unit tangential vector at (r,s)
  */
  template <DRT::Element::DiscretizationType DISTYPE, class M1, class M2>
  inline void computeNormalAndTangentToSurfaceElementViaShape(const M1& deriv, const M2& xyze,
      LINALG::Matrix<3, 1>& unrm, LINALG::Matrix<3, 1>& utang1, LINALG::Matrix<3, 1>& utang2)
  {
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // compute dXYZ / drs
    static LINALG::Matrix<3, 2> dxyzdrs;
    dxyzdrs.Clear();
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 2; j++)
        for (int k = 0; k < numnode; k++) dxyzdrs(i, j) += xyze(i, k) * deriv(j, k);

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
    static LINALG::Matrix<2, 2> metrictensor;
    metrictensor.MultiplyTN(dxyzdrs, dxyzdrs);

    unrm(0) = dxyzdrs(1, 0) * dxyzdrs(2, 1) - dxyzdrs(2, 0) * dxyzdrs(1, 1);
    unrm(1) = dxyzdrs(2, 0) * dxyzdrs(0, 1) - dxyzdrs(0, 0) * dxyzdrs(2, 1);
    unrm(2) = dxyzdrs(0, 0) * dxyzdrs(1, 1) - dxyzdrs(1, 0) * dxyzdrs(0, 1);

    // scale
    unrm.Scale(1.0 / unrm.Norm2());

    utang1(0) = dxyzdrs(0, 0);
    utang1(1) = dxyzdrs(1, 0);
    utang1(2) = dxyzdrs(2, 0);

    utang2(0) = dxyzdrs(0, 1);
    utang2(1) = dxyzdrs(1, 1);
    utang2(2) = dxyzdrs(2, 1);

    // scale
    utang1.Scale(1.0 / utang1.Norm2());
    utang2.Scale(1.0 / utang2.Norm2());

    return;
  }


  /*!
  \brief calculate element normal in gausspoint (r,s)
  \param xyze_surfaceElement   (in)   nodal position array (3,numnode)
  \param eleCoord              (in)   node in element coordinates (r, s)
  \param unitnormalvec         (out)  unit normal vector
  */
  template <DRT::Element::DiscretizationType DISTYPE, class M, class V>
  inline void computeNormalToSurfaceElementT(
      const M& xyze_surfaceElement, const V& eleCoord, LINALG::Matrix<3, 1>& unitnormalvec)
  {
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    static LINALG::Matrix<2, numnode> deriv_boundary;
    DRT::UTILS::shape_function_2D_deriv1(deriv_boundary, eleCoord(0), eleCoord(1), DISTYPE);

    computeNormalToSurfaceElementViaShape<DISTYPE>(
        deriv_boundary, xyze_surfaceElement, unitnormalvec);
    return;
  }


  /*!
  \brief calculate element normal in gausspoint (r,s)
  \param xyze_surfaceElement   (in)   nodal position array (3,numnode)
  \param eleCoord              (in)   node in element coordinates (r, s)
  \param unitnormalvec         (out)  unit normal vector
  \param unittang1vec          (out):  unit tangential vector
  \param unittang2vec          (out):  unit tangential vector
  */
  template <DRT::Element::DiscretizationType DISTYPE, class M>
  inline void computeNormalAndTangentToSurfaceElementT(const M& xyze_surfaceElement,
      const LINALG::Matrix<2, 1>& eleCoord, LINALG::Matrix<3, 1>& unitnormalvec,
      LINALG::Matrix<3, 1>& unittang1vec, LINALG::Matrix<3, 1>& unittang2vec)
  {
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    static LINALG::Matrix<2, numnode> deriv_boundary;
    DRT::UTILS::shape_function_2D_deriv1(deriv_boundary, eleCoord(0), eleCoord(1), DISTYPE);

    computeNormalAndTangentToSurfaceElementViaShape<DISTYPE>(
        deriv_boundary, xyze_surfaceElement, unitnormalvec, unittang1vec, unittang2vec);
    return;
  }


  /*!
  \brief calculate element normal in gausspoint (r,s)
  \param surfaceElement        (in) :  surface element
  \param xyze_surfaceElement   (in) :  nodal position array (3,numnode)
  \param eleCoord              (in) :  node in element coordinates (r, s)
  \param unitnormalvec         (out):  unit normal vector
  */
  template <class M, class V>
  inline void computeNormalToSurfaceElement(const DRT::Element::DiscretizationType distype,
      const M& xyze_surfaceElement, const V& eleCoord, LINALG::Matrix<3, 1>& unitnormalvec)
  {
    switch (distype)
    {
      case DRT::Element::quad4:
        computeNormalToSurfaceElementT<DRT::Element::quad4>(
            xyze_surfaceElement, eleCoord, unitnormalvec);
        break;
      case DRT::Element::quad8:
        computeNormalToSurfaceElementT<DRT::Element::quad8>(
            xyze_surfaceElement, eleCoord, unitnormalvec);
        break;
      case DRT::Element::quad9:
        computeNormalToSurfaceElementT<DRT::Element::quad9>(
            xyze_surfaceElement, eleCoord, unitnormalvec);
        break;
      case DRT::Element::tri3:
        computeNormalToSurfaceElementT<DRT::Element::tri3>(
            xyze_surfaceElement, eleCoord, unitnormalvec);
        break;
      case DRT::Element::tri6:
        computeNormalToSurfaceElementT<DRT::Element::tri6>(
            xyze_surfaceElement, eleCoord, unitnormalvec);
        break;
      default:
        std::cout << DRT::DistypeToString(distype) << " not found!" << std::endl;
        dserror("please add your distype!");
    }
    return;
  }


  /*!
  \brief calculate element normal in gausspoint (r,s)
  \param surfaceElement        (in) :  surface element
  \param xyze_surfaceElement   (in) :  nodal position array (3,numnode)
  \param eleCoord              (in) :  node in element coordinates (r, s)
  \param unitnormalvec         (out):  unit normal vector
  \param unittang1vec          (out):  unit tangential vector
  \param unittang2vec          (out):  unit tangential vector
  */
  template <class M>
  inline void computeNormalAndTangentToSurfaceElement(
      const DRT::Element::DiscretizationType distype, const M& xyze_surfaceElement,
      const LINALG::Matrix<2, 1>& eleCoord, LINALG::Matrix<3, 1>& unitnormalvec,
      LINALG::Matrix<3, 1>& unittang1vec, LINALG::Matrix<3, 1>& unittang2vec)
  {
    switch (distype)
    {
      case DRT::Element::quad4:
        computeNormalAndTangentToSurfaceElementT<DRT::Element::quad4>(
            xyze_surfaceElement, eleCoord, unitnormalvec, unittang1vec, unittang2vec);
        break;
      case DRT::Element::tri3:
        computeNormalAndTangentToSurfaceElementT<DRT::Element::tri3>(
            xyze_surfaceElement, eleCoord, unitnormalvec, unittang1vec, unittang2vec);
        break;
      default:
        std::cout << DRT::DistypeToString(distype) << " not found!" << std::endl;
        dserror("please add your distype!");
    }
    return;
  }

}  // namespace GEO


#endif /*ELEMENT_NORMALS_H_NORMALS_H_*/
