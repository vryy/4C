/*!----------------------------------------------------------------------

\brief collection of templated service methods for intersection computations

\level 1

\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/


#ifndef INTERSECTION_SERVICE_TEMPLATES_H_
#define INTERSECTION_SERVICE_TEMPLATES_H_


#include "intersection_math.H"
#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "geo_utils.H"


namespace GEO
{
  //! checks if a given surface element lies in a plane
  template <class M>
  static inline bool SurfaceElementIsPlane(const M& xyze)
  {
    static LINALG::Matrix<3, 1> elenorm;
    // xyze 3 x 4

    if (xyze.N() == 3) return true;

    elenorm(0) = -(xyze(2, 0) * xyze(1, 1)) + xyze(1, 0) * xyze(2, 1) + xyze(2, 0) * xyze(1, 2) -
                 xyze(2, 1) * xyze(1, 2) - xyze(1, 0) * xyze(2, 2) + xyze(1, 1) * xyze(2, 2);
    elenorm(1) = xyze(2, 0) * xyze(0, 1) - xyze(0, 0) * xyze(2, 1) - xyze(2, 0) * xyze(0, 2) +
                 xyze(2, 1) * xyze(0, 2) + xyze(0, 0) * xyze(2, 2) - xyze(0, 1) * xyze(2, 2);
    elenorm(2) = -(xyze(1, 0) * xyze(0, 1)) + xyze(0, 0) * xyze(1, 1) + xyze(1, 0) * xyze(0, 2) -
                 xyze(1, 1) * xyze(0, 2) - xyze(0, 0) * xyze(1, 2) + xyze(0, 1) * xyze(1, 2);

    for (int i = 3; i < xyze.N(); i++)
    {
      double distance =
          (-(elenorm(0) * (xyze(0, 0) - xyze(0, i))) + elenorm(1) * (-xyze(1, 0) + xyze(1, i)) -
              elenorm(2) * (xyze(2, 0) - xyze(2, i))) /
          (-elenorm.Norm2());

      if (fabs(distance) >= GEO::TOL7) return false;
    }
    return true;
  }



  //! checks if a given point lies in a plane surface element
  template <class M1, class M2>
  static inline bool pointInPlaneSurfaceElement(const M1& xyze, const M2& xyze_points)
  {
    static LINALG::Matrix<3, 1> elenorm;
    // xyze 3 x ...

    elenorm(0) = -(xyze(2, 0) * xyze(1, 1)) + xyze(1, 0) * xyze(2, 1) + xyze(2, 0) * xyze(1, 2) -
                 xyze(2, 1) * xyze(1, 2) - xyze(1, 0) * xyze(2, 2) + xyze(1, 1) * xyze(2, 2);
    elenorm(1) = xyze(2, 0) * xyze(0, 1) - xyze(0, 0) * xyze(2, 1) - xyze(2, 0) * xyze(0, 2) +
                 xyze(2, 1) * xyze(0, 2) + xyze(0, 0) * xyze(2, 2) - xyze(0, 1) * xyze(2, 2);
    elenorm(2) = -(xyze(1, 0) * xyze(0, 1)) + xyze(0, 0) * xyze(1, 1) + xyze(1, 0) * xyze(0, 2) -
                 xyze(1, 1) * xyze(0, 2) - xyze(0, 0) * xyze(1, 2) + xyze(0, 1) * xyze(1, 2);

    for (int i = 0; i < 2; i++)
    {
      double distance = (-(elenorm(0) * (xyze(0, 0) - xyze_points(0, i))) +
                            elenorm(1) * (-xyze(1, 0) + xyze_points(1, i)) -
                            elenorm(2) * (xyze(2, 0) - xyze_points(2, i))) /
                        (-elenorm.Norm2());

      if (fabs(distance) >= GEO::TOL7) return false;
    }
    return true;
  }


  //! checks if a some test points lie in a plane defined by 3 points
  template <class M1, class M2>
  static inline bool pointsInPlaneSurfaceElement(
      const M1& xyze,        ///< (3,3)             plane made of 3 points
      const M2& xyze_points  ///< (3,numtestpoints) test points
  )
  {
    // normal vector on plane xyze
    static LINALG::Matrix<3, 1> elenorm;
    // xyze 3 x ...

    // cross product -> normal vector
    elenorm(0) = -(xyze(2, 0) * xyze(1, 1)) + xyze(1, 0) * xyze(2, 1) + xyze(2, 0) * xyze(1, 2) -
                 xyze(2, 1) * xyze(1, 2) - xyze(1, 0) * xyze(2, 2) + xyze(1, 1) * xyze(2, 2);
    elenorm(1) = xyze(2, 0) * xyze(0, 1) - xyze(0, 0) * xyze(2, 1) - xyze(2, 0) * xyze(0, 2) +
                 xyze(2, 1) * xyze(0, 2) + xyze(0, 0) * xyze(2, 2) - xyze(0, 1) * xyze(2, 2);
    elenorm(2) = -(xyze(1, 0) * xyze(0, 1)) + xyze(0, 0) * xyze(1, 1) + xyze(1, 0) * xyze(0, 2) -
                 xyze(1, 1) * xyze(0, 2) - xyze(0, 0) * xyze(1, 2) + xyze(0, 1) * xyze(1, 2);

    const int numtestpoints = xyze_points.Cols();
    for (int i = 0; i < numtestpoints; i++)
    {
      const double distance = (elenorm(0) * (xyze_points(0, i) - xyze(0, 0)) +
                                  elenorm(1) * (xyze_points(1, i) - xyze(1, 0)) +
                                  elenorm(2) * (xyze_points(2, i) - xyze(2, 0))) /
                              elenorm.Norm2();

      if (fabs(distance) >= GEO::TOL7)
      {
        return false;
      }
    }
    return true;
  }


  //! checks for a coplanar 4-node tetraeder and not corectly numbered tetrahedra
  template <class M1, class M2>
  inline bool checkDegenerateTet(const int numTetNodes,
      M1& xyze,      ///< tet made of 4 points, xyze(nsd = 3, numnode = 4)
      M2& xyze_phys  ///< tet made of 4 points, xyze(nsd = 3, numnode = 4)
  )
  {
    // create planes consisting of 3 nodes each
    LINALG::Matrix<3, 1> v01;
    LINALG::Matrix<3, 1> v02;
    LINALG::Matrix<3, 1> v03;
    //    LINALG::Matrix<3,1> v12;
    //    LINALG::Matrix<3,1> v13;
    for (int isd = 0; isd < 3; ++isd) v01(isd, 0) = xyze(isd, 1) - xyze(isd, 0);
    for (int isd = 0; isd < 3; ++isd) v02(isd, 0) = xyze(isd, 2) - xyze(isd, 0);
    for (int isd = 0; isd < 3; ++isd) v03(isd, 0) = xyze(isd, 3) - xyze(isd, 0);
    /*    for (int isd = 0; isd < 3; ++isd)  v12(isd,0) = xyze(isd,2)-xyze(isd,1);
        for (int isd = 0; isd < 3; ++isd)  v13(isd,0) = xyze(isd,3)-xyze(isd,1);
    */
    // create 4 normal vectors to each tet surface plane
    const LINALG::Matrix<3, 1> nplane012 = GEO::computeCrossProduct(v01, v02);
    /*    const LINALG::Matrix<3,1> nplane013 = GEO::computeCrossProduct(v01,v03);
        const LINALG::Matrix<3,1> nplane023 = GEO::computeCrossProduct(v02,v03);
        const LINALG::Matrix<3,1> nplane123 = GEO::computeCrossProduct(v12,v13);
    */
    // compute norm (area) of each plane
    const double norm012 = nplane012.Norm2();
    /*    const double norm013 = nplane013.Norm2();
        const double norm023 = nplane023.Norm2();
        const double norm123 = nplane123.Norm2();
    */

    const double vol_tet =
        (std::abs(nplane012(0) * v03(0) + nplane012(1) * v03(1) + nplane012(2) * v03(2))) / 6.0;
    // smallest area
    if (vol_tet < 1e-10) return true;


    // compute normal distance of point to plane of the three remaining points
    double distance = 0.0;
    distance = (nplane012(0) * v03(0) + nplane012(1) * v03(1) + nplane012(2) * v03(2)) / norm012;

    if (fabs(distance) < GEO::TOL7)
    {
      return true;
    }

    // tet numbering wrong exchange 1 with 3
    if (distance <= (-1.0) * GEO::TOL7)
    {
      for (int i = 0; i < 3; i++)
      {
        double point_coord = xyze(i, 1);
        xyze(i, 1) = xyze(i, 3);
        xyze(i, 3) = point_coord;

        point_coord = xyze_phys(i, 1);
        xyze_phys(i, 1) = xyze_phys(i, 3);
        xyze_phys(i, 3) = point_coord;
      }
      return false;
    }
    return false;
  }



  //! template of extended aligned bounding boxes
  template <int ndim, DRT::Element::DiscretizationType distype, class M>
  static inline LINALG::Matrix<3, 2> computeFastXAABBT(
      const M& xyze, const GEO::EleGeoType eleGeoType)
  {
    LINALG::Matrix<3, 2> XAABB;

    // first node
    for (int dim = 0; dim < 3; ++dim)
    {
      XAABB(dim, 0) = xyze(dim, 0) - TOL7;
      XAABB(dim, 1) = xyze(dim, 0) + TOL7;
    }
    // remaining nodes
    const int numNodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
    for (int i = 1; i < numNodes; ++i)
      for (int dim = 0; dim < ndim; dim++)
      {
        XAABB(dim, 0) = std::min(XAABB(dim, 0), xyze(dim, i) - TOL7);
        XAABB(dim, 1) = std::max(XAABB(dim, 1), xyze(dim, i) + TOL7);
      }

    /*

    if(eleGeoType == HIGHERORDER)
    {
      double maxDistance = fabs(XAABB(0,1) - XAABB(0,0));
      for(int dim = 1; dim < ndim; ++dim)
        maxDistance = std::max(maxDistance, fabs(XAABB(dim,1)-XAABB(dim,0)) );

      // subtracts a certain percent of maximal distance to minX, minY, minZ
      // adds a certain percent of the maximal distance to maxX, maxY, maxZ
      const double extension = 0.2*maxDistance;
      for(int dim=0; dim < ndim; ++dim)
      {
          XAABB(dim, 0) -= extension;
          XAABB(dim, 1) += extension;
      }
    }

    */
    /*
    printf("\n");
    printf("axis-aligned bounding box:\n minX = %f\n minY = %f\n minZ = %f\n maxX = %f\n maxY = %f\n
    maxZ = %f\n", XAABB(0,0), XAABB(1,0), XAABB(2,0), XAABB(0,1), XAABB(1,1), XAABB(2,1));
    printf("\n");
    */

    return XAABB;
  }



  /*!
  \brief Computes a rough overestimating extended
         axis-aligned bounding box for an element (XAABB)
  \param distype        (in)  distype of element
  \param xyze           (in)  nodal position array (3,numnode)
  \param eleGeoType     (in)  element geometric type CARTESIAN LINEAR or HIGHERORDER
  \return extended axis-aligned bounding box  (XAABB) for an element
   */
  template <class M>
  LINALG::Matrix<3, 2> computeFastXAABB(
      DRT::Element::DiscretizationType distype, const M& xyze, const GEO::EleGeoType eleGeoType)
  {
    switch (distype)
    {
      case DRT::Element::hex8:
        return computeFastXAABBT<3, DRT::Element::hex8>(xyze, eleGeoType);
      case DRT::Element::quad4:
        return computeFastXAABBT<3, DRT::Element::quad4>(xyze, eleGeoType);
      case DRT::Element::hex20:
        return computeFastXAABBT<3, DRT::Element::hex20>(xyze, eleGeoType);
      case DRT::Element::hex27:
        return computeFastXAABBT<3, DRT::Element::hex27>(xyze, eleGeoType);
      case DRT::Element::tet4:
        return computeFastXAABBT<3, DRT::Element::tet4>(xyze, eleGeoType);
      case DRT::Element::tet10:
        return computeFastXAABBT<3, DRT::Element::tet10>(xyze, eleGeoType);
      case DRT::Element::line2:
        return computeFastXAABBT<3, DRT::Element::line2>(xyze, eleGeoType);
      case DRT::Element::line3:
        return computeFastXAABBT<3, DRT::Element::line3>(xyze, eleGeoType);
      case DRT::Element::tri3:
        return computeFastXAABBT<3, DRT::Element::tri3>(xyze, eleGeoType);
      case DRT::Element::tri6:
        return computeFastXAABBT<3, DRT::Element::tri6>(xyze, eleGeoType);
      case DRT::Element::quad8:
        return computeFastXAABBT<3, DRT::Element::quad8>(xyze, eleGeoType);
      case DRT::Element::quad9:
        return computeFastXAABBT<3, DRT::Element::quad9>(xyze, eleGeoType);
      case DRT::Element::pyramid5:
        return computeFastXAABBT<3, DRT::Element::pyramid5>(xyze, eleGeoType);
      default:
        std::cout << DistypeToString(distype) << std::endl;
        dserror("add your distype to this switch!");
    }
    return LINALG::Matrix<3, 2>(true);
  }


  //! checks if two XAABB's intersect
  template <int dim>
  static inline bool intersectionOfXAABB(
      const LINALG::Matrix<3, 2>& cutterXAABB, const LINALG::Matrix<3, 2>& xfemXAABB)
  {
    // check XAABB for intersection
    for (int i = 0; i < dim; i++)
    {
      if (!(((cutterXAABB(i, 0) > (xfemXAABB(i, 0) - GEO::TOL7)) &&
                (cutterXAABB(i, 0) < (xfemXAABB(i, 1) + GEO::TOL7))) ||
              ((cutterXAABB(i, 1) > (xfemXAABB(i, 0) - GEO::TOL7)) &&
                  (cutterXAABB(i, 1) < (xfemXAABB(i, 1) + GEO::TOL7))) ||
              ((xfemXAABB(i, 0) > (cutterXAABB(i, 0) - GEO::TOL7)) &&
                  (xfemXAABB(i, 0) < (cutterXAABB(i, 1) + GEO::TOL7))) ||
              ((xfemXAABB(i, 1) > (cutterXAABB(i, 0) - GEO::TOL7)) &&
                  (xfemXAABB(i, 1) < (cutterXAABB(i, 1) + GEO::TOL7)))))
        return false;
    }
    return true;
  }


}  // namespace GEO


#endif /*INTERSECTION_SERVICE_TEMPLATES_H_*/
