/*----------------------------------------------------------------------*/
/*! \file

\brief collection of templated service methods for intersection computations

\level 1


*----------------------------------------------------------------------*/


#ifndef BACI_DISCRETIZATION_GEOMETRY_INTERSECTION_SERVICE_TEMPLATES_H
#define BACI_DISCRETIZATION_GEOMETRY_INTERSECTION_SERVICE_TEMPLATES_H


#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_discretization_geometry_geo_utils.H"
#include "baci_discretization_geometry_intersection_math.H"
#include "baci_lib_element.H"


namespace CORE::GEO
{
  //! template of extended aligned bounding boxes
  template <int ndim, CORE::FE::CellType distype, class M>
  static inline CORE::LINALG::Matrix<3, 2> computeFastXAABBT(
      const M& xyze, const CORE::GEO::EleGeoType eleGeoType)
  {
    CORE::LINALG::Matrix<3, 2> XAABB;

    // first node
    for (int dim = 0; dim < 3; ++dim)
    {
      XAABB(dim, 0) = xyze(dim, 0) - TOL7;
      XAABB(dim, 1) = xyze(dim, 0) + TOL7;
    }
    // remaining nodes
    const int numNodes = CORE::FE::num_nodes<distype>;
    for (int i = 1; i < numNodes; ++i)
      for (int dim = 0; dim < ndim; dim++)
      {
        XAABB(dim, 0) = std::min(XAABB(dim, 0), xyze(dim, i) - TOL7);
        XAABB(dim, 1) = std::max(XAABB(dim, 1), xyze(dim, i) + TOL7);
      }

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
  CORE::LINALG::Matrix<3, 2> computeFastXAABB(
      CORE::FE::CellType distype, const M& xyze, const CORE::GEO::EleGeoType eleGeoType)
  {
    switch (distype)
    {
      case CORE::FE::CellType::hex8:
        return computeFastXAABBT<3, CORE::FE::CellType::hex8>(xyze, eleGeoType);
      case CORE::FE::CellType::quad4:
        return computeFastXAABBT<3, CORE::FE::CellType::quad4>(xyze, eleGeoType);
      case CORE::FE::CellType::hex20:
        return computeFastXAABBT<3, CORE::FE::CellType::hex20>(xyze, eleGeoType);
      case CORE::FE::CellType::hex27:
        return computeFastXAABBT<3, CORE::FE::CellType::hex27>(xyze, eleGeoType);
      case CORE::FE::CellType::tet4:
        return computeFastXAABBT<3, CORE::FE::CellType::tet4>(xyze, eleGeoType);
      case CORE::FE::CellType::tet10:
        return computeFastXAABBT<3, CORE::FE::CellType::tet10>(xyze, eleGeoType);
      case CORE::FE::CellType::line2:
        return computeFastXAABBT<3, CORE::FE::CellType::line2>(xyze, eleGeoType);
      case CORE::FE::CellType::line3:
        return computeFastXAABBT<3, CORE::FE::CellType::line3>(xyze, eleGeoType);
      case CORE::FE::CellType::tri3:
        return computeFastXAABBT<3, CORE::FE::CellType::tri3>(xyze, eleGeoType);
      case CORE::FE::CellType::tri6:
        return computeFastXAABBT<3, CORE::FE::CellType::tri6>(xyze, eleGeoType);
      case CORE::FE::CellType::quad8:
        return computeFastXAABBT<3, CORE::FE::CellType::quad8>(xyze, eleGeoType);
      case CORE::FE::CellType::quad9:
        return computeFastXAABBT<3, CORE::FE::CellType::quad9>(xyze, eleGeoType);
      case CORE::FE::CellType::pyramid5:
        return computeFastXAABBT<3, CORE::FE::CellType::pyramid5>(xyze, eleGeoType);
      default:
        std::cout << CORE::FE::CellTypeToString(distype) << std::endl;
        dserror("add your distype to this switch!");
    }
    return CORE::LINALG::Matrix<3, 2>(true);
  }


}  // namespace CORE::GEO


#endif  // DISCRETIZATION_GEOMETRY_INTERSECTION_SERVICE_TEMPLATES_H
