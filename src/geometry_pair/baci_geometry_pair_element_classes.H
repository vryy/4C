/*----------------------------------------------------------------------*/
/*! \file

\brief Utility classes for the geometry pairs.

\level 1
*/
// End doxygen header.


#ifndef BACI_GEOMETRY_PAIR_ELEMENT_CLASSES_H
#define BACI_GEOMETRY_PAIR_ELEMENT_CLASSES_H


#include "baci_config.H"

#include "baci_geometry_pair_constants.H"
#include "baci_lib_discret.H"
#include "baci_linalg_utils_densematrix_inverse.H"
#include "baci_utils_exceptions.H"

#include <vector>

BACI_NAMESPACE_OPEN

namespace GEOMETRYPAIR
{
  /**
   * \brief Geometry discretization type of element.
   */
  enum class DiscretizationTypeGeometry
  {
    //! none
    none,
    //! 1D curve
    line,
    //! triangle
    triangle,
    //! quadrilateral
    quad,
    //! hexahedron
    hexahedron,
    //! tetraeder
    tetraeder
  };

  /**
   * \brief This structure "converts" the DRT discretization type to a geometry type.
   *
   * For some geometry pairs we need to know if a geometry is a triangle / a quad / tetraeder or
   * hexahedron (linear, quadratic, ...) this structure "returns" the correct type depending on the
   * DRT discretization type of the element.
   */
  template <CORE::FE::CellType discretization>
  struct ElementDiscretizationToGeometryType
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::none;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::line2>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::line;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::tri3>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::triangle;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::tri6>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::triangle;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::quad4>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::quad;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::quad8>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::quad;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::quad9>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::quad;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::nurbs9>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::quad;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::hex8>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::hexahedron;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::hex20>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::hexahedron;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::hex27>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::hexahedron;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::tet4>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::tetraeder;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::tet10>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::tetraeder;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::nurbs27>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::hexahedron;
  };

}  // namespace GEOMETRYPAIR


BACI_NAMESPACE_CLOSE

#endif
