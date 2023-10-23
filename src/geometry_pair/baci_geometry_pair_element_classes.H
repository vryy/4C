/*----------------------------------------------------------------------*/
/*! \file

\brief Utility classes for the geometry pairs.

\level 1
*/
// End doxygen header.


#ifndef BACI_GEOMETRY_PAIR_ELEMENT_CLASSES_H
#define BACI_GEOMETRY_PAIR_ELEMENT_CLASSES_H


#include "baci_geometry_pair_constants.H"
#include "baci_lib_discret.H"
#include "baci_linalg_utils_densematrix_inverse.H"
#include "baci_utils_exceptions.H"

#include <vector>


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
  template <DRT::Element::DiscretizationType discretization>
  struct ElementDiscretizationToGeometryType
  {
    static const GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::none;
  };
  template <>
  struct ElementDiscretizationToGeometryType<DRT::Element::DiscretizationType::line2>
  {
    static const GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::line;
  };
  template <>
  struct ElementDiscretizationToGeometryType<DRT::Element::DiscretizationType::tri3>
  {
    static const GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::triangle;
  };
  template <>
  struct ElementDiscretizationToGeometryType<DRT::Element::DiscretizationType::tri6>
  {
    static const GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::triangle;
  };
  template <>
  struct ElementDiscretizationToGeometryType<DRT::Element::DiscretizationType::quad4>
  {
    static const GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::quad;
  };
  template <>
  struct ElementDiscretizationToGeometryType<DRT::Element::DiscretizationType::quad8>
  {
    static const GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::quad;
  };
  template <>
  struct ElementDiscretizationToGeometryType<DRT::Element::DiscretizationType::quad9>
  {
    static const GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::quad;
  };
  template <>
  struct ElementDiscretizationToGeometryType<DRT::Element::DiscretizationType::nurbs9>
  {
    static const GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::quad;
  };
  template <>
  struct ElementDiscretizationToGeometryType<DRT::Element::DiscretizationType::hex8>
  {
    static const GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::hexahedron;
  };
  template <>
  struct ElementDiscretizationToGeometryType<DRT::Element::DiscretizationType::hex20>
  {
    static const GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::hexahedron;
  };
  template <>
  struct ElementDiscretizationToGeometryType<DRT::Element::DiscretizationType::hex27>
  {
    static const GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::hexahedron;
  };
  template <>
  struct ElementDiscretizationToGeometryType<DRT::Element::DiscretizationType::tet4>
  {
    static const GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::tetraeder;
  };
  template <>
  struct ElementDiscretizationToGeometryType<DRT::Element::DiscretizationType::tet10>
  {
    static const GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::tetraeder;
  };
  template <>
  struct ElementDiscretizationToGeometryType<DRT::Element::DiscretizationType::nurbs27>
  {
    static const GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::hexahedron;
  };

}  // namespace GEOMETRYPAIR


#endif
