/*----------------------------------------------------------------------*/
/*! \file

\brief Utility classes for the geometry pairs.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_ELEMENT_CLASSES_H
#define GEOMETRY_PAIR_ELEMENT_CLASSES_H


#include "geometry_pair_constants.H"

#include "../linalg/linalg_utils_densematrix_inverse.H"
#include "../drt_lib/drt_dserror.H"

#include <vector>


namespace GEOMETRYPAIR
{
  /**
   * \brief Geometry discretization type of surface element.
   */
  enum class DiscretizationTypeSurface
  {
    //! none
    none,
    //! triangle
    triangle,
    //! quad
    quad
  };

  /**
   * \brief This structure "converts" the DRT discretization type to a surface type.
   *
   * For some geometry pairs we need to know if a surface is a triangle or a quad (linear,
   * quadratic, ...) this structure "returns" the correct type depending on the DRT discretization
   * type of the element.
   */
  template <DRT::Element::DiscretizationType discretization>
  struct ElementDiscretizationToSurfaceType
  {
    static const GEOMETRYPAIR::DiscretizationTypeSurface surface_type_ =
        GEOMETRYPAIR::DiscretizationTypeSurface::none;
  };
  template <>
  struct ElementDiscretizationToSurfaceType<DRT::Element::tri3>
  {
    static const GEOMETRYPAIR::DiscretizationTypeSurface surface_type_ =
        GEOMETRYPAIR::DiscretizationTypeSurface::triangle;
  };
  template <>
  struct ElementDiscretizationToSurfaceType<DRT::Element::tri6>
  {
    static const GEOMETRYPAIR::DiscretizationTypeSurface surface_type_ =
        GEOMETRYPAIR::DiscretizationTypeSurface::triangle;
  };
  template <>
  struct ElementDiscretizationToSurfaceType<DRT::Element::quad4>
  {
    static const GEOMETRYPAIR::DiscretizationTypeSurface surface_type_ =
        GEOMETRYPAIR::DiscretizationTypeSurface::quad;
  };
  template <>
  struct ElementDiscretizationToSurfaceType<DRT::Element::quad8>
  {
    static const GEOMETRYPAIR::DiscretizationTypeSurface surface_type_ =
        GEOMETRYPAIR::DiscretizationTypeSurface::quad;
  };
  template <>
  struct ElementDiscretizationToSurfaceType<DRT::Element::quad9>
  {
    static const GEOMETRYPAIR::DiscretizationTypeSurface surface_type_ =
        GEOMETRYPAIR::DiscretizationTypeSurface::quad;
  };


  /**
   * \brief Geometry discretization type of volume element.
   */
  enum class DiscretizationTypeVolume
  {
    //! none
    none,
    //! hexaeder
    hexaeder,
    //! tetraeder
    tetraeder
  };

  /**
   * \brief This structure "converts" the DRT discretization type to a volume type.
   *
   * For some geometry pairs we need to know if a volume is a hexaeder or a tetraeder (linear,
   * quadratic, ...) this structure "returns" the correct type depending on the DRT discretization
   * type of the element.
   */
  template <DRT::Element::DiscretizationType discretization>
  struct ElementDiscretizationToVolumeType
  {
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        GEOMETRYPAIR::DiscretizationTypeVolume::none;
  };
  template <>
  struct ElementDiscretizationToVolumeType<DRT::Element::hex8>
  {
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        GEOMETRYPAIR::DiscretizationTypeVolume::hexaeder;
  };
  template <>
  struct ElementDiscretizationToVolumeType<DRT::Element::hex20>
  {
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        GEOMETRYPAIR::DiscretizationTypeVolume::hexaeder;
  };
  template <>
  struct ElementDiscretizationToVolumeType<DRT::Element::hex27>
  {
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        GEOMETRYPAIR::DiscretizationTypeVolume::hexaeder;
  };
  template <>
  struct ElementDiscretizationToVolumeType<DRT::Element::tet4>
  {
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        GEOMETRYPAIR::DiscretizationTypeVolume::tetraeder;
  };
  template <>
  struct ElementDiscretizationToVolumeType<DRT::Element::tet10>
  {
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        GEOMETRYPAIR::DiscretizationTypeVolume::tetraeder;
  };
  template <>
  struct ElementDiscretizationToVolumeType<DRT::Element::nurbs27>
  {
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        GEOMETRYPAIR::DiscretizationTypeVolume::hexaeder;
  };

}  // namespace GEOMETRYPAIR


#endif
