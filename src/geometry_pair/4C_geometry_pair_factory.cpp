/*----------------------------------------------------------------------*/
/*! \file

\brief functions to create geometry pairs.

\level 1
*/


#include "4C_geometry_pair_factory.hpp"

#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_line_to_3D_evaluation_data.hpp"
#include "4C_geometry_pair_line_to_surface_evaluation_data.hpp"
#include "4C_geometry_pair_line_to_surface_gauss_point_projection.hpp"
#include "4C_geometry_pair_line_to_surface_segmentation.hpp"
#include "4C_geometry_pair_line_to_volume_gauss_point_projection.hpp"
#include "4C_geometry_pair_line_to_volume_segmentation.hpp"
#include "4C_geometry_pair_scalar_types.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
template <typename ScalarType, typename Line, typename Volume>
Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToVolumeFactory(
    const Core::Elements::Element* element1, const Core::Elements::Element* element2,
    const Teuchos::RCP<GeometryEvaluationDataBase>& geometry_evaluation_data)
{
  // Cast the geometry evaluation data to the correct format.
  auto line_to_3d_evaluation_data =
      Teuchos::rcp_dynamic_cast<LineTo3DEvaluationData>(geometry_evaluation_data, true);

  // Get the strategy for line to volume interaction.
  Inpar::GEOMETRYPAIR::LineTo3DStrategy strategy = line_to_3d_evaluation_data->get_strategy();

  // Create the class depending on the strategy.
  switch (strategy)
  {
    case Inpar::GEOMETRYPAIR::LineTo3DStrategy::
        gauss_point_projection_without_boundary_segmentation:
    case Inpar::GEOMETRYPAIR::LineTo3DStrategy::gauss_point_projection_boundary_segmentation:
      return Teuchos::rcp(
          new GeometryPairLineToVolumeGaussPointProjection<ScalarType, Line, Volume>(
              element1, element2, line_to_3d_evaluation_data));
    case Inpar::GEOMETRYPAIR::LineTo3DStrategy::segmentation:
      return Teuchos::rcp(new GeometryPairLineToVolumeSegmentation<ScalarType, Line, Volume>(
          element1, element2, line_to_3d_evaluation_data));
    default:
    {
      FOUR_C_THROW(
          "The given geometry pair strategy is not universally valid. You might want to create "
          "your pair directly if you need certain features (for example cross section "
          "projection)!");
      return Teuchos::null;
    }
  }
}


/**
 * Explicit template initialization of factory function.
 */
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToVolumeFactory<
    double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex8>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToVolumeFactory<
    double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex20>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToVolumeFactory<
    double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex27>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToVolumeFactory<
    double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet4>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToVolumeFactory<
    double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet10>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToVolumeFactory<
    double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_nurbs27>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);


/**
 *
 */
template <typename ScalarType, typename Line, typename Surface>
Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactory(
    const Core::Elements::Element* element1, const Core::Elements::Element* element2,
    const Teuchos::RCP<GeometryEvaluationDataBase>& geometry_evaluation_data)
{
  // Cast the geometry evaluation data to the correct format.
  auto line_to_surface_evaluation_data =
      Teuchos::rcp_dynamic_cast<LineToSurfaceEvaluationData>(geometry_evaluation_data, true);

  // Get the strategy for line to volume interaction.
  Inpar::GEOMETRYPAIR::LineTo3DStrategy strategy = line_to_surface_evaluation_data->get_strategy();

  // Create the class depending on the strategy.
  switch (strategy)
  {
    case Inpar::GEOMETRYPAIR::LineTo3DStrategy::
        gauss_point_projection_without_boundary_segmentation:
    case Inpar::GEOMETRYPAIR::LineTo3DStrategy::gauss_point_projection_boundary_segmentation:
      return Teuchos::rcp(
          new GeometryPairLineToSurfaceGaussPointProjection<ScalarType, Line, Surface>(
              element1, element2, line_to_surface_evaluation_data));
    case Inpar::GEOMETRYPAIR::LineTo3DStrategy::segmentation:
      return Teuchos::rcp(new GeometryPairLineToSurfaceSegmentation<ScalarType, Line, Surface>(
          element1, element2, line_to_surface_evaluation_data));
    default:
    {
      FOUR_C_THROW("The given geometry pair strategy is not valid.");
      return Teuchos::null;
    }
  }
}

/**
 *
 */
template <typename ScalarType, typename Line, typename Surface>
Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD(
    const Core::Elements::Element* element1, const Core::Elements::Element* element2,
    const Teuchos::RCP<GeometryEvaluationDataBase>& geometry_evaluation_data)
{
  // Create the internal pair as double.
  auto internal_geometry_pair_double =
      Teuchos::rcp_dynamic_cast<GeometryPairLineToSurface<double, Line, Surface>>(
          GeometryPairLineToSurfaceFactory<double, Line, Surface>(
              element1, element2, geometry_evaluation_data),
          true);

  // Create the wrapper.
  return Teuchos::rcp(new GeometryPairLineToSurfaceFADWrapper<ScalarType, Line, Surface>(
      element1, element2, internal_geometry_pair_double));
}


/**
 * Explicit template initialization of factory function.
 */
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair>
GEOMETRYPAIR::GeometryPairLineToSurfaceFactory<double, GEOMETRYPAIR::t_line2, GEOMETRYPAIR::t_tri3>(
    const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair>
GEOMETRYPAIR::GeometryPairLineToSurfaceFactory<double, GEOMETRYPAIR::t_line2, GEOMETRYPAIR::t_tri6>(
    const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactory<
    double, GEOMETRYPAIR::t_line2, GEOMETRYPAIR::t_quad4>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactory<
    double, GEOMETRYPAIR::t_line2, GEOMETRYPAIR::t_quad8>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactory<
    double, GEOMETRYPAIR::t_line2, GEOMETRYPAIR::t_quad9>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactory<
    double, GEOMETRYPAIR::t_line2, GEOMETRYPAIR::t_nurbs9>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);

template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type_1st_order, GEOMETRYPAIR::t_line2,
    GEOMETRYPAIR::t_tri3>(const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type_1st_order, GEOMETRYPAIR::t_line2,
    GEOMETRYPAIR::t_tri6>(const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type_1st_order, GEOMETRYPAIR::t_line2,
    GEOMETRYPAIR::t_quad4>(const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type_1st_order, GEOMETRYPAIR::t_line2,
    GEOMETRYPAIR::t_quad8>(const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type_1st_order, GEOMETRYPAIR::t_line2,
    GEOMETRYPAIR::t_quad9>(const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type_fixed_size_1st_order<GEOMETRYPAIR::t_line2,
        GEOMETRYPAIR::t_nurbs9>,
    GEOMETRYPAIR::t_line2, GEOMETRYPAIR::t_nurbs9>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);

template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type, GEOMETRYPAIR::t_line2, GEOMETRYPAIR::t_tri3>(
    const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type, GEOMETRYPAIR::t_line2, GEOMETRYPAIR::t_tri6>(
    const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type, GEOMETRYPAIR::t_line2, GEOMETRYPAIR::t_quad4>(
    const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type, GEOMETRYPAIR::t_line2, GEOMETRYPAIR::t_quad8>(
    const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type, GEOMETRYPAIR::t_line2, GEOMETRYPAIR::t_quad9>(
    const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type_fixed_size<GEOMETRYPAIR::t_line2,
        GEOMETRYPAIR::t_nurbs9>,
    GEOMETRYPAIR::t_line2, GEOMETRYPAIR::t_nurbs9>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);

template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactory<
    double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri3>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactory<
    double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri6>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactory<
    double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad4>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactory<
    double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad8>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactory<
    double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad9>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactory<
    double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_nurbs9>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);

template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type_1st_order, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tri3>(const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type_1st_order, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tri6>(const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type_1st_order, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_quad4>(const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type_1st_order, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_quad8>(const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type_1st_order, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_quad9>(const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type_fixed_size_1st_order<GEOMETRYPAIR::t_hermite,
        GEOMETRYPAIR::t_nurbs9>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_nurbs9>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);

template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri3>(
    const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri6>(
    const Core::Elements::Element*, const Core::Elements::Element*,
    const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair>
GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<GEOMETRYPAIR::line_to_surface_patch_scalar_type,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad4>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair>
GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<GEOMETRYPAIR::line_to_surface_patch_scalar_type,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad8>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair>
GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<GEOMETRYPAIR::line_to_surface_patch_scalar_type,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad9>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToSurfaceFactoryFAD<
    GEOMETRYPAIR::line_to_surface_patch_scalar_type_fixed_size<GEOMETRYPAIR::t_hermite,
        GEOMETRYPAIR::t_nurbs9>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_nurbs9>(const Core::Elements::Element*,
    const Core::Elements::Element*, const Teuchos::RCP<GeometryEvaluationDataBase>&);

FOUR_C_NAMESPACE_CLOSE
