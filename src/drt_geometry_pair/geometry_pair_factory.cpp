/*----------------------------------------------------------------------*/
/*! \file

\brief functions to create geometry pairs.

\level 1
\maintainer Ivo Steinbrecher
*/


#include "geometry_pair_factory.H"

#include "geometry_pair_element.H"
#include "geometry_pair_line_to_volume.H"
#include "geometry_pair_line_to_volume_gauss_point_projection.H"
#include "geometry_pair_line_to_volume_segmentation.H"
#include "geometry_pair_evaluation_data_global.H"
#include "geometry_pair_line_to_3D_evaluation_data.H"


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToVolumeFactory(
    const Teuchos::RCP<GeometryEvaluationDataGlobal> geometry_evaluation_data_ptr)
{
  // Get the strategy for line to volume interaction.
  INPAR::GEOMETRYPAIR::LineTo3DStrategy strategy =
      geometry_evaluation_data_ptr->LineTo3DEvaluationData()->GetStrategy();

  // Create the class depending on the strategy.
  switch (strategy)
  {
    case INPAR::GEOMETRYPAIR::LineTo3DStrategy::gauss_point_projection:
      return Teuchos::rcp(
          new GeometryPairLineToVolumeGaussPointProjection<scalar_type, line, volume>());
    case INPAR::GEOMETRYPAIR::LineTo3DStrategy::segmentation:
      return Teuchos::rcp(new GeometryPairLineToVolumeSegmentation<scalar_type, line, volume>());
    default:
    {
      dserror(
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
    double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex8>(
    const Teuchos::RCP<GeometryEvaluationDataGlobal>);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToVolumeFactory<
    double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex20>(
    const Teuchos::RCP<GeometryEvaluationDataGlobal>);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToVolumeFactory<
    double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex27>(
    const Teuchos::RCP<GeometryEvaluationDataGlobal>);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToVolumeFactory<
    double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet4>(
    const Teuchos::RCP<GeometryEvaluationDataGlobal>);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToVolumeFactory<
    double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet10>(
    const Teuchos::RCP<GeometryEvaluationDataGlobal>);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToVolumeFactory<
    double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_nurbs27>(
    const Teuchos::RCP<GeometryEvaluationDataGlobal>);
