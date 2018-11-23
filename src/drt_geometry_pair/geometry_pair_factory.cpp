/*!
\file geometry_pair_factory.cpp

\brief functions to create geometry pairs.

<pre>
\level 3
\maintainer Ivo Steinbrecher
            ivo.steinbrecher@unibw.de
            +49 89 6004-4403
</pre>
*/


#include "geometry_pair_factory.H"
#include "geometry_pair_line_to_volume.H"
#include "geometry_pair_line_to_volume_gauss_point_projection.H"
#include "geometry_pair_evaluation_data_global.H"
#include "geometry_pair_line_to_volume_evaluation_data.H"


/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
Teuchos::RCP<GEOMETRYPAIR::GeometryPair> GEOMETRYPAIR::GeometryPairLineToVolumeFactory(
    const Teuchos::RCP<GeometryEvaluationDataGlobal> geometry_evaluation_data_ptr)
{
  // Get the strategy for line to volume interaction.
  INPAR::GEOMETRYPAIR::LineToVolumeStrategy strategy =
      geometry_evaluation_data_ptr->LineToVolumeEvaluationData()->GetStrategy();

  // Create the class depending on the strategy.
  switch (strategy)
  {
    case INPAR::GEOMETRYPAIR::LineToVolumeStrategy::gauss_point_projection:
      return Teuchos::rcp(
          new GeometryPairLineToVolumeGaussPointProjection<scalar_type, n_nodes_element_1,
              n_nodal_values_element_1, n_nodes_element_2, n_nodal_values_element_2>());
    case INPAR::GEOMETRYPAIR::LineToVolumeStrategy::segmentation:
      return Teuchos::rcp(new GeometryPairLineToVolume<scalar_type, n_nodes_element_1,
          n_nodal_values_element_1, n_nodes_element_2, n_nodal_values_element_2>());
    default:
      return Teuchos::null;
  }
}


/**
 * Explicit template initialization of factory function.
 */
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair>
GEOMETRYPAIR::GeometryPairLineToVolumeFactory<double, 2, 2, 8, 1>(
    const Teuchos::RCP<GeometryEvaluationDataGlobal>);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair>
GEOMETRYPAIR::GeometryPairLineToVolumeFactory<double, 2, 2, 20, 1>(
    const Teuchos::RCP<GeometryEvaluationDataGlobal>);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair>
GEOMETRYPAIR::GeometryPairLineToVolumeFactory<double, 2, 2, 27, 1>(
    const Teuchos::RCP<GeometryEvaluationDataGlobal>);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair>
GEOMETRYPAIR::GeometryPairLineToVolumeFactory<double, 2, 2, 4, 1>(
    const Teuchos::RCP<GeometryEvaluationDataGlobal>);
template Teuchos::RCP<GEOMETRYPAIR::GeometryPair>
GEOMETRYPAIR::GeometryPairLineToVolumeFactory<double, 2, 2, 10, 1>(
    const Teuchos::RCP<GeometryEvaluationDataGlobal>);
