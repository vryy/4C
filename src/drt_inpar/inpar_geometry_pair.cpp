/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameter for geometry pairs.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "inpar_geometry_pair.H"
#include "drt_validparameters.H"


/**
 *
 */
void INPAR::GEOMETRYPAIR::SetValidParametersLineTo3D(Teuchos::ParameterList& list)
{
  // Add the input parameters for line to 3D coupling.

  // Segmentation strategy.
  Teuchos::setStringToIntegralParameter<LineTo3DStrategy>("GEOMETRY_PAIR_STRATEGY", "segmentation",
      "Type of employed segmentation strategy",
      Teuchos::tuple<std::string>(
          "segmentation", "gauss_point_projection", "gauss_point_projection_cross_section"),
      Teuchos::tuple<LineTo3DStrategy>(LineTo3DStrategy::segmentation,
          LineTo3DStrategy::gauss_point_projection,
          LineTo3DStrategy::gauss_point_projection_cross_section),
      &list);

  // Number of search points for segmentation.
  DRT::INPUT::IntParameter(
      "GEOMETRY_PAIR_SEARCH_POINTS", 6, "Number of search points for segmentation", &list);
}
