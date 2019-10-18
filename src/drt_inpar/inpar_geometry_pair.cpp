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

  // Number of integration points on the line.
  DRT::INPUT::IntParameter(
      "GAUSS_POINTS", 6, "Number of Gauss Points for the integral evaluations", &list);

  // Number of integration along the circumfence in cross section coupling.
  DRT::INPUT::IntParameter("INTEGRATION_POINTS_CIRCUMFENCE", 6,
      "Number of Integration points along the circumfencial direction of the beam. This is "
      "parameter is only used in beam to cylinder meshtying. No gauss integration is "
      "used along the circumfencial direction, equally spaced integration points are used.",
      &list);
}
