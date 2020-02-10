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

  // Number of integration along the circumference in cross section coupling.
  DRT::INPUT::IntParameter("INTEGRATION_POINTS_CIRCUMFERENCE", 6,
      "Number of Integration points along the circumferencial direction of the beam. This is "
      "parameter is only used in beam to cylinder meshtying. No gauss integration is "
      "used along the circumferencial direction, equally spaced integration points are used.",
      &list);
}

/**
 *
 */
DRT::UTILS::GaussRule1D INPAR::GEOMETRYPAIR::IntToGaussRule1D(const int n_gauss_points)
{
  switch (n_gauss_points)
  {
    case 1:
      return DRT::UTILS::GaussRule1D::intrule_line_1point;
    case 2:
      return DRT::UTILS::GaussRule1D::intrule_line_2point;
    case 3:
      return DRT::UTILS::GaussRule1D::intrule_line_3point;
    case 4:
      return DRT::UTILS::GaussRule1D::intrule_line_4point;
    case 5:
      return DRT::UTILS::GaussRule1D::intrule_line_5point;
    case 6:
      return DRT::UTILS::GaussRule1D::intrule_line_6point;
    case 7:
      return DRT::UTILS::GaussRule1D::intrule_line_7point;
    case 8:
      return DRT::UTILS::GaussRule1D::intrule_line_8point;
    case 9:
      return DRT::UTILS::GaussRule1D::intrule_line_9point;
    case 10:
      return DRT::UTILS::GaussRule1D::intrule_line_10point;
    case 20:
      return DRT::UTILS::GaussRule1D::intrule_line_20point;
    case 32:
      return DRT::UTILS::GaussRule1D::intrule_line_32point;
    case 50:
      return DRT::UTILS::GaussRule1D::intrule_line_50point;
    default:
    {
      dserror("No Gauss rule defined for %d points", n_gauss_points);
      return DRT::UTILS::GaussRule1D::intrule1D_undefined;
    }
  }
};
