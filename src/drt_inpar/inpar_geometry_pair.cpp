/*!
\file inpar_geometry_pair.cpp

\brief Input parameter for geometry pairs.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "inpar_geometry_pair.H"
#include "drt_validparameters.H"


/**
 *
 */
void INPAR::GEOMETRYPAIR::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  // Parent section for general geometry pair options.
  Teuchos::ParameterList& geometry_pair = list->sublist("GEOMETRY PAIR", false, "");

  // Input parameters for line to volume pairs.
  {
    Teuchos::ParameterList& line_to_volume = geometry_pair.sublist("LINE TO VOLUME", false, "");

    // Segmentation strategy.
    Teuchos::setStringToIntegralParameter<LineToVolumeStrategy>("STRATEGY", "segmentation",
        "Type of employed segmentation strategy",
        Teuchos::tuple<std::string>(
            "segmentation", "gauss_point_projection", "gauss_point_projection_cylinder"),
        Teuchos::tuple<LineToVolumeStrategy>(LineToVolumeStrategy::segmentation,
            LineToVolumeStrategy::gauss_point_projection,
            LineToVolumeStrategy::gauss_point_projection_cylinder),
        &line_to_volume);

    // Number of search points for segmentation.
    DRT::INPUT::IntParameter(
        "SEARCH_POINTS", 6, "Number of search points for segmentation", &line_to_volume);
  }
}
