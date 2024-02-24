/*-----------------------------------------------------------*/
/*! \file

\brief input parameter for geometric search strategy

\level 2

*/
/*-----------------------------------------------------------*/

#include "baci_inpar_geometric_search.hpp"

#include "baci_utils_parameter_list.hpp"

BACI_NAMESPACE_OPEN

void INPAR::GEOMETRICSEARCH::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  Teuchos::ParameterList& boundingvolumestrategy =
      list->sublist("BOUNDINGVOLUME STRATEGY", false, "");

  CORE::UTILS::DoubleParameter("BEAM_RADIUS_EXTENSION_FACTOR", 2.0,
      "Beams radius is multiplied with the factor and then the bounding volume only depending on "
      "the beam centerline is extended in all directions (+ and -) by that value.",
      &boundingvolumestrategy);

  CORE::UTILS::DoubleParameter("SPHERE_RADIUS_EXTENSION_FACTOR", 2.0,
      "Bounding volume of the sphere is the sphere center extended by this factor times the sphere "
      "radius in all directions (+ and -).",
      &boundingvolumestrategy);

  CORE::UTILS::BoolParameter("WRITE_GEOMETRIC_SEARCH_VISUALIZATION", "no",
      "If visualization output for the geometric search should be written",
      &boundingvolumestrategy);
}

BACI_NAMESPACE_CLOSE
