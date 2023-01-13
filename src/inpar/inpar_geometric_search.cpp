/*-----------------------------------------------------------*/
/*! \file

\brief input parameter for geometric search strategy

\level 2

*/
/*-----------------------------------------------------------*/

#include "inpar_validparameters.H"
#include "inpar_geometric_search.H"

void INPAR::GEOMETRICSEARCH::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;

  Teuchos::ParameterList& boundingvolumestrategy =
      list->sublist("BOUNDINGVOLUME STRATEGY", false, "");

  DoubleParameter("BEAM_RADIUS_EXTENSION_FACTOR", 2.0,
      "Beams radius is multiplied with the factor and then the bounding volume only depending on "
      "the beam centerline is extended in all directions (+ and -) by that value.",
      &boundingvolumestrategy);

  DoubleParameter("SPHERE_RADIUS_EXTENSION_FACTOR", 2.0,
      "Bounding volume of the sphere is the sphere center extended by this factor times the sphere "
      "radius in all directions (+ and -).",
      &boundingvolumestrategy);
}
