/*!
\file geometry_pair_factory.H

\brief functions to create geometry pairs.

<pre>
\level 3
\maintainer Ivo Steinbrecher
            ivo.steinbrecher@unibw.de
            +49 89 6004-4403
</pre>
*/


#ifndef GEOMETRY_PAIR_FACTORY_H
#define GEOMETRY_PAIR_FACTORY_H


#include <Teuchos_RCP.hpp>


// Forward declarations.
namespace GEOMETRYPAIR
{
  class GeometryPair;
  class GeometryEvaluationDataGlobal;
}  // namespace GEOMETRYPAIR


namespace GEOMETRYPAIR
{
  /**
   * \brief Create the correct geometry pair for line to volume coupling.
   * @return RCP to created geometry pair.
   */
  template <typename scalar_type, unsigned int n_nodes_element_1,
      unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
      unsigned int n_nodal_values_element_2>
  Teuchos::RCP<GeometryPair> GeometryPairLineToVolumeFactory(
      const Teuchos::RCP<GeometryEvaluationDataGlobal> geometry_evaluation_data_ptr);
}  // namespace GEOMETRYPAIR

#endif
