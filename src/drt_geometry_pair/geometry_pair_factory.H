/*----------------------------------------------------------------------*/
/*! \file

\brief functions to create geometry pairs.

\level 1

\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_FACTORY_H
#define GEOMETRY_PAIR_FACTORY_H


#include <Teuchos_RCP.hpp>


// Forward declarations.
namespace GEOMETRYPAIR
{
  class GeometryPair;
  class LineTo3DEvaluationData;
}  // namespace GEOMETRYPAIR


namespace GEOMETRYPAIR
{
  /**
   * \brief Create the correct geometry pair for line to volume coupling.
   * @return RCP to created geometry pair.
   */
  template <typename scalar_type, typename line, typename volume>
  Teuchos::RCP<GeometryPair> GeometryPairLineToVolumeFactory(
      const Teuchos::RCP<LineTo3DEvaluationData>& geometry_evaluation_data_ptr);
}  // namespace GEOMETRYPAIR

#endif
