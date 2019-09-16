/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameter for geometry pairs.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef INPAR_GEOMETRYPAIR_H
#define INPAR_GEOMETRYPAIR_H


#include <Teuchos_RCP.hpp>


// Forward declaration.
namespace Teuchos
{
  class ParameterList;

}  // namespace Teuchos


namespace INPAR
{
  namespace GEOMETRYPAIR
  {
    /**
     * Method to be used for line to volume pairs.
     */
    enum class LineToVolumeStrategy
    {
      none,                    ///< Default value.
      gauss_point_projection,  ///< Project the Gauss points, segment lines that poke out of the
                               ///< volume.
      gauss_point_projection_cross_section,  ///< Project the Gauss points on the cross section
                                             ///< surface.
      segmentation                           ///< Segment the line with the volume.
    };

    /**
     * Set valid input parameters for geometry pairs.
     * @param list
     */
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace GEOMETRYPAIR

}  // namespace INPAR


#endif
