/*!
\file inpar_geometry_pair.H

\brief Input parameter for geometry pairs.

<pre>
\level 3
\maintainer Ivo Steinbrecher
            ivo.steinbrecher@unibw.de
            +49 89 6004-4403
</pre>
*/


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
      segmentation             ///< Segment the line with the volume.
    };

    /**
     * Set valid input parameters for geometry pairs.
     * @param list
     */
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace GEOMETRYPAIR

}  // namespace INPAR


#endif
