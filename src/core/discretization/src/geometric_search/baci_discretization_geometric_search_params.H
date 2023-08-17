/*-----------------------------------------------------------*/
/*! \file

\brief Data container class for geometric search parameters

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_DISCRETIZATION_GEOMETRIC_SEARCH_PARAMS_H
#define BACI_DISCRETIZATION_GEOMETRIC_SEARCH_PARAMS_H

#include "baci_io_pstream.H"

namespace CORE::GEOMETRICSEARCH
{
  //! Data container for geometric search parameters
  class GeometricSearchParams
  {
   public:
    GeometricSearchParams();

    /*! \brief Returns beam bounding volume scaling factor
     *
     *  Beams radius is multiplied with the factor and then the bounding box only depending on
     *  the beam centerline is extended in all directions (+ and -) by that value.
     */
    double GetBeamBoundingVolumeScaling() const { return beam_radius_extension_factor_; }

    /*! \brief Returns sphere bounding volume extension factor
     *
     * Center of the sphere is extended by the factor times the sphere radius.
     */
    double GetSphereBoundingVolumeScaling() const { return sphere_radius_extension_factor_; }

    /*! \brief verbosity level of the geometric search algorithm
     */
    IO::verbositylevel verbosity_;

   private:
    double beam_radius_extension_factor_;
    double sphere_radius_extension_factor_;
  };
}  // namespace CORE::GEOMETRICSEARCH

#endif
