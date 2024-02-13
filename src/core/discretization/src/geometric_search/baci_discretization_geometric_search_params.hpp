/*-----------------------------------------------------------*/
/*! \file

\brief Data container class for geometric search parameters

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_DISCRETIZATION_GEOMETRIC_SEARCH_PARAMS_HPP
#define BACI_DISCRETIZATION_GEOMETRIC_SEARCH_PARAMS_HPP

#include "baci_config.hpp"

#include "baci_io_pstream.hpp"

#include <Teuchos_ParameterList.hpp>

BACI_NAMESPACE_OPEN

namespace CORE::GEOMETRICSEARCH
{
  //! Data container for geometric search parameters
  class GeometricSearchParams
  {
   public:
    GeometricSearchParams(const Teuchos::ParameterList& geometric_search_params,
        const Teuchos::ParameterList& io_params);

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

    /*!
     * \brief Returns a flag indicating if visualization output should be written or not
     */
    bool GetWriteVisualizationFlag() const { return write_visualization_; }

    /*! \brief verbosity level of the geometric search algorithm
     */
    IO::verbositylevel verbosity_;

   private:
    double beam_radius_extension_factor_;
    double sphere_radius_extension_factor_;
    bool write_visualization_;
  };
}  // namespace CORE::GEOMETRICSEARCH

BACI_NAMESPACE_CLOSE

#endif
