/*-----------------------------------------------------------*/
/*! \file

\brief Data container class for geometric search parameters

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FEM_GEOMETRIC_SEARCH_PARAMS_HPP
#define FOUR_C_FEM_GEOMETRIC_SEARCH_PARAMS_HPP

#include "4C_config.hpp"

#include "4C_io_pstream.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::GeometricSearch
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
    double get_beam_bounding_volume_scaling() const { return beam_radius_extension_factor_; }

    /*! \brief Returns sphere bounding volume extension factor
     *
     * Center of the sphere is extended by the factor times the sphere radius.
     */
    double get_sphere_bounding_volume_scaling() const { return sphere_radius_extension_factor_; }

    /*!
     * \brief Returns a flag indicating if visualization output should be written or not
     */
    bool get_write_visualization_flag() const { return write_visualization_; }

    /*! \brief verbosity level of the geometric search algorithm
     */
    Core::IO::Verbositylevel verbosity_;

   private:
    double beam_radius_extension_factor_;
    double sphere_radius_extension_factor_;
    bool write_visualization_;
  };
}  // namespace Core::GeometricSearch

FOUR_C_NAMESPACE_CLOSE

#endif
