/*-----------------------------------------------------------*/
/*! \file

\brief Contains a baci-specific implementation of a bounding
       volume defined as k-DOP (discrete oriented polytope).

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_DISCRETIZATION_GEOMETRIC_SEARCH_BOUNDING_VOLUME_HPP
#define FOUR_C_DISCRETIZATION_GEOMETRIC_SEARCH_BOUNDING_VOLUME_HPP

#include "baci_config.hpp"

#include "baci_linalg_fixedsizematrix.hpp"

#ifdef BACI_WITH_ARBORX
#include <ArborX.hpp>
#include <ArborX_KDOP.hpp>
#endif

#ifdef BACI_WITH_ARBORX
namespace ArborX::Details
{
  /*! \brief  This struct helps to get access to the protected static member directions from ArborX.
   */
  template <int n_directions>
  struct GetKDOPDirections : private KDOP_Directions<2 * n_directions>
  {
    static KOKKOS_FUNCTION Kokkos::Array<Direction, n_directions> const &directions()
    {
      return KDOP_Directions<2 * n_directions>::directions();
    }
  };
}  // namespace ArborX::Details

#endif

FOUR_C_NAMESPACE_OPEN

namespace CORE::GEOMETRICSEARCH
{
#ifdef BACI_WITH_ARBORX
  //! kdop_directions is the number of directions defining the possible faces of the k-DOP. In this
  //! case 13 means that we have the (3) Cartesian basis vectors, all vectors pointing to the edges
  //! of a unit cube (6) and all vectors pointing to corners of the unit cube (4). For more details
  //! have a look in ArborX_KDOP.hpp.
  constexpr int kdop_directions = 13;

  //! kdop_k is the number of bounding slabs for the k-DOP, i.e., the resulting polytope can have a
  //! maximum of kdop_k faces
  constexpr int kdop_k = 2 * kdop_directions;
#endif

  /*!
   * \brief Class representing a bounding volume of a set of points.
   *
   * Internal this class uses k-DOPs implemented in ArborX. A k-DOP (discrete oriented polytope) is
   * a convex polytope containing the object.
   */
  struct BoundingVolume
  {
#ifndef BACI_WITH_ARBORX
    /*! \brief This class can not be used without ArborX, add empy methods and a controlled error.
     */
    BoundingVolume()
    {
      FOUR_C_THROW(
          "The struct 'CORE::GEOMETRICSEARCH::BoundingVolume' can only be used with ArborX."
          "To use it, enable ArborX during the configure process.");
    }
    inline void ExtendBoundaries(const double offset) {}
    inline void AddPoint(const CORE::LINALG::Matrix<3, 1, double> &point) {}
#else
    /*! \brief Constructor initializing the bounding volume corners with numerical limit values.
     */
    BoundingVolume() : bounding_volume_{} {}

    /*! \brief Adds a point to the bounding volume.
     *
     * @param point Point to add to the bounding volume
     */
    inline void AddPoint(const CORE::LINALG::Matrix<3, 1, double>& point)
    {
      bounding_volume_ += ArborX::Point{
          static_cast<float>(point(0)), static_cast<float>(point(1)), static_cast<float>(point(2))};
    }

    /*! \brief Extends the bounding volume based on a scalar value.
     *
     * @param offset Value by which to expand the bounding volume
     */
    inline void ExtendBoundaries(const double offset)
    {
      // Loop over directions.
      for (int i_dir = 0; i_dir < kdop_directions; i_dir++)
      {
        const auto& direction =
            ArborX::Details::GetKDOPDirections<kdop_directions>::directions()[i_dir];
        const auto direction_scaling_factor = std::sqrt(direction._data[0] * direction._data[0] +
                                                        direction._data[1] * direction._data[1] +
                                                        direction._data[2] * direction._data[2]);
        const auto scaled_offset = offset * direction_scaling_factor;
        bounding_volume_._min_values[i_dir] -= static_cast<float>(scaled_offset);
        bounding_volume_._max_values[i_dir] += static_cast<float>(scaled_offset);
      }
    }

    //! Use an ArborX geometry as internal storage.
    ArborX::Experimental::KDOP<kdop_k> bounding_volume_;
#endif
  };
}  // namespace CORE::GEOMETRICSEARCH

FOUR_C_NAMESPACE_CLOSE

#endif
