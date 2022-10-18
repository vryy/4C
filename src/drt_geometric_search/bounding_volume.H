/*-----------------------------------------------------------*/
/*! \file

\brief Contains a baci-specific implementation of a bounding
       volume defined by it's lower left and upper right
       corner.

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BOUNDING_VOLUME_H
#define BOUNDING_VOLUME_H

#ifdef HAVE_ARBORX
#include <ArborX.hpp>
#endif

#include "../linalg/linalg_fixedsizematrix.H"

namespace GEOMETRICSEARCH
{
  struct BoundingVolume
  {
#ifndef HAVE_ARBORX
    /*! \brief This class can not be used without ArborX, add empy methods and a controlled error.
     */
    BoundingVolume()
    {
      dserror(
          "The struct 'GEOMETRICSEARCH::BoundingVolume' can only be used with ArborX."
          "To use it, enable ArborX during the configure process.");
    }
    inline void ExtendBoundaries(const double offset) {}
    inline void AddPoint(const LINALG::Matrix<3, 1, double> &point) {}
#else
    /*! \brief Constructor initializing the bounding volume corners with numerical limit values.
     */
    BoundingVolume() : bounding_volume_{} {}

    /*! \brief Adds a point to the bounding volume.
     *
     * @param point Point to add to the bounding volume
     */
    inline void AddPoint(const LINALG::Matrix<3, 1, double>& point)
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
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      {
        bounding_volume_.minCorner()[i_dim] -= static_cast<float>(offset);
        bounding_volume_.maxCorner()[i_dim] += static_cast<float>(offset);
      }
    }

    //! Data container to store geometry information of bounding volume
    ArborX::Box bounding_volume_;
#endif
  };
}  // namespace GEOMETRICSEARCH

#endif
