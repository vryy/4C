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

#include "../linalg/linalg_fixedsizematrix.H"

namespace GEOMETRICSEARCH
{
  struct BoundingVolume
  {
    /*! \brief Constructor initializing the bounding volume corners with numerical limit values.
     */
    BoundingVolume()
    {
      min_corner_.PutScalar(std::numeric_limits<double>::max());
      max_corner_.PutScalar(std::numeric_limits<double>::lowest());
    }

    /*! \brief Adds a point to the bounding volume.
     *
     * @param point Point to add to the bounding volume
     */
    inline void AddPoint(const LINALG::Matrix<3, 1, double>& point)
    {
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      {
        min_corner_(i_dim) = std::min(min_corner_(i_dim), point(i_dim));
        max_corner_(i_dim) = std::max(max_corner_(i_dim), point(i_dim));
      }
    }

    /*! \brief Extends the bounding volume based on a scalar value.
     *
     * @param offset Safety value by which to expand the bounding volume
     */
    inline void ExtendBoundaries(const double offset)
    {
      dsassert(!std::signbit(offset), "Offset parameter needs to be positiv!");

      LINALG::Matrix<3, 1, double> offset_vector;
      offset_vector.PutScalar(offset);

      min_corner_ -= offset_vector;
      max_corner_ += offset_vector;
    }

    //! Corners of the bounding box
    LINALG::Matrix<3, 1, double> min_corner_;
    LINALG::Matrix<3, 1, double> max_corner_;
  };

  /*! \brief Checks if two bounding volumes intersect with each other.
   *
   * @param box_1 Fist bounding volume to be considered for intersection
   * @param box_2 Second bounding volume to be considered for intersection
   * @return Boolean is true if bounding volumes intersect, false otherwise
   */
  inline bool BoxesIntersect(const BoundingVolume& box_1, const BoundingVolume& box_2)
  {
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
    {
      const bool dimension_overlaps = (box_1.min_corner_(i_dim) < box_2.max_corner_(i_dim)) and
                                      (box_2.min_corner_(i_dim) < box_1.max_corner_(i_dim));
      if (!dimension_overlaps) return false;
    }
    return true;
  }
}  // namespace GEOMETRICSEARCH

#endif
