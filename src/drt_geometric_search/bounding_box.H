/*----------------------------------------------------------------------*/
/*! \file

\brief TODO

\level 3

*/
// End doxygen header.



#ifndef BOUNDING_BOX_H_
#define BOUNDING_BOX_H_

#include "../linalg/linalg_fixedsizematrix.H"


// Todo: rename this to BoundingVolume
// Todo: put this and the collision search in an own namespace, and
// add the namespace to the namespace header
struct BoundingBox
{
  BoundingBox()
  {
    min_corner_.PutScalar(1e100);
    max_corner_.PutScalar(-1e100);
  }

  BoundingBox(const LINALG::Matrix<3, 1, double>& min_corner,
      const LINALG::Matrix<3, 1, double>& max_corner)
      : min_corner_(min_corner), max_corner_(max_corner)
  {
  }

  inline void AddPoint(const LINALG::Matrix<3, 1, double>& point)
  {
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
    {
      min_corner_(i_dim) = std::min(min_corner_(i_dim), point(i_dim));
      max_corner_(i_dim) = std::max(max_corner_(i_dim), point(i_dim));
    }
  }

  inline void ExtendBoundaries(const double offset)
  {
    // Todo: debug check if offset is positive.
    LINALG::Matrix<3, 1, double> offset_vector;
    offset_vector.PutScalar(offset);
    min_corner_ -= offset_vector;
    max_corner_ += offset_vector;
  }

  //! Corners of the bounding box
  LINALG::Matrix<3, 1, double> min_corner_;
  LINALG::Matrix<3, 1, double> max_corner_;
};


inline bool BoxesIntersect(const BoundingBox& box_1, const BoundingBox& box_2)
{
  for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
  {
    const bool dimension_overlaps = (box_1.min_corner_(i_dim) < box_2.max_corner_(i_dim)) and
                                    (box_2.min_corner_(i_dim) < box_1.max_corner_(i_dim));
    if (!dimension_overlaps) return false;
  }
  return true;
}


#endif
