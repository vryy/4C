/*-----------------------------------------------------------*/
/*! \file

\brief Contains all possible intersections of the k-dop planes.

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_BACI_DISCRETIZATION_GEOMETRIC_SEARCH_KDOP_INTERSECTIONS_H
#define BACI_BACI_DISCRETIZATION_GEOMETRIC_SEARCH_KDOP_INTERSECTIONS_H

#include <array>

// The k-dop visualization is based on polygons defined which are defined by intersections
// of the k-dop planes. This pre-computed array (see the Mathematica scrip in the scrips/
// sub-directory) contains data that describes all possible intersection points between all
// k-dop planes.
//
// The first array index is the k-dop direction.
// The second array index specifies if the min (0) or max (1) value along a direction is to be
// considered. The element addressed by these two indices is a vector of pairs. Each element of
// the vector is a possible intersection plane to the described by the first two indices. This
// is stored in a tuple where the first index is the direction id and the second index is the
// min/max flag.

namespace CORE::GEOMETRICSEARCH
{
  struct KdopIntersection
  {
    /*! \brief Constructor initializing the pre-computed intersection pairs of the kdop
     */
    KdopIntersection()
    {
      pairs_[0][0].emplace_back(std::make_pair(11, false));
      pairs_[0][0].emplace_back(std::make_pair(7, false));
      pairs_[0][0].emplace_back(std::make_pair(12, false));
      pairs_[0][0].emplace_back(std::make_pair(6, false));
      pairs_[0][0].emplace_back(std::make_pair(10, false));
      pairs_[0][0].emplace_back(std::make_pair(4, false));
      pairs_[0][0].emplace_back(std::make_pair(9, false));
      pairs_[0][0].emplace_back(std::make_pair(3, false));
      pairs_[0][1].emplace_back(std::make_pair(12, true));
      pairs_[0][1].emplace_back(std::make_pair(7, true));
      pairs_[0][1].emplace_back(std::make_pair(11, true));
      pairs_[0][1].emplace_back(std::make_pair(3, true));
      pairs_[0][1].emplace_back(std::make_pair(9, true));
      pairs_[0][1].emplace_back(std::make_pair(4, true));
      pairs_[0][1].emplace_back(std::make_pair(10, true));
      pairs_[0][1].emplace_back(std::make_pair(6, true));
      pairs_[1][0].emplace_back(std::make_pair(9, false));
      pairs_[1][0].emplace_back(std::make_pair(5, false));
      pairs_[1][0].emplace_back(std::make_pair(12, true));
      pairs_[1][0].emplace_back(std::make_pair(6, true));
      pairs_[1][0].emplace_back(std::make_pair(10, true));
      pairs_[1][0].emplace_back(std::make_pair(8, false));
      pairs_[1][0].emplace_back(std::make_pair(11, false));
      pairs_[1][0].emplace_back(std::make_pair(3, false));
      pairs_[1][1].emplace_back(std::make_pair(12, false));
      pairs_[1][1].emplace_back(std::make_pair(5, true));
      pairs_[1][1].emplace_back(std::make_pair(9, true));
      pairs_[1][1].emplace_back(std::make_pair(3, true));
      pairs_[1][1].emplace_back(std::make_pair(11, true));
      pairs_[1][1].emplace_back(std::make_pair(8, true));
      pairs_[1][1].emplace_back(std::make_pair(10, false));
      pairs_[1][1].emplace_back(std::make_pair(6, false));
      pairs_[2][0].emplace_back(std::make_pair(10, false));
      pairs_[2][0].emplace_back(std::make_pair(8, true));
      pairs_[2][0].emplace_back(std::make_pair(11, true));
      pairs_[2][0].emplace_back(std::make_pair(7, true));
      pairs_[2][0].emplace_back(std::make_pair(12, true));
      pairs_[2][0].emplace_back(std::make_pair(5, false));
      pairs_[2][0].emplace_back(std::make_pair(9, false));
      pairs_[2][0].emplace_back(std::make_pair(4, false));
      pairs_[2][1].emplace_back(std::make_pair(11, false));
      pairs_[2][1].emplace_back(std::make_pair(8, false));
      pairs_[2][1].emplace_back(std::make_pair(10, true));
      pairs_[2][1].emplace_back(std::make_pair(4, true));
      pairs_[2][1].emplace_back(std::make_pair(9, true));
      pairs_[2][1].emplace_back(std::make_pair(5, true));
      pairs_[2][1].emplace_back(std::make_pair(12, false));
      pairs_[2][1].emplace_back(std::make_pair(7, false));
      pairs_[3][0].emplace_back(std::make_pair(10, false));
      pairs_[3][0].emplace_back(std::make_pair(4, false));
      pairs_[3][0].emplace_back(std::make_pair(9, false));
      pairs_[3][0].emplace_back(std::make_pair(5, false));
      pairs_[3][0].emplace_back(std::make_pair(12, true));
      pairs_[3][0].emplace_back(std::make_pair(1, false));
      pairs_[3][0].emplace_back(std::make_pair(10, true));
      pairs_[3][0].emplace_back(std::make_pair(8, false));
      pairs_[3][0].emplace_back(std::make_pair(11, false));
      pairs_[3][0].emplace_back(std::make_pair(7, false));
      pairs_[3][0].emplace_back(std::make_pair(12, false));
      pairs_[3][0].emplace_back(std::make_pair(0, false));
      pairs_[3][1].emplace_back(std::make_pair(12, false));
      pairs_[3][1].emplace_back(std::make_pair(5, true));
      pairs_[3][1].emplace_back(std::make_pair(9, true));
      pairs_[3][1].emplace_back(std::make_pair(4, true));
      pairs_[3][1].emplace_back(std::make_pair(10, true));
      pairs_[3][1].emplace_back(std::make_pair(0, true));
      pairs_[3][1].emplace_back(std::make_pair(12, true));
      pairs_[3][1].emplace_back(std::make_pair(7, true));
      pairs_[3][1].emplace_back(std::make_pair(11, true));
      pairs_[3][1].emplace_back(std::make_pair(8, true));
      pairs_[3][1].emplace_back(std::make_pair(10, false));
      pairs_[3][1].emplace_back(std::make_pair(1, true));
      pairs_[4][0].emplace_back(std::make_pair(12, false));
      pairs_[4][0].emplace_back(std::make_pair(6, false));
      pairs_[4][0].emplace_back(std::make_pair(10, false));
      pairs_[4][0].emplace_back(std::make_pair(8, true));
      pairs_[4][0].emplace_back(std::make_pair(11, true));
      pairs_[4][0].emplace_back(std::make_pair(2, false));
      pairs_[4][0].emplace_back(std::make_pair(12, true));
      pairs_[4][0].emplace_back(std::make_pair(5, false));
      pairs_[4][0].emplace_back(std::make_pair(9, false));
      pairs_[4][0].emplace_back(std::make_pair(3, false));
      pairs_[4][0].emplace_back(std::make_pair(11, false));
      pairs_[4][0].emplace_back(std::make_pair(0, false));
      pairs_[4][1].emplace_back(std::make_pair(11, false));
      pairs_[4][1].emplace_back(std::make_pair(8, false));
      pairs_[4][1].emplace_back(std::make_pair(10, true));
      pairs_[4][1].emplace_back(std::make_pair(6, true));
      pairs_[4][1].emplace_back(std::make_pair(12, true));
      pairs_[4][1].emplace_back(std::make_pair(0, true));
      pairs_[4][1].emplace_back(std::make_pair(11, true));
      pairs_[4][1].emplace_back(std::make_pair(3, true));
      pairs_[4][1].emplace_back(std::make_pair(9, true));
      pairs_[4][1].emplace_back(std::make_pair(5, true));
      pairs_[4][1].emplace_back(std::make_pair(12, false));
      pairs_[4][1].emplace_back(std::make_pair(2, true));
      pairs_[5][0].emplace_back(std::make_pair(4, false));
      pairs_[5][0].emplace_back(std::make_pair(10, false));
      pairs_[5][0].emplace_back(std::make_pair(2, false));
      pairs_[5][0].emplace_back(std::make_pair(11, true));
      pairs_[5][0].emplace_back(std::make_pair(7, true));
      pairs_[5][0].emplace_back(std::make_pair(12, true));
      pairs_[5][0].emplace_back(std::make_pair(6, true));
      pairs_[5][0].emplace_back(std::make_pair(10, true));
      pairs_[5][0].emplace_back(std::make_pair(1, false));
      pairs_[5][0].emplace_back(std::make_pair(11, false));
      pairs_[5][0].emplace_back(std::make_pair(3, false));
      pairs_[5][0].emplace_back(std::make_pair(9, false));
      pairs_[5][1].emplace_back(std::make_pair(7, false));
      pairs_[5][1].emplace_back(std::make_pair(11, false));
      pairs_[5][1].emplace_back(std::make_pair(2, true));
      pairs_[5][1].emplace_back(std::make_pair(10, true));
      pairs_[5][1].emplace_back(std::make_pair(4, true));
      pairs_[5][1].emplace_back(std::make_pair(9, true));
      pairs_[5][1].emplace_back(std::make_pair(3, true));
      pairs_[5][1].emplace_back(std::make_pair(11, true));
      pairs_[5][1].emplace_back(std::make_pair(1, true));
      pairs_[5][1].emplace_back(std::make_pair(10, false));
      pairs_[5][1].emplace_back(std::make_pair(6, false));
      pairs_[5][1].emplace_back(std::make_pair(12, false));
      pairs_[6][0].emplace_back(std::make_pair(11, false));
      pairs_[6][0].emplace_back(std::make_pair(7, false));
      pairs_[6][0].emplace_back(std::make_pair(12, false));
      pairs_[6][0].emplace_back(std::make_pair(5, true));
      pairs_[6][0].emplace_back(std::make_pair(9, true));
      pairs_[6][0].emplace_back(std::make_pair(1, true));
      pairs_[6][0].emplace_back(std::make_pair(11, true));
      pairs_[6][0].emplace_back(std::make_pair(8, true));
      pairs_[6][0].emplace_back(std::make_pair(10, false));
      pairs_[6][0].emplace_back(std::make_pair(4, false));
      pairs_[6][0].emplace_back(std::make_pair(9, false));
      pairs_[6][0].emplace_back(std::make_pair(0, false));
      pairs_[6][1].emplace_back(std::make_pair(9, false));
      pairs_[6][1].emplace_back(std::make_pair(5, false));
      pairs_[6][1].emplace_back(std::make_pair(12, true));
      pairs_[6][1].emplace_back(std::make_pair(7, true));
      pairs_[6][1].emplace_back(std::make_pair(11, true));
      pairs_[6][1].emplace_back(std::make_pair(0, true));
      pairs_[6][1].emplace_back(std::make_pair(9, true));
      pairs_[6][1].emplace_back(std::make_pair(4, true));
      pairs_[6][1].emplace_back(std::make_pair(10, true));
      pairs_[6][1].emplace_back(std::make_pair(8, false));
      pairs_[6][1].emplace_back(std::make_pair(11, false));
      pairs_[6][1].emplace_back(std::make_pair(1, false));
      pairs_[7][0].emplace_back(std::make_pair(9, false));
      pairs_[7][0].emplace_back(std::make_pair(3, false));
      pairs_[7][0].emplace_back(std::make_pair(11, false));
      pairs_[7][0].emplace_back(std::make_pair(8, false));
      pairs_[7][0].emplace_back(std::make_pair(10, true));
      pairs_[7][0].emplace_back(std::make_pair(2, true));
      pairs_[7][0].emplace_back(std::make_pair(9, true));
      pairs_[7][0].emplace_back(std::make_pair(5, true));
      pairs_[7][0].emplace_back(std::make_pair(12, false));
      pairs_[7][0].emplace_back(std::make_pair(6, false));
      pairs_[7][0].emplace_back(std::make_pair(10, false));
      pairs_[7][0].emplace_back(std::make_pair(0, false));
      pairs_[7][1].emplace_back(std::make_pair(10, false));
      pairs_[7][1].emplace_back(std::make_pair(8, true));
      pairs_[7][1].emplace_back(std::make_pair(11, true));
      pairs_[7][1].emplace_back(std::make_pair(3, true));
      pairs_[7][1].emplace_back(std::make_pair(9, true));
      pairs_[7][1].emplace_back(std::make_pair(0, true));
      pairs_[7][1].emplace_back(std::make_pair(10, true));
      pairs_[7][1].emplace_back(std::make_pair(6, true));
      pairs_[7][1].emplace_back(std::make_pair(12, true));
      pairs_[7][1].emplace_back(std::make_pair(5, false));
      pairs_[7][1].emplace_back(std::make_pair(9, false));
      pairs_[7][1].emplace_back(std::make_pair(2, false));
      pairs_[8][0].emplace_back(std::make_pair(3, false));
      pairs_[8][0].emplace_back(std::make_pair(9, false));
      pairs_[8][0].emplace_back(std::make_pair(1, false));
      pairs_[8][0].emplace_back(std::make_pair(12, true));
      pairs_[8][0].emplace_back(std::make_pair(6, true));
      pairs_[8][0].emplace_back(std::make_pair(10, true));
      pairs_[8][0].emplace_back(std::make_pair(4, true));
      pairs_[8][0].emplace_back(std::make_pair(9, true));
      pairs_[8][0].emplace_back(std::make_pair(2, true));
      pairs_[8][0].emplace_back(std::make_pair(12, false));
      pairs_[8][0].emplace_back(std::make_pair(7, false));
      pairs_[8][0].emplace_back(std::make_pair(11, false));
      pairs_[8][1].emplace_back(std::make_pair(6, false));
      pairs_[8][1].emplace_back(std::make_pair(12, false));
      pairs_[8][1].emplace_back(std::make_pair(1, true));
      pairs_[8][1].emplace_back(std::make_pair(9, true));
      pairs_[8][1].emplace_back(std::make_pair(3, true));
      pairs_[8][1].emplace_back(std::make_pair(11, true));
      pairs_[8][1].emplace_back(std::make_pair(7, true));
      pairs_[8][1].emplace_back(std::make_pair(12, true));
      pairs_[8][1].emplace_back(std::make_pair(2, false));
      pairs_[8][1].emplace_back(std::make_pair(9, false));
      pairs_[8][1].emplace_back(std::make_pair(4, false));
      pairs_[8][1].emplace_back(std::make_pair(10, false));
      pairs_[9][0].emplace_back(std::make_pair(6, false));
      pairs_[9][0].emplace_back(std::make_pair(4, false));
      pairs_[9][0].emplace_back(std::make_pair(8, true));
      pairs_[9][0].emplace_back(std::make_pair(2, false));
      pairs_[9][0].emplace_back(std::make_pair(7, true));
      pairs_[9][0].emplace_back(std::make_pair(5, false));
      pairs_[9][0].emplace_back(std::make_pair(6, true));
      pairs_[9][0].emplace_back(std::make_pair(1, false));
      pairs_[9][0].emplace_back(std::make_pair(8, false));
      pairs_[9][0].emplace_back(std::make_pair(3, false));
      pairs_[9][0].emplace_back(std::make_pair(7, false));
      pairs_[9][0].emplace_back(std::make_pair(0, false));
      pairs_[9][1].emplace_back(std::make_pair(7, false));
      pairs_[9][1].emplace_back(std::make_pair(2, true));
      pairs_[9][1].emplace_back(std::make_pair(8, false));
      pairs_[9][1].emplace_back(std::make_pair(4, true));
      pairs_[9][1].emplace_back(std::make_pair(6, true));
      pairs_[9][1].emplace_back(std::make_pair(0, true));
      pairs_[9][1].emplace_back(std::make_pair(7, true));
      pairs_[9][1].emplace_back(std::make_pair(3, true));
      pairs_[9][1].emplace_back(std::make_pair(8, true));
      pairs_[9][1].emplace_back(std::make_pair(1, true));
      pairs_[9][1].emplace_back(std::make_pair(6, false));
      pairs_[9][1].emplace_back(std::make_pair(5, true));
      pairs_[10][0].emplace_back(std::make_pair(7, false));
      pairs_[10][0].emplace_back(std::make_pair(6, false));
      pairs_[10][0].emplace_back(std::make_pair(5, true));
      pairs_[10][0].emplace_back(std::make_pair(1, true));
      pairs_[10][0].emplace_back(std::make_pair(3, true));
      pairs_[10][0].emplace_back(std::make_pair(8, true));
      pairs_[10][0].emplace_back(std::make_pair(7, true));
      pairs_[10][0].emplace_back(std::make_pair(2, false));
      pairs_[10][0].emplace_back(std::make_pair(5, false));
      pairs_[10][0].emplace_back(std::make_pair(4, false));
      pairs_[10][0].emplace_back(std::make_pair(3, false));
      pairs_[10][0].emplace_back(std::make_pair(0, false));
      pairs_[10][1].emplace_back(std::make_pair(3, false));
      pairs_[10][1].emplace_back(std::make_pair(1, false));
      pairs_[10][1].emplace_back(std::make_pair(5, false));
      pairs_[10][1].emplace_back(std::make_pair(6, true));
      pairs_[10][1].emplace_back(std::make_pair(7, true));
      pairs_[10][1].emplace_back(std::make_pair(0, true));
      pairs_[10][1].emplace_back(std::make_pair(3, true));
      pairs_[10][1].emplace_back(std::make_pair(4, true));
      pairs_[10][1].emplace_back(std::make_pair(5, true));
      pairs_[10][1].emplace_back(std::make_pair(2, true));
      pairs_[10][1].emplace_back(std::make_pair(7, false));
      pairs_[10][1].emplace_back(std::make_pair(8, false));
      pairs_[11][0].emplace_back(std::make_pair(4, false));
      pairs_[11][0].emplace_back(std::make_pair(3, false));
      pairs_[11][0].emplace_back(std::make_pair(5, false));
      pairs_[11][0].emplace_back(std::make_pair(1, false));
      pairs_[11][0].emplace_back(std::make_pair(6, true));
      pairs_[11][0].emplace_back(std::make_pair(8, false));
      pairs_[11][0].emplace_back(std::make_pair(4, true));
      pairs_[11][0].emplace_back(std::make_pair(2, true));
      pairs_[11][0].emplace_back(std::make_pair(5, true));
      pairs_[11][0].emplace_back(std::make_pair(7, false));
      pairs_[11][0].emplace_back(std::make_pair(6, false));
      pairs_[11][0].emplace_back(std::make_pair(0, false));
      pairs_[11][1].emplace_back(std::make_pair(6, false));
      pairs_[11][1].emplace_back(std::make_pair(1, true));
      pairs_[11][1].emplace_back(std::make_pair(5, true));
      pairs_[11][1].emplace_back(std::make_pair(3, true));
      pairs_[11][1].emplace_back(std::make_pair(4, true));
      pairs_[11][1].emplace_back(std::make_pair(0, true));
      pairs_[11][1].emplace_back(std::make_pair(6, true));
      pairs_[11][1].emplace_back(std::make_pair(7, true));
      pairs_[11][1].emplace_back(std::make_pair(5, false));
      pairs_[11][1].emplace_back(std::make_pair(2, false));
      pairs_[11][1].emplace_back(std::make_pair(4, false));
      pairs_[11][1].emplace_back(std::make_pair(8, true));
      pairs_[12][0].emplace_back(std::make_pair(3, false));
      pairs_[12][0].emplace_back(std::make_pair(7, false));
      pairs_[12][0].emplace_back(std::make_pair(8, false));
      pairs_[12][0].emplace_back(std::make_pair(2, true));
      pairs_[12][0].emplace_back(std::make_pair(4, true));
      pairs_[12][0].emplace_back(std::make_pair(5, true));
      pairs_[12][0].emplace_back(std::make_pair(3, true));
      pairs_[12][0].emplace_back(std::make_pair(1, true));
      pairs_[12][0].emplace_back(std::make_pair(8, true));
      pairs_[12][0].emplace_back(std::make_pair(6, false));
      pairs_[12][0].emplace_back(std::make_pair(4, false));
      pairs_[12][0].emplace_back(std::make_pair(0, false));
      pairs_[12][1].emplace_back(std::make_pair(4, false));
      pairs_[12][1].emplace_back(std::make_pair(2, false));
      pairs_[12][1].emplace_back(std::make_pair(8, true));
      pairs_[12][1].emplace_back(std::make_pair(7, true));
      pairs_[12][1].emplace_back(std::make_pair(3, true));
      pairs_[12][1].emplace_back(std::make_pair(0, true));
      pairs_[12][1].emplace_back(std::make_pair(4, true));
      pairs_[12][1].emplace_back(std::make_pair(6, true));
      pairs_[12][1].emplace_back(std::make_pair(8, false));
      pairs_[12][1].emplace_back(std::make_pair(1, false));
      pairs_[12][1].emplace_back(std::make_pair(3, false));
      pairs_[12][1].emplace_back(std::make_pair(5, false));
    }

    std::array<std::array<std::vector<std::pair<int, bool>>, 2>, kdop_directions> pairs_;
  };
}  // namespace CORE::GEOMETRICSEARCH

#endif
