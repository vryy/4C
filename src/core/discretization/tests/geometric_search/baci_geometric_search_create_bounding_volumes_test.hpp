/*----------------------------------------------------------------------*/
/*! \file

\brief Create bounding volumes for the geometric search unit tests

\level 3

*-----------------------------------------------------------------------*/

#ifndef FOUR_C_GEOMETRIC_SEARCH_CREATE_BOUNDING_VOLUMES_TEST_HPP
#define FOUR_C_GEOMETRIC_SEARCH_CREATE_BOUNDING_VOLUMES_TEST_HPP

#ifdef FOUR_C_WITH_ARBORX

#include "baci_config.hpp"

#include "baci_discretization_geometric_search_bounding_volume.hpp"
#include "baci_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  /**
   * Create a set of three different bounding volumes that is used for some unit tests.
   * The geometry is based on https://github.com/arborx/ArborX/issues/867
   */
  std::vector<CORE::GEOMETRICSEARCH::BoundingVolume> CreateKDOPBoundingVolumes()
  {
    std::vector<CORE::GEOMETRICSEARCH::BoundingVolume> volumes(3);

    CORE::LINALG::Matrix<3, 1, double> point(true);

    // setting up bounding volume 1
    {
      point(0) = 0.0;
      point(1) = 0.0;
      point(2) = 0.0;
      volumes[0].AddPoint(point);

      point(0) = 1.0;
      point(1) = 1.0;
      point(2) = 0.0;
      volumes[0].AddPoint(point);

      point(0) = 0.5;
      point(1) = 0.0;
      point(2) = 0.0;
      volumes[0].AddPoint(point);
    }

    // setting up bounding volume 2
    {
      point(0) = 1.0;
      point(1) = -1.0;
      point(2) = 0.0;
      volumes[1].AddPoint(point);

      point(0) = 1.0;
      point(1) = 0.25;
      point(2) = 0.0;
      volumes[1].AddPoint(point);

      point(0) = 0.75;
      point(1) = 0.0;
      point(2) = 0.0;
      volumes[1].AddPoint(point);
    }

    // setting up bounding volume 3
    {
      point(0) = 0.5;
      point(1) = 0.25;
      point(2) = 0.0;
      volumes[2].AddPoint(point);

      point(0) = 0.8;
      point(1) = 0.25;
      point(2) = 0.0;
      volumes[2].AddPoint(point);

      point(0) = 0.75;
      point(1) = 0.125;
      point(2) = 0.0;
      volumes[2].AddPoint(point);
    }

    return volumes;
  }
}  // namespace

FOUR_C_NAMESPACE_CLOSE

#endif

#endif
