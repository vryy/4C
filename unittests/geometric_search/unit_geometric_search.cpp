/*----------------------------------------------------------------------*/
/*! \file

\brief Unittests for the geometric search class

\level 3

*-----------------------------------------------------------------------*/

#include <gtest/gtest.h>
#include <Epetra_SerialComm.h>

#ifdef BACI_WITH_ARBORX

#include "baci_discretization_geometric_search.H"
#include "baci_discretization_geometric_search_utils.H"

const double testTolerance = 1e-14;

namespace
{
  class GeometricSearch : public ::testing::Test
  {
   public:
    GeometricSearch()
    {
      comm_ = Epetra_SerialComm();
      verbosity_ = IO::minimal;
    }

   protected:
    std::vector<std::pair<int, CORE::GEOMETRICSEARCH::BoundingVolume>> primitives_, predicates_;
    Epetra_SerialComm comm_;
    IO::verbositylevel verbosity_;
  };

  /**
   * Test setup is based on https://github.com/arborx/ArborX/issues/867
   * Checking collision search of kdop primitives against kdop predicates
   */
  TEST_F(GeometricSearch, CollisionSearchKDOPVsKDOP)
  {
    Kokkos::ScopeGuard kokkos_guard;

    CORE::GEOMETRICSEARCH::BoundingVolume boundingVolume1, boundingVolume2, boundingVolume3;

    // setting up bounding volume 1
    {
      CORE::LINALG::Matrix<3, 1, double> coord1(true);
      coord1(0) = 0.0;
      coord1(1) = 0.0;
      coord1(2) = 0.0;
      boundingVolume1.AddPoint(coord1);

      CORE::LINALG::Matrix<3, 1, double> coord2(true);
      coord2(0) = 1.0;
      coord2(1) = 1.0;
      coord2(2) = 0.0;
      boundingVolume1.AddPoint(coord2);

      CORE::LINALG::Matrix<3, 1, double> coord3(true);
      coord3(0) = 0.5;
      coord3(1) = 0.0;
      coord3(2) = 0.0;
      boundingVolume1.AddPoint(coord3);

      primitives_.emplace_back(std::pair{0, boundingVolume1});
    }

    // setting up bounding volume 2
    {
      CORE::LINALG::Matrix<3, 1, double> coord1(true);
      coord1(0) = 1.0;
      coord1(1) = -1.0;
      coord1(2) = 0.0;
      boundingVolume2.AddPoint(coord1);

      CORE::LINALG::Matrix<3, 1, double> coord2(true);
      coord2(0) = 1.0;
      coord2(1) = 0.25;
      coord2(2) = 0.0;
      boundingVolume2.AddPoint(coord2);

      CORE::LINALG::Matrix<3, 1, double> coord3(true);
      coord3(0) = 0.75;
      coord3(1) = 0.0;
      coord3(2) = 0.0;
      boundingVolume2.AddPoint(coord3);

      primitives_.emplace_back(std::pair{1, boundingVolume2});
    }

    // setting up bounding volume 3
    {
      CORE::LINALG::Matrix<3, 1, double> coord1(true);
      coord1(0) = 0.5;
      coord1(1) = 0.25;
      coord1(2) = 0.0;
      boundingVolume3.AddPoint(coord1);

      CORE::LINALG::Matrix<3, 1, double> coord2(true);
      coord2(0) = 0.8;
      coord2(1) = 0.25;
      coord2(2) = 0.0;
      boundingVolume3.AddPoint(coord2);

      CORE::LINALG::Matrix<3, 1, double> coord3(true);
      coord3(0) = 0.75;
      coord3(1) = 0.125;
      coord3(2) = 0.0;
      boundingVolume3.AddPoint(coord3);

      predicates_.emplace_back(std::pair{0, boundingVolume3});
    }

    EXPECT_EQ(primitives_.size(), 2);
    EXPECT_EQ(predicates_.size(), 1);

    const auto &[indices, offsets] =
        CORE::GEOMETRICSEARCH::CollisionSearch(primitives_, predicates_, comm_, verbosity_);
    const auto pairs = CORE::GEOMETRICSEARCH::GetPairs(indices, offsets);

    EXPECT_EQ(pairs.size(), 1);
    EXPECT_EQ(pairs[0].first, 0);
    EXPECT_EQ(pairs[0].second, 0);
  }
}  // namespace

#endif
