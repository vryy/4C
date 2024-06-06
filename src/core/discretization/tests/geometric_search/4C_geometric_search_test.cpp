/*----------------------------------------------------------------------*/
/*! \file

\brief Unittests for the geometric search class

\level 3

*-----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "4C_config.hpp"

#ifdef FOUR_C_WITH_ARBORX

#include "4C_discretization_geometric_search_bvh.hpp"
#include "4C_discretization_geometric_search_utils.hpp"
#include "4C_geometric_search_create_bounding_volumes_test.hpp"

#include <Epetra_SerialComm.h>

namespace
{
  using namespace FourC;

  class GeometricSearch : public ::testing::Test
  {
   public:
    static void SetUpTestSuite() { Kokkos::initialize(); }

    static void TearDownTestSuite() { Kokkos::finalize(); }

    GeometricSearch()
    {
      comm_ = Epetra_SerialComm();
      verbosity_ = Core::IO::minimal;
    }

   protected:
    std::vector<std::pair<int, Core::GeometricSearch::BoundingVolume>> primitives_, predicates_;
    Epetra_SerialComm comm_;
    Core::IO::Verbositylevel verbosity_;
  };

  /**
   * Test setup is based on https://github.com/arborx/ArborX/issues/867
   * Checking collision search of kdop primitives against kdop predicates
   */
  TEST_F(GeometricSearch, CollisionSearchKDOPVsKDOP)
  {
    const auto volumes = CreateKDOPBoundingVolumes();

    // Add the volumes to the primitives and predicates
    primitives_.emplace_back(std::pair{0, volumes[0]});
    primitives_.emplace_back(std::pair{1, volumes[1]});
    predicates_.emplace_back(std::pair{2, volumes[2]});

    EXPECT_EQ(primitives_.size(), 2);
    EXPECT_EQ(predicates_.size(), 1);

    const auto &[indices, offsets] =
        Core::GeometricSearch::CollisionSearch(primitives_, predicates_, comm_, verbosity_);

    const auto pairs = Core::GeometricSearch::GetPairs(indices, offsets);

    EXPECT_EQ(pairs.size(), 1);
    EXPECT_EQ(pairs[0].first, 0);
    EXPECT_EQ(pairs[0].second, 0);
  }

  /**
   * Checking collision search with no primitives
   */
  TEST_F(GeometricSearch, CollisionSearchNoPrimitives)
  {
    const auto volumes = CreateKDOPBoundingVolumes();

    // Add the volumes to the primitives and predicates
    predicates_.emplace_back(std::pair{0, volumes[0]});
    predicates_.emplace_back(std::pair{1, volumes[1]});
    predicates_.emplace_back(std::pair{2, volumes[2]});

    EXPECT_EQ(primitives_.size(), 0);
    EXPECT_EQ(predicates_.size(), 3);

    const auto &[indices, offsets] =
        Core::GeometricSearch::CollisionSearch(primitives_, predicates_, comm_, verbosity_);

    const auto pairs = Core::GeometricSearch::GetPairs(indices, offsets);

    EXPECT_EQ(pairs.size(), 0);
  }

  /**
   * Checking collision search with no predicates
   */
  TEST_F(GeometricSearch, CollisionSearchNoPredicates)
  {
    const auto volumes = CreateKDOPBoundingVolumes();

    // Add the volumes to the primitives and predicates
    primitives_.emplace_back(std::pair{0, volumes[0]});
    primitives_.emplace_back(std::pair{1, volumes[1]});
    primitives_.emplace_back(std::pair{2, volumes[2]});

    EXPECT_EQ(primitives_.size(), 3);
    EXPECT_EQ(predicates_.size(), 0);

    const auto &[indices, offsets] =
        Core::GeometricSearch::CollisionSearch(primitives_, predicates_, comm_, verbosity_);

    const auto pairs = Core::GeometricSearch::GetPairs(indices, offsets);

    EXPECT_EQ(pairs.size(), 0);
  }

  /**
   * Checking collision search with no primitives and no predicates
   */
  TEST_F(GeometricSearch, CollisionSearchNoPrimitivesNoPredicates)
  {
    EXPECT_EQ(primitives_.size(), 0);
    EXPECT_EQ(predicates_.size(), 0);

    const auto &[indices, offsets] =
        Core::GeometricSearch::CollisionSearch(primitives_, predicates_, comm_, verbosity_);

    const auto pairs = Core::GeometricSearch::GetPairs(indices, offsets);

    EXPECT_EQ(pairs.size(), 0);
  }
}  // namespace

#endif
