/*----------------------------------------------------------------------*/
/*! \file

\brief Unittests for the geometric search class

\level 3

*-----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "baci_utils_exceptions.H"

#include <array>

#ifdef BACI_WITH_ARBORX

#include "baci_discretization_geometric_search_distributed_tree.H"
#include "baci_geometric_search_create_bounding_volumes_test.H"

#include <Epetra_MpiComm.h>

namespace
{
  class GeometricSearchDistributed : public ::testing::Test
  {
   public:
    static void SetUpTestSuite() { Kokkos::initialize(); }

    static void TearDownTestSuite() { Kokkos::finalize(); }

    GeometricSearchDistributed()
    {
      comm_ = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
      verbosity_ = IO::minimal;
      my_rank_ = comm_->MyPID();
    }

   protected:
    std::vector<std::pair<int, CORE::GEOMETRICSEARCH::BoundingVolume>> primitives_, predicates_;
    Teuchos::RCP<Epetra_Comm> comm_;
    int my_rank_;
    IO::verbositylevel verbosity_;
  };

  /**
   * Test setup is based on https://github.com/arborx/ArborX/issues/867
   * Checking collision search of kdop primitives against kdop predicates
   */
  TEST_F(GeometricSearchDistributed, CollisionSearchKDOPVsKDOP)
  {
    const auto volumes = CreateKDOPBoundingVolumes();

    // Add the volumes to the primitives
    if (my_rank_ == 0)
    {
      primitives_.emplace_back(std::pair{0, volumes[0]});
      primitives_.emplace_back(std::pair{1, volumes[2]});
      primitives_.emplace_back(std::pair{2, volumes[1]});
      predicates_.emplace_back(std::pair{3, volumes[0]});
      predicates_.emplace_back(std::pair{4, volumes[2]});
      predicates_.emplace_back(std::pair{5, volumes[1]});
    }
    else if (my_rank_ == 1)
    {
      primitives_.emplace_back(std::pair{6, volumes[0]});
      primitives_.emplace_back(std::pair{7, volumes[1]});
      predicates_.emplace_back(std::pair{8, volumes[2]});
    }
    else if (my_rank_ == 2)
    {
      primitives_.emplace_back(std::pair{9, volumes[0]});
      primitives_.emplace_back(std::pair{10, volumes[2]});
      predicates_.emplace_back(std::pair{11, volumes[1]});
    }

    const auto pairs =
        CORE::GEOMETRICSEARCH::GlobalCollisionSearch(primitives_, predicates_, *comm_, verbosity_);

    // The order of the results we get are is not deterministic. Therefore, we save the reference
    // results in a map, with the colliding GIDs being the keys. Thus the ordering of the pairs in
    // the collision search does not matter.
    auto compare_results = [](const auto& pairs_from_collision_search, const auto& reference_map)
    {
      EXPECT_EQ(pairs_from_collision_search.size(), reference_map.size());
      for (const auto& pair : pairs_from_collision_search)
      {
        if (reference_map.count({std::get<1>(pair), std::get<3>(pair)}))
        {
          const auto& reference_pair = reference_map.at({std::get<1>(pair), std::get<3>(pair)});
          EXPECT_EQ(std::get<0>(pair), std::get<0>(reference_pair));
          EXPECT_EQ(std::get<2>(pair), std::get<1>(reference_pair));
          EXPECT_EQ(std::get<4>(pair), std::get<2>(reference_pair));
        }
        else
        {
          dserror("Pair {%d, %d} not found in reference map", std::get<1>(pair), std::get<3>(pair));
        }
      }
    };

    if (my_rank_ == 0)
    {
      std::map<std::pair<int, int>, std::array<int, 4>> reference_pairs{
          {{3, 0}, {0, 0, 0}},   //
          {{3, 1}, {0, 1, 0}},   //
          {{3, 6}, {0, 0, 1}},   //
          {{3, 9}, {0, 0, 2}},   //
          {{3, 10}, {0, 1, 2}},  //
          {{4, 0}, {1, 0, 0}},   //
          {{4, 1}, {1, 1, 0}},   //
          {{4, 6}, {1, 0, 1}},   //
          {{4, 9}, {1, 0, 2}},   //
          {{4, 10}, {1, 1, 2}},  //
          {{5, 2}, {2, 2, 0}},   //
          {{5, 7}, {2, 1, 1}}    //
      };
      compare_results(pairs, reference_pairs);
    }
    else if (my_rank_ == 1)
    {
      std::map<std::pair<int, int>, std::array<int, 4>> reference_pairs{
          {{8, 0}, {0, 0, 0}},   //
          {{8, 1}, {0, 1, 0}},   //
          {{8, 6}, {0, 0, 1}},   //
          {{8, 9}, {0, 0, 2}},   //
          {{8, 10}, {0, 1, 2}},  //
      };
      compare_results(pairs, reference_pairs);
    }
    else if (my_rank_ == 2)
    {
      std::map<std::pair<int, int>, std::array<int, 4>> reference_pairs{
          {{11, 2}, {0, 2, 0}},  //
          {{11, 7}, {0, 1, 1}}   //
      };
      compare_results(pairs, reference_pairs);
    }
  }

  /**
   * Checking collision search with no predicates on rank 0 and 1 and no primitives on rank 1 and 2
   */
  TEST_F(GeometricSearchDistributed, CollisionSearchNoPrimitivesNoPredicatesMixed)
  {
    const auto volumes = CreateKDOPBoundingVolumes();

    // Add the volumes to the primitives and predicates for rank 0 and 1
    if (my_rank_ == 0)
    {
      primitives_.emplace_back(std::pair{0, volumes[0]});
      primitives_.emplace_back(std::pair{1, volumes[2]});
      primitives_.emplace_back(std::pair{2, volumes[1]});
    }
    if (my_rank_ == 1)
    {
      predicates_.emplace_back(std::pair{3, volumes[2]});
    }

    const auto pairs =
        CORE::GEOMETRICSEARCH::GlobalCollisionSearch(primitives_, predicates_, *comm_, verbosity_);

    if (my_rank_ == 1)
    {
      EXPECT_EQ(pairs.size(), 2);

      EXPECT_EQ(std::get<0>(pairs[0]), 0);
      EXPECT_EQ(std::get<1>(pairs[0]), 3);
      EXPECT_EQ(std::get<2>(pairs[0]), 1);
      EXPECT_EQ(std::get<3>(pairs[0]), 1);
      EXPECT_EQ(std::get<4>(pairs[0]), 0);

      EXPECT_EQ(std::get<0>(pairs[1]), 0);
      EXPECT_EQ(std::get<1>(pairs[1]), 3);
      EXPECT_EQ(std::get<2>(pairs[1]), 0);
      EXPECT_EQ(std::get<3>(pairs[1]), 0);
      EXPECT_EQ(std::get<4>(pairs[1]), 0);
    }
    else
    {
      EXPECT_EQ(pairs.size(), 0);
    }
  }

  /**
   * Checking collision search with no primitives and no predicates
   */
  TEST_F(GeometricSearchDistributed, CollisionSearchNoPrimitivesNoPredicatesAll)
  {
    EXPECT_EQ(primitives_.size(), 0);
    EXPECT_EQ(predicates_.size(), 0);

    const auto pairs =
        CORE::GEOMETRICSEARCH::GlobalCollisionSearch(primitives_, predicates_, *comm_, verbosity_);

    EXPECT_EQ(pairs.size(), 0);
  }
}  // namespace

#endif
