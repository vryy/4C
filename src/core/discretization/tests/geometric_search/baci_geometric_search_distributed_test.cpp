/*----------------------------------------------------------------------*/
/*! \file

\brief Unittests for the geometric search class

\level 3

*-----------------------------------------------------------------------*/

#include <gtest/gtest.h>

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

    if (my_rank_ == 0)
    {
      EXPECT_EQ(pairs.size(), 12);

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

      EXPECT_EQ(std::get<0>(pairs[2]), 0);
      EXPECT_EQ(std::get<1>(pairs[2]), 3);
      EXPECT_EQ(std::get<2>(pairs[2]), 0);
      EXPECT_EQ(std::get<3>(pairs[2]), 6);
      EXPECT_EQ(std::get<4>(pairs[2]), 1);

      EXPECT_EQ(std::get<0>(pairs[3]), 0);
      EXPECT_EQ(std::get<1>(pairs[3]), 3);
      EXPECT_EQ(std::get<2>(pairs[3]), 1);
      EXPECT_EQ(std::get<3>(pairs[3]), 10);
      EXPECT_EQ(std::get<4>(pairs[3]), 2);

      EXPECT_EQ(std::get<0>(pairs[4]), 0);
      EXPECT_EQ(std::get<1>(pairs[4]), 3);
      EXPECT_EQ(std::get<2>(pairs[4]), 0);
      EXPECT_EQ(std::get<3>(pairs[4]), 9);
      EXPECT_EQ(std::get<4>(pairs[4]), 2);

      EXPECT_EQ(std::get<0>(pairs[5]), 1);
      EXPECT_EQ(std::get<1>(pairs[5]), 4);
      EXPECT_EQ(std::get<2>(pairs[5]), 1);
      EXPECT_EQ(std::get<3>(pairs[5]), 1);
      EXPECT_EQ(std::get<4>(pairs[5]), 0);

      EXPECT_EQ(std::get<0>(pairs[6]), 1);
      EXPECT_EQ(std::get<1>(pairs[6]), 4);
      EXPECT_EQ(std::get<2>(pairs[6]), 0);
      EXPECT_EQ(std::get<3>(pairs[6]), 0);
      EXPECT_EQ(std::get<4>(pairs[6]), 0);

      EXPECT_EQ(std::get<0>(pairs[7]), 1);
      EXPECT_EQ(std::get<1>(pairs[7]), 4);
      EXPECT_EQ(std::get<2>(pairs[7]), 0);
      EXPECT_EQ(std::get<3>(pairs[7]), 6);
      EXPECT_EQ(std::get<4>(pairs[7]), 1);

      EXPECT_EQ(std::get<0>(pairs[8]), 1);
      EXPECT_EQ(std::get<1>(pairs[8]), 4);
      EXPECT_EQ(std::get<2>(pairs[8]), 1);
      EXPECT_EQ(std::get<3>(pairs[8]), 10);
      EXPECT_EQ(std::get<4>(pairs[8]), 2);

      EXPECT_EQ(std::get<0>(pairs[9]), 1);
      EXPECT_EQ(std::get<1>(pairs[9]), 4);
      EXPECT_EQ(std::get<2>(pairs[9]), 0);
      EXPECT_EQ(std::get<3>(pairs[9]), 9);
      EXPECT_EQ(std::get<4>(pairs[9]), 2);

      EXPECT_EQ(std::get<0>(pairs[10]), 2);
      EXPECT_EQ(std::get<1>(pairs[10]), 5);
      EXPECT_EQ(std::get<2>(pairs[10]), 2);
      EXPECT_EQ(std::get<3>(pairs[10]), 2);
      EXPECT_EQ(std::get<4>(pairs[10]), 0);

      EXPECT_EQ(std::get<0>(pairs[11]), 2);
      EXPECT_EQ(std::get<1>(pairs[11]), 5);
      EXPECT_EQ(std::get<2>(pairs[11]), 1);
      EXPECT_EQ(std::get<3>(pairs[11]), 7);
      EXPECT_EQ(std::get<4>(pairs[11]), 1);
    }
    else if (my_rank_ == 1)
    {
      EXPECT_EQ(pairs.size(), 5);

      EXPECT_EQ(std::get<0>(pairs[0]), 0);
      EXPECT_EQ(std::get<1>(pairs[0]), 8);
      EXPECT_EQ(std::get<2>(pairs[0]), 1);
      EXPECT_EQ(std::get<3>(pairs[0]), 1);
      EXPECT_EQ(std::get<4>(pairs[0]), 0);

      EXPECT_EQ(std::get<0>(pairs[1]), 0);
      EXPECT_EQ(std::get<1>(pairs[1]), 8);
      EXPECT_EQ(std::get<2>(pairs[1]), 0);
      EXPECT_EQ(std::get<3>(pairs[1]), 0);
      EXPECT_EQ(std::get<4>(pairs[1]), 0);

      EXPECT_EQ(std::get<0>(pairs[2]), 0);
      EXPECT_EQ(std::get<1>(pairs[2]), 8);
      EXPECT_EQ(std::get<2>(pairs[2]), 0);
      EXPECT_EQ(std::get<3>(pairs[2]), 6);
      EXPECT_EQ(std::get<4>(pairs[2]), 1);

      EXPECT_EQ(std::get<0>(pairs[3]), 0);
      EXPECT_EQ(std::get<1>(pairs[3]), 8);
      EXPECT_EQ(std::get<2>(pairs[3]), 1);
      EXPECT_EQ(std::get<3>(pairs[3]), 10);
      EXPECT_EQ(std::get<4>(pairs[3]), 2);

      EXPECT_EQ(std::get<0>(pairs[4]), 0);
      EXPECT_EQ(std::get<1>(pairs[4]), 8);
      EXPECT_EQ(std::get<2>(pairs[4]), 0);
      EXPECT_EQ(std::get<3>(pairs[4]), 9);
      EXPECT_EQ(std::get<4>(pairs[4]), 2);
    }
    else if (my_rank_ == 2)
    {
      EXPECT_EQ(pairs.size(), 2);

      EXPECT_EQ(std::get<0>(pairs[0]), 0);
      EXPECT_EQ(std::get<1>(pairs[0]), 11);
      EXPECT_EQ(std::get<2>(pairs[0]), 2);
      EXPECT_EQ(std::get<3>(pairs[0]), 2);
      EXPECT_EQ(std::get<4>(pairs[0]), 0);

      EXPECT_EQ(std::get<0>(pairs[1]), 0);
      EXPECT_EQ(std::get<1>(pairs[1]), 11);
      EXPECT_EQ(std::get<2>(pairs[1]), 1);
      EXPECT_EQ(std::get<3>(pairs[1]), 7);
      EXPECT_EQ(std::get<4>(pairs[1]), 1);
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
