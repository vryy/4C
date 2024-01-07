/*----------------------------------------------------------------------*/
/*! \file
\brief Testcases for the helpfer functions for adaptive history integation of full constrained
mixture fibers
\level 3
*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include <gmock/gmock.h>

#include "baci_mixture_full_constrained_mixture_fiber_adaptive_history.H"

#include "baci_mixture_growth_evolution_linear_cauchy_poisson_turnover.H"

#include <Sacado.hpp>

#include <bitset>
#include <memory>

namespace
{
  using namespace BACI;

  TEST(TimestepAdaptivityTest, EmplaceBack)
  {
    MIXTURE::TimestepAdaptivityInfo level{};

    level.EmplaceBack(3, 5);
    EXPECT_EQ(level.GetNumberOfLevels(), 1);
    EXPECT_EQ(level.GetTotalNumberOfSimpsonIntervals(), 5);

    level.EmplaceBack(2, 3);
    EXPECT_EQ(level.GetNumberOfLevels(), 2);
    EXPECT_EQ(level.GetTotalNumberOfSimpsonIntervals(), 8);

    level.EmplaceBack(2, 6);
    EXPECT_EQ(level.GetNumberOfLevels(), 2);
    EXPECT_EQ(level.GetTotalNumberOfSimpsonIntervals(), 14);
  }

  TEST(TimestepAdaptivityTest, EmplaceBackInWrongOrder)
  {
#ifndef BACI_DEBUG
    GTEST_SKIP() << "Skip debug assertion tests in release mode.";
#endif
    MIXTURE::TimestepAdaptivityInfo level{};

    level.EmplaceBack(1, 5);

    EXPECT_ANY_THROW(level.EmplaceBack(2, 1));
  }

  TEST(TimestepAdaptivityTest, SplitLevel)
  {
    MIXTURE::TimestepAdaptivityInfo level{};
    level.EmplaceBack(3, 6);
    level.EmplaceBack(2, 12);
    level.EmplaceBack(1, 8);

    level.SplitLevel(3, 2);
    EXPECT_EQ(level[0], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(4, 2));
    EXPECT_EQ(level[1], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(3, 2));
    EXPECT_EQ(level[2], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(2, 12));
    EXPECT_EQ(level[3], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(1, 8));
    EXPECT_EQ(level.GetNumberOfLevels(), 4);

    level.SplitLevel(4, 1);
    EXPECT_EQ(level[0], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(5, 1));
    EXPECT_EQ(level[1], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(3, 2));
    EXPECT_EQ(level[2], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(2, 12));
    EXPECT_EQ(level[3], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(1, 8));
    EXPECT_EQ(level.GetNumberOfLevels(), 4);

    level.SplitLevel(1, 1);
    EXPECT_EQ(level[0], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(5, 1));
    EXPECT_EQ(level[1], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(3, 2));
    EXPECT_EQ(level[2], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(2, 13));
    EXPECT_EQ(level[3], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(1, 6));
    EXPECT_EQ(level.GetNumberOfLevels(), 4);

    level.SplitLevel(1, 3);
    EXPECT_EQ(level[0], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(5, 1));
    EXPECT_EQ(level[1], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(3, 2));
    EXPECT_EQ(level[2], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(2, 16));
    EXPECT_EQ(level.GetNumberOfLevels(), 3);

    EXPECT_ANY_THROW(level.SplitLevel(1, 0));
  }

  TEST(TimestepAdaptivityTest, GetTotalNumberOfSimpsonIntervals)
  {
    MIXTURE::TimestepAdaptivityInfo level{};
    level.EmplaceBack(3, 2);
    level.EmplaceBack(2, 1);

    EXPECT_EQ(level.GetTotalNumberOfSimpsonIntervals(), 3);
  }

  TEST(TimestepAdaptivityTest, MaxLevel)
  {
    MIXTURE::TimestepAdaptivityInfo level{};
    level.EmplaceBack(3, 2);
    level.EmplaceBack(2, 1);

    EXPECT_EQ(level.MaxLevel(), 3);
  }

  TEST(TimestepAdaptivityTest, GetNumberOfSimpsonIntervals)
  {
    MIXTURE::TimestepAdaptivityInfo level{};
    level.EmplaceBack(3, 2);
    level.EmplaceBack(2, 1);

    EXPECT_EQ(level.GetNumberOfSimpsonIntervals(3), 2);
    EXPECT_EQ(level.GetNumberOfSimpsonIntervals(2), 1);
    EXPECT_EQ(level.GetNumberOfSimpsonIntervals(1), 0);
  }

  TEST(TimestepAdaptivityTest, GetBeginIndex)
  {
    MIXTURE::TimestepAdaptivityInfo level{};
    level.EmplaceBack(3, 2);
    level.EmplaceBack(2, 1);

    EXPECT_EQ(level.GetBeginIndex(3), 0);
    EXPECT_EQ(level.GetBeginIndex(2), 4);
    EXPECT_EQ(level.GetBeginIndex(1), 6);
    EXPECT_EQ(level.GetBeginIndex(0), 6);
  }

  TEST(TimestepAdaptivityTest, GetBeginTime)
  {
    MIXTURE::TimestepAdaptivityInfo level{};
    level.EmplaceBack(3, 2);
    level.EmplaceBack(2, 1);

    EXPECT_NEAR(level.GetBeginTime(3, 1.0, 0.1), 1.0, 1e-15);
    EXPECT_NEAR(level.GetBeginTime(2, 1.0, 0.1), 4.2, 1e-15);
    EXPECT_NEAR(level.GetBeginTime(1, 1.0, 0.1), 5.0, 1e-15);
    EXPECT_NEAR(level.GetBeginTime(0, 1.0, 0.1), 5.0, 1e-15);
  }

  TEST(TimestepAdaptivityTest, GetIndexTime)
  {
    MIXTURE::TimestepAdaptivityInfo level{};
    level.EmplaceBack(3, 2);
    level.EmplaceBack(2, 1);

    EXPECT_NEAR(level.GetIndexTime(0, 1.0, 0.1), 1.0, 1e-15);
    EXPECT_NEAR(level.GetIndexTime(1, 1.0, 0.1), 1.8, 1e-15);
    EXPECT_NEAR(level.GetIndexTime(2, 1.0, 0.1), 2.6, 1e-15);
    EXPECT_NEAR(level.GetIndexTime(3, 1.0, 0.1), 3.4, 1e-15);
    EXPECT_NEAR(level.GetIndexTime(4, 1.0, 0.1), 4.2, 1e-15);
    EXPECT_NEAR(level.GetIndexTime(5, 1.0, 0.1), 4.6, 1e-15);
    EXPECT_NEAR(level.GetIndexTime(6, 1.0, 0.1), 5.0, 1e-15);
    EXPECT_NEAR(level.GetIndexTime(7, 1.0, 0.1), 5.1, 1e-15);
    EXPECT_NEAR(level.GetIndexTime(8, 1.0, 0.1), 5.2, 1e-15);
  }

  TEST(TimestepAdaptivityTest, GetNumberOfLevels)
  {
    MIXTURE::TimestepAdaptivityInfo level{};
    level.EmplaceBack(3, 2);
    level.EmplaceBack(2, 1);

    EXPECT_EQ(level.GetNumberOfLevels(), 2);
  }

  TEST(TimestepAdaptivityTest, OptimizeHistoryFromUnoptimizedHistory)
  {
    MIXTURE::TimestepAdaptivityInfo level{};
    const MIXTURE::LinearCauchyGrowthWithPoissonTurnoverGrowthEvolution<double> growth_evolution(
        0.0, 1.0);
    const double time = 30;
    const double dt = 1.0;
    const double tolerance = 1e-8 / 24;

    const auto [remove_items, new_level] = MIXTURE::OptimizeHistoryIntegration(level, 25,
        [&](const std::array<std::optional<unsigned int>, 5>& indices)
        {
          const double begin_time = level.GetIndexTime(indices[0].value(), 0.0, dt);
          const double end_time = level.GetIndexTime(indices[4].value(), 0.0, dt);
          return MIXTURE::IsModelEquationSimpsonRuleIntegrationBelowTolerance(
              growth_evolution, time, begin_time, end_time, tolerance * (end_time - begin_time));
        });

    level = new_level;

    EXPECT_EQ(level[0], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(2, 1));
    EXPECT_EQ(level[1], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(1, 1));

    std::array comparison = {false, true, true, true, false, true, true, true, false, true, false,
        true, false, false, false, false, false, false, false, false, false, false, false, false,
        false};

    EXPECT_THAT(remove_items, testing::ElementsAreArray(comparison));
  }

  TEST(TimestepAdaptivityTest, OptimizeHistoryFromPartiallyOptimizedHistory)
  {
    MIXTURE::TimestepAdaptivityInfo level{};
    level.EmplaceBack(1, 2);

    const MIXTURE::LinearCauchyGrowthWithPoissonTurnoverGrowthEvolution<double> growth_evolution(
        0.0, 1.0);
    const double time = 30;
    const double dt = 1.0;
    const double tolerance = 1e-8 / 24;

    const auto [remove_items, new_level] = MIXTURE::OptimizeHistoryIntegration(level, 21,
        [&](const std::array<std::optional<unsigned int>, 5>& indices)
        {
          const double begin_time = level.GetIndexTime(indices[0].value(), 0.0, dt);
          const double end_time = level.GetIndexTime(indices[4].value(), 0.0, dt);
          return MIXTURE::IsModelEquationSimpsonRuleIntegrationBelowTolerance(
              growth_evolution, time, begin_time, end_time, tolerance * (end_time - begin_time));
        });


    EXPECT_EQ(new_level[0], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(2, 1));
    EXPECT_EQ(new_level[1], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(1, 1));

    std::array comparison = {false, true, false, true, false, true, false, true, false, false,
        false, false, false, false, false, false, false, false, false, false, false};

    EXPECT_THAT(remove_items, testing::ElementsAreArray(comparison));
  }

  TEST(TimestepAdaptivityTest, OptimizeHistoryFromOptimizedHistory)
  {
    MIXTURE::TimestepAdaptivityInfo level{};
    level.EmplaceBack(2, 1);
    level.EmplaceBack(1, 1);
    const MIXTURE::LinearCauchyGrowthWithPoissonTurnoverGrowthEvolution<double> growth_evolution(
        0.0, 1.0);
    const double time = 30;
    const double dt = 1.0;
    const double tolerance = 1e-8 / 24;

    const auto [remove_items, new_level] = MIXTURE::OptimizeHistoryIntegration(level, 17,
        [&](const std::array<std::optional<unsigned int>, 5>& indices)
        {
          const double begin_time = level.GetIndexTime(indices[0].value(), 0.0, dt);
          const double end_time = level.GetIndexTime(indices[4].value(), 0.0, dt);
          return MIXTURE::IsModelEquationSimpsonRuleIntegrationBelowTolerance(
              growth_evolution, time, begin_time, end_time, tolerance * (end_time - begin_time));
        });

    EXPECT_EQ(level[0], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(2, 1));
    EXPECT_EQ(level[1], MIXTURE::TimestepAdaptivityInfo::TimestepAdaptivityInfoItem(1, 1));

    std::array comparison = {false, false, false, false, false, false, false, false, false, false,
        false, false, false, false, false, false, false};

    EXPECT_THAT(remove_items, testing::ElementsAreArray(comparison));
  }

  TEST(TimestepAdaptivityTest, GetBaseTimestep)
  {
    MIXTURE::TimestepAdaptivityInfo level{};
    level.EmplaceBack(3, 1);
    level.EmplaceBack(2, 2);
    level.EmplaceBack(1, 1);

    EXPECT_EQ(level.GetBaseIndex(0), 0);
    EXPECT_EQ(level.GetBaseIndex(1), 8);
    EXPECT_EQ(level.GetBaseIndex(2), 16);
    EXPECT_EQ(level.GetBaseIndex(3), 16 + 4);
    EXPECT_EQ(level.GetBaseIndex(4), 16 + 8);
    EXPECT_EQ(level.GetBaseIndex(5), 16 + 12);
    EXPECT_EQ(level.GetBaseIndex(6), 16 + 16);
    EXPECT_EQ(level.GetBaseIndex(7), 16 + 16 + 2);
    EXPECT_EQ(level.GetBaseIndex(8), 16 + 16 + 4);
    EXPECT_EQ(level.GetBaseIndex(9), 16 + 16 + 4 + 1);
    EXPECT_EQ(level.GetBaseIndex(10), 16 + 16 + 4 + 2);
  }

  TEST(TimestepAdaptivityTest, GetBaseIndices)
  {
    MIXTURE::TimestepAdaptivityInfo level{};
    level.EmplaceBack(3, 1);
    level.EmplaceBack(2, 2);
    level.EmplaceBack(1, 1);

    auto base_indices = level.GetBaseIndices<unsigned int, 11>({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10});

    EXPECT_THAT(base_indices, testing::ElementsAre(0, 8, 16, 20, 24, 28, 32, 34, 36, 37, 38));
  }

  TEST(TimestepAdaptivityTest, GetIndexFromBase)
  {
    MIXTURE::TimestepAdaptivityInfo level{};
    level.EmplaceBack(2, 2);
    level.EmplaceBack(1, 1);

    EXPECT_EQ(level.GetIndexFromBase(0), 0);
    EXPECT_FALSE(level.GetIndexFromBase(1).has_value());
    EXPECT_FALSE(level.GetIndexFromBase(2).has_value());
    EXPECT_FALSE(level.GetIndexFromBase(3).has_value());
    EXPECT_EQ(level.GetIndexFromBase(4), 1);
    EXPECT_FALSE(level.GetIndexFromBase(5).has_value());
    EXPECT_FALSE(level.GetIndexFromBase(6).has_value());
    EXPECT_FALSE(level.GetIndexFromBase(7).has_value());
    EXPECT_EQ(level.GetIndexFromBase(8), 2);
    EXPECT_FALSE(level.GetIndexFromBase(9).has_value());
    EXPECT_FALSE(level.GetIndexFromBase(10).has_value());
    EXPECT_FALSE(level.GetIndexFromBase(11).has_value());
    EXPECT_EQ(level.GetIndexFromBase(12), 3);
    EXPECT_FALSE(level.GetIndexFromBase(13).has_value());
    EXPECT_FALSE(level.GetIndexFromBase(14).has_value());
    EXPECT_FALSE(level.GetIndexFromBase(15).has_value());
    EXPECT_EQ(level.GetIndexFromBase(16), 4);
    EXPECT_FALSE(level.GetIndexFromBase(17).has_value());
    EXPECT_EQ(level.GetIndexFromBase(18), 5);
    EXPECT_FALSE(level.GetIndexFromBase(19).has_value());
    EXPECT_EQ(level.GetIndexFromBase(20), 6);
    EXPECT_EQ(level.GetIndexFromBase(21), 7);
    EXPECT_EQ(level.GetIndexFromBase(22), 8);
  }

  TEST(TimestepAdaptivityTest, GetBaseIndexWithBase)
  {
    MIXTURE::TimestepAdaptivityInfo level{};
    level.EmplaceBack(3, 1);
    level.EmplaceBack(1, 1);

    MIXTURE::TimestepAdaptivityInfo base{};
    base.EmplaceBack(1, 4);

    EXPECT_EQ(level.GetBaseIndex(base, 0), 0);
    EXPECT_EQ(level.GetBaseIndex(base, 1), 4);
    EXPECT_EQ(level.GetBaseIndex(base, 2), 8);
    EXPECT_EQ(level.GetBaseIndex(base, 3), 10);
    EXPECT_EQ(level.GetBaseIndex(base, 4), 12);
    EXPECT_EQ(level.GetBaseIndex(base, 5), 13);
  }

  TEST(TimestepAdaptivityTest, GetBaseIndicesWithBase)
  {
    MIXTURE::TimestepAdaptivityInfo level{};
    level.EmplaceBack(3, 1);
    level.EmplaceBack(1, 1);

    MIXTURE::TimestepAdaptivityInfo base{};
    base.EmplaceBack(1, 4);

    auto base_indices = level.GetBaseIndices<unsigned int, 6>(base, {0, 1, 2, 3, 4, 5});

    EXPECT_THAT(base_indices, testing::ElementsAre(0, 4, 8, 10, 12, 13));
  }

  TEST(TimestepAdaptivityTest, GetIndicesFromBase)
  {
    MIXTURE::TimestepAdaptivityInfo level{};
    level.EmplaceBack(2, 2);
    level.EmplaceBack(1, 1);

    auto indices = level.GetIndicesFromBase<unsigned int, 23>(
        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22});
    EXPECT_THAT(indices,
        testing::ElementsAre(0, std::nullopt, std::nullopt, std::nullopt, 1, std::nullopt,
            std::nullopt, std::nullopt, 2, std::nullopt, std::nullopt, std::nullopt, 3,
            std::nullopt, std::nullopt, std::nullopt, 4, std::nullopt, 5, std::nullopt, 6, 7, 8));
  }

  TEST(TimestepAdaptivityTest, GetIndicesFromBaseWrongOrder)
  {
#ifndef BACI_DEBUG
    GTEST_SKIP() << "Skip debug assertion tests in release mode.";
#endif
    MIXTURE::TimestepAdaptivityInfo level{};
    auto get_indices_from_base = [&]() {
      return level.GetIndicesFromBase<unsigned int, 3>({0, 3, 1});
    };
    EXPECT_ANY_THROW(get_indices_from_base());
  }

  TEST(TimestepAdaptivityTest, GetBaseIndicesWrongOrder)
  {
#ifndef BACI_DEBUG
    GTEST_SKIP() << "Skip debug assertion tests in release mode.";
#endif
    MIXTURE::TimestepAdaptivityInfo level{};

    auto get_base_indices = [&]() { return level.GetBaseIndices<unsigned int, 3>({0, 3, 1}); };
    EXPECT_ANY_THROW(get_base_indices());
  }
}  // namespace