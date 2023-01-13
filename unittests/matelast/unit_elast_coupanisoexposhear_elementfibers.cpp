/*----------------------------------------------------------------------*/
/*! \file

\brief Testcases for the CoupAnisoExpoShear summand with element fibers

\level 2


*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "mat_anisotropy.H"
#include "matelast_coupanisoexpo.H"
#include "lib_voigt_notation.H"
#include "matelast_coupanisoexposhear.H"
#include "unittests_assertions.h"

namespace
{


  void SetupSingleStructuralTensor(const LINALG::Matrix<3, 1>& fiber1,
      const LINALG::Matrix<3, 1>& fiber2, LINALG::Matrix<3, 3>& structuralTensor)
  {
    LINALG::Matrix<3, 3> fiber1fiber2T(false);

    fiber1fiber2T.MultiplyNT(fiber1, fiber2);

    structuralTensor.Update(0.5, fiber1fiber2T);
    structuralTensor.UpdateT(0.5, fiber1fiber2T, 1.0);
  }

  class CoupAnisoExpoShearElementFibersTest
      : public ::testing::TestWithParam<std::tuple<std::array<int, 2>, int>>
  {
   public:
    CoupAnisoExpoShearElementFibersTest()
        : anisotropy_(),
          eleFibers_(3),
          eleTensors_(false),
          eleTensors_stress_(false),
          eleScalarProducts_(0.0)
    {
      // setup fibers fibers
      eleFibers_[0](0) = 0.858753861115007;
      eleFibers_[0](1) = 0.449823451060242;
      eleFibers_[0](2) = 0.245358246032859;

      eleFibers_[1](0) = 0.103448275862069;
      eleFibers_[1](1) = 0.137931034482759;
      eleFibers_[1](2) = 0.068965517241379;

      eleFibers_[2](0) = 0.872502871778232;
      eleFibers_[2](1) = 0.134231211042805;
      eleFibers_[2](2) = 0.469809238649817;

      // setup structural tensor
      SetupSingleStructuralTensor(
          eleFibers_[GetFiberIds()[0]], eleFibers_[GetFiberIds()[1]], eleTensors_);

      // Setup structural tensors in stress like Voigt notation
      UTILS::VOIGT::Stresses::MatrixToVector(eleTensors_, eleTensors_stress_);

      // setup scalar product
      eleScalarProducts_ = eleFibers_[GetFiberIds()[0]].Dot(eleFibers_[GetFiberIds()[1]]);

      SetupAnisotropyExtension(GetFiberIds());
    }

    void SetupAnisotropyExtension(std::array<int, 2> fiber_ids)
    {
      anisotropyExtension_ =
          std::make_unique<MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension>(1, fiber_ids);

      anisotropy_.RegisterAnisotropyExtension(*anisotropyExtension_);

      anisotropy_.SetNumberOfGaussPoints(2);

      // Setup element fibers
      anisotropy_.SetElementFibers(eleFibers_);
    }

    [[nodiscard]] int GetGaussPoint() const { return std::get<1>(GetParam()); }

    [[nodiscard]] std::array<int, 2> GetFiberIds() const { return std::get<0>(GetParam()); }

    MAT::Anisotropy anisotropy_;
    std::unique_ptr<MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension> anisotropyExtension_;

    std::vector<LINALG::Matrix<3, 1>> eleFibers_;
    LINALG::Matrix<3, 3> eleTensors_;
    LINALG::Matrix<6, 1> eleTensors_stress_;
    double eleScalarProducts_;
  };

  TEST_P(CoupAnisoExpoShearElementFibersTest, GetScalarProduct)
  {
    EXPECT_NEAR(anisotropyExtension_->GetScalarProduct(GetGaussPoint()), eleScalarProducts_, 1e-10);
  }

  TEST_P(CoupAnisoExpoShearElementFibersTest, GetStructuralTensor)
  {
    BACI_EXPECT_NEAR(
        anisotropyExtension_->GetStructuralTensor(GetGaussPoint()), eleTensors_, 1e-10);
  }

  TEST_P(CoupAnisoExpoShearElementFibersTest, GetStructuralTensorStress)
  {
    BACI_EXPECT_NEAR(anisotropyExtension_->GetStructuralTensor_stress(GetGaussPoint()),
        eleTensors_stress_, 1e-10);
  }

  INSTANTIATE_TEST_SUITE_P(GaussPoints, CoupAnisoExpoShearElementFibersTest,
      ::testing::Combine(::testing::Values(std::array<int, 2>({0, 1}), std::array<int, 2>({0, 2}),
                             std::array<int, 2>({1, 2})),
          ::testing::Values(0, 1)));
}  // namespace
