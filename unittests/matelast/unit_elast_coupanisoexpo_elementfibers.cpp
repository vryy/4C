/*----------------------------------------------------------------------*/
/*! \file

\brief Testcases for the CoupAnisoExpoBase summand with element fibers

\level 2


*/
/*----------------------------------------------------------------------*/


#include <gtest/gtest.h>
#include "unittest_utils_assertions.h"

#include <Teuchos_RCPDecl.hpp>
#include <boost/mpl/protect.hpp>
#include <tuple>

#include "linalg_fixedsizematrix.H"
#include "matelast_coupanisoexpo.H"
#include "mat_anisotropy.H"
#include "matelast_aniso_structuraltensor_strategy.H"
#include "lib_voigt_notation.H"



namespace
{
  class CoupAnisoExpoAnisotropyExtensionElementFiberTest
      : public ::testing::TestWithParam<std::tuple<int, int>>
  {
   public:
    CoupAnisoExpoAnisotropyExtensionElementFiberTest()
        : anisotropy_(), eleFibers_(2), eleTensors_(2), eleTensors_stress_(2)
    {
      /// initialize dummy fibers
      // Element fibers
      eleFibers_[0](0) = 0.858753861115007;
      eleFibers_[0](1) = 0.449823451060242;
      eleFibers_[0](2) = 0.245358246032859;

      eleFibers_[1](0) = 0.103448275862069;
      eleFibers_[1](1) = 0.137931034482759;
      eleFibers_[1](2) = 0.068965517241379;
      for (std::size_t i = 0; i < 2; ++i)
      {
        eleTensors_[i].MultiplyNT(eleFibers_[i], eleFibers_[i]);
        UTILS::VOIGT::Stresses::MatrixToVector(eleTensors_[i], eleTensors_stress_[i]);
      }

      SetupAnisotropyExtension();
    }

    void SetupAnisotropyExtension()
    {
      int fiber_id = std::get<0>(GetParam());
      auto strategy = Teuchos::rcp(new MAT::ELASTIC::StructuralTensorStrategyStandard(nullptr));
      anisotropyExtension_ = std::make_unique<MAT::ELASTIC::CoupAnisoExpoAnisotropyExtension>(
          1, 0.0, false, strategy, fiber_id);
      anisotropyExtension_->RegisterNeededTensors(
          MAT::FiberAnisotropyExtension<1>::FIBER_VECTORS |
          MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR_STRESS |
          MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR);
      anisotropy_.RegisterAnisotropyExtension(*anisotropyExtension_);
      anisotropy_.SetNumberOfGaussPoints(2);

      // Setup element fibers
      anisotropy_.SetElementFibers(eleFibers_);
    }

    [[nodiscard]] int GetGaussPoint() const { return std::get<1>(GetParam()); }

    [[nodiscard]] int GetFiberId() const { return std::get<0>(GetParam()); }

    MAT::Anisotropy anisotropy_;
    std::unique_ptr<MAT::ELASTIC::CoupAnisoExpoAnisotropyExtension> anisotropyExtension_;

    std::vector<CORE::LINALG::Matrix<3, 1>> eleFibers_;
    std::vector<CORE::LINALG::Matrix<3, 3>> eleTensors_;
    std::vector<CORE::LINALG::Matrix<6, 1>> eleTensors_stress_;
  };

  TEST_P(CoupAnisoExpoAnisotropyExtensionElementFiberTest, GetScalarProduct)
  {
    EXPECT_NEAR(anisotropyExtension_->GetScalarProduct(GetGaussPoint()), 1.0, 1e-10);
  }

  TEST_P(CoupAnisoExpoAnisotropyExtensionElementFiberTest, GetFiber)
  {
    BACI_EXPECT_NEAR(
        anisotropyExtension_->GetFiber(GetGaussPoint()), eleFibers_.at(GetFiberId() - 1), 1e-10);
  }

  TEST_P(CoupAnisoExpoAnisotropyExtensionElementFiberTest, GetStructuralTensor)
  {
    BACI_EXPECT_NEAR(anisotropyExtension_->GetStructuralTensor(GetGaussPoint()),
        eleTensors_.at(GetFiberId() - 1), 1e-10);
  }

  TEST_P(CoupAnisoExpoAnisotropyExtensionElementFiberTest, GetStructuralTensorStress)
  {
    BACI_EXPECT_NEAR(anisotropyExtension_->GetStructuralTensor_stress(GetGaussPoint()),
        eleTensors_stress_.at(GetFiberId() - 1), 1e-10);
  }

  INSTANTIATE_TEST_SUITE_P(GaussPoints, CoupAnisoExpoAnisotropyExtensionElementFiberTest,
      ::testing::Combine(::testing::Values(1, 2), ::testing::Values(0, 1)));
}  // namespace
