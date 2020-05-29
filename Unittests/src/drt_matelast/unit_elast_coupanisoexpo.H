/*----------------------------------------------------------------------*/
/*! \file

\brief Testcases for the CoupAnisoExpoBase summand

\level 2

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_UNIT_ELAST_COUPANISOEXPO_H
#define BACI_UNIT_ELAST_COUPANISOEXPO_H

#include <cxxtest/TestSuite.h>
#include <Teuchos_RCPDecl.hpp>
#include "src/linalg/linalg_fixedsizematrix.H"
#include "src/drt_matelast/elast_coupanisoexpo.H"
#include "src/drt_mat/anisotropy.H"
#include "src/drt_matelast/elast_aniso_structuraltensor_strategy.H"

#include "../common/special_assertions.H"

namespace MAT
{
  namespace ELASTIC
  {
    class CoupAnisoExpoAnisotropyExtension_TestSuite;
  }  // namespace ELASTIC
}  // namespace MAT

// class implementation
class MAT::ELASTIC::CoupAnisoExpoAnisotropyExtension_TestSuite : public CxxTest::TestSuite
{
 public:
  /// Setup Testsuite
  void setUp() override
  {
    SetupFibers();
    SetupStructuralTensors();
    SetupStructuralTensors_stress();
  }

  /// Tear down test suite
  void tearDown() override {}

  void SetupAnisotropyExtension(int init_mode, int fiber_id)
  {
    auto strategy = Teuchos::rcp(new StructuralTensorStrategyStandard(nullptr));
    anisotropyExtension_ = Teuchos::rcp(new MAT::ELASTIC::CoupAnisoExpoAnisotropyExtension(
        init_mode, 0.0, false, strategy, fiber_id));
    anisotropyExtension_->RegisterNeededTensors(
        FiberAnisotropyExtension<1>::FIBER_VECTORS |
        FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR_STRESS |
        FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR);

    anisotropy_ = Teuchos::rcp(new MAT::Anisotropy());
    anisotropy_->RegisterAnisotropyExtension(*anisotropyExtension_);
    anisotropy_->SetNumberOfGaussPoints(2);

    // Setup element fibers
    anisotropy_->SetElementFibers(*fibersele_);

    // Setup Gauss point fibers
    std::vector<std::vector<LINALG::Matrix<3, 1>>> gpfibers(0);
    gpfibers.emplace_back(*fibersgp0_);
    gpfibers.emplace_back(*fibersgp1_);
    anisotropy_->SetGaussPointFibers(gpfibers);
  }

  void test_GetScalarProduct()
  {
    SetupAnisotropyExtension(1, 1);
    TS_ASSERT_DELTA(anisotropyExtension_->GetScalarProduct(0), 1.0, 1e-10);

    SetupAnisotropyExtension(1, 2);
    TS_ASSERT_DELTA(anisotropyExtension_->GetScalarProduct(0), 1.0, 1e-10);

    SetupAnisotropyExtension(3, 1);
    TS_ASSERT_DELTA(anisotropyExtension_->GetScalarProduct(0), 1.0, 1e-10);

    SetupAnisotropyExtension(3, 2);
    TS_ASSERT_DELTA(anisotropyExtension_->GetScalarProduct(0), 1.0, 1e-10);
  }

  void test_GetFiber()
  {
    // Test element fibers
    SetupAnisotropyExtension(1, 1);
    TESTING::AssertDelta(anisotropyExtension_->GetFiber(0), fibersele_->at(0), 1e-10);

    SetupAnisotropyExtension(1, 2);
    TESTING::AssertDelta(anisotropyExtension_->GetFiber(0), fibersele_->at(1), 1e-10);

    SetupAnisotropyExtension(1, 1);
    TESTING::AssertDelta(anisotropyExtension_->GetFiber(1), fibersele_->at(0), 1e-10);

    SetupAnisotropyExtension(1, 2);
    TESTING::AssertDelta(anisotropyExtension_->GetFiber(1), fibersele_->at(1), 1e-10);

    // Test GP fibers
    SetupAnisotropyExtension(3, 1);
    TESTING::AssertDelta(anisotropyExtension_->GetFiber(0), fibersgp0_->at(0), 1e-10);

    SetupAnisotropyExtension(3, 2);
    TESTING::AssertDelta(anisotropyExtension_->GetFiber(0), fibersgp0_->at(1), 1e-10);

    SetupAnisotropyExtension(3, 1);
    TESTING::AssertDelta(anisotropyExtension_->GetFiber(1), fibersgp1_->at(0), 1e-10);

    SetupAnisotropyExtension(3, 2);
    TESTING::AssertDelta(anisotropyExtension_->GetFiber(1), fibersgp1_->at(1), 1e-10);
  }

  void test_GetStructuralTensor()
  {
    // Test element fibers
    SetupAnisotropyExtension(1, 1);
    TESTING::AssertDelta(anisotropyExtension_->GetStructuralTensor(0), tensorsele_->at(0), 1e-10);

    SetupAnisotropyExtension(1, 2);
    TESTING::AssertDelta(anisotropyExtension_->GetStructuralTensor(0), tensorsele_->at(1), 1e-10);

    SetupAnisotropyExtension(1, 1);
    TESTING::AssertDelta(anisotropyExtension_->GetStructuralTensor(1), tensorsele_->at(0), 1e-10);

    SetupAnisotropyExtension(1, 2);
    TESTING::AssertDelta(anisotropyExtension_->GetStructuralTensor(1), tensorsele_->at(1), 1e-10);

    // Test GP fibers
    SetupAnisotropyExtension(3, 1);
    TESTING::AssertDelta(anisotropyExtension_->GetStructuralTensor(0), tensorsgp0_->at(0), 1e-10);

    SetupAnisotropyExtension(3, 2);
    TESTING::AssertDelta(anisotropyExtension_->GetStructuralTensor(0), tensorsgp0_->at(1), 1e-10);

    SetupAnisotropyExtension(3, 1);
    TESTING::AssertDelta(anisotropyExtension_->GetStructuralTensor(1), tensorsgp1_->at(0), 1e-10);

    SetupAnisotropyExtension(3, 2);
    TESTING::AssertDelta(anisotropyExtension_->GetStructuralTensor(1), tensorsgp1_->at(1), 1e-10);
  }

  void test_StructuralTensors_stress()
  {
    // Test element fibers
    SetupAnisotropyExtension(1, 1);
    TESTING::AssertDelta(
        anisotropyExtension_->GetStructuralTensor_stress(0), tensorsele_stress_->at(0), 1e-10);

    SetupAnisotropyExtension(1, 2);
    TESTING::AssertDelta(
        anisotropyExtension_->GetStructuralTensor_stress(0), tensorsele_stress_->at(1), 1e-10);

    SetupAnisotropyExtension(1, 1);
    TESTING::AssertDelta(
        anisotropyExtension_->GetStructuralTensor_stress(1), tensorsele_stress_->at(0), 1e-10);

    SetupAnisotropyExtension(1, 2);
    TESTING::AssertDelta(
        anisotropyExtension_->GetStructuralTensor_stress(1), tensorsele_stress_->at(1), 1e-10);

    // Test GP fibers
    SetupAnisotropyExtension(3, 1);
    TESTING::AssertDelta(
        anisotropyExtension_->GetStructuralTensor_stress(0), tensorsgp0_stress_->at(0), 1e-10);

    SetupAnisotropyExtension(3, 2);
    TESTING::AssertDelta(
        anisotropyExtension_->GetStructuralTensor_stress(0), tensorsgp0_stress_->at(1), 1e-10);

    SetupAnisotropyExtension(3, 1);
    TESTING::AssertDelta(
        anisotropyExtension_->GetStructuralTensor_stress(1), tensorsgp1_stress_->at(0), 1e-10);

    SetupAnisotropyExtension(3, 2);
    TESTING::AssertDelta(
        anisotropyExtension_->GetStructuralTensor_stress(1), tensorsgp1_stress_->at(1), 1e-10);
  }

 private:
  void SetupFibers()
  {
    fibersgp0_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 1>>(2));
    fibersgp1_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 1>>(2));
    fibersele_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 1>>(2));

    /// initialize dummy fibers
    // gp 0
    fibersgp0_->at(0)(0) = 0.469809238649817;
    fibersgp0_->at(0)(1) = 0.872502871778232;
    fibersgp0_->at(0)(2) = 0.134231211042805;

    fibersgp0_->at(0)(0) = 0.071428571428571;
    fibersgp0_->at(0)(1) = 0.142857142857143;
    fibersgp0_->at(0)(2) = 0.214285714285714;

    // gp 1
    fibersgp1_->at(1)(0) = 0.245358246032859;
    fibersgp1_->at(1)(1) = 0.858753861115007;
    fibersgp1_->at(1)(2) = 0.449823451060242;

    fibersgp1_->at(1)(0) = 0.068965517241379;
    fibersgp1_->at(1)(1) = 0.103448275862069;
    fibersgp1_->at(1)(2) = 0.137931034482759;

    // Element fibers
    fibersele_->at(0)(0) = 0.858753861115007;
    fibersele_->at(0)(1) = 0.449823451060242;
    fibersele_->at(0)(2) = 0.245358246032859;

    fibersele_->at(1)(0) = 0.103448275862069;
    fibersele_->at(1)(1) = 0.137931034482759;
    fibersele_->at(1)(2) = 0.068965517241379;
  }

  void SetupStructuralTensors()
  {
    tensorsgp0_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>(2));
    tensorsgp1_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>(2));
    tensorsele_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>(2));

    for (int i = 0; i < 2; ++i)
    {
      tensorsgp0_->at(i).MultiplyNT(fibersgp0_->at(i), fibersgp0_->at(i));
      tensorsgp1_->at(i).MultiplyNT(fibersgp1_->at(i), fibersgp1_->at(i));
      tensorsele_->at(i).MultiplyNT(fibersele_->at(i), fibersele_->at(i));
    }
  }
  void SetupStructuralTensors_stress()
  {
    tensorsgp0_stress_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>(2));
    tensorsgp1_stress_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>(2));
    tensorsele_stress_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>(2));

    for (int i = 0; i < 2; ++i)
    {
      UTILS::VOIGT::Stresses::MatrixToVector(tensorsgp0_->at(i), tensorsgp0_stress_->at(i));
      UTILS::VOIGT::Stresses::MatrixToVector(tensorsgp1_->at(i), tensorsgp1_stress_->at(i));
      UTILS::VOIGT::Stresses::MatrixToVector(tensorsele_->at(i), tensorsele_stress_->at(i));
    }
  }

  Teuchos::RCP<MAT::Anisotropy> anisotropy_;
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpoAnisotropyExtension> anisotropyExtension_;

  Teuchos::RCP<std::vector<LINALG::Matrix<3, 1>>> fibersgp0_;
  Teuchos::RCP<std::vector<LINALG::Matrix<3, 1>>> fibersgp1_;
  Teuchos::RCP<std::vector<LINALG::Matrix<3, 1>>> fibersele_;
  Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> tensorsgp0_;
  Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> tensorsgp1_;
  Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> tensorsele_;
  Teuchos::RCP<std::vector<LINALG::Matrix<6, 1>>> tensorsgp0_stress_;
  Teuchos::RCP<std::vector<LINALG::Matrix<6, 1>>> tensorsgp1_stress_;
  Teuchos::RCP<std::vector<LINALG::Matrix<6, 1>>> tensorsele_stress_;
};

#endif  // BACI_UNIT_ELAST_COUPANISOEXPO_H
