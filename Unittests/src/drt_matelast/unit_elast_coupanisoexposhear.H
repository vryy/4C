/*----------------------------------------------------------------------*/
/*! \file

\brief Testcases for the CoupAnisoExpoShear summand

\level 2

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_UNIT_ELAST_COUPANISOEXPOSHEAR_H
#define BACI_UNIT_ELAST_COUPANISOEXPOSHEAR_H

#include <cxxtest/TestSuite.h>
#include <Teuchos_RCPDecl.hpp>
#include "src/linalg/linalg_fixedsizematrix.H"
#include "src/drt_matelast/elast_coupanisoexposhear.H"
#include "src/drt_mat/anisotropy.H"
#include "src/drt_matelast/elast_aniso_structuraltensor_strategy.H"

namespace MAT
{
  namespace ELASTIC
  {
    class CoupAnisoExpoShearAnisotropyExtension_TestSuite;
  }  // namespace ELASTIC
}  // namespace MAT

// class implementation
class MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension_TestSuite : public CxxTest::TestSuite
{
 public:
  /// Setup Testsuite
  void setUp() override
  {
    SetupFibers();
    SetupStructuralTensors();
    SetupStructuralTensors_stress();
    SetupScalarProducts();

    anisotropy_ = Teuchos::rcp(new MAT::Anisotropy());
    anisotropy_->numgp_ = 2;
    anisotropy_->element_fibers_initialized_ = true;
    anisotropy_->gp_fibers_initialized_ = true;
    anisotropy_->elementCylinderCoordinateSystemManager_ = boost::none;

    anisotropy_->elementFibers_ = *fibersele_;
    anisotropy_->gpFibers_.emplace_back(*fibersgp0_);
    anisotropy_->gpFibers_.emplace_back(*fibersgp1_);
  }

  /// Tear down test suite
  void tearDown() override {}

  void SetupAnisotropyExtension(int init_mode, std::array<int, 2> fiber_ids)
  {
    anisotropyExtension_ =
        Teuchos::rcp(new MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension(init_mode, fiber_ids));

    anisotropy_->extensions_.clear();
    anisotropy_->RegisterAnisotropyExtension(*anisotropyExtension_);

    anisotropy_->OnElementFibersInitialized();
    anisotropy_->OnGPFibersInitialized();
  }

  void test_GetScalarProduct()
  {
    // element fibers
    SetupAnisotropyExtension(1, {0, 1});
    TS_ASSERT_DELTA(anisotropyExtension_->GetScalarProduct(0), scalarProductsele_->at(0), 1e-10);

    SetupAnisotropyExtension(1, {0, 2});
    TS_ASSERT_DELTA(anisotropyExtension_->GetScalarProduct(0), scalarProductsele_->at(1), 1e-10);

    SetupAnisotropyExtension(1, {1, 2});
    TS_ASSERT_DELTA(anisotropyExtension_->GetScalarProduct(0), scalarProductsele_->at(2), 1e-10);

    SetupAnisotropyExtension(1, {0, 1});
    TS_ASSERT_DELTA(anisotropyExtension_->GetScalarProduct(1), scalarProductsele_->at(0), 1e-10);

    SetupAnisotropyExtension(1, {0, 2});
    TS_ASSERT_DELTA(anisotropyExtension_->GetScalarProduct(1), scalarProductsele_->at(1), 1e-10);

    SetupAnisotropyExtension(1, {1, 2});
    TS_ASSERT_DELTA(anisotropyExtension_->GetScalarProduct(1), scalarProductsele_->at(2), 1e-10);

    // Gauss point fibers
    SetupAnisotropyExtension(3, {0, 1});
    TS_ASSERT_DELTA(anisotropyExtension_->GetScalarProduct(0), scalarProductsgp0_->at(0), 1e-10);

    SetupAnisotropyExtension(3, {0, 2});
    TS_ASSERT_DELTA(anisotropyExtension_->GetScalarProduct(0), scalarProductsgp0_->at(1), 1e-10);

    SetupAnisotropyExtension(3, {1, 2});
    TS_ASSERT_DELTA(anisotropyExtension_->GetScalarProduct(0), scalarProductsgp0_->at(2), 1e-10);

    SetupAnisotropyExtension(3, {0, 1});
    TS_ASSERT_DELTA(anisotropyExtension_->GetScalarProduct(1), scalarProductsgp1_->at(0), 1e-10);

    SetupAnisotropyExtension(3, {0, 2});
    TS_ASSERT_DELTA(anisotropyExtension_->GetScalarProduct(1), scalarProductsgp1_->at(1), 1e-10);

    SetupAnisotropyExtension(3, {1, 2});
    TS_ASSERT_DELTA(anisotropyExtension_->GetScalarProduct(1), scalarProductsgp1_->at(2), 1e-10);
  }

  void test_GetStructuralTensor()
  {
    // element fibers
    SetupAnisotropyExtension(1, {0, 1});
    AssertMatrixEqual(anisotropyExtension_->GetStructuralTensor(0), tensorsele_->at(0), 1e-10);

    SetupAnisotropyExtension(1, {0, 2});
    AssertMatrixEqual(anisotropyExtension_->GetStructuralTensor(0), tensorsele_->at(1), 1e-10);

    SetupAnisotropyExtension(1, {1, 2});
    AssertMatrixEqual(anisotropyExtension_->GetStructuralTensor(0), tensorsele_->at(2), 1e-10);

    SetupAnisotropyExtension(1, {0, 1});
    AssertMatrixEqual(anisotropyExtension_->GetStructuralTensor(1), tensorsele_->at(0), 1e-10);

    SetupAnisotropyExtension(1, {0, 2});
    AssertMatrixEqual(anisotropyExtension_->GetStructuralTensor(1), tensorsele_->at(1), 1e-10);

    SetupAnisotropyExtension(1, {1, 2});
    AssertMatrixEqual(anisotropyExtension_->GetStructuralTensor(1), tensorsele_->at(2), 1e-10);

    // Gauss point fibers
    SetupAnisotropyExtension(3, {0, 1});
    AssertMatrixEqual(anisotropyExtension_->GetStructuralTensor(0), tensorsgp0_->at(0), 1e-10);

    SetupAnisotropyExtension(3, {0, 2});
    AssertMatrixEqual(anisotropyExtension_->GetStructuralTensor(0), tensorsgp0_->at(1), 1e-10);

    SetupAnisotropyExtension(3, {1, 2});
    AssertMatrixEqual(anisotropyExtension_->GetStructuralTensor(0), tensorsgp0_->at(2), 1e-10);

    SetupAnisotropyExtension(3, {0, 1});
    AssertMatrixEqual(anisotropyExtension_->GetStructuralTensor(1), tensorsgp1_->at(0), 1e-10);

    SetupAnisotropyExtension(3, {0, 2});
    AssertMatrixEqual(anisotropyExtension_->GetStructuralTensor(1), tensorsgp1_->at(1), 1e-10);

    SetupAnisotropyExtension(3, {1, 2});
    AssertMatrixEqual(anisotropyExtension_->GetStructuralTensor(1), tensorsgp1_->at(2), 1e-10);
  }

  void test_StructuralTensors_stress() {}

 private:
  template <unsigned int N, unsigned int M>
  void AssertMatrixEqual(
      const LINALG::Matrix<N, M>& mat1, const LINALG::Matrix<N, M>& mat2, double delta = 1e-10)
  {
    for (unsigned int i = 0; i < N; ++i)
    {
      for (unsigned int j = 0; j < M; ++j)
      {
        TS_ASSERT_DELTA(mat1(i, j), mat2(i, j), delta);
      }
    }
  }

 private:
  void SetupFibers()
  {
    fibersgp0_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 1>>(3));
    fibersgp1_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 1>>(3));
    fibersele_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 1>>(3));

    /// initialize dummy fibers
    // gp 0
    fibersgp0_->at(0)(0) = 0.469809238649817;
    fibersgp0_->at(0)(1) = 0.872502871778232;
    fibersgp0_->at(0)(2) = 0.134231211042805;

    fibersgp0_->at(1)(0) = 0.071428571428571;
    fibersgp0_->at(1)(1) = 0.142857142857143;
    fibersgp0_->at(1)(2) = 0.214285714285714;

    fibersgp0_->at(2)(0) = 0.142857142857143;
    fibersgp0_->at(2)(1) = 0.214285714285714;
    fibersgp0_->at(2)(2) = 0.071428571428571;

    // gp 1
    fibersgp1_->at(0)(0) = 0.245358246032859;
    fibersgp1_->at(0)(1) = 0.858753861115007;
    fibersgp1_->at(0)(2) = 0.449823451060242;

    fibersgp1_->at(1)(0) = 0.068965517241379;
    fibersgp1_->at(1)(1) = 0.103448275862069;
    fibersgp1_->at(1)(2) = 0.137931034482759;

    fibersgp1_->at(2)(0) = 0.137931034482759;
    fibersgp1_->at(2)(1) = 0.068965517241379;
    fibersgp1_->at(2)(2) = 0.103448275862069;

    // Element fibers
    fibersele_->at(0)(0) = 0.858753861115007;
    fibersele_->at(0)(1) = 0.449823451060242;
    fibersele_->at(0)(2) = 0.245358246032859;

    fibersele_->at(1)(0) = 0.103448275862069;
    fibersele_->at(1)(1) = 0.137931034482759;
    fibersele_->at(1)(2) = 0.068965517241379;

    fibersele_->at(2)(0) = 0.872502871778232;
    fibersele_->at(2)(1) = 0.134231211042805;
    fibersele_->at(2)(2) = 0.469809238649817;
  }

  void SetupStructuralTensors()
  {
    tensorsgp0_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>(3));
    tensorsgp1_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>(3));
    tensorsele_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>(3));

    // 1 - 2
    SetupSingleStructuralTensor(fibersgp0_->at(0), fibersgp0_->at(1), tensorsgp0_->at(0));
    SetupSingleStructuralTensor(fibersgp1_->at(0), fibersgp1_->at(1), tensorsgp1_->at(0));
    SetupSingleStructuralTensor(fibersele_->at(0), fibersele_->at(1), tensorsele_->at(0));

    // 1 - 3
    SetupSingleStructuralTensor(fibersgp0_->at(0), fibersgp0_->at(2), tensorsgp0_->at(1));
    SetupSingleStructuralTensor(fibersgp1_->at(0), fibersgp1_->at(2), tensorsgp1_->at(1));
    SetupSingleStructuralTensor(fibersele_->at(0), fibersele_->at(2), tensorsele_->at(1));

    // 2 - 3
    SetupSingleStructuralTensor(fibersgp0_->at(1), fibersgp0_->at(2), tensorsgp0_->at(2));
    SetupSingleStructuralTensor(fibersgp1_->at(1), fibersgp1_->at(2), tensorsgp1_->at(2));
    SetupSingleStructuralTensor(fibersele_->at(1), fibersele_->at(2), tensorsele_->at(2));
  }

  void SetupSingleStructuralTensor(const LINALG::Matrix<3, 1>& fiber1,
      const LINALG::Matrix<3, 1>& fiber2, LINALG::Matrix<3, 3>& structuralTensor)
  {
    LINALG::Matrix<3, 3> fiber1fiber2T(false);

    fiber1fiber2T.MultiplyNT(fiber1, fiber2);

    structuralTensor.Update(0.5, fiber1fiber2T);
    structuralTensor.UpdateT(0.5, fiber1fiber2T, 1.0);
  }

  void SetupScalarProducts()
  {
    scalarProductsgp0_ = Teuchos::rcp(new std::vector<double>(3));
    scalarProductsgp1_ = Teuchos::rcp(new std::vector<double>(3));
    scalarProductsele_ = Teuchos::rcp(new std::vector<double>(3));

    // 1 - 2
    scalarProductsgp0_->at(0) = fibersgp0_->at(0).Dot(fibersgp0_->at(1));
    scalarProductsgp1_->at(0) = fibersgp1_->at(0).Dot(fibersgp1_->at(1));
    scalarProductsele_->at(0) = fibersele_->at(0).Dot(fibersele_->at(1));

    // 1 - 3
    scalarProductsgp0_->at(1) = fibersgp0_->at(0).Dot(fibersgp0_->at(2));
    scalarProductsgp1_->at(1) = fibersgp1_->at(0).Dot(fibersgp1_->at(2));
    scalarProductsele_->at(1) = fibersele_->at(0).Dot(fibersele_->at(2));

    // 2 - 3
    scalarProductsgp0_->at(2) = fibersgp0_->at(1).Dot(fibersgp0_->at(2));
    scalarProductsgp1_->at(2) = fibersgp1_->at(1).Dot(fibersgp1_->at(2));
    scalarProductsele_->at(2) = fibersele_->at(1).Dot(fibersele_->at(2));
  }

  void SetupStructuralTensors_stress()
  {
    tensorsgp0_stress_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>(3));
    tensorsgp1_stress_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>(3));
    tensorsele_stress_ = Teuchos::rcp(new std::vector<LINALG::Matrix<6, 1>>(3));

    for (int i = 0; i < 3; ++i)
    {
      UTILS::VOIGT::Stresses::MatrixToVector(tensorsgp0_->at(i), tensorsgp0_stress_->at(i));
      UTILS::VOIGT::Stresses::MatrixToVector(tensorsgp1_->at(i), tensorsgp1_stress_->at(i));
      UTILS::VOIGT::Stresses::MatrixToVector(tensorsele_->at(i), tensorsele_stress_->at(i));
    }
  }

  Teuchos::RCP<MAT::Anisotropy> anisotropy_;
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension> anisotropyExtension_;

  Teuchos::RCP<std::vector<LINALG::Matrix<3, 1>>> fibersgp0_;
  Teuchos::RCP<std::vector<LINALG::Matrix<3, 1>>> fibersgp1_;
  Teuchos::RCP<std::vector<LINALG::Matrix<3, 1>>> fibersele_;

  Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> tensorsgp0_;
  Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> tensorsgp1_;
  Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> tensorsele_;

  Teuchos::RCP<std::vector<LINALG::Matrix<6, 1>>> tensorsgp0_stress_;
  Teuchos::RCP<std::vector<LINALG::Matrix<6, 1>>> tensorsgp1_stress_;
  Teuchos::RCP<std::vector<LINALG::Matrix<6, 1>>> tensorsele_stress_;

  Teuchos::RCP<std::vector<double>> scalarProductsgp0_;
  Teuchos::RCP<std::vector<double>> scalarProductsgp1_;
  Teuchos::RCP<std::vector<double>> scalarProductsele_;
};

#endif  // BACI_UNIT_ELAST_COUPANISOEXPOSHEAR_H
