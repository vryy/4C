/*----------------------------------------------------------------------*/
/*! \file

\brief Testcases for the material_service functions

\level 2

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_UNIT_MATERIAL_SERVICE_H
#define BACI_UNIT_MATERIAL_SERVICE_H

#include <cxxtest/TestSuite.h>
#include <src/linalg/linalg_fixedsizematrix.H>
#include <Teuchos_RCPDecl.hpp>
#include "src/drt_mat/material_service.H"

namespace MAT
{
  class MaterialService_TestSuite;
}

// class implementation
class MAT::MaterialService_TestSuite : public CxxTest::TestSuite
{
 public:
  /// Setup Testsuite
  void setUp() override
  {
    tens = Teuchos::rcp(new LINALG::Matrix<3, 3>);
    itens = Teuchos::rcp(new LINALG::Matrix<3, 3>);
    tens_strain = Teuchos::rcp(new LINALG::Matrix<6, 1>);
    itens_strain = Teuchos::rcp(new LINALG::Matrix<6, 1>);
    tens_stress = Teuchos::rcp(new LINALG::Matrix<6, 1>);
    itens_stress = Teuchos::rcp(new LINALG::Matrix<6, 1>);

    (*tens)(0, 0) = 1.1;
    (*tens)(1, 1) = 1.2;
    (*tens)(2, 2) = 1.3;
    (*tens)(0, 1) = (*tens)(1, 0) = 0.01;
    (*tens)(1, 2) = (*tens)(2, 1) = 0.02;
    (*tens)(0, 2) = (*tens)(2, 0) = 0.03;

    (*itens)(0, 0) = 0.90972618;
    (*itens)(1, 1) = 0.83360457;
    (*itens)(2, 2) = 0.76990741;
    (*itens)(0, 1) = (*itens)(1, 0) = -0.00723301;
    (*itens)(1, 2) = (*itens)(2, 1) = -0.01265777;
    (*itens)(0, 2) = (*itens)(2, 0) = -0.0208824;

    (*tens_strain)(0) = (*tens_stress)(0) = 1.1;
    (*tens_strain)(1) = (*tens_stress)(1) = 1.2;
    (*tens_strain)(2) = (*tens_stress)(2) = 1.3;

    (*tens_stress)(3) = 0.01;
    (*tens_stress)(4) = 0.02;
    (*tens_stress)(5) = 0.03;

    (*itens_stress)(0) = 0.90972618;
    (*itens_stress)(1) = 0.83360457;
    (*itens_stress)(2) = 0.76990741;
    (*itens_stress)(3) = -0.00723301;
    (*itens_stress)(4) = -0.01265777;
    (*itens_stress)(5) = -0.0208824;

    (*tens_strain)(3) = 2 * 0.01;
    (*tens_strain)(4) = 2 * 0.02;
    (*tens_strain)(5) = 2 * 0.03;

    (*itens_strain)(0) = 0.90972618;
    (*itens_strain)(1) = 0.83360457;
    (*itens_strain)(2) = 0.76990741;
    (*itens_strain)(3) = 2 * -0.00723301;
    (*itens_strain)(4) = 2 * -0.01265777;
    (*itens_strain)(5) = 2 * -0.0208824;
  }

  /// Tear down test suite
  void tearDown() override {}

  void test_conversion_stress()
  {
    LINALG::Matrix<6, 1> cmp_stress(false);

    MAT::MatrixtoStressLikeVoigtNotation(*tens, cmp_stress);

    for (int i = 0; i < 6; ++i) TS_ASSERT_DELTA(cmp_stress(i), (*tens_stress)(i), 1e-10);
  }

  void test_conversion_strain()
  {
    LINALG::Matrix<6, 1> cmp_strain(false);

    MAT::MatrixtoStrainLikeVoigtNotation(*tens, cmp_strain);

    for (int i = 0; i < 6; ++i) TS_ASSERT_DELTA(cmp_strain(i), (*tens_strain)(i), 1e-10);
  }

  void test_determinant_stress()
  {
    TS_ASSERT_DELTA(MAT::VStressUtils::Determinant(*tens_stress), 1.7143620000000002, 1e-10);
  }

  void test_determinant_strain()
  {
    TS_ASSERT_DELTA(MAT::VStrainUtils::Determinant(*tens_strain), 1.7143620000000002, 1e-10);
  }

  void test_invariants_principal()
  {
    LINALG::Matrix<3, 1> prinv(false);
    MAT::InvariantsPrincipal(prinv, *tens);
    TS_ASSERT_DELTA(prinv(0), 3.5999999999999996, 1e-10);
    TS_ASSERT_DELTA(prinv(1), 4.3085999999999984, 1e-10);
    TS_ASSERT_DELTA(prinv(2), 1.7143620000000002, 1e-10);
  }

  void test_invariants_principal_stress()
  {
    LINALG::Matrix<3, 1> prinv(false);
    MAT::InvariantsPrincipal<MAT::VoigtNotation::stress>(prinv, *tens_stress);
    TS_ASSERT_DELTA(prinv(0), 3.5999999999999996, 1e-10);
    TS_ASSERT_DELTA(prinv(1), 4.3085999999999984, 1e-10);
    TS_ASSERT_DELTA(prinv(2), 1.7143620000000002, 1e-10);
  }

  void test_invariants_principal_strain()
  {
    LINALG::Matrix<3, 1> prinv(false);
    MAT::InvariantsPrincipal<MAT::VoigtNotation::strain>(prinv, *tens_strain);
    TS_ASSERT_DELTA(prinv(0), 3.5999999999999996, 1e-10);
    TS_ASSERT_DELTA(prinv(1), 4.3085999999999984, 1e-10);
    TS_ASSERT_DELTA(prinv(2), 1.7143620000000002, 1e-10);
  }

  void test_inverse_stress()
  {
    LINALG::Matrix<6, 1> itens_stress_result(false);
    VStressUtils::InverseTensor(*tens_stress, itens_stress_result);

    for (unsigned i = 0; i < 6; ++i)
      TS_ASSERT_DELTA(itens_stress_result(i), (*itens_stress)(i), 1e-5);
  }

  void test_inverse_strain()
  {
    LINALG::Matrix<6, 1> itens_strain_result(false);
    VStrainUtils::InverseTensor(*tens_strain, itens_strain_result);

    for (unsigned i = 0; i < 6; ++i)
      TS_ASSERT_DELTA(itens_strain_result(i), (*itens_strain)(i), 1e-5);
  }

  void test_to_stress_like()
  {
    LINALG::Matrix<6, 1> strain_to_stress(false);
    LINALG::Matrix<6, 1> stress_to_stress(false);
    VStrainUtils::ToStressLike(*tens_strain, strain_to_stress);
    VStressUtils::ToStressLike(*tens_stress, stress_to_stress);

    for (unsigned i = 0; i < 6; ++i)
    {
      TS_ASSERT_DELTA(strain_to_stress(i), (*tens_stress)(i), 1e-5);
      TS_ASSERT_DELTA(stress_to_stress(i), (*tens_stress)(i), 1e-5);
    }
  }

  void test_to_strain_like()
  {
    LINALG::Matrix<6, 1> strain_to_strain(false);
    LINALG::Matrix<6, 1> stress_to_strain(false);
    VStrainUtils::ToStrainLike(*tens_strain, strain_to_strain);
    VStressUtils::ToStrainLike(*tens_stress, stress_to_strain);

    for (unsigned i = 0; i < 6; ++i)
    {
      TS_ASSERT_DELTA(strain_to_strain(i), (*tens_strain)(i), 1e-5);
      TS_ASSERT_DELTA(stress_to_strain(i), (*tens_strain)(i), 1e-5);
    }
  }

  void test_identity_matrix()
  {
    LINALG::Matrix<6, 1> id(false);

    IdentityMatrix(id);

    TS_ASSERT_DELTA(id(0), 1.0, 1e-10);
    TS_ASSERT_DELTA(id(1), 1.0, 1e-10);
    TS_ASSERT_DELTA(id(2), 1.0, 1e-10);
    TS_ASSERT_DELTA(id(3), 0.0, 1e-10);
    TS_ASSERT_DELTA(id(4), 0.0, 1e-10);
    TS_ASSERT_DELTA(id(5), 0.0, 1e-10);
  }

 private:
  Teuchos::RCP<LINALG::Matrix<3, 3>> tens;
  Teuchos::RCP<LINALG::Matrix<3, 3>> itens;
  Teuchos::RCP<LINALG::Matrix<6, 1>> tens_strain;
  Teuchos::RCP<LINALG::Matrix<6, 1>> itens_strain;
  Teuchos::RCP<LINALG::Matrix<6, 1>> tens_stress;
  Teuchos::RCP<LINALG::Matrix<6, 1>> itens_stress;
};

#endif  // BACI_UNIT_MATERIAL_SERVICE_H
