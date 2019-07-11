/*----------------------------------------------------------------------*/
/*!

\brief Testcases for the ElastHyper service functions

\level 2

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_UNIT_ELASTHYPER_SERVICE_H
#define BACI_UNIT_ELASTHYPER_SERVICE_H



#include <cxxtest/TestSuite.h>
#include <src/linalg/linalg_fixedsizematrix.H>
#include <Teuchos_RCPDecl.hpp>
#include <src/drt_mat/elasthyper_service.H>
#include <src/drt_matelast/elast_isoneohooke.H>
#include <src/drt_matelast/elast_coupneohooke.H>
#include <src/drt_matelast/elast_coupanisoexpo.H>
#include <src/drt_matelast/elast_coupanisoneohooke.H>
namespace MAT
{
  class ElastHyperService_TestSuite;
}

// class implementation
class MAT::ElastHyperService_TestSuite : public CxxTest::TestSuite
{
 public:
  /// Setup Testsuite
  void setUp() override
  {
    E_strain_ = Teuchos::rcp(new LINALG::Matrix<6, 1>);

    // create dummy Right Cauchy-Green tensor for testing of the other methods
    (*E_strain_)(0, 0) = 0.1;
    (*E_strain_)(1, 0) = 0.2;
    (*E_strain_)(2, 0) = 0.3;
    (*E_strain_)(3, 0) = 0.01;
    (*E_strain_)(4, 0) = 0.02;
    (*E_strain_)(5, 0) = 0.03;

    // Create potsum
    // ToDo: Here we need more summands to actually test everything. Unfortunately this is not
    //  easy since for anisotropic materials, the global problem is needed, to create the
    //  structural Tensor strategy.
    potsum_.clear();
    Teuchos::RCP<MAT::PAR::Material> isoNeoHookeParams = Teuchos::rcp(new MAT::PAR::Material());
    isoNeoHookeParams->Add("MUE", 1.3);
    isoNeoHookeParams->SetParameter(new MAT::ELASTIC::PAR::IsoNeoHooke(isoNeoHookeParams));
    auto* isoNeoHookeParams2 =
        dynamic_cast<MAT::ELASTIC::PAR::IsoNeoHooke*>(isoNeoHookeParams->Parameter());
    potsum_.emplace_back(Teuchos::rcp(new MAT::ELASTIC::IsoNeoHooke(isoNeoHookeParams2)));

    MAT::ElastHyperProperties(potsum_, properties);
  }

  void test_EvaluateRightCauchyGreenStrainLikeVoigt()
  {
    LINALG::Matrix<6, 1> C_strain;

    MAT::EvaluateRightCauchyGreenStrainLikeVoigt(*E_strain_, C_strain);

    TS_ASSERT_DELTA(C_strain(0), 1.2, 1e-10)
    TS_ASSERT_DELTA(C_strain(1), 1.4, 1e-10)
    TS_ASSERT_DELTA(C_strain(2), 1.6, 1e-10)
    TS_ASSERT_DELTA(C_strain(3), 0.02, 1e-10)
    TS_ASSERT_DELTA(C_strain(4), 0.04, 1e-10)
    TS_ASSERT_DELTA(C_strain(5), 0.06, 1e-10)
  }

  void test_EvaluateInvariantDerivatives()
  {
    LINALG::Matrix<3, 1> prinv_new(true);
    LINALG::Matrix<3, 1> dPI_new(true);
    LINALG::Matrix<6, 1> ddPII_new(true);

    prinv_new(0) = 3.1;
    prinv_new(1) = 3.2;
    prinv_new(2) = 1.05;

    ElastHyperEvaluateInvariantDerivatives(prinv_new, dPI_new, ddPII_new, potsum_, properties, 0);

    TS_ASSERT_DELTA(dPI_new(0), 0.639514, 1e-4);
    TS_ASSERT_DELTA(dPI_new(1), 0., 1e-4);
    TS_ASSERT_DELTA(dPI_new(2), -0.629363, 1e-4);
    TS_ASSERT_DELTA(ddPII_new(0), 0., 1e-4);
    TS_ASSERT_DELTA(ddPII_new(1), 0., 1e-4);
    TS_ASSERT_DELTA(ddPII_new(2), 0.799191, 1e-4);
    TS_ASSERT_DELTA(ddPII_new(3), 0., 1e-4);
    TS_ASSERT_DELTA(ddPII_new(4), -0.20302, 1e-4);
    TS_ASSERT_DELTA(ddPII_new(5), 0., 1e-4);
  }

  void test_Comparison()
  {
    const int eleGID = 0;

    LINALG::Matrix<6, 1> id2_new(false);
    id2_new.Clear();
    LINALG::Matrix<6, 1> C_strain_new(false);
    C_strain_new.Clear();
    LINALG::Matrix<6, 1> iC_strain_new(false);
    iC_strain_new.Clear();
    LINALG::Matrix<6, 1> iC_stress_new(false);
    iC_stress_new.Clear();
    LINALG::Matrix<3, 1> prinv_new(false);
    prinv_new.Clear();
    LINALG::Matrix<3, 1> dPI_new(false);
    dPI_new.Clear();
    LINALG::Matrix<6, 1> ddPII_new(false);
    ddPII_new.Clear();

    LINALG::Matrix<NUM_STRESS_3D, 1> stressiso_new(true);
    stressiso_new.Clear();
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatiso_new(true);
    cmatiso_new.Clear();

    // Evaluate identity tensor
    IdentityMatrix(id2_new);

    // Evalutate Right Cauchy-Green strain tensor in strain-like Voigt notation
    EvaluateRightCauchyGreenStrainLikeVoigt(*E_strain_, C_strain_new);

    // Invert Right Cauchy Green Strain tensor
    VStrainUtils::InverseTensor(C_strain_new, iC_strain_new);
    VStrainUtils::ToStressLike(iC_strain_new, iC_stress_new);

    // Evaluate principle invariants
    InvariantsPrincipal<MAT::VoigtNotation::strain>(prinv_new, C_strain_new);

    // Evaluate derivatives of potsum w.r.t the principal invariants
    ElastHyperEvaluateInvariantDerivatives(
        prinv_new, dPI_new, ddPII_new, potsum_, properties, eleGID);

    ElastHyperAddIsotropicStressCmat(
        stressiso_new, cmatiso_new, C_strain_new, iC_strain_new, prinv_new, dPI_new, ddPII_new);
  }

 private:
  Teuchos::RCP<LINALG::Matrix<6, 1>> E_strain_;
  std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>> potsum_;
  SummandProperties properties;
};

#endif  // BACI_UNIT_ELASTHYPER_SERVICE_H
