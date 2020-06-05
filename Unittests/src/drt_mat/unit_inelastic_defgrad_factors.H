/*----------------------------------------------------------------------*/
/*! \file

\brief test cases for the inelastic_defgrad_factors methods

\level 3

\maintainer Christoph Schmidt

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_UNIT_INELASTIC_DEFGRAD_FACTORS_H
#define BACI_UNIT_INELASTIC_DEFGRAD_FACTORS_H

#include "cxxtest/TestSuite.h"
#include "../common/special_assertions.H"

#include "src/drt_inpar/inpar_material.H"

#include "src/drt_lib/drt_globalproblem.H"

#include "src/drt_mat/inelastic_defgrad_factors.H"
#include "src/drt_mat/matpar_bundle.H"

#include "src/linalg/linalg_fixedsizematrix.H"

// forward declaration
namespace MAT
{
  class InelasticDefgradFactors_TestSuite;
}

// class implementation
class MAT::InelasticDefgradFactors_TestSuite : public CxxTest::TestSuite
{
 public:
  /// set up test suite
  void setUp() override
  {
    // clang-format off
    // set up the deformation gradient
    FM_(0, 0) = 1.1;  FM_(0, 1) = 0.01; FM_(0, 2) = 0.03;
    FM_(1, 0) = 0.04; FM_(1, 1) = 1.2;  FM_(1, 2) = 0.02;
    FM_(2, 0) = 0.06; FM_(2, 1) = 0.05; FM_(2, 2) = 1.3;

    // set up the derivative of the second Piola-Kirchhoff stress tensor
    dSdiFin_(0, 0) =  9.83e+03; dSdiFin_(0, 1) =  2.15e+03; dSdiFin_(0, 2) =  2.14e+03; dSdiFin_(0, 3) = -5.49e+00; dSdiFin_(0, 4) =  6.38e-01; dSdiFin_(0, 5) = -3.29e+00; dSdiFin_(0, 6) =  2.13e+00; dSdiFin_(0, 7) =  6.38e-01; dSdiFin_(0, 8) =  1.28e+00;
    dSdiFin_(1, 0) =  1.80e+03; dSdiFin_(1, 1) =  9.49e+03; dSdiFin_(1, 2) =  1.80e+03; dSdiFin_(1, 3) =  1.79e+00; dSdiFin_(1, 4) = -1.75e+00; dSdiFin_(1, 5) =  1.07e+00; dSdiFin_(1, 6) = -5.83e+00; dSdiFin_(1, 7) =  5.36e-01; dSdiFin_(1, 8) =  1.07e+00;
    dSdiFin_(2, 0) =  1.55e+03; dSdiFin_(2, 1) =  1.54e+03; dSdiFin_(2, 2) =  9.23e+03; dSdiFin_(2, 3) =  1.53e+00; dSdiFin_(2, 4) =  4.59e-01; dSdiFin_(2, 5) =  9.18e-01; dSdiFin_(2, 6) =  1.53e+00; dSdiFin_(2, 7) = -1.82e+00; dSdiFin_(2, 8) = -3.65e+00;
    dSdiFin_(3, 0) = -8.76e+01; dSdiFin_(3, 1) = -8.75e+01; dSdiFin_(3, 2) = -8.37e+01; dSdiFin_(3, 3) =  3.84e+03; dSdiFin_(3, 4) = -2.31e+00; dSdiFin_(3, 5) = -1.19e+00; dSdiFin_(3, 6) =  3.84e+03; dSdiFin_(3, 7) = -2.49e-02; dSdiFin_(3, 8) = -4.98e-02;
    dSdiFin_(4, 0) = -8.97e+01; dSdiFin_(4, 1) = -9.08e+01; dSdiFin_(4, 2) = -9.07e+01; dSdiFin_(4, 3) = -8.89e-02; dSdiFin_(4, 4) =  3.85e+03; dSdiFin_(4, 5) = -5.33e-02; dSdiFin_(4, 6) = -2.37e+00; dSdiFin_(4, 7) =  3.84e+03; dSdiFin_(4, 8) = -3.86e+00;
    dSdiFin_(5, 0) = -1.40e+02; dSdiFin_(5, 1) = -1.37e+02; dSdiFin_(5, 2) = -1.40e+02; dSdiFin_(5, 3) = -1.28e+00; dSdiFin_(5, 4) = -4.09e-02; dSdiFin_(5, 5) =  3.85e+03; dSdiFin_(5, 6) = -1.36e-01; dSdiFin_(5, 7) = -3.85e+00; dSdiFin_(5, 8) =  3.84e+03;

    // set up reference solution of inverse inelastic deformation gradient of InelasticDefgradLinScalarIso
    iFin_lin_scalar_iso_solution_.Clear();
    iFin_lin_scalar_iso_solution_(0, 0) = iFin_lin_scalar_iso_solution_(1, 1) = iFin_lin_scalar_iso_solution_(2, 2) = 0.9948803003804159;

    // set up reference solution of inverse inelastic deformation gradient of InelasticDefgradLinScalarAniso
    iFin_lin_scalar_aniso_solution_(0, 0) =  0.9885965678505413; iFin_lin_scalar_aniso_solution_(0, 1) = -0.0057017160747293; iFin_lin_scalar_aniso_solution_(0, 2) = -0.0034210296448376;
    iFin_lin_scalar_aniso_solution_(1, 0) = -0.0057017160747293; iFin_lin_scalar_aniso_solution_(1, 1) =  0.9971491419626352; iFin_lin_scalar_aniso_solution_(1, 2) = -0.0017105148224188;
    iFin_lin_scalar_aniso_solution_(2, 0) = -0.0034210296448376; iFin_lin_scalar_aniso_solution_(2, 1) = -0.0017105148224188; iFin_lin_scalar_aniso_solution_(2, 2) =  0.9989736911065487;

    // set up reference solution of inverse inelastic deformation gradient of InelasticDefgradPolyIntercalFracIso
    iFin_poly_intercal_frac_iso_solution_.Clear();
    iFin_poly_intercal_frac_iso_solution_(0, 0) = iFin_poly_intercal_frac_iso_solution_(1, 1) = iFin_poly_intercal_frac_iso_solution_(2, 2) = 0.9991151119224016;

    // set up reference solution of inverse inelastic deformation gradient of InelasticDefgradPolyIntercalFracAniso
    iFin_poly_intercal_frac_aniso_solution_(0, 0) =  0.9980206598171963; iFin_poly_intercal_frac_aniso_solution_(0, 1) = -0.0009896700914018; iFin_poly_intercal_frac_aniso_solution_(0, 2) = -0.0005938020548410;
    iFin_poly_intercal_frac_aniso_solution_(1, 0) = -0.0009896700914018; iFin_poly_intercal_frac_aniso_solution_(1, 1) =  0.9995051649542991; iFin_poly_intercal_frac_aniso_solution_(1, 2) = -0.0002969010274205;
    iFin_poly_intercal_frac_aniso_solution_(2, 0) = -0.0005938020548410; iFin_poly_intercal_frac_aniso_solution_(2, 1) = -0.0002969010274205; iFin_poly_intercal_frac_aniso_solution_(2, 2) =  0.9998218593835476;
    // clang-format on

    // prepare variables needed to instantiate linear shape evaluation object
    const double GrowthFac(5.27e-7);
    const double RefConc(46456.0);

    // create linear shape evaluation object
    linear_shape_ = Teuchos::rcp(new InelasticDefgradLinearShape(GrowthFac, RefConc));

    // prepare variables needed to instantiate polynomial shape evaluation object
    std::vector<double> PolyCoeffs{0.1051717305, -3.9012322937, 31.9658107225, -122.8624633232,
        258.6769103514, -306.7800791732, 192.5096604774, -49.7490196448};
    const double Xmin(0.152);
    const double Xmax(0.887);

    // create polynomial shape evaluation object
    polynomial_shape_ = Teuchos::rcp(new InelasticDefgradPolynomialShape(PolyCoeffs, Xmin, Xmax));

    // parameter list to be passed to PreEvaluate
    Teuchos::ParameterList params_lin;
    // set up a dummy concentration vector and store it to the parameter list
    auto gpconc_lin =
        Teuchos::rcp(new std::vector<std::vector<double>>(1, std::vector<double>(1, 44327.362)));
    params_lin.set<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc", gpconc_lin);

    // create InelasticDefgradLinScalarIso object initialize container for material parameters
    const Teuchos::RCP<MAT::PAR::Material> InelasticDefgradScalarData =
        Teuchos::rcp(new MAT::PAR::Material());

    InelasticDefgradScalarData->Add("SCALAR1", 1);
    InelasticDefgradScalarData->Add("SCALAR1_MolarGrowthFac", GrowthFac);
    InelasticDefgradScalarData->Add("SCALAR1_RefConc", RefConc);

    InelasticDefgradScalarData->SetParameter(
        new MAT::PAR::InelasticDefgradScalar(InelasticDefgradScalarData));

    // get pointer to parameter class
    params_inelastic_defgrad_scalar_ =
        Teuchos::rcp(new MAT::PAR::InelasticDefgradScalar(InelasticDefgradScalarData));

    // setup pointer to InelasticDefgradLinScalarIso object
    lin_scalar_iso_ = Teuchos::rcp(new MAT::InelasticDefgradLinScalarIso(
        params_inelastic_defgrad_scalar_.get(), linear_shape_));

    // call PreEvaluate to set the concentration value
    lin_scalar_iso_->PreEvaluate(params_lin, 0);

    // create InelasticDefgradLinScalarAniso object initialize container for material parameters
    const Teuchos::RCP<MAT::PAR::Material> InelasticDefgradLinScalarAnisoData =
        Teuchos::rcp(new MAT::PAR::Material());

    InelasticDefgradLinScalarAnisoData->Add("SCALAR1", 1);
    InelasticDefgradLinScalarAnisoData->Add("SCALAR1_MolarGrowthFac", GrowthFac);
    InelasticDefgradLinScalarAnisoData->Add("SCALAR1_RefConc", RefConc);

    // vector to instantiate the deformation direction object
    const std::vector<double> growthdir{1.0, 0.5, 0.3};
    InelasticDefgradLinScalarAnisoData->Add("GrowthDirection", growthdir);

    InelasticDefgradLinScalarAnisoData->SetParameter(
        new MAT::PAR::InelasticDefgradLinScalarAniso(InelasticDefgradLinScalarAnisoData));

    // get pointer to parameter class
    params_inelastic_defgrad_lin_scalar_aniso_ = Teuchos::rcp(
        new MAT::PAR::InelasticDefgradLinScalarAniso(InelasticDefgradLinScalarAnisoData));

    // set up pointer to InelasticDefgradLinScalarAniso object
    lin_scalar_aniso_ = Teuchos::rcp(new MAT::InelasticDefgradLinScalarAniso(
        params_inelastic_defgrad_lin_scalar_aniso_.get(), linear_shape_));

    // call PreEvaluate to set the concentration value
    lin_scalar_aniso_->PreEvaluate(params_lin, 0);

    // InelasticDefgradPolyIntercalFracIso object initialize container for required electrode
    // material parameters
    const int problemid(0);
    DRT::Problem& problem = (*DRT::Problem::Instance());
    problem.Materials()->SetReadFromProblem(problemid);
    // set up material to be added to problem instance
    const int matid(1);
    const Teuchos::RCP<MAT::PAR::Material> electrode_mat =
        Teuchos::rcp(new MAT::PAR::Material(matid, INPAR::MAT::m_electrode, "MAT_electrode"));

    // add dummy parameters to electrode material
    electrode_mat->Add("DIFFCOEF", 0);
    electrode_mat->Add("DIFF_PARA_NUM", 0);
    electrode_mat->Add("DIFF_PARA", std::vector<double>(0, 0.0));
    electrode_mat->Add("COND", 0);
    electrode_mat->Add("COND_PARA_NUM", 0);
    electrode_mat->Add("COND_PARA", std::vector<double>(0, 0.0));
    electrode_mat->Add("OCP_MODEL", "Polynomial");
    electrode_mat->Add("X_MIN", -1.0);
    electrode_mat->Add("X_MAX", -1.0);
    electrode_mat->Add("OCP_PARA_NUM", 1);
    std::vector<double> ocp_para(1, 0.0);
    electrode_mat->Add("OCP_PARA", ocp_para);

    // add actually required parameters to electrode material
    const double c_max(4.91375e4);
    const double chi_max(1.0);
    electrode_mat->Add("C_MAX", c_max);
    electrode_mat->Add("CHI_MAX", chi_max);

    // add material to problem instance
    problem.Materials()->Insert(matid, electrode_mat);

    // parameter list to be passed to PreEvaluate
    Teuchos::ParameterList params_poly;
    // set up a dummy concentration vector and store it to the parameter list
    auto gpconc_poly =
        Teuchos::rcp(new std::vector<std::vector<double>>(1, std::vector<double>(1, 22641.893)));
    params_poly.set<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc", gpconc_poly);

    // initialize container for material parameters
    const Teuchos::RCP<MAT::PAR::Material> InelasticDefgradPolyIntercalFracData =
        Teuchos::rcp(new MAT::PAR::Material());

    InelasticDefgradPolyIntercalFracData->Add("MATID", matid);
    InelasticDefgradPolyIntercalFracData->Add("SCALAR1", 1);
    InelasticDefgradPolyIntercalFracData->Add("SCALAR1_RefConc", RefConc);

    InelasticDefgradPolyIntercalFracData->SetParameter(
        new MAT::PAR::InelasticDefgradPolyIntercalFrac(InelasticDefgradPolyIntercalFracData));

    // get pointer to parameter class
    params_inelastic_defgrad_poly_intercal_frac_ = Teuchos::rcp(
        new MAT::PAR::InelasticDefgradPolyIntercalFrac(InelasticDefgradPolyIntercalFracData));

    // get reference intercalation fraction
    const double X_ref = MAT::Electrode::ComputeIntercalationFraction(RefConc, chi_max, c_max, 1.0);

    // set the value of the reference polynomial
    params_inelastic_defgrad_poly_intercal_frac_->SetPolynomReferenceValue(
        polynomial_shape_->ComputePolynomial(X_ref));

    // set up pointer to InelasticDefgradPolyIntercalFracIso object
    poly_intercal_frac_iso_ = Teuchos::rcp(new MAT::InelasticDefgradPolyIntercalFracIso(
        params_inelastic_defgrad_poly_intercal_frac_.get(), polynomial_shape_));

    // call PreEvaluate to set the concentration value
    poly_intercal_frac_iso_->PreEvaluate(params_poly, 0);

    // create InelasticDefgradPolyIntercalFracAniso object initialize container for material
    // parameters
    const Teuchos::RCP<MAT::PAR::Material> InelasticDefgradPolyIntercalFracAnisoData =
        Teuchos::rcp(new MAT::PAR::Material());

    InelasticDefgradPolyIntercalFracAnisoData->Add("MATID", matid);
    InelasticDefgradPolyIntercalFracAnisoData->Add("SCALAR1", 1);
    InelasticDefgradPolyIntercalFracAnisoData->Add("SCALAR1_RefConc", RefConc);
    InelasticDefgradPolyIntercalFracAnisoData->Add("GrowthDirection", growthdir);

    InelasticDefgradPolyIntercalFracAnisoData->SetParameter(
        new MAT::PAR::InelasticDefgradPolyIntercalFracAniso(
            InelasticDefgradPolyIntercalFracAnisoData));

    // get pointer to parameter class
    params_inelastic_defgrad_poly_intercal_frac_aniso_ =
        Teuchos::rcp(new MAT::PAR::InelasticDefgradPolyIntercalFracAniso(
            InelasticDefgradPolyIntercalFracAnisoData));

    // set the value of the reference polynomial
    params_inelastic_defgrad_poly_intercal_frac_aniso_->SetPolynomReferenceValue(
        polynomial_shape_->ComputePolynomial(X_ref));

    // set up pointer to InelasticDefgradPolyIntercalFracIso object
    poly_intercal_frac_aniso_ = Teuchos::rcp(new MAT::InelasticDefgradPolyIntercalFracAniso(
        params_inelastic_defgrad_poly_intercal_frac_aniso_.get(), polynomial_shape_));

    // call PreEvaluate to set the concentration value
    poly_intercal_frac_aniso_->PreEvaluate(params_poly, 0);
  }

  /// tear down test suite
  void tearDown() override
  {
    linear_shape_ = Teuchos::null;
    polynomial_shape_ = Teuchos::null;
    lin_scalar_iso_ = Teuchos::null;
    params_inelastic_defgrad_scalar_ = Teuchos::null;
    lin_scalar_aniso_ = Teuchos::null;
    params_inelastic_defgrad_lin_scalar_aniso_ = Teuchos::null;
    poly_intercal_frac_iso_ = Teuchos::null;
    params_inelastic_defgrad_poly_intercal_frac_ = Teuchos::null;
    poly_intercal_frac_aniso_ = Teuchos::null;
    params_inelastic_defgrad_poly_intercal_frac_aniso_ = Teuchos::null;

    // We need to make sure the DRT::Problem instance created in setUp is deleted again. If this is
    // not done, some troubles arise where unit tests influence each other on some configurations.
    // We suspect that missing singleton destruction might be the reason for that.
    DRT::Problem::Done();
  }

  /*--------------------------------------------------------------------*
   *--------------------------------------------------------------------*/
  void test_EvaluateLinearGrowth()
  {
    const std::vector<double> TestValues{952834.0233, 44327.362, 12823.902};
    const std::vector<double> Results{0.4776612182790999, -0.0011217922259999, -0.0177241156459999};

    // loop over test values and check whether the result is correct
    for (auto i = 0U; i < TestValues.size(); ++i)
      TS_ASSERT_DELTA(linear_shape_->EvaluateLinearGrowth(TestValues[i]), Results[i], 1.0e-12)
  }

  /*--------------------------------------------------------------------*
   *--------------------------------------------------------------------*/
  void test_InelasticDeformationDirection()
  {
    // growth directions to be tested
    const std::vector<std::vector<double>> growth_directions{
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.5, 0.3}};

    // result of first growth direction
    LINALG::Matrix<3, 3> growth_dir_1(true);
    growth_dir_1(0, 0) = 1.0;

    // result of second growth direction
    LINALG::Matrix<3, 3> growth_dir_2(true);
    growth_dir_2(1, 1) = 1.0;

    // result of third growth direction
    LINALG::Matrix<3, 3> growth_dir_3(true);
    growth_dir_3(2, 2) = 1.0;

    // result of fourth growth direction
    LINALG::Matrix<3, 3> growth_dir_4(true);
    // clang-format off
    growth_dir_4(0, 0) = 0.74626865671641791044776119402985; growth_dir_4(0, 1) = 0.37313432835820895522388059701493; growth_dir_4(0, 2) = 0.22388059701492537313432835820896;
    growth_dir_4(1, 0) = growth_dir_4(0,1); growth_dir_4(1, 1) = 0.18656716417910447761194029850746; growth_dir_4(1, 2) = 0.11194029850746268656716417910448;
    growth_dir_4(2, 0) = growth_dir_4(0, 2); growth_dir_4(2, 1) = growth_dir_4(1, 2); growth_dir_4(2, 2) = 0.067164179104477611940298507462687;
    // clang-format on

    // put all results together
    const std::vector<LINALG::Matrix<3, 3>> growth_direction_solutions{
        growth_dir_1, growth_dir_2, growth_dir_3, growth_dir_4};

    // loop over all growth directions to be tested and perform the test
    for (auto i = 0U; i < growth_directions.size(); ++i)
    {
      // create object
      auto growth_direction = MAT::PAR::InelasticDeformationDirection(growth_directions[i]);
      // check the results
      TESTING::AssertDelta(growth_direction.Growthdirmat(), growth_direction_solutions[i], 1.0e-12);
    }
  }

  /*--------------------------------------------------------------------*
   *--------------------------------------------------------------------*/
  void test_ComputePolynomial()
  {
    const std::vector<double> TestValues{0.215, 0.462, 0.675, 0.802};
    const std::vector<double> Results{
        -0.0472760218320868, -0.0136018625081755, -0.0054533973886232, -0.0027181447027916};

    // loop over test values and check whether the result is correct
    for (auto i = 0U; i < TestValues.size(); ++i)
      TS_ASSERT_DELTA(polynomial_shape_->ComputePolynomial(TestValues[i]), Results[i], 1.0e-12)
  }

  /*--------------------------------------------------------------------*
   *--------------------------------------------------------------------*/
  void test_ComputePolynomialDerivative()
  {
    const std::vector<double> TestValues{0.215, 0.462, 0.675, 0.802};
    const std::vector<double> Results{
        0.3081029423168080, 0.0393764326640364, 0.0224486538440280, 0.0329592765802879};

    // loop over test values and check whether the result is correct
    for (auto i = 0U; i < TestValues.size(); ++i)
      TS_ASSERT_DELTA(
          polynomial_shape_->ComputePolynomialDerivative(TestValues[i]), Results[i], 1.0e-12)
  }

  /*--------------------------------------------------------------------*
   *--------------------------------------------------------------------*/
  void test_EvaluateInverseInelasticDefGrad()
  {
    // create matrix to be filled by tested methods
    LINALG::Matrix<3, 3> iFin(true);

    // test InelasticDefgradLinScalarIso evaluate the method
    lin_scalar_iso_->EvaluateInverseInelasticDefGrad(&FM_, iFin);

    // compare the results
    TESTING::AssertDelta(iFin, iFin_lin_scalar_iso_solution_, 1.0e-10);

    // clear matrix to be filled again
    iFin.Clear();

    // test InelasticDefgradLinScalarAniso evaluate the method
    lin_scalar_aniso_->EvaluateInverseInelasticDefGrad(&FM_, iFin);

    // compare the results
    TESTING::AssertDelta(iFin, iFin_lin_scalar_aniso_solution_, 1.0e-10);

    // clear matrix to be filled again
    iFin.Clear();

    // test InelasticDefgradPolyIntercalFracIso evaluate the method
    poly_intercal_frac_iso_->EvaluateInverseInelasticDefGrad(&FM_, iFin);

    // compare the results
    TESTING::AssertDelta(iFin, iFin_poly_intercal_frac_iso_solution_, 1.0e-10);

    // clear matrix to be filled again
    iFin.Clear();

    // test InelasticDefgradPolyIntercalFracAniso evaluate the method
    poly_intercal_frac_aniso_->EvaluateInverseInelasticDefGrad(&FM_, iFin);

    // compare the results
    TESTING::AssertDelta(iFin, iFin_poly_intercal_frac_aniso_solution_, 1.0e-10);
  }

  /*--------------------------------------------------------------------*
   *--------------------------------------------------------------------*/
  void test_EvaluateAdditionalCmat()
  {
    // calculate the inverse Cauchy-Green tensor in voigt notation
    LINALG::Matrix<3, 3> CM;
    LINALG::Matrix<3, 3> iCM;
    LINALG::Matrix<6, 1> iCV;
    CM.MultiplyTN(1.0, FM_, FM_, 0.0);
    iCM.Invert(CM);
    UTILS::VOIGT::Stresses::MatrixToVector(iCM, iCV);

    // matrix to be filled by the methods
    LINALG::Matrix<6, 6> CMatAdd(true);
    LINALG::Matrix<6, 6> CMatAdd_ref_solution(true);

    // test InelasticDefgradLinScalarIso set up reference solution
    // clang-format off
    CMatAdd_ref_solution(0, 0) = -1.5298408106321838e+02; CMatAdd_ref_solution(0, 1) = -1.2850516942246250e+02; CMatAdd_ref_solution(0, 2) = -1.1003777104023889e+02; CMatAdd_ref_solution(0, 3) =  5.9684525433019955e+00; CMatAdd_ref_solution(0, 4) =  6.3904917135249439e+00; CMatAdd_ref_solution(0, 5) =  9.7975743219718243e+00;
    CMatAdd_ref_solution(1, 0) = -1.4182447741625558e+02; CMatAdd_ref_solution(1, 1) = -1.1913120876345852e+02; CMatAdd_ref_solution(1, 2) = -1.0201093646719031e+02; CMatAdd_ref_solution(1, 3) =  5.5330767557948386e+00; CMatAdd_ref_solution(1, 4) =  5.9243297825808430e+00; CMatAdd_ref_solution(1, 5) =  9.0828787446608494e+00;
    CMatAdd_ref_solution(2, 0) = -1.3348186109765231e+02; CMatAdd_ref_solution(2, 1) = -1.1212349060090213e+02; CMatAdd_ref_solution(2, 2) = -9.6010293145590879e+01; CMatAdd_ref_solution(2, 3) =  5.2076016525127891e+00; CMatAdd_ref_solution(2, 4) =  5.5758397953702055e+00; CMatAdd_ref_solution(2, 5) =  8.5485917596808001e+00;
    CMatAdd_ref_solution(3, 0) =  2.8039858483825015e+00; CMatAdd_ref_solution(3, 1) =  2.3553213772332362e+00; CMatAdd_ref_solution(3, 2) =  2.0168395995193928e+00; CMatAdd_ref_solution(3, 3) = -1.0939345029791475e-01; CMatAdd_ref_solution(3, 4) = -1.1712884245469231e-01; CMatAdd_ref_solution(3, 5) = -1.7957593728939861e-01;
    CMatAdd_ref_solution(4, 0) =  2.9383344748119566e+00; CMatAdd_ref_solution(4, 1) =  2.4681729424484300e+00; CMatAdd_ref_solution(4, 2) =  2.1134733361269680e+00; CMatAdd_ref_solution(4, 3) = -1.1463486754557373e-01; CMatAdd_ref_solution(4, 4) = -1.2274088900198049e-01; CMatAdd_ref_solution(4, 5) = -1.8818003938518124e-01;
    CMatAdd_ref_solution(5, 0) =  4.5180142920228095e+00; CMatAdd_ref_solution(5, 1) =  3.7950889269948203e+00; CMatAdd_ref_solution(5, 2) =  3.2496990455934576e+00; CMatAdd_ref_solution(5, 3) = -1.7626378969949943e-01; CMatAdd_ref_solution(5, 4) = -1.8872769437251424e-01; CMatAdd_ref_solution(5, 5) = -2.8934762693075433e-01;
    // clang-format on

    // evaluate the method
    lin_scalar_iso_->EvaluateAdditionalCmat(
        &FM_, iFin_lin_scalar_iso_solution_, iCV, dSdiFin_, CMatAdd);

    // compare the results
    TESTING::AssertDelta(CMatAdd, CMatAdd_ref_solution, 1.0e-10);

    // clear matrix before test of next method
    CMatAdd.Clear();

    // test InelasticDefgradLinScalarAniso set up reference solution
    // clang-format off
    CMatAdd_ref_solution(0, 0) = -2.5348465151755443e+02; CMatAdd_ref_solution(0, 1) = -2.1292469035256323e+02; CMatAdd_ref_solution(0, 2) = -1.8232541485396175e+02; CMatAdd_ref_solution(0, 3) =  9.8893368677541016e+00; CMatAdd_ref_solution(0, 4) =  1.0588628266226561e+01; CMatAdd_ref_solution(0, 5) =  1.6233942090328238e+01;
    CMatAdd_ref_solution(1, 0) = -1.0402849846169919e+02; CMatAdd_ref_solution(1, 1) = -8.7382946818243298e+01; CMatAdd_ref_solution(1, 2) = -7.4825197601167204e+01; CMatAdd_ref_solution(1, 3) =  4.0585213304843331e+00; CMatAdd_ref_solution(1, 4) =  4.3455060995216570e+00; CMatAdd_ref_solution(1, 5) =  6.6623072034563453e+00;
    CMatAdd_ref_solution(2, 0) = -6.6413171948881100e+01; CMatAdd_ref_solution(2, 1) = -5.5786431201605694e+01; CMatAdd_ref_solution(2, 2) = -4.7769397692738316e+01; CMatAdd_ref_solution(2, 3) =  2.5910137987706867e+00; CMatAdd_ref_solution(2, 4) =  2.7742286782952910e+00; CMatAdd_ref_solution(2, 5) =  4.2533052040765602e+00;
    CMatAdd_ref_solution(3, 0) = -8.9367069403395547e+01; CMatAdd_ref_solution(3, 1) = -7.5067486202872757e+01; CMatAdd_ref_solution(3, 2) = -6.4279584210361918e+01; CMatAdd_ref_solution(3, 3) =  3.4865268919563563e+00; CMatAdd_ref_solution(3, 4) =  3.7330649863393925e+00; CMatAdd_ref_solution(3, 5) =  5.7233438821308598e+00;
    CMatAdd_ref_solution(4, 0) = -2.4741829181888583e+01; CMatAdd_ref_solution(4, 1) = -2.0782900604712967e+01; CMatAdd_ref_solution(4, 2) = -1.7796202818699221e+01; CMatAdd_ref_solution(4, 3) =  9.6526666225856417e-01; CMatAdd_ref_solution(4, 4) =  1.0335222675813636e+00; CMatAdd_ref_solution(4, 5) =  1.5845433628542860e+00;
    CMatAdd_ref_solution(5, 0) = -5.0871154285791278e+01; CMatAdd_ref_solution(5, 1) = -4.2731284554439576e+01; CMatAdd_ref_solution(5, 2) = -3.6590398092069350e+01; CMatAdd_ref_solution(5, 3) =  1.9846644701043878e+00; CMatAdd_ref_solution(5, 4) =  2.1250033837602929e+00; CMatAdd_ref_solution(5, 5) =  3.2579462614386228e+00;
    // clang-format on

    // evaluate the method
    lin_scalar_aniso_->EvaluateAdditionalCmat(
        &FM_, iFin_lin_scalar_aniso_solution_, iCV, dSdiFin_, CMatAdd);

    // compare the results
    TESTING::AssertDelta(CMatAdd, CMatAdd_ref_solution, 1.0e-10);

    // clear matrix before test of next method
    CMatAdd.Clear();

    // test InelasticDefgradPolyIntercalFracIso set up reference solution
    // clang-format off
    CMatAdd_ref_solution(0, 0) = -8.9980093707125818e+01; CMatAdd_ref_solution(0, 1) = -7.5582420772949916e+01; CMatAdd_ref_solution(0, 2) = -6.4720517852000711e+01; CMatAdd_ref_solution(0, 3) =  3.5104431480744864e+00; CMatAdd_ref_solution(0, 4) =  3.7586724005615206e+00; CMatAdd_ref_solution(0, 5) =  5.7626038569938007e+00;
    CMatAdd_ref_solution(1, 0) = -8.3416389987696661e+01; CMatAdd_ref_solution(1, 1) = -7.0068972232146919e+01; CMatAdd_ref_solution(1, 2) = -5.9999403589425583e+01; CMatAdd_ref_solution(1, 3) =  3.2543697456299596e+00; CMatAdd_ref_solution(1, 4) =  3.4844916234667358e+00; CMatAdd_ref_solution(1, 5) =  5.3422439439128073e+00;
    CMatAdd_ref_solution(2, 0) = -7.8509543517832142e+01; CMatAdd_ref_solution(2, 1) = -6.5947267983197094e+01; CMatAdd_ref_solution(2, 2) = -5.6470026907694674e+01; CMatAdd_ref_solution(2, 3) =  3.0629362311811383e+00; CMatAdd_ref_solution(2, 4) =  3.2795215279686927e+00; CMatAdd_ref_solution(2, 5) =  5.0279943001532299e+00;
    CMatAdd_ref_solution(3, 0) =  1.6492102161051103e+00; CMatAdd_ref_solution(3, 1) =  1.3853208566600168e+00; CMatAdd_ref_solution(3, 2) =  1.1862372535480017e+00; CMatAdd_ref_solution(3, 3) = -6.4341550051110280e-02; CMatAdd_ref_solution(3, 4) = -6.8891247681680007e-02; CMatAdd_ref_solution(3, 5) = -1.0562052961685520e-01;
    CMatAdd_ref_solution(4, 0) =  1.7282295618535779e+00; CMatAdd_ref_solution(4, 1) =  1.4516963536560918e+00; CMatAdd_ref_solution(4, 2) =  1.2430739689421100e+00; CMatAdd_ref_solution(4, 3) = -6.7424375478597795e-02; CMatAdd_ref_solution(4, 4) = -7.2192064803986156e-02; CMatAdd_ref_solution(4, 5) = -1.1068117323064577e-01;
    CMatAdd_ref_solution(5, 0) =  2.6573441271863643e+00; CMatAdd_ref_solution(5, 1) =  2.2321437296260704e+00; CMatAdd_ref_solution(5, 2) =  1.9113637354308992e+00; CMatAdd_ref_solution(5, 3) = -1.0367243574695899e-01; CMatAdd_ref_solution(5, 4) = -1.1100328548400526e-01; CMatAdd_ref_solution(5, 5) = -1.7018454733473193e-01;
    // clang-format on

    // evaluate the method
    poly_intercal_frac_iso_->EvaluateAdditionalCmat(
        &FM_, iFin_poly_intercal_frac_iso_solution_, iCV, dSdiFin_, CMatAdd);

    // compare the results
    TESTING::AssertDelta(CMatAdd, CMatAdd_ref_solution, 1.0e-10);

    // clear matrix before test of next method
    CMatAdd.Clear();

    // test InelasticDefgradPolyIntercalFracAniso set up reference solution
    // clang-format off
    CMatAdd_ref_solution(0, 0) = -1.5036309611574805e+02; CMatAdd_ref_solution(0, 1) = -1.2630356705712084e+02; CMatAdd_ref_solution(0, 2) = -1.0815255958852929e+02; CMatAdd_ref_solution(0, 3) =  5.8661986083372124e+00; CMatAdd_ref_solution(0, 4) =  6.2810072333641545e+00; CMatAdd_ref_solution(0, 5) =  9.6297183291054971e+00;
    CMatAdd_ref_solution(1, 0) = -6.1708064055666021e+01; CMatAdd_ref_solution(1, 1) = -5.1834185433507010e+01; CMatAdd_ref_solution(1, 2) = -4.4385126718431692e+01; CMatAdd_ref_solution(1, 3) =  2.4074508229590710e+00; CMatAdd_ref_solution(1, 4) =  2.5776856602645020e+00; CMatAdd_ref_solution(1, 5) =  3.9519755235222713e+00;
    CMatAdd_ref_solution(2, 0) = -3.9395245815937521e+01; CMatAdd_ref_solution(2, 1) = -3.3091630860106335e+01; CMatAdd_ref_solution(2, 2) = -2.8336053065395109e+01; CMatAdd_ref_solution(2, 3) =  1.5369485076488196e+00; CMatAdd_ref_solution(2, 4) =  1.6456286836471103e+00; CMatAdd_ref_solution(2, 5) =  2.5229935437171265e+00;
    CMatAdd_ref_solution(3, 0) = -5.3011135648009507e+01; CMatAdd_ref_solution(3, 1) = -4.4528848494436232e+01; CMatAdd_ref_solution(3, 2) = -3.8129635230532394e+01; CMatAdd_ref_solution(3, 3) =  2.0681527462386371e+00; CMatAdd_ref_solution(3, 4) =  2.2143952542562975e+00; CMatAdd_ref_solution(3, 5) =  3.3949972951033844e+00;
    CMatAdd_ref_solution(4, 0) = -1.4676462724995005e+01; CMatAdd_ref_solution(4, 1) = -1.2328088752048501e+01; CMatAdd_ref_solution(4, 2) = -1.0556426745775207e+01; CMatAdd_ref_solution(4, 3) =  5.7258095527910324e-01; CMatAdd_ref_solution(4, 4) =  6.1306910350483479e-01; CMatAdd_ref_solution(4, 5) =  9.3992612389760388e-01;
    CMatAdd_ref_solution(5, 0) = -3.0175966140749807e+01; CMatAdd_ref_solution(5, 1) = -2.5347523836817480e+01; CMatAdd_ref_solution(5, 2) = -2.1704846870581804e+01; CMatAdd_ref_solution(5, 3) =  1.1772716521069135e+00; CMatAdd_ref_solution(5, 4) =  1.2605184815953689e+00; CMatAdd_ref_solution(5, 5) =  1.9325623224754220e+00;
    // clang-format on

    // evaluate the method
    poly_intercal_frac_aniso_->EvaluateAdditionalCmat(
        &FM_, iFin_poly_intercal_frac_aniso_solution_, iCV, dSdiFin_, CMatAdd);

    // compare the results
    TESTING::AssertDelta(CMatAdd, CMatAdd_ref_solution, 1.0e-10);
  }

  /*--------------------------------------------------------------------*
   *--------------------------------------------------------------------*/
  void test_EvaluateODStiffMat()
  {
    LINALG::Matrix<6, 1> dSdc(true);
    LINALG::Matrix<6, 1> dSdc_ref_solution(true);

    // test InelasticDefgradLinScalarIso set up reference solution
    // clang-format off
    dSdc_ref_solution(0) = -4.1608801904598307e-03; dSdc_ref_solution(1) = -3.8573598932803953e-03; dSdc_ref_solution(2) = -3.6304563701462541e-03; dSdc_ref_solution(3) = 7.6263158165085275e-05; dSdc_ref_solution(4) = 7.9917188927245458e-05; dSdc_ref_solution(5) = 1.2288151837264513e-04;
    // clang-format on

    // evaluate the method
    lin_scalar_iso_->EvaluateODStiffMat(&FM_, iFin_lin_scalar_iso_solution_, dSdiFin_, dSdc);

    // compare the results
    TESTING::AssertDelta(dSdc, dSdc_ref_solution, 1.0e-10);

    // clear the derivative to be refilled again by next method to be tested
    dSdc.Clear();

    // test InelasticDefgradLinScalarAniso set up reference solution
    // clang-format off
    dSdc_ref_solution(0) = -6.8943072884109979e-03; dSdc_ref_solution(1) = -2.8293801255942189e-03; dSdc_ref_solution(2) = -1.8063137656362403e-03; dSdc_ref_solution(3) = -2.4306167424464095e-03; dSdc_ref_solution(4) = -6.7293136778145881e-04; dSdc_ref_solution(5) = -1.3836000233652130e-03;
    // clang-format on

    // evaluate the method
    lin_scalar_aniso_->EvaluateODStiffMat(&FM_, iFin_lin_scalar_aniso_solution_, dSdiFin_, dSdc);

    // compare the results
    TESTING::AssertDelta(dSdc, dSdc_ref_solution, 1.0e-10);

    // clear the derivative to be refilled again by next method to be tested
    dSdc.Clear();

    // test InelasticDefgradPolyIntercalFracIso set up reference solution
    // clang-format off
    dSdc_ref_solution(0) = -4.7912028948510080e-03; dSdc_ref_solution(1) = -4.4417029669688171e-03; dSdc_ref_solution(2) = -4.1804263218530044e-03; dSdc_ref_solution(3) = 8.7816098384379664e-05; dSdc_ref_solution(4) = 9.2023670331699243e-05; dSdc_ref_solution(5) = 1.4149657274453755e-04;
    // clang-format on

    // evaluate the method
    poly_intercal_frac_iso_->EvaluateODStiffMat(
        &FM_, iFin_poly_intercal_frac_iso_solution_, dSdiFin_, dSdc);

    // compare the results
    TESTING::AssertDelta(dSdc, dSdc_ref_solution, 1.0e-10);

    // clear the derivative to be refilled again by next method to be tested
    dSdc.Clear();

    // test InelasticDefgradPolyIntercalFracAniso set up reference solution
    // clang-format off
    dSdc_ref_solution(0) = -8.0064386655720986e-03; dSdc_ref_solution(1) = -3.2857918119254472e-03; dSdc_ref_solution(2) = -2.0976930343176458e-03; dSdc_ref_solution(3) = -2.8227032903830871e-03; dSdc_ref_solution(4) = -7.8148296803340908e-04; dSdc_ref_solution(5) = -1.6067906841603652e-03;
    // clang-format on

    // evaluate the method
    poly_intercal_frac_aniso_->EvaluateODStiffMat(
        &FM_, iFin_poly_intercal_frac_aniso_solution_, dSdiFin_, dSdc);

    // compare the results
    TESTING::AssertDelta(dSdc, dSdc_ref_solution, 1.0e-10);
  }

 private:
  // deformation gradient
  LINALG::Matrix<3, 3> FM_;
  // derivative of second Piola-Kirchhoff stress tensor w.r.t. inverse inelastic deformation
  // gradient
  LINALG::Matrix<6, 9> dSdiFin_;
  // reference solution of inverse inelastic deformation gradient using InelasticDefgradLinScalarIso
  LINALG::Matrix<3, 3> iFin_lin_scalar_iso_solution_;
  // reference solution of inverse inelastic deformation gradient using
  // InelasticDefgradLinScalarAniso
  LINALG::Matrix<3, 3> iFin_lin_scalar_aniso_solution_;
  // reference solution of inverse inelastic deformation gradient using
  // InelasticDefgradPolyIntercalFracIso
  LINALG::Matrix<3, 3> iFin_poly_intercal_frac_iso_solution_;
  // reference solution of inverse inelastic deformation gradient using
  // InelasticDefgradPolyIntercalFracAniso
  LINALG::Matrix<3, 3> iFin_poly_intercal_frac_aniso_solution_;
  // pointer to object that evaluates a linear shape
  Teuchos::RCP<MAT::InelasticDefgradLinearShape> linear_shape_;
  // pointer to object that evaluates a polynomial shape
  Teuchos::RCP<MAT::InelasticDefgradPolynomialShape> polynomial_shape_;
  // pointer to InelasticDefgradLinScalarIso object
  Teuchos::RCP<MAT::InelasticDefgradLinScalarIso> lin_scalar_iso_;
  // pointer to parameters of InelasticDefgradScalar
  Teuchos::RCP<MAT::PAR::InelasticDefgradScalar> params_inelastic_defgrad_scalar_;
  // pointer to InelasticDefgradLinScalarAniso object
  Teuchos::RCP<MAT::InelasticDefgradLinScalarAniso> lin_scalar_aniso_;
  // pointer to parameters of InelasticDefgradScalar
  Teuchos::RCP<MAT::PAR::InelasticDefgradLinScalarAniso> params_inelastic_defgrad_lin_scalar_aniso_;
  // pointer to InelasticDefgradPolyIntercalFracIso object
  Teuchos::RCP<MAT::InelasticDefgradPolyIntercalFracIso> poly_intercal_frac_iso_;
  // pointer to parameters of InelasticDefgradIntercalFrac
  Teuchos::RCP<MAT::PAR::InelasticDefgradPolyIntercalFrac>
      params_inelastic_defgrad_poly_intercal_frac_;
  // pointer to InelasticDefgradPolyIntercalFracAniso object
  Teuchos::RCP<MAT::InelasticDefgradPolyIntercalFracAniso> poly_intercal_frac_aniso_;
  // pointer to parameters of InelasticDefgradPolyIntercalFracAniso
  Teuchos::RCP<MAT::PAR::InelasticDefgradPolyIntercalFracAniso>
      params_inelastic_defgrad_poly_intercal_frac_aniso_;
};

#endif  // BACI_UNIT_INELASTIC_DEFGRAD_FACTORS_H