/*----------------------------------------------------------------------*/
/*! \file

\brief unit testing functionality for electrode material

\level 2

\maintainer Christoph Schmidt

*/
/*----------------------------------------------------------------------*/
#ifndef UNIT_ELECTRODE_H
#define UNIT_ELECTRODE_H

#include <cxxtest/TestSuite.h>

#include "src/drt_inpar/inpar_elch.H"

#include "src/drt_mat/electrode.H"
#include "src/drt_mat/matpar_material.H"

// forward declaration
namespace MAT
{
  class Electrode_TestSuite;
}

// class implementation
class MAT::Electrode_TestSuite : public CxxTest::TestSuite
{
 public:
  //! set up unit testing for electrode material
  void setUp()
  {
    // initialize container for material parameters
    const Teuchos::RCP<MAT::PAR::Material> container = Teuchos::rcp(new MAT::PAR::Material());

    // add dummy parameters to container
    container->Add("DIFFCOEF", 0);
    container->Add("DIFF_PARA_NUM", 0);
    container->Add("DIFF_PARA", std::vector<double>(0, 0.));
    container->Add("COND", 0);
    container->Add("COND_PARA_NUM", 0);
    container->Add("COND_PARA", std::vector<double>(0, 0.));

    // obtain half-cell open-circuit equilibrium potential from cubic spline interpolation of *.csv
    // data points
    container->Add("OCP_MODEL", "csv");

    // add cathode parameters to container according to master thesis by Alexander Rupp (2017)
    container->Add("C_MAX", 4793.3);
    container->Add("CHI_MAX", 1.0);
    container->Add("OCP_PARA_NUM", 0);
    container->Add("OCP_PARA", std::vector<double>(0, 0.));
    std::string ocpcsv(__FILE__);
    ocpcsv.replace(ocpcsv.end() - 1, ocpcsv.end(), "csv");
    container->Add("OCP_CSV", ocpcsv);
    container->Add("X_MIN", -1.);
    container->Add("X_MAX", -1.);

    // initialize parameter class for cathode material
    parameters_cathode_csv_ = Teuchos::rcp(new MAT::PAR::Electrode(container));

    // initialize cathode material
    cathode_csv_ = Teuchos::rcp(new MAT::Electrode(parameters_cathode_csv_.get()));

    // define sample concentration values for cathode material
    // cf. master thesis by Alexander Rupp (2017)
    concentrations_cathode_csv_.resize(3, 0.);
    concentrations_cathode_csv_[0] = 1677.6;  // cathode concentration at 100% state of charge
    concentrations_cathode_csv_[1] = 3115.6;  // cathode concentration at 50% intercalation fraction
    concentrations_cathode_csv_[2] = 4553.6;  // cathode concentration at 0% state of charge

    // choose semi-empirical Redlich-Kister expansion to model half-cell open-circuit equilibrium
    // potential
    container->Add("OCP_MODEL", "Redlich-Kister");

    // add anode parameters to container according to Goldin et al., Electrochimica Acta 64 (2012)
    // 118-129
    container->Add("C_MAX", 16.1);
    container->Add("OCP_PARA_NUM", 16);
    std::vector<double> ocp_para(16, 0.);
    ocp_para[0] = 1.1652e4;
    ocp_para[1] = -3.268e3;
    ocp_para[2] = 3.955e3;
    ocp_para[3] = -4.573e3;
    ocp_para[4] = 6.147e3;
    ocp_para[5] = -3.339e3;
    ocp_para[6] = 1.117e4;
    ocp_para[7] = 2.997e2;
    ocp_para[8] = -4.866e4;
    ocp_para[9] = 1.362e2;
    ocp_para[10] = 1.373e5;
    ocp_para[11] = -2.129e4;
    ocp_para[12] = -1.722e5;
    ocp_para[13] = 3.956e4;
    ocp_para[14] = 9.302e4;
    ocp_para[15] = -3.280e4;
    container->Add("OCP_PARA", ocp_para);
    container->Add("OCP_CSV", "");
    container->Add("X_MIN", -1.);
    container->Add("X_MAX", -1.);

    // initialize parameter class for anode material
    parameters_anode_redlichkister_ = Teuchos::rcp(new MAT::PAR::Electrode(container));

    // initialize anode material
    anode_redlichkister_ = Teuchos::rcp(new MAT::Electrode(parameters_anode_redlichkister_.get()));

    // define sample concentration values for anode material
    // cf. Goldin et al., Electrochimica Acta 64 (2012) 118-129
    concentrations_anode_redlichkister_.resize(3, 0.);
    concentrations_anode_redlichkister_[0] = 2.029;  // anode concentration at 0% state of charge
    concentrations_anode_redlichkister_[1] =
        8.05;  // anode concentration at 50% intercalation fraction (constitutes singularity in
               // original form of Redlich-Kister expansion)
    concentrations_anode_redlichkister_[2] = 10.88;  // anode concentration at 100% state of charge

    // add cathode parameters to container according to Goldin et al., Electrochimica Acta 64 (2012)
    // 118-129
    container->Add("C_MAX", 23.9);
    container->Add("OCP_PARA_NUM", 21);
    ocp_para.resize(21, 0.);
    ocp_para[0] = 3.954616e5;
    ocp_para[1] = -7.676e4;
    ocp_para[2] = 3.799e4;
    ocp_para[3] = -2.873e4;
    ocp_para[4] = 1.169e4;
    ocp_para[5] = 1.451e4;
    ocp_para[6] = -8.938e4;
    ocp_para[7] = 1.671e5;
    ocp_para[8] = -7.236e4;
    ocp_para[9] = -1.746e5;
    ocp_para[10] = -4.067e5;
    ocp_para[11] = 9.534e5;
    ocp_para[12] = 5.897e5;
    ocp_para[13] = -7.455e5;
    ocp_para[14] = -1.102e6;
    ocp_para[15] = -2.927e5;
    ocp_para[16] = 7.214e5;
    ocp_para[17] = 9.029e5;
    ocp_para[18] = -1.599e5;
    ocp_para[19] = 6.658e5;
    ocp_para[20] = -1.084e6;
    container->Add("OCP_PARA", ocp_para);
    container->Add("OCP_CSV", "");
    container->Add("X_MIN", -1.);
    container->Add("X_MAX", -1.);

    // initialize parameter class for cathode material
    parameters_cathode_redlichkister_ = Teuchos::rcp(new MAT::PAR::Electrode(container));

    // initialize cathode material
    cathode_redlichkister_ =
        Teuchos::rcp(new MAT::Electrode(parameters_cathode_redlichkister_.get()));

    // define sample concentration values for cathode material
    // cf. Goldin et al., Electrochimica Acta 64 (2012) 118-129
    concentrations_cathode_redlichkister_.resize(3, 0.);
    concentrations_cathode_redlichkister_[0] =
        10.56;  // cathode concentration at 100% state of charge
    concentrations_cathode_redlichkister_[1] =
        11.95;  // cathode concentration at 50% intercalation fraction (constitutes singularity in
                // original form of Redlich-Kister expansion)
    concentrations_cathode_redlichkister_[2] =
        22.37;  // cathode concentration at 0% state of charge

    // choose half-cell open-circuit equilibrium potential according to Taralov, Taralova, Popov,
    // Iliev, Latz, and Zausch (2012)
    container->Add("OCP_MODEL", "Taralov");

    // add anode parameters to container according to Taralov, Taralova, Popov, Iliev, Latz, and
    // Zausch (2012)
    container->Add("C_MAX", 24.681);
    container->Add("OCP_PARA_NUM", 13);
    ocp_para.resize(13, 0.);
    std::fill(ocp_para.begin(), ocp_para.end(), 0.);
    ocp_para[0] = -0.132;
    ocp_para[10] = 1.41;
    ocp_para[11] = -3.52;
    container->Add("OCP_PARA", ocp_para);
    container->Add("OCP_CSV", "");
    container->Add("X_MIN", -1.);
    container->Add("X_MAX", -1.);

    // initialize parameter class for anode material
    parameters_anode_taralov_ = Teuchos::rcp(new MAT::PAR::Electrode(container));

    // initialize anode material
    anode_taralov_ = Teuchos::rcp(new MAT::Electrode(parameters_anode_taralov_.get()));

    // define sample concentration values for anode material
    // cf. Taralov, Taralova, Popov, Iliev, Latz, and Zausch (2012)
    concentrations_anode_taralov_.resize(3, 0.);
    concentrations_anode_taralov_[0] = 2.4681;   // anode concentration at 10% state of charge
    concentrations_anode_taralov_[1] = 12.3405;  // anode concentration at 50% state of charge
    concentrations_anode_taralov_[2] = 22.2129;  // anode concentration at 90% state of charge

    // add cathode parameters to container according to Taralov, Taralova, Popov, Iliev, Latz, and
    // Zausch (2012)
    container->Add("C_MAX", 23.671);
    ocp_para[0] = 4.06279;
    ocp_para[1] = 0.0677504;
    ocp_para[2] = -21.8502;
    ocp_para[3] = 12.8268;
    ocp_para[4] = -0.045;
    ocp_para[5] = -71.69;
    ocp_para[6] = -0.105734;
    ocp_para[7] = 1.00167;
    ocp_para[8] = 0.379571;
    ocp_para[9] = -1.576;
    ocp_para[10] = 0.01;
    ocp_para[11] = -200.;
    ocp_para[12] = -0.19;
    container->Add("OCP_PARA", ocp_para);
    container->Add("OCP_CSV", "");
    container->Add("X_MIN", -1.);
    container->Add("X_MAX", -1.);

    // initialize parameter class for cathode material
    parameters_cathode_taralov_ = Teuchos::rcp(new MAT::PAR::Electrode(container));

    // initialize cathode material
    cathode_taralov_ = Teuchos::rcp(new MAT::Electrode(parameters_cathode_taralov_.get()));

    // define sample concentration values for cathode material
    // cf. Taralov, Taralova, Popov, Iliev, Latz, and Zausch (2012)
    concentrations_cathode_taralov_.resize(3, 0.);
    concentrations_cathode_taralov_[0] = 4.02407;  // cathode concentration at 90% state of charge
    concentrations_cathode_taralov_[1] = 11.8355;  // cathode concentration at 50% state of charge
    concentrations_cathode_taralov_[2] = 21.3039;  // cathode concentration at 17% state of charge

    // choose polynomial half-cell open-circuit equilibrium potential
    container->Add("OCP_MODEL", "Polynomial");

    // add parameters for lithium metal anode to container
    container->Add("C_MAX", 100000.0);
    container->Add("OCP_PARA_NUM", 1);
    ocp_para.resize(1, 0.);
    ocp_para[0] = 0.;
    container->Add("OCP_PARA", ocp_para);
    container->Add("OCP_CSV", "");
    container->Add("X_MIN", -1.);
    container->Add("X_MAX", -1.);

    // initialize parameter class for anode material
    parameters_anode_polynomial_ = Teuchos::rcp(new MAT::PAR::Electrode(container));

    // initialize anode material
    anode_polynomial_ = Teuchos::rcp(new MAT::Electrode(parameters_anode_polynomial_.get()));

    // define fictitious sample concentration values for anode material
    concentrations_anode_polynomial_.resize(3, 0.);
    concentrations_anode_polynomial_[0] = 500.;   // anode concentration at 10% state of charge
    concentrations_anode_polynomial_[1] = 2500.;  // anode concentration at 50% state of charge
    concentrations_anode_polynomial_[2] = 4500.;  // anode concentration at 90% state of charge

    // add cathode parameters to container according to Ji et al., Journal of The Electrochemical
    // Society 160 (4) (2013) A636-A649
    container->Add("C_MAX", 4793.3);
    container->Add("OCP_PARA_NUM", 5);
    ocp_para.resize(5, 0.);
    ocp_para[0] = 4.563;
    ocp_para[1] = 2.595;
    ocp_para[2] = -16.77;
    ocp_para[3] = 23.88;
    ocp_para[4] = -10.72;
    container->Add("OCP_PARA", ocp_para);
    container->Add("OCP_CSV", "");
    container->Add("X_MIN", 0.3);
    container->Add("X_MAX", 1.0);

    // initialize parameter class for cathode material
    parameters_cathode_polynomial_ = Teuchos::rcp(new MAT::PAR::Electrode(container));

    // initialize cathode material
    cathode_polynomial_ = Teuchos::rcp(new MAT::Electrode(parameters_cathode_polynomial_.get()));

    // define sample concentration values for cathode material
    // cf. Ji et al., Journal of The Electrochemical Society 160 (4) (2013) A636-A649
    concentrations_cathode_polynomial_.resize(3, 0.);
    concentrations_cathode_polynomial_[0] =
        1677.6;  // cathode concentration at 100% state of charge
    concentrations_cathode_polynomial_[1] = 3115.6;  // cathode concentration at 50% state of charge
    concentrations_cathode_polynomial_[2] = 4553.6;  // cathode concentration at 0% state of charge

    // as we do not have access to data defined in the input file here, we define the constants here
    // using the SI unit system
    faraday_ = 9.64853399e4;
    const double gasconstant = 8.314472;

    // calculate factor F/(RT) for temperature 298.0 K
    frt_ = faraday_ / (gasconstant * 298.);

    return;
  }

  //! tear down unit testing for electrode material
  void tearDown()
  {
    // free memory by resetting pointers and clearing vectors
    anode_redlichkister_ = Teuchos::null;
    cathode_redlichkister_ = Teuchos::null;
    anode_taralov_ = Teuchos::null;
    cathode_taralov_ = Teuchos::null;
    anode_polynomial_ = Teuchos::null;
    cathode_polynomial_ = Teuchos::null;
    parameters_anode_redlichkister_ = Teuchos::null;
    parameters_cathode_redlichkister_ = Teuchos::null;
    parameters_anode_taralov_ = Teuchos::null;
    parameters_cathode_taralov_ = Teuchos::null;
    parameters_anode_polynomial_ = Teuchos::null;
    parameters_cathode_polynomial_ = Teuchos::null;
    concentrations_anode_redlichkister_.clear();
    concentrations_cathode_redlichkister_.clear();
    concentrations_anode_taralov_.clear();
    concentrations_cathode_taralov_.clear();
    concentrations_anode_polynomial_.clear();
    concentrations_cathode_polynomial_.clear();
    frt_ = 0.;

    return;
  }

  //! test member function ComputeOpenCircuitPotential of electrode material
  void test_ComputeOpenCircuitPotential()
  {
    // define results and tolerances for anode materials
    const double results_anode_redlichkister[] = {
        1.78149102529067321354e-01, 1.00269118707846299765e-01, 8.08307883078566852264e-02};
    const double tolerances_anode_redlichkister[] = {1.e-17, 1.e-17, 1.e-18};
    const double results_anode_taralov[] = {
        8.59624971986640673549e-01, 1.10583257990501226953e-01, -7.26563582809017227682e-02};
    const double tolerances_anode_taralov[] = {1.e-17, 1.e-17, 1.e-18};
    const double results_anode_polynomial[] = {0., 0., 0.};
    const double tolerances_anode_polynomial[] = {1.e-20, 1.e-20, 1.e-20};

    // test member function using sample concentration values for anode materials
    for (unsigned i = 0; i < concentrations_anode_redlichkister_.size(); ++i)
      TS_ASSERT_DELTA(anode_redlichkister_->ComputeOpenCircuitPotential(
                          concentrations_anode_redlichkister_[i], faraday_, frt_),
          results_anode_redlichkister[i], tolerances_anode_redlichkister[i]);
    for (unsigned i = 0; i < concentrations_anode_taralov_.size(); ++i)
      TS_ASSERT_DELTA(anode_taralov_->ComputeOpenCircuitPotential(
                          concentrations_anode_taralov_[i], faraday_, frt_),
          results_anode_taralov[i], tolerances_anode_taralov[i]);
    for (unsigned i = 0; i < concentrations_anode_polynomial_.size(); ++i)
      TS_ASSERT_DELTA(anode_polynomial_->ComputeOpenCircuitPotential(
                          concentrations_anode_polynomial_[i], faraday_, frt_),
          results_anode_polynomial[i], tolerances_anode_polynomial[i]);

    // define results and tolerances for cathode materials
    const double results_cathode_csv[] = {
        4.26568197738244947459e+00, 3.85982096426498033637e+00, 3.49646582304368758187e+00};
    const double tolerances_cathode_csv[] = {1.e-16, 1.e-16, 1.e-16};
    const double results_cathode_redlichkister[] = {
        3.97322289583641552468e+00, 3.90180104449214848472e+00, 3.55923291247763096123e+00};
    const double tolerances_cathode_redlichkister[] = {1.e-16, 1.e-16, 1.e-16};
    const double results_cathode_taralov[] = {
        4.68476462172226959524e+00, 4.12283197556701885844e+00, 3.90987668507063856893e+00};
    const double tolerances_cathode_taralov[] = {1.e-16, 1.e-16, 1.e-16};
    const double results_cathode_polynomial[] = {
        4.27993831912628497349e+00, 3.80888970352517386431e+00, 3.63594305532142847426e+00};
    const double tolerances_cathode_polynomial[] = {1.e-16, 1.e-16, 1.e-16};

    // test member function using sample concentration values for cathode materials
    for (unsigned i = 0; i < concentrations_cathode_csv_.size(); ++i)
      TS_ASSERT_DELTA(
          cathode_csv_->ComputeOpenCircuitPotential(concentrations_cathode_csv_[i], faraday_, frt_),
          results_cathode_csv[i], tolerances_cathode_csv[i]);
    for (unsigned i = 0; i < concentrations_cathode_redlichkister_.size(); ++i)
      TS_ASSERT_DELTA(cathode_redlichkister_->ComputeOpenCircuitPotential(
                          concentrations_cathode_redlichkister_[i], faraday_, frt_),
          results_cathode_redlichkister[i], tolerances_cathode_redlichkister[i]);
    for (unsigned i = 0; i < concentrations_cathode_taralov_.size(); ++i)
      TS_ASSERT_DELTA(cathode_taralov_->ComputeOpenCircuitPotential(
                          concentrations_cathode_taralov_[i], faraday_, frt_),
          results_cathode_taralov[i], tolerances_cathode_taralov[i]);
    for (unsigned i = 0; i < concentrations_cathode_polynomial_.size(); ++i)
      TS_ASSERT_DELTA(cathode_polynomial_->ComputeOpenCircuitPotential(
                          concentrations_cathode_polynomial_[i], faraday_, frt_),
          results_cathode_polynomial[i], tolerances_cathode_polynomial[i]);

    return;
  }

  //! test member function ComputeFirstDerivOpenCircuitPotential of electrode material
  void test_ComputeFirstDerivOpenCircuitPotential()
  {
    // define results and tolerances for anode materials
    const double results_anode_redlichkister[] = {
        -4.67265915814377499893e-02, -4.69987444027383099998e-03, -8.94809228199507436519e-03};
    const double tolerances_anode_redlichkister[] = {1.e-18, 1.e-19, 1.e-19};
    const double results_anode_taralov[] = {
        -1.41425383954984579260e-01, -3.45971827773009360518e-02, -8.46358003529945959742e-03};
    const double tolerances_anode_taralov[] = {1.e-17, 1.e-18, 1.e-19};
    const double results_anode_polynomial[] = {0., 0., 0.};
    const double tolerances_anode_polynomial[] = {1.e-20, 1.e-20, 1.e-20};

    // test member function using sample concentration values for anode materials
    for (unsigned i = 0; i < concentrations_anode_redlichkister_.size(); ++i)
      TS_ASSERT_DELTA(anode_redlichkister_->ComputeFirstDerivOpenCircuitPotential(
                          concentrations_anode_redlichkister_[i], faraday_, frt_),
          results_anode_redlichkister[i], tolerances_anode_redlichkister[i]);
    for (unsigned i = 0; i < concentrations_anode_taralov_.size(); ++i)
      TS_ASSERT_DELTA(anode_taralov_->ComputeFirstDerivOpenCircuitPotential(
                          concentrations_anode_taralov_[i], faraday_, frt_),
          results_anode_taralov[i], tolerances_anode_taralov[i]);
    for (unsigned i = 0; i < concentrations_anode_polynomial_.size(); ++i)
      TS_ASSERT_DELTA(anode_polynomial_->ComputeFirstDerivOpenCircuitPotential(
                          concentrations_anode_polynomial_[i], faraday_, frt_),
          results_anode_polynomial[i], tolerances_anode_polynomial[i]);

    // define results and tolerances for cathode materials
    const double results_cathode_csv[] = {
        -3.97916196311376825587e-04, -1.82254204470957764054e-04, -1.19179495232888516508e-03};
    const double tolerances_cathode_csv[] = {1.e-16, 1.e-16, 1.e-16};
    const double results_cathode_redlichkister[] = {
        -5.87374545823974372749e-02, -4.59544093325152322449e-02, -2.40034452721541273490e-02};
    const double tolerances_cathode_redlichkister[] = {1.e-18, 1.e-18, 1.e-18};
    const double results_cathode_taralov[] = {
        -4.61526542535566619563e+00, -3.29033139167962890639e-03, -3.97138702703261764482e-02};
    const double tolerances_cathode_taralov[] = {1.e-16, 1.e-19, 1.e-18};
    const double results_cathode_polynomial[] = {
        -4.60348767050483748987e-04, -1.48972359281208444927e-04, -2.87284599311867232425e-04};
    const double tolerances_cathode_polynomial[] = {1.e-20, 1.e-20, 1.e-20};

    // test member function using sample concentration values for cathode materials
    for (unsigned i = 0; i < concentrations_cathode_csv_.size(); ++i)
      TS_ASSERT_DELTA(cathode_csv_->ComputeFirstDerivOpenCircuitPotential(
                          concentrations_cathode_csv_[i], faraday_, frt_),
          results_cathode_csv[i], tolerances_cathode_csv[i]);
    for (unsigned i = 0; i < concentrations_cathode_redlichkister_.size(); ++i)
      TS_ASSERT_DELTA(cathode_redlichkister_->ComputeFirstDerivOpenCircuitPotential(
                          concentrations_cathode_redlichkister_[i], faraday_, frt_),
          results_cathode_redlichkister[i], tolerances_cathode_redlichkister[i]);
    for (unsigned i = 0; i < concentrations_cathode_taralov_.size(); ++i)
      TS_ASSERT_DELTA(cathode_taralov_->ComputeFirstDerivOpenCircuitPotential(
                          concentrations_cathode_taralov_[i], faraday_, frt_),
          results_cathode_taralov[i], tolerances_cathode_taralov[i]);
    for (unsigned i = 0; i < concentrations_cathode_polynomial_.size(); ++i)
      TS_ASSERT_DELTA(cathode_polynomial_->ComputeFirstDerivOpenCircuitPotential(
                          concentrations_cathode_polynomial_[i], faraday_, frt_),
          results_cathode_polynomial[i], tolerances_cathode_polynomial[i]);

    return;
  }

  //! test member function ComputeSecondDerivOpenCircuitPotential of electrode material
  void test_ComputeSecondDerivOpenCircuitPotential()
  {
    // define results and tolerances for anode materials
    const double results_anode_redlichkister[] = {
        3.04927774052812586292e-02, -1.05174077799887128251e-03, -7.35165984480485711766e-04};
    const double tolerances_anode_redlichkister[] = {1.e-18, 1.e-19, 1.e-20};
    const double results_anode_taralov[] = {
        2.01700640785035384406e-02, 4.93424429221260389677e-03, 1.20707433751687906627e-03};
    const double tolerances_anode_taralov[] = {1.e-18, 1.e-19, 1.e-19};
    const double results_anode_polynomial[] = {0., 0., 0.};
    const double tolerances_anode_polynomial[] = {1.e-20, 1.e-20, 1.e-20};

    // test member function using sample concentration values for anode materials
    for (unsigned i = 0; i < concentrations_anode_redlichkister_.size(); ++i)
      TS_ASSERT_DELTA(anode_redlichkister_->ComputeSecondDerivOpenCircuitPotential(
                          concentrations_anode_redlichkister_[i], faraday_, frt_),
          results_anode_redlichkister[i], tolerances_anode_redlichkister[i]);
    for (unsigned i = 0; i < concentrations_anode_taralov_.size(); ++i)
      TS_ASSERT_DELTA(anode_taralov_->ComputeSecondDerivOpenCircuitPotential(
                          concentrations_anode_taralov_[i], faraday_, frt_),
          results_anode_taralov[i], tolerances_anode_taralov[i]);
    for (unsigned i = 0; i < concentrations_anode_polynomial_.size(); ++i)
      TS_ASSERT_DELTA(anode_polynomial_->ComputeSecondDerivOpenCircuitPotential(
                          concentrations_anode_polynomial_[i], faraday_, frt_),
          results_anode_polynomial[i], tolerances_anode_polynomial[i]);

    // define results and tolerances for cathode materials
    const double results_cathode_csv[] = {
        1.63562986543888766237e-06, 1.43228888999912318160e-07, -1.10872570647637550212e-05};
    const double tolerances_cathode_csv[] = {1.e-16, 1.e-16, 1.e-16};
    const double results_cathode_redlichkister[] = {
        1.41256487214221678611e-02, 5.72637600860917590773e-03, -7.56499159049272597299e-03};
    const double tolerances_cathode_redlichkister[] = {1.e-18, 1.e-19, 1.e-19};
    const double results_cathode_taralov[] = {
        3.89765224056106092121e+01, -7.34373323389212427637e-03, -2.27647861985425976894e-02};
    const double tolerances_cathode_taralov[] = {1.e-15, 1.e-19, 1.e-18};
    const double results_cathode_polynomial[] = {
        3.69515732169165666976e-08, 2.28146222921396267486e-07, -5.88484654599405069383e-07};
    const double tolerances_cathode_polynomial[] = {1.e-24, 1.e-23, 1.e-23};

    // test member function using sample concentration values for cathode materials
    for (unsigned i = 0; i < concentrations_cathode_csv_.size(); ++i)
      TS_ASSERT_DELTA(cathode_csv_->ComputeSecondDerivOpenCircuitPotential(
                          concentrations_cathode_csv_[i], faraday_, frt_),
          results_cathode_csv[i], tolerances_cathode_csv[i]);
    for (unsigned i = 0; i < concentrations_cathode_redlichkister_.size(); ++i)
      TS_ASSERT_DELTA(cathode_redlichkister_->ComputeSecondDerivOpenCircuitPotential(
                          concentrations_cathode_redlichkister_[i], faraday_, frt_),
          results_cathode_redlichkister[i], tolerances_cathode_redlichkister[i]);
    for (unsigned i = 0; i < concentrations_cathode_taralov_.size(); ++i)
      TS_ASSERT_DELTA(cathode_taralov_->ComputeSecondDerivOpenCircuitPotential(
                          concentrations_cathode_taralov_[i], faraday_, frt_),
          results_cathode_taralov[i], tolerances_cathode_taralov[i]);
    for (unsigned i = 0; i < concentrations_cathode_polynomial_.size(); ++i)
      TS_ASSERT_DELTA(cathode_polynomial_->ComputeSecondDerivOpenCircuitPotential(
                          concentrations_cathode_polynomial_[i], faraday_, frt_),
          results_cathode_polynomial[i], tolerances_cathode_polynomial[i]);

    return;
  }

 private:
  //! cathode material based on half cell open circuit potential obtained from cubic spline
  //! interpolation of *.csv data points
  Teuchos::RCP<const MAT::Electrode> cathode_csv_;

  //! anode material based on half cell open circuit potential according to Redlich-Kister expansion
  Teuchos::RCP<const MAT::Electrode> anode_redlichkister_;

  //! cathode material based on half cell open circuit potential according to Redlich-Kister
  //! expansion
  Teuchos::RCP<const MAT::Electrode> cathode_redlichkister_;

  //! anode material based on half cell open circuit potential according to Taralov, Taralova,
  //! Popov, Iliev, Latz, and Zausch (2012)
  Teuchos::RCP<const MAT::Electrode> anode_taralov_;

  //! cathode material based on half cell open circuit potential according to Taralov, Taralova,
  //! Popov, Iliev, Latz, and Zausch (2012)
  Teuchos::RCP<const MAT::Electrode> cathode_taralov_;

  //! anode material based on polynomial half cell open circuit potential
  Teuchos::RCP<const MAT::Electrode> anode_polynomial_;

  //! cathode material based on polynomial half cell open circuit potential
  Teuchos::RCP<const MAT::Electrode> cathode_polynomial_;

  //! parameters for cathode material based on half cell open circuit potential obtained from cubic
  //! spline interpolation of *.csv data points
  Teuchos::RCP<MAT::PAR::Electrode> parameters_cathode_csv_;

  //! parameters for anode material based on half cell open circuit potential according to
  //! Redlich-Kister expansion
  Teuchos::RCP<MAT::PAR::Electrode> parameters_anode_redlichkister_;

  //! parameters for cathode material based on half cell open circuit potential according to
  //! Redlich-Kister expansion
  Teuchos::RCP<MAT::PAR::Electrode> parameters_cathode_redlichkister_;

  //! parameters for anode material based on half cell open circuit potential according to Taralov,
  //! Taralova, Popov, Iliev, Latz, and Zausch (2012)
  Teuchos::RCP<MAT::PAR::Electrode> parameters_anode_taralov_;

  //! parameters for cathode material based on half cell open circuit potential according to
  //! Taralov, Taralova, Popov, Iliev, Latz, and Zausch (2012)
  Teuchos::RCP<MAT::PAR::Electrode> parameters_cathode_taralov_;

  //! parameters for anode material based on polynomial half cell open circuit potential
  Teuchos::RCP<MAT::PAR::Electrode> parameters_anode_polynomial_;

  //! parameters for cathode material based on polynomial half cell open circuit potential
  Teuchos::RCP<MAT::PAR::Electrode> parameters_cathode_polynomial_;

  //! sample concentration values for cathode material based on half cell open circuit potential
  //! obtained from cubic spline interpolation of *.csv data points
  std::vector<double> concentrations_cathode_csv_;

  //! sample concentration values for anode material based on half cell open circuit potential
  //! according to Redlich-Kister expansion
  std::vector<double> concentrations_anode_redlichkister_;

  //! sample concentration values for cathode material based on half cell open circuit potential
  //! according to Redlich-Kister expansion
  std::vector<double> concentrations_cathode_redlichkister_;

  //! sample concentration values for anode material based on half cell open circuit potential
  //! according to Taralov, Taralova, Popov, Iliev, Latz, and Zausch (2012)
  std::vector<double> concentrations_anode_taralov_;

  //! sample concentration values for cathode material based on half cell open circuit potential
  //! according to Taralov, Taralova, Popov, Iliev, Latz, and Zausch (2012)
  std::vector<double> concentrations_cathode_taralov_;

  //! sample concentration values for anode material based on polynomial half cell open circuit
  //! potential
  std::vector<double> concentrations_anode_polynomial_;

  //! sample concentration values for cathode material based on polynomial half cell open circuit
  //! potential
  std::vector<double> concentrations_cathode_polynomial_;

  //! faraday constant
  double faraday_;

  //! factor F/(RT)
  double frt_;
};
#endif
