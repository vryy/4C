/*----------------------------------------------------------------------*/
/*! \file
\brief Testcases for the multiplicative_split_defgrad_elasthyper functions
\level 3

*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mat_elasthyper_service.hpp"
#include "4C_mat_multiplicative_split_defgrad_elasthyper.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_par_material.hpp"
#include "4C_unittest_utils_assertions_test.hpp"

namespace
{
  using namespace FourC;

  class MultiplicativeSplitDefgradElastHyperTest : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      // set up the deformation gradient
      FM_(0, 0) = 1.1;
      FM_(0, 1) = 0.01;
      FM_(0, 2) = 0.03;
      FM_(1, 0) = 0.04;
      FM_(1, 1) = 1.2;
      FM_(1, 2) = 0.02;
      FM_(2, 0) = 0.06;
      FM_(2, 1) = 0.05;
      FM_(2, 2) = 1.3;

      // set up the inverse inelastic deformation gradient
      iFinM_(0, 0) = 0.988;
      iFinM_(0, 1) = 0.006;
      iFinM_(0, 2) = 0.008;
      iFinM_(1, 0) = -0.005;
      iFinM_(1, 1) = 0.997;
      iFinM_(1, 2) = -0.004;
      iFinM_(2, 0) = -0.003;
      iFinM_(2, 1) = -0.002;
      iFinM_(2, 2) = 0.999;

      // set corresponding determinant of inelastic deformation gradient
      detFin_ = 1.016159878774565;

      SetRefValuesEvaluateKinQuantElast();

      SetRefValuesEvaluateInvariantsDerivative();

      SetRefValuesEvaluatedSdiFin();

      // this method is tested in unit_elasthyper_service.H
      MAT::CalculateGammaDelta(gamma_ref_, delta_ref_, prinv_ref_, dPIe_ref_, ddPIIe_ref_);

      SetUpMultiplicativeSplitObjectWithRequirements();
    }

    void TearDown() override
    {
      // We need to make sure the GLOBAL::Problem instance created in SetUp is deleted again. If
      // this is not done, some troubles arise where unit tests influence each other on some
      // configurations. We suspect that missing singleton destruction might be the reason for that.
      GLOBAL::Problem::Done();
    }

    // note: requirement of the generation of an object of the multiplicative split material is that
    // the global problem is equipped with the corresponding elastic and inelastic material. This is
    // set up within this method.
    void SetUpMultiplicativeSplitObjectWithRequirements()
    {
      // do problem instance specific stuff
      const int problemid(0);
      GLOBAL::Problem& problem = (*GLOBAL::Problem::Instance());
      problem.Materials()->SetReadFromProblem(problemid);

      // set up elastic material to be added to problem instance
      const int matid_elastic(1);
      const auto mat_elastic_neo_hooke = Teuchos::rcp(new MAT::PAR::Material(
          matid_elastic, CORE::Materials::mes_coupneohooke, "ELAST_CoupNeoHooke"));

      // add actually required parameters to electrode material
      mat_elastic_neo_hooke->Add("YOUNG", 1.5e2);
      mat_elastic_neo_hooke->Add("NUE", 0.3);

      // add elastic material to problem instance
      problem.Materials()->Insert(matid_elastic, mat_elastic_neo_hooke);

      // set up inelastic material to be added to problem instance
      const int inelastic_defgrad_id(2);
      const auto mat_inelastic = Teuchos::rcp(new MAT::PAR::Material(inelastic_defgrad_id,
          CORE::Materials::mfi_lin_scalar_iso, "MAT_InelasticDefgradLinScalarIso"));

      mat_inelastic->Add("SCALAR1", 1);
      mat_inelastic->Add("SCALAR1_MolarGrowthFac", 1.1);
      mat_inelastic->Add("SCALAR1_RefConc", 1.2);

      // add inelastic material to problem instance
      problem.Materials()->Insert(inelastic_defgrad_id, mat_inelastic);

      // set parameter list
      auto parameter_list_pointer = Teuchos::rcp(new Teuchos::ParameterList());
      auto structural_dynamic_params = parameter_list_pointer->sublist("STRUCTURAL DYNAMIC", false);

      structural_dynamic_params.set("MASSLIN", "No");

      // set the parameter list in the global problem
      problem.setParameterList(parameter_list_pointer);

      // create MultiplicativeSplitDefgrad_ElastHyper object;
      // initialize container for material parameters first
      auto multiplicativeSplitDefgradData = Teuchos::rcp(new MAT::PAR::Material());

      multiplicativeSplitDefgradData->Add("NUMMATEL", 1);
      std::vector<int> matids_elastic = {matid_elastic};
      multiplicativeSplitDefgradData->Add("MATIDSEL", matids_elastic);
      multiplicativeSplitDefgradData->Add("NUMFACINEL", 1);
      std::vector<int> inelastic_defgrad_factor_ids = {inelastic_defgrad_id};
      multiplicativeSplitDefgradData->Add("INELDEFGRADFACIDS", inelastic_defgrad_factor_ids);
      multiplicativeSplitDefgradData->Add("DENS", 1.32e1);

      // get pointer to parameter class
      parameters_multiplicative_split_defgrad_ = Teuchos::rcp(
          new MAT::PAR::MultiplicativeSplitDefgradElastHyper(multiplicativeSplitDefgradData));

      // setup pointer to MultiplicativeSplitDefgrad_ElastHyper object
      multiplicative_split_defgrad_ = Teuchos::rcp(new MAT::MultiplicativeSplitDefgradElastHyper(
          parameters_multiplicative_split_defgrad_.get()));
    }

    void SetRefValuesEvaluateKinQuantElast()
    {
      iCinV_ref_(0) = 0.976244;
      iCinV_ref_(1) = 0.99405;
      iCinV_ref_(2) = 0.998014;
      iCinV_ref_(3) = 0.00101;
      iCinV_ref_(4) = -0.005975;
      iCinV_ref_(5) = 0.005016;

      iCinCM_ref_(0, 0) = 1.1869551176;
      iCinCM_ref_(0, 1) = 0.0624320828;
      iCinCM_ref_(0, 2) = 0.117717833;
      iCinCM_ref_(1, 0) = 0.062190447;
      iCinCM_ref_(1, 1) = 1.4335455825;
      iCinCM_ref_(1, 2) = 0.0787760655;
      iCinCM_ref_(2, 0) = 0.1173029584;
      iCinCM_ref_(2, 1) = 0.0808141072;
      iCinCM_ref_(2, 2) = 1.6879682995;

      iCinCiCinV_ref_(0) = 1.159411340880251;
      iCinCiCinV_ref_(1) = 1.424608111644232;
      iCinCiCinV_ref_(2) = 1.684721521806007;
      iCinCiCinV_ref_(3) = 0.062556072523941;
      iCinCiCinV_ref_(4) = 0.070366128660632;
      iCinCiCinV_ref_(5) = 0.123064780558814;

      iCV_ref_(0) = 0.82944796303475;
      iCV_ref_(1) = 0.696728380339461;
      iCV_ref_(2) = 0.596601975917312;
      iCV_ref_(3) = -0.032359712004713;
      iCV_ref_(4) = -0.034647920866899;
      iCV_ref_(5) = -0.053120416239144;

      iFinCeM_ref_(0, 0) = 1.1720463422758;
      iFinCeM_ref_(0, 1) = 0.0691310815912;
      iFinCeM_ref_(0, 2) = 0.1268460277766;
      iFinCeM_ref_(1, 0) = 0.054040105527;
      iFinCeM_ref_(1, 1) = 1.4294605363035;
      iFinCeM_ref_(1, 2) = 0.0734606306805;
      iFinCeM_ref_(2, 0) = 0.1104273474647;
      iFinCeM_ref_(2, 1) = 0.0778995460298;
      iFinCeM_ref_(2, 2) = 1.6868954984389;

      CiFin9x1_ref_(0) = 1.1999722;
      CiFin9x1_ref_(1) = 1.4384656;
      CiFin9x1_ref_(2) = 1.6901459;
      CiFin9x1_ref_(3) = 0.0688816;
      CiFin9x1_ref_(4) = 0.0839363;
      CiFin9x1_ref_(5) = 0.1211618;
      CiFin9x1_ref_(6) = 0.0537751;
      CiFin9x1_ref_(7) = 0.0863203;
      CiFin9x1_ref_(8) = 0.104938;

      CiFinCe9x1_ref_(0) = 1.43996697912278;
      CiFinCe9x1_ref_(1) = 2.073382326190544;
      CiFinCe9x1_ref_(2) = 2.873787776734904;
      CiFinCe9x1_ref_(3) = 0.181343812846575;
      CiFinCe9x1_ref_(4) = 0.264478527552432;
      CiFinCe9x1_ref_(5) = 0.347292768781784;
      CiFinCe9x1_ref_(6) = 0.160486291582948;
      CiFinCe9x1_ref_(7) = 0.267131183013999;
      CiFinCe9x1_ref_(8) = 0.322626335257043;

      CiFiniCe9x1_ref_(0) = 1.012090158460073;
      CiFiniCe9x1_ref_(1) = 1.002986382106132;
      CiFiniCe9x1_ref_(2) = 1.000984547144946;
      CiFiniCe9x1_ref_(3) = 0.005087912513024;
      CiFiniCe9x1_ref_(4) = 0.001989641042641;
      CiFiniCe9x1_ref_(5) = 0.003049495796202;
      CiFiniCe9x1_ref_(6) = -0.006107120871435;
      CiFiniCe9x1_ref_(7) = 0.003975217445766;
      CiFiniCe9x1_ref_(8) = -0.008129279030197;

      prinv_ref_(0) = 4.3084689996;
      prinv_ref_(1) = 6.100818829441683;
      prinv_ref_(2) = 2.839432625034153;
    }

    void SetRefValuesEvaluateInvariantsDerivative()
    {
      dPIe_ref_.Clear();
      dPIe_ref_(0) = 28.846153846153847;
      dPIe_ref_(2) = -4.644432669839725;

      ddPIIe_ref_.Clear();
      ddPIIe_ref_(2) = 2.862458189907485;
    }

    void SetConcentrationToInelasticMaterial(const double concentration) const
    {
      Teuchos::ParameterList params;

      // set up a concentration vector and store it to the parameter list
      auto gpconc_lin = Teuchos::rcp(
          new std::vector<std::vector<double>>(1, std::vector<double>(1, concentration)));
      params.set<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc", gpconc_lin);

      // call pre evaluate to have concentration during actual call
      multiplicative_split_defgrad_->PreEvaluate(params, 0);
    }

    void SetRefValuesEvaluatedSdiFin()
    {
      dSdiFin_ref_(0, 0) = 114.1660435509188;
      dSdiFin_ref_(0, 1) = -1.661105324599284;
      dSdiFin_ref_(0, 2) = -1.657789966811478;
      dSdiFin_ref_(0, 3) = 0.6950689065562556;
      dSdiFin_ref_(0, 4) = -0.00329516271498935;
      dSdiFin_ref_(0, 5) = 0.9329432831258889;
      dSdiFin_ref_(0, 6) = 0.01011436563696388;
      dSdiFin_ref_(0, 7) = -0.006583593739068971;
      dSdiFin_ref_(0, 8) = 0.013463381879464;
      dSdiFin_ref_(1, 0) = -11.73283847873412;
      dSdiFin_ref_(1, 1) = 105.2701678289236;
      dSdiFin_ref_(1, 2) = -11.6040946680382;
      dSdiFin_ref_(1, 3) = -0.05898254736522643;
      dSdiFin_ref_(1, 4) = -0.4920621412911008;
      dSdiFin_ref_(1, 5) = -0.03535183236330931;
      dSdiFin_ref_(1, 6) = -0.5154481790422178;
      dSdiFin_ref_(1, 7) = -0.04608342825898787;
      dSdiFin_ref_(1, 8) = 0.0942401395889132;
      dSdiFin_ref_(2, 0) = -18.75794747350803;
      dSdiFin_ref_(2, 1) = -18.58921926562027;
      dSdiFin_ref_(2, 2) = 98.5798500483063;
      dSdiFin_ref_(2, 3) = -0.09429870933073022;
      dSdiFin_ref_(2, 4) = -0.03687574852596564;
      dSdiFin_ref_(2, 5) = -0.05651895879799482;
      dSdiFin_ref_(2, 6) = 0.1131885845970979;
      dSdiFin_ref_(2, 7) = -0.3081745971147502;
      dSdiFin_ref_(2, 8) = -0.201080649217473;
      dSdiFin_ref_(3, 0) = -2.547471706027534;
      dSdiFin_ref_(3, 1) = -1.882323090966196;
      dSdiFin_ref_(3, 2) = -2.22961181644894;
      dSdiFin_ref_(3, 3) = 58.4374016535939;
      dSdiFin_ref_(3, 4) = 0.4645651032280521;
      dSdiFin_ref_(3, 5) = -0.241290937884049;
      dSdiFin_ref_(3, 6) = 57.93471620610067;
      dSdiFin_ref_(3, 7) = -0.008854474142796649;
      dSdiFin_ref_(3, 8) = 0.01810730908547461;
      dSdiFin_ref_(4, 0) = -1.995076152570434;
      dSdiFin_ref_(4, 1) = -2.094379609332357;
      dSdiFin_ref_(4, 2) = -2.207682721946282;
      dSdiFin_ref_(4, 3) = -0.01002951450149832;
      dSdiFin_ref_(4, 4) = 58.56206171546444;
      dSdiFin_ref_(4, 5) = -0.006011298785499145;
      dSdiFin_ref_(4, 6) = -0.1638352028130229;
      dSdiFin_ref_(4, 7) = 58.44089844409201;
      dSdiFin_ref_(4, 8) = -0.2770982534557008;
      dSdiFin_ref_(5, 0) = -4.075765888689157;
      dSdiFin_ref_(5, 1) = -3.864812436613845;
      dSdiFin_ref_(5, 2) = -3.388101888830967;
      dSdiFin_ref_(5, 3) = -0.1368544954849143;
      dSdiFin_ref_(5, 4) = -0.007666693768910916;
      dSdiFin_ref_(5, 5) = 58.5542331450641;
      dSdiFin_ref_(5, 6) = 0.02353259936218022;
      dSdiFin_ref_(5, 7) = 0.3364299251029703;
      dSdiFin_ref_(5, 8) = 57.95243768164728;
    }

    // defined input quantities
    CORE::LINALG::Matrix<3, 3> FM_;
    CORE::LINALG::Matrix<3, 3> iFinM_;
    double detFin_;
    // reference solutions
    CORE::LINALG::Matrix<6, 1> iCinV_ref_;
    CORE::LINALG::Matrix<6, 1> iCinCiCinV_ref_;
    CORE::LINALG::Matrix<6, 1> iCV_ref_;
    CORE::LINALG::Matrix<3, 3> iCinCM_ref_;
    CORE::LINALG::Matrix<3, 3> iFinCeM_ref_;
    CORE::LINALG::Matrix<9, 1> CiFin9x1_ref_;
    CORE::LINALG::Matrix<9, 1> CiFinCe9x1_ref_;
    CORE::LINALG::Matrix<9, 1> CiFiniCe9x1_ref_;
    CORE::LINALG::Matrix<3, 1> prinv_ref_;
    CORE::LINALG::Matrix<3, 1> dPIe_ref_;
    CORE::LINALG::Matrix<6, 1> ddPIIe_ref_;
    CORE::LINALG::Matrix<6, 9> dSdiFin_ref_;
    CORE::LINALG::Matrix<3, 1> gamma_ref_;
    CORE::LINALG::Matrix<8, 1> delta_ref_;

    // pointer to material parameters
    Teuchos::RCP<MAT::PAR::MultiplicativeSplitDefgradElastHyper>
        parameters_multiplicative_split_defgrad_;

    // pointer to material
    Teuchos::RCP<MAT::MultiplicativeSplitDefgradElastHyper> multiplicative_split_defgrad_;
  };

  TEST_F(MultiplicativeSplitDefgradElastHyperTest, TestEvaluateAdditionalCmat)
  {
    const double concentration(44327.362);
    SetConcentrationToInelasticMaterial(concentration);

    // actual call that is tested
    CORE::LINALG::Matrix<6, 6> cMatAdd(true);
    // reference solution
    CORE::LINALG::Matrix<6, 6> cMatAdd_ref;
    cMatAdd_ref(0, 0) = -0.701208493301168;
    cMatAdd_ref(0, 1) = -0.5890084484992926;
    cMatAdd_ref(0, 2) = -0.5043624088277512;
    cMatAdd_ref(0, 3) = 0.02735663466513797;
    cMatAdd_ref(0, 4) = 0.029291067637571;
    cMatAdd_ref(0, 5) = 0.04490756345738416;
    cMatAdd_ref(1, 0) = -0.5183018322515506;
    cMatAdd_ref(1, 1) = -0.4353685971936842;
    cMatAdd_ref(1, 2) = -0.3728020454851408;
    cMatAdd_ref(1, 3) = 0.02022079596387239;
    cMatAdd_ref(1, 4) = 0.02165064195626716;
    cMatAdd_ref(1, 5) = 0.03319365444696364;
    cMatAdd_ref(2, 0) = -0.3873521175618817;
    cMatAdd_ref(2, 1) = -0.3253720854320112;
    cMatAdd_ref(2, 2) = -0.2786130643659098;
    cMatAdd_ref(2, 3) = 0.01511198233926229;
    cMatAdd_ref(2, 4) = 0.01618057565396367;
    cMatAdd_ref(2, 5) = 0.02480722918495967;
    cMatAdd_ref(3, 0) = 0.04212677142365208;
    cMatAdd_ref(3, 1) = 0.03538608632607147;
    cMatAdd_ref(3, 2) = 0.03030077375609256;
    cMatAdd_ref(3, 3) = -0.00164351502651245;
    cMatAdd_ref(3, 4) = -0.001759730697661004;
    cMatAdd_ref(3, 5) = -0.00269792890279474;
    cMatAdd_ref(4, 0) = 0.0398350977693272;
    cMatAdd_ref(4, 1) = 0.03346110230705895;
    cMatAdd_ref(4, 2) = 0.02865242799932495;
    cMatAdd_ref(4, 3) = -0.001554108695111724;
    cMatAdd_ref(4, 4) = -0.001664002296403264;
    cMatAdd_ref(4, 5) = -0.002551163024973279;
    cMatAdd_ref(5, 0) = 0.07166415111957784;
    cMatAdd_ref(5, 1) = 0.06019720363802254;
    cMatAdd_ref(5, 2) = 0.05154630075158284;
    cMatAdd_ref(5, 3) = -0.002795873152556784;
    cMatAdd_ref(5, 4) = -0.002993573976479339;
    cMatAdd_ref(5, 5) = -0.004589594171729136;

    multiplicative_split_defgrad_->EvaluateAdditionalCmat(&FM_, iCV_ref_, dSdiFin_ref_, cMatAdd);

    FOUR_C_EXPECT_NEAR(cMatAdd, cMatAdd_ref, 1.0e-10);
  }

  TEST_F(MultiplicativeSplitDefgradElastHyperTest, TestEvaluateCauchyNDirAndDerivatives)
  {
    // we need to first set this dummy concentration. otherwise the call below results in a segfault
    const double dummy_conc(-1.0);
    SetConcentrationToInelasticMaterial(dummy_conc);

    // input variables
    CORE::LINALG::Matrix<3, 1> n, dir;
    n(0) = 1.0 / std::sqrt(2.0);
    n(1) = 0.0;
    n(2) = -1.0 / std::sqrt(2.0);
    dir(0) = 1.0 / std::sqrt(3.0);
    dir(1) = -1.0 / std::sqrt(3.0);
    dir(2) = -1.0 / std::sqrt(3.0);
    const double concentration(1.0);

    // output variables
    double cauchy_n_dir(0.0);
    CORE::LINALG::Matrix<3, 1> d_cauchyndir_dn(true), d_cauchyndir_ddir(true);
    CORE::LINALG::Matrix<9, 1> d_cauchyndir_dF(true);

    multiplicative_split_defgrad_->EvaluateCauchyNDirAndDerivatives(FM_, n, dir, cauchy_n_dir,
        &d_cauchyndir_dn, &d_cauchyndir_ddir, &d_cauchyndir_dF, nullptr, nullptr, nullptr, 0, 0,
        &concentration, nullptr, nullptr, nullptr);

    const double cauchy_n_dir_ref(6.019860168755);
    CORE::LINALG::Matrix<3, 1> d_cauchyndir_dn_ref(true), d_cauchyndir_ddir_ref(true);
    CORE::LINALG::Matrix<9, 1> d_cauchyndir_dF_ref(true);
    d_cauchyndir_dn_ref(0) = -2.856437080521;
    d_cauchyndir_dn_ref(1) = -6.736850094992;
    d_cauchyndir_dn_ref(2) = -1.136980497476e+01;

    d_cauchyndir_ddir_ref(0) = -1.933304201727;
    d_cauchyndir_ddir_ref(1) = -8.793331859785e-01;
    d_cauchyndir_ddir_ref(2) = -1.148067468249e+01;

    d_cauchyndir_dF_ref(0) = 7.834365817988e+01;
    d_cauchyndir_dF_ref(1) = 4.257334930649e+01;
    d_cauchyndir_dF_ref(2) = 7.961614320095e+01;
    d_cauchyndir_dF_ref(3) = -2.192663678809e+01;
    d_cauchyndir_dF_ref(4) = 1.867810957230e+01;
    d_cauchyndir_dF_ref(5) = -4.292543557935e+01;
    d_cauchyndir_dF_ref(6) = -1.694113997401e+01;
    d_cauchyndir_dF_ref(7) = 1.982432031426e+01;
    d_cauchyndir_dF_ref(8) = -3.361951777558e+01;

    EXPECT_NEAR(cauchy_n_dir, cauchy_n_dir_ref, 1.0e-10);
    FOUR_C_EXPECT_NEAR(d_cauchyndir_dn, d_cauchyndir_dn_ref, 1.0e-10);
    FOUR_C_EXPECT_NEAR(d_cauchyndir_ddir, d_cauchyndir_ddir_ref, 1.0e-10);
    FOUR_C_EXPECT_NEAR(d_cauchyndir_dF, d_cauchyndir_dF_ref, 1.0e-10);
  }

  TEST_F(MultiplicativeSplitDefgradElastHyperTest, TestEvaluatedSdiFin)
  {
    CORE::LINALG::Matrix<6, 9> dSdiFin(true);
    multiplicative_split_defgrad_->EvaluatedSdiFin(gamma_ref_, delta_ref_, iFinM_, iCinCM_ref_,
        iCinV_ref_, CiFin9x1_ref_, CiFinCe9x1_ref_, iCinCiCinV_ref_, CiFiniCe9x1_ref_, iCV_ref_,
        iFinCeM_ref_, detFin_, dSdiFin);

    FOUR_C_EXPECT_NEAR(dSdiFin, dSdiFin_ref_, 1.0e-10);
  }

  TEST_F(MultiplicativeSplitDefgradElastHyperTest, TestEvaluateInvariantDerivatives)
  {
    // derivatives of principle invariants
    const int gp(0);
    const int eleGID(0);
    CORE::LINALG::Matrix<3, 1> dPIe(true);
    CORE::LINALG::Matrix<6, 1> ddPIIe(true);
    multiplicative_split_defgrad_->EvaluateInvariantDerivatives(
        prinv_ref_, gp, eleGID, dPIe, ddPIIe);

    FOUR_C_EXPECT_NEAR(dPIe, dPIe_ref_, 1.0e-10);
    FOUR_C_EXPECT_NEAR(ddPIIe, ddPIIe_ref_, 1.0e-10);
  }

  TEST_F(MultiplicativeSplitDefgradElastHyperTest, TestEvaluateKinQuantElast)
  {
    // variables of kinematic quantities
    CORE::LINALG::Matrix<6, 1> iCinV(true);
    CORE::LINALG::Matrix<6, 1> iCinCiCinV(true);
    CORE::LINALG::Matrix<6, 1> iCV(true);
    CORE::LINALG::Matrix<3, 3> iCinCM(true);
    CORE::LINALG::Matrix<3, 3> iFinCeM(true);
    CORE::LINALG::Matrix<9, 1> CiFin9x1(true);
    CORE::LINALG::Matrix<9, 1> CiFinCe9x1(true);
    CORE::LINALG::Matrix<9, 1> CiFiniCe9x1(true);
    CORE::LINALG::Matrix<3, 1> prinv(true);

    multiplicative_split_defgrad_->EvaluateKinQuantElast(&FM_, iFinM_, iCinV, iCinCiCinV, iCV,
        iCinCM, iFinCeM, CiFin9x1, CiFinCe9x1, CiFiniCe9x1, prinv);

    FOUR_C_EXPECT_NEAR(iCinV, iCinV_ref_, 1.0e-10);
    FOUR_C_EXPECT_NEAR(iCinCiCinV, iCinCiCinV_ref_, 1.0e-10);
    FOUR_C_EXPECT_NEAR(iCV, iCV_ref_, 1.0e-10);
    FOUR_C_EXPECT_NEAR(iCinCM, iCinCM_ref_, 1.0e-10);
    FOUR_C_EXPECT_NEAR(iFinCeM, iFinCeM_ref_, 1.0e-10);
    FOUR_C_EXPECT_NEAR(CiFin9x1, CiFin9x1_ref_, 1.0e-10);
    FOUR_C_EXPECT_NEAR(CiFinCe9x1, CiFinCe9x1_ref_, 1.0e-10);
    FOUR_C_EXPECT_NEAR(CiFiniCe9x1, CiFiniCe9x1_ref_, 1.0e-10);
    FOUR_C_EXPECT_NEAR(prinv, prinv_ref_, 1.0e-10);
  }

  TEST_F(MultiplicativeSplitDefgradElastHyperTest, TestEvaluateLinearizationOD)
  {
    // we need to first set this dummy concentration. otherwise the call below results in a segfault
    const double dummy_conc(-1.0);
    SetConcentrationToInelasticMaterial(dummy_conc);

    // actual material call
    const double concentration(1.0);
    CORE::LINALG::Matrix<9, 1> DFDx(true);
    multiplicative_split_defgrad_->EvaluateLinearizationOD(FM_, concentration, &DFDx);

    // define the reference solution
    CORE::LINALG::Matrix<9, 1> DFDx_ref;
    DFDx_ref(0) = 4.417109534556e-01;
    DFDx_ref(1) = 4.818664946788e-01;
    DFDx_ref(2) = 5.220220359020e-01;
    DFDx_ref(3) = 4.015554122323e-03;
    DFDx_ref(4) = 8.031108244647e-03;
    DFDx_ref(5) = 1.204666236697e-02;
    DFDx_ref(6) = 1.606221648929e-02;
    DFDx_ref(7) = 2.007777061162e-02;
    DFDx_ref(8) = 2.409332473394e-02;

    FOUR_C_EXPECT_NEAR(DFDx, DFDx_ref, 1.0e-10);
  }

  TEST_F(MultiplicativeSplitDefgradElastHyperTest, TestEvaluateOdStiffMat)
  {
    const double concentration(44327.362);
    SetConcentrationToInelasticMaterial(concentration);

    // do the actual call that is tested
    auto source(MAT::PAR::InelasticSource::concentration);
    CORE::LINALG::Matrix<6, 1> dSdx(true);
    // reference solution
    CORE::LINALG::Matrix<6, 1> dSdx_ref;
    dSdx_ref(0) = -1.907155639254611e-05;
    dSdx_ref(1) = -1.409683812529051e-05;
    dSdx_ref(2) = -1.05352513901749e-05;
    dSdx_ref(3) = 1.14576920347856e-06;
    dSdx_ref(4) = 1.08343997650926e-06;
    dSdx_ref(5) = 1.949130554546719e-06;

    multiplicative_split_defgrad_->EvaluateODStiffMat(source, &FM_, dSdiFin_ref_, dSdx);

    FOUR_C_EXPECT_NEAR(dSdx, dSdx_ref, 1.0e-10);
  }

  TEST_F(MultiplicativeSplitDefgradElastHyperTest, TestEvaluateStressCmatIso)
  {
    // second Piola-Kirchhoff stress
    CORE::LINALG::Matrix<6, 1> S(true);
    // reference solution
    CORE::LINALG::Matrix<6, 1> S_ref(true);
    S_ref(0) = 35.001617076265632;
    S_ref(1) = 39.602547633321855;
    S_ref(2) = 42.518455970246585;
    S_ref(3) = 0.926494039729434;
    S_ref(4) = 0.578328181405601;
    S_ref(5) = 1.717758619623368;

    // elasticity tensor: 2 \partial S / \partial C
    CORE::LINALG::Matrix<6, 6> cMatIso(true);
    // reference solution
    CORE::LINALG::Matrix<6, 6> cMatIso_ref;
    cMatIso_ref(0, 0) = 64.536084541141619;
    cMatIso_ref(0, 1) = 23.288856786802732;
    cMatIso_ref(0, 2) = 20.045220754159409;
    cMatIso_ref(0, 3) = -2.517781949843259;
    cMatIso_ref(0, 4) = -1.063209891248748;
    cMatIso_ref(0, 5) = -4.133090713403081;
    cMatIso_ref(1, 0) = 23.288856786802729;
    cMatIso_ref(1, 1) = 45.535627490549523;
    cMatIso_ref(1, 2) = 16.775090529489887;
    cMatIso_ref(1, 3) = -2.114912831353509;
    cMatIso_ref(1, 4) = -2.264461822480193;
    cMatIso_ref(1, 5) = -1.427796649278471;
    cMatIso_ref(2, 0) = 20.045220754159409;
    cMatIso_ref(2, 1) = 16.775090529489891;
    cMatIso_ref(2, 2) = 33.388253156810187;
    cMatIso_ref(2, 3) = -0.677477970551063;
    cMatIso_ref(2, 4) = -1.939037415158507;
    cMatIso_ref(2, 5) = -2.972832770894956;
    cMatIso_ref(3, 0) = -2.517781949843259;
    cMatIso_ref(3, 1) = -2.114912831353509;
    cMatIso_ref(3, 2) = -0.677477970551063;
    cMatIso_ref(3, 3) = 15.558647064320079;
    cMatIso_ref(3, 4) = -0.916806613818434;
    cMatIso_ref(3, 5) = -0.655057743068236;
    cMatIso_ref(4, 0) = -1.063209891248748;
    cMatIso_ref(4, 1) = -2.264461822480193;
    cMatIso_ref(4, 2) = -1.939037415158507;
    cMatIso_ref(4, 3) = -0.916806613818434;
    cMatIso_ref(4, 4) = 11.220930501841977;
    cMatIso_ref(4, 5) = -0.394102458936374;
    cMatIso_ref(5, 0) = -4.133090713403081;
    cMatIso_ref(5, 1) = -1.427796649278471;
    cMatIso_ref(5, 2) = -2.972832770894956;
    cMatIso_ref(5, 3) = -0.655057743068236;
    cMatIso_ref(5, 4) = -0.394102458936374;
    cMatIso_ref(5, 5) = 13.451712479073098;

    multiplicative_split_defgrad_->EvaluateStressCmatIso(
        iCV_ref_, iCinV_ref_, iCinCiCinV_ref_, gamma_ref_, delta_ref_, detFin_, S, cMatIso);

    FOUR_C_EXPECT_NEAR(S, S_ref, 1.0e-10);
    FOUR_C_EXPECT_NEAR(cMatIso, cMatIso_ref, 1.0e-10);
  }
}  // namespace
