/*----------------------------------------------------------------------*/
/*! \file

\brief Testcases for the IsoOgden material

\level 2


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_UNIT_ELAST_ISOOGDEN_H
#define BACI_UNIT_ELAST_ISOOGDEN_H

#include "src/common/unit_cxx_test_wrapper.H"
#include "src/drt_matelast/elast_isoogden.H"

// declaration
namespace MAT
{
  namespace ELASTIC
  {
    class IsoOgden_TestSuite;
  }  // namespace ELASTIC
}  // namespace MAT

// class implementation
class MAT::ELASTIC::IsoOgden_TestSuite : public BACICxxTestWrapper
{
 public:
  //! set up unit test for IsoOgden material
  void Setup() override
  {
    // initialize container for material parameters
    const Teuchos::RCP<MAT::PAR::Material> container = Teuchos::rcp(new MAT::PAR::Material());

    // add material parameters to container
    container->Add("ALPHA", -25.0);
    container->Add("MUE", 0.8);

    // initialize parameter class for IsoOgden material with container
    parameters_isoogden_ = Teuchos::rcp(new MAT::ELASTIC::PAR::IsoOgden(container));

    // initialize IsoOgden material with parameter class
    isoogden_ = Teuchos::rcp(new MAT::ELASTIC::IsoOgden(parameters_isoogden_.get()));
  }

  //! tear down unit testing for IsoOgden material
  void TearDown() override {}

  //! AddCoefficientsStretchesModified function unit testing for IsoOgden material
  void TestAddCoefficientsStretchesModified()
  {
    // define tolerances for test
    const double tolerances_modgamma[] = {1.e-9, 1.e-9, 1.e-9};
    const double tolerances_moddelta[] = {1.e-9, 1.e-9, 1.e-9, 1.e-9, 1.e-9, 1.e-9};

    // define correct reference values
    const double ref_modgamma[] = {-0.990546182, -21.175823681, -681.759084419};
    const double ref_moddelta[] = {28.615778605, 688.214269644, 25322.480278435, 0, 0, 0};

    // define modified principal strains
    LINALG::Matrix<3, 1> modstr(true);
    modstr(0) = 0.9;
    modstr(1) = 0.8;
    modstr(2) = 0.7;

    // initialize resulting coefficients
    LINALG::Matrix<3, 1> modgamma(true);
    LINALG::Matrix<6, 1> moddelta(true);

    // call AddCoefficientsStretchesModified function with test modified principal strains
    isoogden_->AddCoefficientsStretchesModified(modgamma, moddelta, modstr);

    // test member function results using reference stress values
    for (int i = 0; i < 3; ++i)
      TS_ASSERT_DELTA(modgamma(i), ref_modgamma[i], tolerances_modgamma[i]);

    for (int i = 0; i < 6; ++i)
      TS_ASSERT_DELTA(moddelta(i), ref_moddelta[i], tolerances_moddelta[i]);
  }

  //! SpecifyFormulation function unit testing for IsoOgden material
  void TestSpecifyFormulation()
  {
    // initizalize function inputs
    bool isoprinc = false;
    bool isomod = false;
    bool anisoprinc = false;
    bool anisomod = false;
    bool viscogeneral = false;

    // call SpecifyFormulation function
    isoogden_->SpecifyFormulation(isoprinc, isomod, anisoprinc, anisomod, viscogeneral);

    // test if function correctly sets isomod to true
    TS_ASSERT_EQUALS(isomod, true);
  }

 private:
  //! material parameters
  Teuchos::RCP<MAT::ELASTIC::PAR::IsoOgden> parameters_isoogden_;

  //! material class
  Teuchos::RCP<MAT::ELASTIC::IsoOgden> isoogden_;
};

#endif