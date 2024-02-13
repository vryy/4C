/*----------------------------------------------------------------------*/
/*! \file

\brief Testcases for the IsoOgden material

\level 2


*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "baci_mat_par_material.hpp"
#include "baci_matelast_isoogden.hpp"
#include "baci_unittest_utils_assertions_test.hpp"

namespace
{
  using namespace BACI;

  class IsoOgdenTest : public ::testing::Test
  {
   protected:
    IsoOgdenTest()
        : parameters_isoogden_(std::invoke(
              []()
              {
                // initialize container for material parameters
                const Teuchos::RCP<MAT::PAR::Material> container =
                    Teuchos::rcp(new MAT::PAR::Material());

                // add material parameters to container
                container->Add("ALPHA", -25.0);
                container->Add("MUE", 0.8);

                return MAT::ELASTIC::PAR::IsoOgden(container);
              })),
          isoogden_(&parameters_isoogden_)
    {
    }

    //! material parameters
    MAT::ELASTIC::PAR::IsoOgden parameters_isoogden_;

    //! material class
    MAT::ELASTIC::IsoOgden isoogden_;
  };

  TEST_F(IsoOgdenTest, AddCoefficientsStretchesModified)
  {
    // define correct reference values
    const std::array<double, 3> ref_modgamma = {-0.990546182, -21.175823681, -681.759084419};
    const std::array<double, 6> ref_moddelta = {
        28.615778605, 688.214269644, 25322.480278435, 0.0, 0.0, 0.0};

    // define modified principal strains
    CORE::LINALG::Matrix<3, 1> modstr(true);
    modstr(0) = 0.9;
    modstr(1) = 0.8;
    modstr(2) = 0.7;

    // initialize resulting coefficients
    CORE::LINALG::Matrix<3, 1> modgamma(true);
    CORE::LINALG::Matrix<6, 1> moddelta(true);

    // call AddCoefficientsStretchesModified function with test modified principal strains
    isoogden_.AddCoefficientsStretchesModified(modgamma, moddelta, modstr);

    // test member function results using reference stress values
    BACI_EXPECT_NEAR(modgamma, ref_modgamma, 1.e-9);

    BACI_EXPECT_NEAR(moddelta, ref_moddelta, 1.e-9);
  }

  TEST_F(IsoOgdenTest, TestSpecifyFormulationAllFalse)
  {
    // initizalize function inputs
    bool isoprinc = false;
    bool isomod = false;
    bool anisoprinc = false;
    bool anisomod = false;
    bool viscogeneral = false;

    // call SpecifyFormulation function
    isoogden_.SpecifyFormulation(isoprinc, isomod, anisoprinc, anisomod, viscogeneral);

    // test if function correctly sets isomod to true
    EXPECT_TRUE(isomod);
    EXPECT_FALSE(isoprinc);
    EXPECT_FALSE(anisoprinc);
    EXPECT_FALSE(anisomod);
    EXPECT_FALSE(viscogeneral);
  }

  TEST_F(IsoOgdenTest, TestSpecifyFormulationAllTrue)
  {
    // initizalize function inputs
    bool isoprinc = true;
    bool isomod = true;
    bool anisoprinc = true;
    bool anisomod = true;
    bool viscogeneral = true;

    // call SpecifyFormulation function
    isoogden_.SpecifyFormulation(isoprinc, isomod, anisoprinc, anisomod, viscogeneral);

    // test if function correctly sets isomod to true
    EXPECT_TRUE(isomod);
    EXPECT_TRUE(isoprinc);
    EXPECT_TRUE(anisoprinc);
    EXPECT_TRUE(anisomod);
    EXPECT_TRUE(viscogeneral);
  }
}  // namespace
