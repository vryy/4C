/*----------------------------------------------------------------------*/
/*! \file

\brief unittests for utils of poromultiphase_scatra-framework


\level 1

*-----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "baci_poromultiphase_scatra_utils.hpp"

namespace
{
  using namespace FourC;

  class TestCalculationOxyPartialPressure : public ::testing::Test
  {
   protected:
    /**
     * helper function to get oxygen concentration from oxygen partial pressure, this is the known
     * forward relation
     *
     * @param Pb [in]: oxygen partial pressure
     * @param CaO2_max [in]: maximum oxygen concentration
     * @param Pb50 [in]: partial pressure at 50% maximum oxygen concentration
     * @param n [in]: exponent in Hill equation
     * @param alpha_eff [in]: effective solubility of oxygen in blood
     * @return CaO2: oxygen concentration
     */
    double GetOxyConcentrationFromPartialPressure(const double& Pb, const double& CaO2_max,
        const double& Pb50, const double& n, const double& alpha_eff)
    {
      // Hill equation for saturation
      const double saturation_hill = std::pow(Pb, n) / (std::pow(Pb, n) + std::pow(Pb50, n));
      // oxygen dissolved in plasma + bound to hemoglobin
      const double CaO2 = alpha_eff * Pb + CaO2_max * saturation_hill;

      return CaO2;
    }
  };

  //! test corner value Pb = 0 mmHg
  TEST_F(TestCalculationOxyPartialPressure, Pb0)
  {
    // Parameters
    const double n = 2.7;
    const double Pb50 = 37.0;
    const double alpha_eff = 3.1e-5;
    const double CaO2_max = 0.225;

    // invert to get partial pressure for this concentration
    double inverted_Pb = 0.0;
    POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<double>(
        inverted_Pb, 0.0, CaO2_max, Pb50, n, alpha_eff);

    // check if partial pressure and numerically inverted partial pressure are equal
    EXPECT_NEAR(0.0, inverted_Pb, 1e-14);
  }

  //! test value Pb = 25 mmHg
  TEST_F(TestCalculationOxyPartialPressure, Pb25)
  {
    // Parameters
    const double n = 2.5;
    const double Pb50 = 45.0;
    const double alpha_eff = 3.1e-5;
    const double CaO2_max = 0.225;
    const double Pb = 25.0;

    // get oxygen concentration for this partial pressure
    const double CaO2 = GetOxyConcentrationFromPartialPressure(Pb, CaO2_max, Pb50, n, alpha_eff);

    // invert to get partial pressure for this concentration
    double inverted_Pb = 0.0;
    POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<double>(
        inverted_Pb, CaO2, CaO2_max, Pb50, n, alpha_eff);

    // check if partial pressure and numerically inverted partial pressure are equal
    EXPECT_NEAR(Pb, inverted_Pb, 1e-14);
  }

  //! test value Pb = 50 mmHg
  TEST_F(TestCalculationOxyPartialPressure, Pb50)
  {
    // Parameters
    const double n = 2.7;
    const double Pb50 = 37.0;
    const double alpha_eff = 3.1e-5;
    const double CaO2_max = 0.225;
    const double Pb = 50.0;

    // get oxygen concentration for this partial pressure
    const double CaO2 = GetOxyConcentrationFromPartialPressure(Pb, CaO2_max, Pb50, n, alpha_eff);

    // invert to get partial pressure for this concentration
    double inverted_Pb = 0.0;
    POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<double>(
        inverted_Pb, CaO2, CaO2_max, Pb50, n, alpha_eff);

    // check if partial pressure and numerically inverted partial pressure are equal
    EXPECT_NEAR(Pb, inverted_Pb, 1e-14);
  }

  //! test value Pb = 100 mmHg
  TEST_F(TestCalculationOxyPartialPressure, Pb100)
  {
    // Parameters
    const double n = 2.2;
    const double Pb50 = 42.0;
    const double alpha_eff = 3.7e-5;
    const double CaO2_max = 0.225;
    const double Pb = 100.0;

    // get oxygen concentration for this partial pressure
    const double CaO2 = GetOxyConcentrationFromPartialPressure(Pb, CaO2_max, Pb50, n, alpha_eff);

    // invert to get partial pressure for this concentration
    double inverted_Pb = 0.0;
    POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<double>(
        inverted_Pb, CaO2, CaO2_max, Pb50, n, alpha_eff);

    // check if partial pressure and numerically inverted partial pressure are equal
    EXPECT_NEAR(Pb, inverted_Pb, 1e-14);
  }

}  // namespace