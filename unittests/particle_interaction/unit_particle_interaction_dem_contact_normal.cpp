/*---------------------------------------------------------------------------*/
/*! \file
\brief unittests for normal contact handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include "baci_particle_interaction_dem_contact_normal.H"
#include "baci_particle_interaction_utils.H"

#include "baci_inpar_validparameters.H"


namespace
{
  class DEMContactNormalLinearSpringTest : public ::testing::Test
  {
   protected:
    std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLinearSpring> contactnormal_;

    const double r_max_ = 1.5;
    const double v_max_ = 1.25;
    const double c_ = 0.05;

    const double dens_max_ = 1.0;

    DEMContactNormalLinearSpringTest()
    {
      // create a parameter list
      Teuchos::ParameterList params_dem;
      params_dem.set("MAX_RADIUS", r_max_);
      params_dem.set("MAX_VELOCITY", v_max_);
      params_dem.set("REL_PENETRATION", c_);
      params_dem.set("NORMAL_STIFF", -1.0);

      // create normal contact handler
      contactnormal_ =
          std::make_unique<PARTICLEINTERACTION::DEMContactNormalLinearSpring>(params_dem);

      // init normal contact handler
      contactnormal_->Init();

      // setup normal contact handler
      contactnormal_->Setup(dens_max_);
    }

    // note: the public functions Init() and Setup() of class DEMContactNormalLinearSpring are
    // called in the constructor and thus implicitly tested by all following unittests
  };

  TEST_F(DEMContactNormalLinearSpringTest, GetNormalContactStiffness)
  {
    const double k_normal = 2.0 / 3.0 * r_max_ * M_PI * dens_max_ *
                            PARTICLEINTERACTION::UTILS::Pow<2>(v_max_) /
                            PARTICLEINTERACTION::UTILS::Pow<2>(c_);

    EXPECT_NEAR(contactnormal_->GetNormalContactStiffness(), k_normal, 1.0e-12);
  }

  TEST_F(DEMContactNormalLinearSpringTest, GetTimeCriticalStiffness)
  {
    const double k_normal = 2.0 / 3.0 * r_max_ * M_PI * dens_max_ *
                            PARTICLEINTERACTION::UTILS::Pow<2>(v_max_) /
                            PARTICLEINTERACTION::UTILS::Pow<2>(c_);

    const double k_normal_crit = k_normal;

    EXPECT_NEAR(contactnormal_->GetCriticalNormalContactStiffness(), k_normal_crit, 1.0e-12);
  }

  TEST_F(DEMContactNormalLinearSpringTest, NormalContactForce)
  {
    const double gap = -0.15;
    const double rad_i = 1.2;
    const double rad_j = 0.8;
    const double v_rel_normal = 0.3;
    const double m_eff = 1.0;

    const double normalcontactforce_ref = contactnormal_->GetNormalContactStiffness() * gap;

    double normalcontactforce = 0.0;
    contactnormal_->NormalContactForce(
        gap, &rad_i, &rad_j, v_rel_normal, m_eff, normalcontactforce);

    EXPECT_NEAR(normalcontactforce, normalcontactforce_ref, 1.0e-12);
  }

  TEST_F(DEMContactNormalLinearSpringTest, NormalPotentialEnergy)
  {
    const double gap = -0.15;

    const double normalpotentialenergy_ref =
        0.5 * contactnormal_->GetNormalContactStiffness() * PARTICLEINTERACTION::UTILS::Pow<2>(gap);

    double normalpotentialenergy = 0.0;
    contactnormal_->NormalPotentialEnergy(gap, normalpotentialenergy);

    EXPECT_NEAR(normalpotentialenergy, normalpotentialenergy_ref, 1.0e-12);
  }


  class DEMContactNormalLinearSpringDampTest : public ::testing::Test
  {
   protected:
    std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp> contactnormal_;
    std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp> contactnormal_ezero_;

    const double r_max_ = 1.5;
    const double v_max_ = 1.25;
    const double c_ = 0.05;
    const double e_ = 0.8;

    const double dens_max_ = 1.0;

    DEMContactNormalLinearSpringDampTest()
    {
      // create a parameter list
      Teuchos::ParameterList params_dem;
      params_dem.set("MAX_RADIUS", r_max_);
      params_dem.set("MAX_VELOCITY", v_max_);
      params_dem.set("REL_PENETRATION", c_);
      params_dem.set("NORMAL_STIFF", -1.0);
      params_dem.set("COEFF_RESTITUTION", e_);
      params_dem.set("DAMP_REG_FAC", -1.0);

      // create normal contact handler
      contactnormal_ =
          std::make_unique<PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp>(params_dem);

      params_dem.set("COEFF_RESTITUTION", 0.0);
      contactnormal_ezero_ =
          std::make_unique<PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp>(params_dem);

      // init normal contact handler
      contactnormal_->Init();
      contactnormal_ezero_->Init();

      // setup normal contact handler
      contactnormal_->Setup(dens_max_);
      contactnormal_ezero_->Setup(dens_max_);
    }
    // note: the public functions Init() and Setup() of class DEMContactNormalLinearSpringDamp are
    // called in the constructor and thus implicitly tested by all following unittests
  };

  TEST_F(DEMContactNormalLinearSpringDampTest, GetNormalContactStiffness)
  {
    const double k_normal = 2.0 / 3.0 * r_max_ * M_PI * dens_max_ *
                            PARTICLEINTERACTION::UTILS::Pow<2>(v_max_) /
                            PARTICLEINTERACTION::UTILS::Pow<2>(c_);

    EXPECT_NEAR(contactnormal_->GetNormalContactStiffness(), k_normal, 1.0e-12);
  }

  TEST_F(DEMContactNormalLinearSpringDampTest, GetTimeCriticalStiffness)
  {
    const double k_normal = 2.0 / 3.0 * r_max_ * M_PI * dens_max_ *
                            PARTICLEINTERACTION::UTILS::Pow<2>(v_max_) /
                            PARTICLEINTERACTION::UTILS::Pow<2>(c_);

    const double k_normal_crit = k_normal;

    EXPECT_NEAR(contactnormal_->GetCriticalNormalContactStiffness(), k_normal_crit, 1.0e-12);
  }

  TEST_F(DEMContactNormalLinearSpringDampTest, NormalContactForce)
  {
    const double gap = -0.15;
    const double rad_i = 1.2;
    const double rad_j = 0.8;
    const double v_rel_normal = 0.3;
    const double m_eff = 1.0;

    const double k_normal = contactnormal_->GetNormalContactStiffness();
    const double lne = std::log(e_);
    const double d_normal_fac = 2.0 * std::abs(lne) *
                                std::sqrt(k_normal / (PARTICLEINTERACTION::UTILS::Pow<2>(lne) +
                                                         PARTICLEINTERACTION::UTILS::Pow<2>(M_PI)));
    const double d_normal = d_normal_fac * std::sqrt(m_eff);
    const double normalcontactforce_ref = k_normal * gap - d_normal * v_rel_normal;

    double normalcontactforce = 0.0;
    contactnormal_->NormalContactForce(
        gap, &rad_i, &rad_j, v_rel_normal, m_eff, normalcontactforce);

    EXPECT_NEAR(normalcontactforce, normalcontactforce_ref, 1.0e-12);



    const double d_ezero_normal = 2.0 * std::sqrt(k_normal) * std::sqrt(m_eff);
    const double normalcontactforce_ezero_ref = k_normal * gap - d_ezero_normal * v_rel_normal;

    double normalcontactforce_ezero = 0.0;
    contactnormal_ezero_->NormalContactForce(
        gap, &rad_i, &rad_j, v_rel_normal, m_eff, normalcontactforce_ezero);

    EXPECT_NEAR(normalcontactforce_ezero, normalcontactforce_ezero_ref, 1.0e-12);
  }

  class DEMContactNormalHertzTest : public ::testing::Test
  {
   protected:
    std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalHertz> contactnormal_;

    const double r_max_ = 1.5;
    const double v_max_ = 1.25;
    const double c_ = 0.05;

    const double dens_max_ = 1.0;

    DEMContactNormalHertzTest()
    {
      // create a parameter list
      Teuchos::ParameterList params_dem;
      params_dem.set("MAX_RADIUS", r_max_);
      params_dem.set("MAX_VELOCITY", v_max_);
      params_dem.set("REL_PENETRATION", c_);
      params_dem.set("NORMAL_STIFF", -1.0);

      // create normal contact handler
      contactnormal_ = std::make_unique<PARTICLEINTERACTION::DEMContactNormalHertz>(params_dem);

      // init normal contact handler
      contactnormal_->Init();

      // setup normal contact handler
      contactnormal_->Setup(dens_max_);
    }

    // note: the public functions Init() and Setup() of class DEMContactNormalHertz are called in
    // Setup() and thus implicitly tested by all following unittests
  };

  TEST_F(DEMContactNormalHertzTest, GetNormalContactStiffness)
  {
    const double k_normal = 10.0 / 3.0 * M_PI * dens_max_ *
                            PARTICLEINTERACTION::UTILS::Pow<2>(v_max_) * std::pow(r_max_, 0.5) /
                            std::pow(2.0 * c_, 2.5);

    EXPECT_NEAR(contactnormal_->GetNormalContactStiffness(), k_normal, 1.0e-12);
  }

  TEST_F(DEMContactNormalHertzTest, GetTimeCriticalStiffness)
  {
    const double k_normal_crit = 2.0 / 3.0 * r_max_ * M_PI * dens_max_ *
                                 PARTICLEINTERACTION::UTILS::Pow<2>(v_max_) /
                                 PARTICLEINTERACTION::UTILS::Pow<2>(c_);

    EXPECT_NEAR(contactnormal_->GetCriticalNormalContactStiffness(), k_normal_crit, 1.0e-12);
  }

  TEST_F(DEMContactNormalHertzTest, NormalContactForce)
  {
    const double gap = -0.15;
    const double rad_i = 1.2;
    const double rad_j = 0.8;
    const double v_rel_normal = 0.3;
    const double m_eff = 1.0;

    const double normalcontactforce_ref =
        -contactnormal_->GetNormalContactStiffness() * std::pow(-gap, 1.5);

    double normalcontactforce = 0.0;
    contactnormal_->NormalContactForce(
        gap, &rad_i, &rad_j, v_rel_normal, m_eff, normalcontactforce);

    EXPECT_NEAR(normalcontactforce, normalcontactforce_ref, 1.0e-12);
  }

  TEST_F(DEMContactNormalHertzTest, NormalPotentialEnergy)
  {
    const double gap = -0.15;

    const double normalpotentialenergy_ref = 0.4 * contactnormal_->GetNormalContactStiffness() *
                                             PARTICLEINTERACTION::UTILS::Pow<2>(gap) *
                                             std::sqrt(-gap);

    double normalpotentialenergy = 0.0;
    contactnormal_->NormalPotentialEnergy(gap, normalpotentialenergy);

    EXPECT_NEAR(normalpotentialenergy, normalpotentialenergy_ref, 1.0e-12);
  }

  class DEMContactNormalLeeHerrmannTest : public ::testing::Test
  {
   protected:
    std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLeeHerrmann> contactnormal_;

    const double r_max_ = 1.5;
    const double v_max_ = 1.25;
    const double c_ = 0.05;
    const double d_normal_ = 6.0;

    const double dens_max_ = 1.0;

    DEMContactNormalLeeHerrmannTest()
    {
      // create a parameter list
      Teuchos::ParameterList params_dem;
      params_dem.set("MAX_RADIUS", r_max_);
      params_dem.set("MAX_VELOCITY", v_max_);
      params_dem.set("REL_PENETRATION", c_);
      params_dem.set("NORMAL_STIFF", -1.0);
      params_dem.set("NORMAL_DAMP", d_normal_);

      // create normal contact handler
      contactnormal_ =
          std::make_unique<PARTICLEINTERACTION::DEMContactNormalLeeHerrmann>(params_dem);

      // init normal contact handler
      contactnormal_->Init();

      // setup normal contact handler
      contactnormal_->Setup(dens_max_);
    }
    // note: the public functions Init() and Setup() of class DEMContactNormalLeeHerrmann are called
    // in Setup() and thus implicitly tested by all following unittests
  };

  TEST_F(DEMContactNormalLeeHerrmannTest, GetNormalContactStiffness)
  {
    const double k_normal = 10.0 / 3.0 * M_PI * dens_max_ *
                            PARTICLEINTERACTION::UTILS::Pow<2>(v_max_) * std::pow(r_max_, 0.5) /
                            std::pow(2.0 * c_, 2.5);

    EXPECT_NEAR(contactnormal_->GetNormalContactStiffness(), k_normal, 1.0e-12);
  }

  TEST_F(DEMContactNormalLeeHerrmannTest, GetTimeCriticalStiffness)
  {
    const double k_normal_crit = 2.0 / 3.0 * r_max_ * M_PI * dens_max_ *
                                 PARTICLEINTERACTION::UTILS::Pow<2>(v_max_) /
                                 PARTICLEINTERACTION::UTILS::Pow<2>(c_);

    EXPECT_NEAR(contactnormal_->GetCriticalNormalContactStiffness(), k_normal_crit, 1.0e-12);
  }

  TEST_F(DEMContactNormalLeeHerrmannTest, NormalContactForce)
  {
    const double gap = -0.15;
    const double rad_i = 1.2;
    const double rad_j = 0.8;
    const double v_rel_normal = 0.3;
    const double m_eff = 1.0;

    const double normalcontactforce_ref =
        -contactnormal_->GetNormalContactStiffness() * std::pow(-gap, 1.5) -
        m_eff * d_normal_ * v_rel_normal;

    double normalcontactforce = 0.0;
    contactnormal_->NormalContactForce(
        gap, &rad_i, &rad_j, v_rel_normal, m_eff, normalcontactforce);

    EXPECT_NEAR(normalcontactforce, normalcontactforce_ref, 1.0e-12);
  }

  class DEMContactNormalKuwabaraKonoTest : public ::testing::Test
  {
   protected:
    std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalKuwabaraKono> contactnormal_;

    const double r_max_ = 1.5;
    const double v_max_ = 1.25;
    const double c_ = 0.05;
    const double d_normal_ = 6.0;

    const double dens_max_ = 1.0;

    DEMContactNormalKuwabaraKonoTest()
    {
      // create a parameter list
      Teuchos::ParameterList params_dem;
      params_dem.set("MAX_RADIUS", r_max_);
      params_dem.set("MAX_VELOCITY", v_max_);
      params_dem.set("REL_PENETRATION", c_);
      params_dem.set("NORMAL_STIFF", -1.0);
      params_dem.set("NORMAL_DAMP", d_normal_);

      // create normal contact handler
      contactnormal_ =
          std::make_unique<PARTICLEINTERACTION::DEMContactNormalKuwabaraKono>(params_dem);

      // init normal contact handler
      contactnormal_->Init();

      // setup normal contact handler
      contactnormal_->Setup(dens_max_);
    }
    // note: the public functions Init() and Setup() of class DEMContactNormalKuwabaraKono are
    // called in the constructor and thus implicitly tested by all following unittests
  };

  TEST_F(DEMContactNormalKuwabaraKonoTest, GetNormalContactStiffness)
  {
    const double k_normal = 10.0 / 3.0 * M_PI * dens_max_ *
                            PARTICLEINTERACTION::UTILS::Pow<2>(v_max_) * std::pow(r_max_, 0.5) /
                            std::pow(2.0 * c_, 2.5);

    EXPECT_NEAR(contactnormal_->GetNormalContactStiffness(), k_normal, 1.0e-12);
  }

  TEST_F(DEMContactNormalKuwabaraKonoTest, GetTimeCriticalStiffness)
  {
    const double k_normal_crit = 2.0 / 3.0 * r_max_ * M_PI * dens_max_ *
                                 PARTICLEINTERACTION::UTILS::Pow<2>(v_max_) /
                                 PARTICLEINTERACTION::UTILS::Pow<2>(c_);

    EXPECT_NEAR(contactnormal_->GetCriticalNormalContactStiffness(), k_normal_crit, 1.0e-12);
  }

  TEST_F(DEMContactNormalKuwabaraKonoTest, NormalContactForce)
  {
    const double gap = -0.15;
    const double rad_i = 1.2;
    const double rad_j = 0.8;
    const double v_rel_normal = 0.3;
    const double m_eff = 1.0;

    const double normalcontactforce_ref =
        -contactnormal_->GetNormalContactStiffness() * std::pow(-gap, 1.5) -
        d_normal_ * v_rel_normal * std::pow(-gap, 0.5);

    double normalcontactforce = 0.0;
    contactnormal_->NormalContactForce(
        gap, &rad_i, &rad_j, v_rel_normal, m_eff, normalcontactforce);

    EXPECT_NEAR(normalcontactforce, normalcontactforce_ref, 1.0e-12);
  }

  class DEMContactNormalTsujiTest : public ::testing::Test
  {
   protected:
    std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalTsuji> contactnormal_;

    const double r_max_ = 1.5;
    const double v_max_ = 1.25;
    const double c_ = 0.05;
    const double d_normal_ = 6.0;

    const double dens_max_ = 1.0;

    DEMContactNormalTsujiTest()
    {
      // create a parameter list
      Teuchos::ParameterList params_dem;
      params_dem.set("MAX_RADIUS", r_max_);
      params_dem.set("MAX_VELOCITY", v_max_);
      params_dem.set("REL_PENETRATION", c_);
      params_dem.set("NORMAL_STIFF", -1.0);
      params_dem.set("NORMAL_DAMP", d_normal_);

      // create normal contact handler
      contactnormal_ = std::make_unique<PARTICLEINTERACTION::DEMContactNormalTsuji>(params_dem);

      // init normal contact handler
      contactnormal_->Init();

      // setup normal contact handler
      contactnormal_->Setup(dens_max_);
    }
    // note: the public functions Init() and Setup() of class DEMContactNormalTsuji are called in
    // Setup() and thus implicitly tested by all following unittests
  };

  TEST_F(DEMContactNormalTsujiTest, GetNormalContactStiffness)
  {
    const double k_normal = 10.0 / 3.0 * M_PI * dens_max_ *
                            PARTICLEINTERACTION::UTILS::Pow<2>(v_max_) * std::pow(r_max_, 0.5) /
                            std::pow(2.0 * c_, 2.5);

    EXPECT_NEAR(contactnormal_->GetNormalContactStiffness(), k_normal, 1.0e-12);
  }

  TEST_F(DEMContactNormalTsujiTest, GetTimeCriticalStiffness)
  {
    const double k_normal_crit = 2.0 / 3.0 * r_max_ * M_PI * dens_max_ *
                                 PARTICLEINTERACTION::UTILS::Pow<2>(v_max_) /
                                 PARTICLEINTERACTION::UTILS::Pow<2>(c_);

    EXPECT_NEAR(contactnormal_->GetCriticalNormalContactStiffness(), k_normal_crit, 1.0e-12);
  }

  TEST_F(DEMContactNormalTsujiTest, NormalContactForce)
  {
    const double gap = -0.15;
    const double rad_i = 1.2;
    const double rad_j = 0.8;
    const double v_rel_normal = 0.3;
    const double m_eff = 1.0;

    const double normalcontactforce_ref =
        -contactnormal_->GetNormalContactStiffness() * std::pow(-gap, 1.5) -
        d_normal_ * v_rel_normal * std::pow(-gap, 0.25);

    double normalcontactforce = 0.0;
    contactnormal_->NormalContactForce(
        gap, &rad_i, &rad_j, v_rel_normal, m_eff, normalcontactforce);

    EXPECT_NEAR(normalcontactforce, normalcontactforce_ref, 1.0e-12);
  }
}  // namespace
