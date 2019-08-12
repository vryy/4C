/*---------------------------------------------------------------------------*/
/*!
\brief unittests for rolling contact handler for discrete element method (DEM) interactions

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
#ifndef UNIT_PARTICLE_INTERACTION_DEM_CONTACT_ROLLING_H
#define UNIT_PARTICLE_INTERACTION_DEM_CONTACT_ROLLING_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
#include <cxxtest/TestSuite.h>
#include "src/drt_particle_interaction/particle_interaction_dem_contact_rolling.H"
#include "src/drt_particle_interaction/particle_interaction_utils.H"

#include "src/drt_inpar/drt_validparameters.H"

/*---------------------------------------------------------------------------*
 | forward declaration                                        sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMContactRollingViscous_TestSuite;
  class DEMContactRollingCoulomb_TestSuite;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | rolling contact handler test suite                         sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
class PARTICLEINTERACTION::DEMContactRollingViscous_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::DEMContactRollingViscous> contactrolling_;

  const double e_ = 0.8;
  const double nue_ = 0.4;
  const double mu_rolling_ = 0.2;

  const double young_ = 200.0e3;
  const double v_max_ = 0.025;

  const double k_normal_ = 4.0;

 public:
  void setUp()
  {
    // create a parameter list
    Teuchos::ParameterList params_dem;
    params_dem.set("COEFF_RESTITUTION", e_);
    params_dem.set("POISSON_RATIO", nue_);
    params_dem.set("FRICT_COEFF_ROLL", mu_rolling_);

    params_dem.set("YOUNG_MODULUS", young_);
    params_dem.set("MAX_VELOCITY", v_max_);

    // create rolling contact handler
    contactrolling_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactRollingViscous>(
        new PARTICLEINTERACTION::DEMContactRollingViscous(params_dem));

    // init rolling contact handler
    contactrolling_->Init();

    // setup rolling contact handler
    contactrolling_->Setup(k_normal_);
  }

  void tearDown() { contactrolling_ = nullptr; }

  // note: the public functions Init() and Setup() of class DEMContactRollingViscous are
  // called in setUp() and thus implicitly tested by all following unittests

  void test_EffectiveRadiusParticle()
  {
    const double rad_i = 1.2;
    const double rad_j = 0.8;
    const double gap = -0.3;

    double r_eff = 0.0;
    contactrolling_->EffectiveRadiusParticle(&rad_i, &rad_j, gap, r_eff);

    const double r_eff_ref = rad_i * rad_j / (rad_i + rad_j);

    TS_ASSERT_DELTA(r_eff, r_eff_ref, 1.0e-12);
  }

  void test_EffectiveRadiusParticle_nullptr_rad_j()
  {
    const double rad_i = 1.2;
    const double gap = -0.3;

    double r_eff = 0.0;
    contactrolling_->EffectiveRadiusParticle(&rad_i, nullptr, gap, r_eff);

    const double r_eff_ref = rad_i;

    TS_ASSERT_DELTA(r_eff, r_eff_ref, 1.0e-12);
  }

  void test_RelativeRollingVelocity()
  {
    const double r_eff = 0.5;

    double e_ji[3] = {0.0};
    e_ji[0] = 1.0 / std::sqrt(21);
    e_ji[1] = 2.0 / std::sqrt(21);
    e_ji[2] = 4.0 / std::sqrt(21);

    double angvel_i[3] = {0.0};
    angvel_i[0] = -0.03;
    angvel_i[1] = 0.1;
    angvel_i[2] = 0.12;

    double angvel_j[3] = {0.0};
    angvel_j[0] = 0.15;
    angvel_j[1] = -0.2;
    angvel_j[2] = 0.0;

    double v_rel_rolling[3] = {0.0};
    contactrolling_->RelativeRollingVelocity(r_eff, e_ji, angvel_i, angvel_j, v_rel_rolling);

    double v_rel_rolling_ref[3] = {0.0};
    UTILS::vec_setcross(v_rel_rolling_ref, angvel_i, e_ji);
    UTILS::vec_addcross(v_rel_rolling_ref, e_ji, angvel_j);

    for (int i = 0; i < 3; ++i) TS_ASSERT_DELTA(v_rel_rolling[i], v_rel_rolling_ref[i], 1.0e-12);
  }

  void test_RelativeRollingVelocity_nullptr_angvel_j()
  {
    const double r_eff = 0.5;

    double e_ji[3] = {0.0};
    e_ji[0] = 1.0 / std::sqrt(21);
    e_ji[1] = 2.0 / std::sqrt(21);
    e_ji[2] = 4.0 / std::sqrt(21);

    double angvel_i[3] = {0.0};
    angvel_i[0] = -0.03;
    angvel_i[1] = 0.1;
    angvel_i[2] = 0.12;

    double v_rel_rolling[3] = {0.0};
    contactrolling_->RelativeRollingVelocity(r_eff, e_ji, angvel_i, nullptr, v_rel_rolling);

    double v_rel_rolling_ref[3] = {0.0};
    UTILS::vec_setcross(v_rel_rolling_ref, angvel_i, e_ji);

    for (int i = 0; i < 3; ++i) TS_ASSERT_DELTA(v_rel_rolling[i], v_rel_rolling_ref[i], 1.0e-12);
  }

  void test_RollingContactMoment()
  {
    double gap_rolling[3] = {0.0};
    gap_rolling[0] = 0.1;
    gap_rolling[1] = 0.05;
    gap_rolling[2] = -0.25;

    bool stick_rolling = false;

    double e_ji[3] = {0.0};
    e_ji[0] = 1.0 / std::sqrt(21);
    e_ji[1] = 2.0 / std::sqrt(21);
    e_ji[2] = 4.0 / std::sqrt(21);

    double v_rel_rolling[3] = {0.0};
    v_rel_rolling[0] = -0.03;
    v_rel_rolling[1] = 0.1;
    v_rel_rolling[2] = 0.12;

    const double m_eff = 2.5;
    const double r_eff = 0.5;
    const double normalcontactforce = 1.5;

    double rollingcontactmoment[3] = {0.0};
    contactrolling_->RollingContactMoment(gap_rolling, stick_rolling, e_ji, v_rel_rolling, m_eff,
        r_eff, mu_rolling_, normalcontactforce, rollingcontactmoment);

    double rollingcontactmoment_ref[3] = {0.0};
    const double fac = young_ / (1.0 - UTILS::pow<2>(nue_));
    const double c_1 = 1.15344;
    const double d_rolling_fac =
        mu_rolling_ * (1.0 - e_) / (c_1 * std::pow(fac, 0.4) * std::pow(v_max_, 0.2));
    const double d_rolling = d_rolling_fac * std::pow(0.5 * r_eff, -0.2);

    double rollingcontactforce[3];
    UTILS::vec_setscale(rollingcontactforce, -(d_rolling * normalcontactforce), v_rel_rolling);

    UTILS::vec_setcross(rollingcontactmoment_ref, rollingcontactforce, e_ji);
    UTILS::vec_scale(rollingcontactmoment_ref, r_eff);

    for (int i = 0; i < 3; ++i)
      TS_ASSERT_DELTA(rollingcontactmoment[i], rollingcontactmoment_ref[i], 1.0e-12);
  }
};

class PARTICLEINTERACTION::DEMContactRollingCoulomb_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::DEMContactRollingCoulomb> contactrolling_;

  const double e_ = 0.8;
  const double nue_ = 0.4;
  const double mu_rolling_ = 0.2;

  const double k_normal_ = 4.0;

 public:
  void setUp()
  {
    // create a parameter list
    Teuchos::ParameterList params_dem;
    params_dem.set("COEFF_RESTITUTION", e_);
    params_dem.set("POISSON_RATIO", nue_);
    params_dem.set("FRICT_COEFF_ROLL", mu_rolling_);

    // create rolling contact handler
    contactrolling_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactRollingCoulomb>(
        new PARTICLEINTERACTION::DEMContactRollingCoulomb(params_dem));

    // init rolling contact handler
    contactrolling_->Init();

    // setup rolling contact handler
    contactrolling_->Setup(k_normal_);
  }

  void tearDown() { contactrolling_ = nullptr; }

  // note: the public functions Init() and Setup() of class DEMContactRollingViscous are
  // called in setUp() and thus implicitly tested by all following unittests

  void test_EffectiveRadiusParticle()
  {
    const double rad_i = 1.2;
    const double rad_j = 0.8;
    const double gap = -0.3;

    double r_eff = 0.0;
    contactrolling_->EffectiveRadiusParticle(&rad_i, &rad_j, gap, r_eff);

    const double r_eff_ref = (rad_i + 0.5 * gap) * (rad_j + 0.5 * gap) / (rad_i + rad_j + gap);

    TS_ASSERT_DELTA(r_eff, r_eff_ref, 1.0e-12);
  }

  void test_EffectiveRadiusParticle_nullptr_rad_j()
  {
    const double rad_i = 1.2;
    const double gap = -0.3;

    double r_eff = 0.0;
    contactrolling_->EffectiveRadiusParticle(&rad_i, nullptr, gap, r_eff);

    const double r_eff_ref = rad_i + gap;

    TS_ASSERT_DELTA(r_eff, r_eff_ref, 1.0e-12);
  }

  void test_RelativeRollingVelocity()
  {
    const double r_eff = 0.5;

    double e_ji[3] = {0.0};
    e_ji[0] = 1.0 / std::sqrt(21);
    e_ji[1] = 2.0 / std::sqrt(21);
    e_ji[2] = 4.0 / std::sqrt(21);

    double angvel_i[3] = {0.0};
    angvel_i[0] = -0.03;
    angvel_i[1] = 0.1;
    angvel_i[2] = 0.12;

    double angvel_j[3] = {0.0};
    angvel_j[0] = 0.15;
    angvel_j[1] = -0.2;
    angvel_j[2] = 0.0;

    double v_rel_rolling[3] = {0.0};
    contactrolling_->RelativeRollingVelocity(r_eff, e_ji, angvel_i, angvel_j, v_rel_rolling);

    double v_rel_rolling_ref[3] = {0.0};
    UTILS::vec_setcross(v_rel_rolling_ref, e_ji, angvel_i);
    UTILS::vec_addcross(v_rel_rolling_ref, angvel_j, e_ji);
    UTILS::vec_scale(v_rel_rolling_ref, r_eff);

    for (int i = 0; i < 3; ++i) TS_ASSERT_DELTA(v_rel_rolling[i], v_rel_rolling_ref[i], 1.0e-12);
  }

  void test_RelativeRollingVelocity_nullptr_angvel_j()
  {
    const double r_eff = 0.5;

    double e_ji[3] = {0.0};
    e_ji[0] = 1.0 / std::sqrt(21);
    e_ji[1] = 2.0 / std::sqrt(21);
    e_ji[2] = 4.0 / std::sqrt(21);

    double angvel_i[3] = {0.0};
    angvel_i[0] = -0.03;
    angvel_i[1] = 0.1;
    angvel_i[2] = 0.12;

    double v_rel_rolling[3] = {0.0};
    contactrolling_->RelativeRollingVelocity(r_eff, e_ji, angvel_i, nullptr, v_rel_rolling);

    double v_rel_rolling_ref[3] = {0.0};
    UTILS::vec_setcross(v_rel_rolling_ref, e_ji, angvel_i);
    UTILS::vec_scale(v_rel_rolling_ref, r_eff);

    for (int i = 0; i < 3; ++i) TS_ASSERT_DELTA(v_rel_rolling[i], v_rel_rolling_ref[i], 1.0e-12);
  }

  void test_RollingContactMoment_stick()
  {
    double gap_rolling[3] = {0.0};
    gap_rolling[0] = 0.1;
    gap_rolling[1] = 0.05;
    gap_rolling[2] = -0.25;

    bool stick_rolling = true;

    double e_ji[3] = {0.0};
    e_ji[0] = 1.0 / std::sqrt(21);
    e_ji[1] = 2.0 / std::sqrt(21);
    e_ji[2] = 4.0 / std::sqrt(21);

    double v_rel_rolling[3] = {0.0};
    v_rel_rolling[0] = -0.03;
    v_rel_rolling[1] = 0.1;
    v_rel_rolling[2] = 0.12;

    const double m_eff = 2.5;
    const double r_eff = 0.5;
    const double normalcontactforce = 1.5;

    double rollingcontactmoment[3] = {0.0};
    contactrolling_->RollingContactMoment(gap_rolling, stick_rolling, e_ji, v_rel_rolling, m_eff,
        r_eff, mu_rolling_, normalcontactforce, rollingcontactmoment);

    double gap_rolling_ref[3] = {0.06858666958093279, 0.05062425428544091, -0.05782673704478762};
    double rollingcontactmoment_ref[3] = {
        -0.1119618097864935, 0.09699534835316952, -0.02050722172996138};

    for (int i = 0; i < 3; ++i) TS_ASSERT_DELTA(gap_rolling[i], gap_rolling_ref[i], 1.0e-12);

    for (int i = 0; i < 3; ++i)
      TS_ASSERT_DELTA(rollingcontactmoment[i], rollingcontactmoment_ref[i], 1.0e-12);

    TS_ASSERT_EQUALS(stick_rolling, false);
  }

  void test_RollingContactMoment_slip()
  {
    double gap_rolling[3] = {0.0};
    gap_rolling[0] = 0.1;
    gap_rolling[1] = 0.05;
    gap_rolling[2] = -0.25;

    bool stick_rolling = false;

    double e_ji[3] = {0.0};
    e_ji[0] = 1.0 / std::sqrt(21);
    e_ji[1] = 2.0 / std::sqrt(21);
    e_ji[2] = 4.0 / std::sqrt(21);

    double v_rel_rolling[3] = {0.0};
    v_rel_rolling[0] = -0.03;
    v_rel_rolling[1] = 0.1;
    v_rel_rolling[2] = 0.12;

    const double m_eff = 2.5;
    const double r_eff = 0.5;
    const double normalcontactforce = 1.5;

    double rollingcontactmoment[3] = {0.0};
    contactrolling_->RollingContactMoment(gap_rolling, stick_rolling, e_ji, v_rel_rolling, m_eff,
        r_eff, mu_rolling_, normalcontactforce, rollingcontactmoment);

    double gap_rolling_ref[3] = {0.06858666958093279, 0.05062425428544091, -0.05782673704478762};
    double rollingcontactmoment_ref[3] = {
        -0.1119618097864935, 0.09699534835316952, -0.02050722172996138};

    for (int i = 0; i < 3; ++i) TS_ASSERT_DELTA(gap_rolling[i], gap_rolling_ref[i], 1.0e-12);

    for (int i = 0; i < 3; ++i)
      TS_ASSERT_DELTA(rollingcontactmoment[i], rollingcontactmoment_ref[i], 1.0e-12);

    TS_ASSERT_EQUALS(stick_rolling, false);
  }
};

/*---------------------------------------------------------------------------*/
#endif
