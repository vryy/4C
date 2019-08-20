/*---------------------------------------------------------------------------*/
/*! \file
\brief unittests for normal contact handler for discrete element method (DEM) interactions

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
#ifndef UNIT_PARTICLE_INTERACTION_DEM_CONTACT_NORMAL_H
#define UNIT_PARTICLE_INTERACTION_DEM_CONTACT_NORMAL_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
#include <cxxtest/TestSuite.h>
#include "src/drt_particle_interaction/particle_interaction_dem_contact_normal.H"
#include "src/drt_particle_interaction/particle_interaction_utils.H"

#include "src/drt_inpar/drt_validparameters.H"

/*---------------------------------------------------------------------------*
 | forward declaration                                        sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMContactNormalLinearSpring_TestSuite;
  class DEMContactNormalLinearSpringDamp_TestSuite;
  class DEMContactNormalHertz_TestSuite;
  class DEMContactNormalLeeHerrmann_TestSuite;
  class DEMContactNormalKuwabaraKono_TestSuite;
  class DEMContactNormalTsuji_TestSuite;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | normal contact handler test suite                          sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
class PARTICLEINTERACTION::DEMContactNormalLinearSpring_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLinearSpring> contactnormal_;

  const double r_max_ = 1.5;
  const double v_max_ = 1.25;
  const double c_ = 0.05;

  const double dens_max_ = 1.0;

 public:
  void setUp()
  {
    // create a parameter list
    Teuchos::ParameterList params_dem;
    params_dem.set("MAX_RADIUS", r_max_);
    params_dem.set("MAX_VELOCITY", v_max_);
    params_dem.set("REL_PENETRATION", c_);
    params_dem.set("NORMAL_STIFF", -1.0);

    // create normal contact handler
    contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLinearSpring>(
        new PARTICLEINTERACTION::DEMContactNormalLinearSpring(params_dem));

    // init normal contact handler
    contactnormal_->Init();

    // setup normal contact handler
    contactnormal_->Setup(dens_max_);
  }

  void tearDown() { contactnormal_ = nullptr; }

  // note: the public functions Init() and Setup() of class DEMContactNormalLinearSpring are called
  // in setUp() and thus implicitly tested by all following unittests

  void test_GetNormalContactStiffness()
  {
    const double k_normal =
        2.0 / 3.0 * r_max_ * M_PI * dens_max_ * UTILS::pow<2>(v_max_) / UTILS::pow<2>(c_);

    TS_ASSERT_DELTA(contactnormal_->GetNormalContactStiffness(), k_normal, 1.0e-12);
  }

  void test_GetTimeCriticalStiffness()
  {
    const double k_normal =
        2.0 / 3.0 * r_max_ * M_PI * dens_max_ * UTILS::pow<2>(v_max_) / UTILS::pow<2>(c_);

    const double k_tcrit = k_normal;

    TS_ASSERT_DELTA(contactnormal_->GetTimeCriticalStiffness(), k_tcrit, 1.0e-12);
  }

  void test_NormalContactForce()
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

    TS_ASSERT_DELTA(normalcontactforce, normalcontactforce_ref, 1.0e-12);
  }
};

class PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp> contactnormal_;
  std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp> contactnormal_ezero_;

  const double r_max_ = 1.5;
  const double v_max_ = 1.25;
  const double c_ = 0.05;
  const double e_ = 0.8;

  const double dens_max_ = 1.0;

 public:
  void setUp()
  {
    // create a parameter list
    Teuchos::ParameterList params_dem;
    params_dem.set("MAX_RADIUS", r_max_);
    params_dem.set("MAX_VELOCITY", v_max_);
    params_dem.set("REL_PENETRATION", c_);
    params_dem.set("NORMAL_STIFF", -1.0);

    params_dem.set("COEFF_RESTITUTION", e_);
    params_dem.set("DAMP_REG_FAC", -1.0);
    DRT::INPUT::BoolParameter("TENSION_CUTOFF", "no", "dummy", &params_dem);

    // create normal contact handler
    contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp>(
        new PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp(params_dem));

    params_dem.set("COEFF_RESTITUTION", 0.0);
    contactnormal_ezero_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp>(
        new PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp(params_dem));

    // init normal contact handler
    contactnormal_->Init();
    contactnormal_ezero_->Init();

    // setup normal contact handler
    contactnormal_->Setup(dens_max_);
    contactnormal_ezero_->Setup(dens_max_);
  }

  void tearDown()
  {
    contactnormal_ = nullptr;
    contactnormal_ezero_ = nullptr;
  }

  // note: the public functions Init() and Setup() of class DEMContactNormalLinearSpringDamp are
  // called in setUp() and thus implicitly tested by all following unittests

  void test_GetNormalContactStiffness()
  {
    const double k_normal =
        2.0 / 3.0 * r_max_ * M_PI * dens_max_ * UTILS::pow<2>(v_max_) / UTILS::pow<2>(c_);

    TS_ASSERT_DELTA(contactnormal_->GetNormalContactStiffness(), k_normal, 1.0e-12);
  }

  void test_GetTimeCriticalStiffness()
  {
    const double k_normal =
        2.0 / 3.0 * r_max_ * M_PI * dens_max_ * UTILS::pow<2>(v_max_) / UTILS::pow<2>(c_);

    const double k_tcrit = k_normal;

    TS_ASSERT_DELTA(contactnormal_->GetTimeCriticalStiffness(), k_tcrit, 1.0e-12);
  }

  void test_NormalContactForce()
  {
    const double gap = -0.15;
    const double rad_i = 1.2;
    const double rad_j = 0.8;
    const double v_rel_normal = 0.3;
    const double m_eff = 1.0;

    const double k_normal = contactnormal_->GetNormalContactStiffness();
    const double lne = std::log(e_);
    const double d_normal_fac =
        2.0 * std::abs(lne) * std::sqrt(k_normal / (UTILS::pow<2>(lne) + UTILS::pow<2>(M_PI)));
    const double d_normal = d_normal_fac * std::sqrt(m_eff);
    const double normalcontactforce_ref = k_normal * gap - d_normal * v_rel_normal;

    double normalcontactforce = 0.0;
    contactnormal_->NormalContactForce(
        gap, &rad_i, &rad_j, v_rel_normal, m_eff, normalcontactforce);

    TS_ASSERT_DELTA(normalcontactforce, normalcontactforce_ref, 1.0e-12);



    const double d_ezero_normal = 2.0 * std::sqrt(k_normal) * std::sqrt(m_eff);
    const double normalcontactforce_ezero_ref = k_normal * gap - d_ezero_normal * v_rel_normal;

    double normalcontactforce_ezero = 0.0;
    contactnormal_ezero_->NormalContactForce(
        gap, &rad_i, &rad_j, v_rel_normal, m_eff, normalcontactforce_ezero);

    TS_ASSERT_DELTA(normalcontactforce_ezero, normalcontactforce_ezero_ref, 1.0e-12);
  }
};

class PARTICLEINTERACTION::DEMContactNormalHertz_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalHertz> contactnormal_;

  const double r_max_ = 1.5;
  const double v_max_ = 1.25;
  const double c_ = 0.05;

  const double dens_max_ = 1.0;

 public:
  void setUp()
  {
    // create a parameter list
    Teuchos::ParameterList params_dem;
    params_dem.set("MAX_RADIUS", r_max_);
    params_dem.set("MAX_VELOCITY", v_max_);
    params_dem.set("REL_PENETRATION", c_);
    params_dem.set("NORMAL_STIFF", -1.0);

    // create normal contact handler
    contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalHertz>(
        new PARTICLEINTERACTION::DEMContactNormalHertz(params_dem));

    // init normal contact handler
    contactnormal_->Init();

    // setup normal contact handler
    contactnormal_->Setup(dens_max_);
  }

  void tearDown() { contactnormal_ = nullptr; }

  // note: the public functions Init() and Setup() of class DEMContactNormalHertz are called in
  // setUp() and thus implicitly tested by all following unittests

  void test_GetNormalContactStiffness()
  {
    const double k_normal = 10.0 / 3.0 * M_PI * dens_max_ * UTILS::pow<2>(v_max_) *
                            std::pow(r_max_, 0.5) / std::pow(2.0 * c_, 2.5);

    TS_ASSERT_DELTA(contactnormal_->GetNormalContactStiffness(), k_normal, 1.0e-12);
  }

  void test_GetTimeCriticalStiffness()
  {
    const double k_tcrit =
        2.0 / 3.0 * r_max_ * M_PI * dens_max_ * UTILS::pow<2>(v_max_) / UTILS::pow<2>(c_);

    TS_ASSERT_DELTA(contactnormal_->GetTimeCriticalStiffness(), k_tcrit, 1.0e-12);
  }

  void test_NormalContactForce()
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

    TS_ASSERT_DELTA(normalcontactforce, normalcontactforce_ref, 1.0e-12);
  }
};

class PARTICLEINTERACTION::DEMContactNormalLeeHerrmann_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLeeHerrmann> contactnormal_;

  const double r_max_ = 1.5;
  const double v_max_ = 1.25;
  const double c_ = 0.05;
  const double d_normal_ = 6.0;

  const double dens_max_ = 1.0;

 public:
  void setUp()
  {
    // create a parameter list
    Teuchos::ParameterList params_dem;
    params_dem.set("MAX_RADIUS", r_max_);
    params_dem.set("MAX_VELOCITY", v_max_);
    params_dem.set("REL_PENETRATION", c_);
    params_dem.set("NORMAL_STIFF", -1.0);

    params_dem.set("NORMAL_DAMP", d_normal_);
    DRT::INPUT::BoolParameter("TENSION_CUTOFF", "no", "dummy", &params_dem);

    // create normal contact handler
    contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLeeHerrmann>(
        new PARTICLEINTERACTION::DEMContactNormalLeeHerrmann(params_dem));

    // init normal contact handler
    contactnormal_->Init();

    // setup normal contact handler
    contactnormal_->Setup(dens_max_);
  }

  void tearDown() { contactnormal_ = nullptr; }

  // note: the public functions Init() and Setup() of class DEMContactNormalLeeHerrmann are called
  // in setUp() and thus implicitly tested by all following unittests

  void test_GetNormalContactStiffness()
  {
    const double k_normal = 10.0 / 3.0 * M_PI * dens_max_ * UTILS::pow<2>(v_max_) *
                            std::pow(r_max_, 0.5) / std::pow(2.0 * c_, 2.5);

    TS_ASSERT_DELTA(contactnormal_->GetNormalContactStiffness(), k_normal, 1.0e-12);
  }

  void test_GetTimeCriticalStiffness()
  {
    const double k_tcrit =
        2.0 / 3.0 * r_max_ * M_PI * dens_max_ * UTILS::pow<2>(v_max_) / UTILS::pow<2>(c_);

    TS_ASSERT_DELTA(contactnormal_->GetTimeCriticalStiffness(), k_tcrit, 1.0e-12);
  }

  void test_NormalContactForce()
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

    TS_ASSERT_DELTA(normalcontactforce, normalcontactforce_ref, 1.0e-12);
  }
};

class PARTICLEINTERACTION::DEMContactNormalKuwabaraKono_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalKuwabaraKono> contactnormal_;

  const double r_max_ = 1.5;
  const double v_max_ = 1.25;
  const double c_ = 0.05;
  const double d_normal_ = 6.0;

  const double dens_max_ = 1.0;

 public:
  void setUp()
  {
    // create a parameter list
    Teuchos::ParameterList params_dem;
    params_dem.set("MAX_RADIUS", r_max_);
    params_dem.set("MAX_VELOCITY", v_max_);
    params_dem.set("REL_PENETRATION", c_);
    params_dem.set("NORMAL_STIFF", -1.0);

    params_dem.set("NORMAL_DAMP", d_normal_);
    DRT::INPUT::BoolParameter("TENSION_CUTOFF", "no", "dummy", &params_dem);

    // create normal contact handler
    contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalKuwabaraKono>(
        new PARTICLEINTERACTION::DEMContactNormalKuwabaraKono(params_dem));

    // init normal contact handler
    contactnormal_->Init();

    // setup normal contact handler
    contactnormal_->Setup(dens_max_);
  }

  void tearDown() { contactnormal_ = nullptr; }

  // note: the public functions Init() and Setup() of class DEMContactNormalKuwabaraKono are called
  // in setUp() and thus implicitly tested by all following unittests

  void test_GetNormalContactStiffness()
  {
    const double k_normal = 10.0 / 3.0 * M_PI * dens_max_ * UTILS::pow<2>(v_max_) *
                            std::pow(r_max_, 0.5) / std::pow(2.0 * c_, 2.5);

    TS_ASSERT_DELTA(contactnormal_->GetNormalContactStiffness(), k_normal, 1.0e-12);
  }

  void test_GetTimeCriticalStiffness()
  {
    const double k_tcrit =
        2.0 / 3.0 * r_max_ * M_PI * dens_max_ * UTILS::pow<2>(v_max_) / UTILS::pow<2>(c_);

    TS_ASSERT_DELTA(contactnormal_->GetTimeCriticalStiffness(), k_tcrit, 1.0e-12);
  }

  void test_NormalContactForce()
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

    TS_ASSERT_DELTA(normalcontactforce, normalcontactforce_ref, 1.0e-12);
  }
};

class PARTICLEINTERACTION::DEMContactNormalTsuji_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalTsuji> contactnormal_;

  const double r_max_ = 1.5;
  const double v_max_ = 1.25;
  const double c_ = 0.05;
  const double d_normal_ = 6.0;

  const double dens_max_ = 1.0;

 public:
  void setUp()
  {
    // create a parameter list
    Teuchos::ParameterList params_dem;
    params_dem.set("MAX_RADIUS", r_max_);
    params_dem.set("MAX_VELOCITY", v_max_);
    params_dem.set("REL_PENETRATION", c_);
    params_dem.set("NORMAL_STIFF", -1.0);

    params_dem.set("NORMAL_DAMP", d_normal_);
    DRT::INPUT::BoolParameter("TENSION_CUTOFF", "no", "dummy", &params_dem);

    // create normal contact handler
    contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalTsuji>(
        new PARTICLEINTERACTION::DEMContactNormalTsuji(params_dem));

    // init normal contact handler
    contactnormal_->Init();

    // setup normal contact handler
    contactnormal_->Setup(dens_max_);
  }

  void tearDown() { contactnormal_ = nullptr; }

  // note: the public functions Init() and Setup() of class DEMContactNormalTsuji are called in
  // setUp() and thus implicitly tested by all following unittests

  void test_GetNormalContactStiffness()
  {
    const double k_normal = 10.0 / 3.0 * M_PI * dens_max_ * UTILS::pow<2>(v_max_) *
                            std::pow(r_max_, 0.5) / std::pow(2.0 * c_, 2.5);

    TS_ASSERT_DELTA(contactnormal_->GetNormalContactStiffness(), k_normal, 1.0e-12);
  }

  void test_GetTimeCriticalStiffness()
  {
    const double k_tcrit =
        2.0 / 3.0 * r_max_ * M_PI * dens_max_ * UTILS::pow<2>(v_max_) / UTILS::pow<2>(c_);

    TS_ASSERT_DELTA(contactnormal_->GetTimeCriticalStiffness(), k_tcrit, 1.0e-12);
  }

  void test_NormalContactForce()
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

    TS_ASSERT_DELTA(normalcontactforce, normalcontactforce_ref, 1.0e-12);
  }
};

/*---------------------------------------------------------------------------*/
#endif
