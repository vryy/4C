/*---------------------------------------------------------------------------*/
/*! \file
\brief unittests for tangential contact handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
#ifndef UNIT_PARTICLE_INTERACTION_DEM_CONTACT_TANGENTIAL_H
#define UNIT_PARTICLE_INTERACTION_DEM_CONTACT_TANGENTIAL_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
#include <cxxtest/TestSuite.h>
#include "src/drt_particle_interaction/particle_interaction_dem_contact_tangential.H"
#include "src/drt_particle_interaction/particle_interaction_utils.H"

#include "src/drt_inpar/drt_validparameters.H"

/*---------------------------------------------------------------------------*
 | forward declaration                                        sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMContactTangentialLinearSpringDamp_TestSuite;
}

/*---------------------------------------------------------------------------*
 | tangential contact handler test suite                      sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
// clang-format off
class PARTICLEINTERACTION::DEMContactTangentialLinearSpringDamp_TestSuite : public CxxTest::TestSuite
// clang-format on
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::DEMContactTangentialLinearSpringDamp> contacttangential_;

  const double e_ = 0.8;
  const double nue_ = 0.4;
  const double mu_tangential_ = 0.2;

  const double k_normal_ = 4.0;

 public:
  void setUp()
  {
    // create a parameter list
    Teuchos::ParameterList params_dem;
    params_dem.set("COEFF_RESTITUTION", e_);
    params_dem.set("POISSON_RATIO", nue_);
    params_dem.set("FRICT_COEFF_TANG", mu_tangential_);

    // create tangential contact handler
    contacttangential_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactTangentialLinearSpringDamp>(
        new PARTICLEINTERACTION::DEMContactTangentialLinearSpringDamp(params_dem));

    // init tangential contact handler
    contacttangential_->Init();

    // setup tangential contact handler
    contacttangential_->Setup(k_normal_);
  }

  void tearDown() { contacttangential_ = nullptr; }

  // note: the public functions Init() and Setup() of class DEMContactTangentialLinearSpringDamp are
  // called in setUp() and thus implicitly tested by all following unittests

  void test_TangentialContactForce_stick()
  {
    double gap_tangential[3] = {0.0};
    gap_tangential[0] = 0.1;
    gap_tangential[1] = 0.05;
    gap_tangential[2] = -0.25;

    bool stick_tangential = true;

    double e_ji[3] = {0.0};
    e_ji[0] = 1.0 / std::sqrt(21);
    e_ji[1] = 2.0 / std::sqrt(21);
    e_ji[2] = 4.0 / std::sqrt(21);

    double vel_rel_tangential[3] = {0.0};
    vel_rel_tangential[0] = -0.03;
    vel_rel_tangential[1] = 0.1;
    vel_rel_tangential[2] = 0.12;

    const double m_eff = 2.5;
    const double normalcontactforce = 2.5e2;

    double tangentialcontactforce[3];
    contacttangential_->TangentialContactForce(gap_tangential, stick_tangential, e_ji,
        vel_rel_tangential, m_eff, mu_tangential_, normalcontactforce, tangentialcontactforce);

    double gap_tangential_ref[3] = {0.17923102017884, 0.16378007016342, -0.12669779012642};
    double tangentialcontactforce_ref[3] = {-0.52425016124507, -0.53614987479507, 0.32632177321349};

    for (int i = 0; i < 3; ++i) TS_ASSERT_DELTA(gap_tangential[i], gap_tangential_ref[i], 1.0e-14);

    for (int i = 0; i < 3; ++i)
      TS_ASSERT_DELTA(tangentialcontactforce[i], tangentialcontactforce_ref[i], 1.0e-14);

    TS_ASSERT_EQUALS(stick_tangential, true);
  }

  void test_TangentialContactForce_slip()
  {
    double gap_tangential[3] = {0.0};
    gap_tangential[0] = 0.1;
    gap_tangential[1] = 0.05;
    gap_tangential[2] = -0.25;

    bool stick_tangential = false;

    double e_ji[3] = {0.0};
    e_ji[0] = 1.0 / std::sqrt(21);
    e_ji[1] = 2.0 / std::sqrt(21);
    e_ji[2] = 4.0 / std::sqrt(21);

    double vel_rel_tangential[3] = {0.0};
    vel_rel_tangential[0] = -0.03;
    vel_rel_tangential[1] = 0.1;
    vel_rel_tangential[2] = 0.12;

    const double m_eff = 2.5;
    const double normalcontactforce = 1.5;

    double tangentialcontactforce[3];
    contacttangential_->TangentialContactForce(gap_tangential, stick_tangential, e_ji,
        vel_rel_tangential, m_eff, mu_tangential_, normalcontactforce, tangentialcontactforce);

    double gap_tangential_ref[3] = {0.068586669580933, 0.050624254285441, -0.057826737044788};
    double tangentialcontactforce_ref[3] = {-0.19231710945136, -0.19668242716113, 0.11970861396859};

    for (int i = 0; i < 3; ++i) TS_ASSERT_DELTA(gap_tangential[i], gap_tangential_ref[i], 1.0e-14);

    for (int i = 0; i < 3; ++i)
      TS_ASSERT_DELTA(tangentialcontactforce[i], tangentialcontactforce_ref[i], 1.0e-14);

    TS_ASSERT_EQUALS(stick_tangential, false);
  }

  void test_TangentialPotentialEnergy()
  {
    double gap_tangential[3] = {0.0};
    gap_tangential[0] = 0.1;
    gap_tangential[1] = 0.05;
    gap_tangential[2] = -0.25;

    const double k_tangential = (1.0 - nue_) / (1.0 - 0.5 * nue_) * k_normal_;
    const double tangentialpotentialenergy_ref =
        0.5 * k_tangential * UTILS::vec_dot(gap_tangential, gap_tangential);

    double tangentialpotentialenergy = 0.0;
    contacttangential_->TangentialPotentialEnergy(gap_tangential, tangentialpotentialenergy);

    TS_ASSERT_DELTA(tangentialpotentialenergy, tangentialpotentialenergy_ref, 1.0e-12);
  }
};

/*---------------------------------------------------------------------------*/
#endif
