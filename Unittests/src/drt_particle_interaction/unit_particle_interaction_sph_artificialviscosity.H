/*---------------------------------------------------------------------------*/
/*!
\file unit_particle_interaction_sph_artificialviscosity.H

\brief unittests for artificial viscosity handler for smoothed particle hydrodynamics (SPH)
interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#ifndef UNIT_PARTICLE_INTERACTION_SPH_ARTIFICIALVISCOSITY_H
#define UNIT_PARTICLE_INTERACTION_SPH_ARTIFICIALVISCOSITY_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include <cxxtest/TestSuite.h>
#include "src/drt_particle_interaction/particle_interaction_sph_artificialviscosity.H"

/*---------------------------------------------------------------------------*
 | forward declaration                                        sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHArtificialViscosity_TestSuite;
}

/*---------------------------------------------------------------------------*
 | artificial viscosity handler test suite                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
class PARTICLEINTERACTION::SPHArtificialViscosity_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::SPHArtificialViscosity> artificialviscosity_;

 public:
  void setUp()
  {
    // create artificial viscosity handler
    artificialviscosity_ = std::unique_ptr<PARTICLEINTERACTION::SPHArtificialViscosity>(
        new PARTICLEINTERACTION::SPHArtificialViscosity());

    // init artificial viscosity handler
    artificialviscosity_->Init();

    // setup artificial viscosity handler
    artificialviscosity_->Setup();
  }

  void tearDown() { artificialviscosity_ = nullptr; }

  // note: the public functions Init() and Setup() of class SPHEquationOfStateGenTait are called in
  // setUp()
  //       and thus implicitly tested by all following unittests

  void test_ArtificialViscosity()
  {
    const double dens_i = 1.01;
    const double dens_j = 0.97;

    double vel_i[3];
    vel_i[0] = 0.2;
    vel_i[1] = 0.3;
    vel_i[2] = -0.12;
    double vel_j[3];
    vel_j[0] = -0.12;
    vel_j[1] = -0.97;
    vel_j[2] = 0.98;

    const double mass_j = 10.32;
    const double artificialviscosity = 0.1;
    const double dWdrij = 0.76;
    const double h_i = 0.2;
    const double h_j = 0.25;
    const double c_i = 10.0;
    const double c_j = 12.5;
    const double abs_rij = 0.3;

    double e_ij[3];
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);

    double acc_i[3];
    acc_i[0] = 0.0;
    acc_i[1] = 0.0;
    acc_i[2] = 0.0;

    artificialviscosity_->ArtificialViscosity(&dens_i, &dens_j, vel_i, vel_j, &mass_j,
        artificialviscosity, dWdrij, h_i, h_j, c_i, c_j, abs_rij, e_ij, acc_i);

    // compute reference solution
    const double h_ij = 0.5 * (h_i + h_j);
    const double c_ij = 0.5 * (c_i + c_j);
    const double dens_ij = 0.5 * (dens_i + dens_j);
    const double e_ij_vrel_ij = ((vel_i[0] - vel_j[0]) * e_ij[0] + (vel_i[1] - vel_j[1]) * e_ij[1] +
                                 (vel_i[2] - vel_j[2]) * e_ij[2]);
    const double epsilon = 0.01;
    const double fac = h_ij * c_ij * e_ij_vrel_ij * abs_rij /
                       (dens_ij * (std::pow(abs_rij, 2) + epsilon * std::pow(h_ij, 2)));

    double acc_i_ref[3];
    for (int i = 0; i < 3; ++i)
      acc_i_ref[i] = artificialviscosity * mass_j * dWdrij * fac * e_ij[i];

    // compare results
    for (int i = 0; i < 3; ++i) TS_ASSERT_EQUALS(acc_i[i], acc_i_ref[i]);
  }
};

/*---------------------------------------------------------------------------*/
#endif
