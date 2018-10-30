/*---------------------------------------------------------------------------*/
/*!
\file unit_particle_interaction_sph_kernel.H

\brief unittests for kernel handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#ifndef UNIT_PARTICLE_INTERACTION_SPH_KERNEL_H
#define UNIT_PARTICLE_INTERACTION_SPH_KERNEL_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#include <cxxtest/TestSuite.h>
#include "src/drt_particle_interaction/particle_interaction_sph_kernel.H"
#include "src/drt_particle_interaction/particle_interaction_utils.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "src/drt_inpar/drt_validparameters.H"

/*---------------------------------------------------------------------------*
 | forward declaration                                        sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHKernelCubicSpline_TestSuite;
  class SPHKernelQuinticSpline_TestSuite;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | kernel handler test suite                                  sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
class PARTICLEINTERACTION::SPHKernelCubicSpline_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::SPHKernelCubicSpline> kernel_1D_;
  std::unique_ptr<PARTICLEINTERACTION::SPHKernelCubicSpline> kernel_2D_;
  std::unique_ptr<PARTICLEINTERACTION::SPHKernelCubicSpline> kernel_3D_;

 public:
  void setUp()
  {
    // create a parameter list
    Teuchos::ParameterList params_1D;
    Teuchos::ParameterList& params_1D_sph = params_1D.sublist("SPH");

    Teuchos::setStringToIntegralParameter<int>("KERNEL_SPACE_DIM", "Kernel1D",
        "kernel space dimension number", Teuchos::tuple<std::string>("Kernel1D"),
        Teuchos::tuple<int>(INPAR::PARTICLE::Kernel1D), &params_1D_sph);

    Teuchos::ParameterList params_2D;
    Teuchos::ParameterList& params_2D_sph = params_2D.sublist("SPH");

    Teuchos::setStringToIntegralParameter<int>("KERNEL_SPACE_DIM", "Kernel2D",
        "kernel space dimension number", Teuchos::tuple<std::string>("Kernel2D"),
        Teuchos::tuple<int>(INPAR::PARTICLE::Kernel2D), &params_2D_sph);

    Teuchos::ParameterList params_3D;
    Teuchos::ParameterList& params_3D_sph = params_3D.sublist("SPH");

    Teuchos::setStringToIntegralParameter<int>("KERNEL_SPACE_DIM", "Kernel3D",
        "kernel space dimension number", Teuchos::tuple<std::string>("Kernel3D"),
        Teuchos::tuple<int>(INPAR::PARTICLE::Kernel3D), &params_3D_sph);

    // create contact law handler
    kernel_1D_ = std::unique_ptr<PARTICLEINTERACTION::SPHKernelCubicSpline>(
        new PARTICLEINTERACTION::SPHKernelCubicSpline(params_1D));
    kernel_2D_ = std::unique_ptr<PARTICLEINTERACTION::SPHKernelCubicSpline>(
        new PARTICLEINTERACTION::SPHKernelCubicSpline(params_2D));
    kernel_3D_ = std::unique_ptr<PARTICLEINTERACTION::SPHKernelCubicSpline>(
        new PARTICLEINTERACTION::SPHKernelCubicSpline(params_3D));

    // init contact law handler
    kernel_1D_->Init();
    kernel_2D_->Init();
    kernel_3D_->Init();

    // setup contact law handler
    kernel_1D_->Setup();
    kernel_2D_->Setup();
    kernel_3D_->Setup();
  }

  void tearDown()
  {
    kernel_1D_ = nullptr;
    kernel_2D_ = nullptr;
    kernel_3D_ = nullptr;
  }

  // note: the public functions Init() and Setup() of class SPHKernelCubicSpline are called in
  // setUp()
  //       and thus implicitly tested by all following unittests

  void test_KernelSpaceDimension()
  {
    int dim = 0;

    kernel_1D_->KernelSpaceDimension(dim);
    TS_ASSERT_EQUALS(dim, 1);

    kernel_2D_->KernelSpaceDimension(dim);
    TS_ASSERT_EQUALS(dim, 2);

    kernel_3D_->KernelSpaceDimension(dim);
    TS_ASSERT_EQUALS(dim, 3);
  }

  void test_SmoothingLength()
  {
    const double support = 0.8;
    const double h = 0.4;

    TS_ASSERT_DELTA(kernel_1D_->SmoothingLength(support), h, 1e-14);
    TS_ASSERT_DELTA(kernel_2D_->SmoothingLength(support), h, 1e-14);
    TS_ASSERT_DELTA(kernel_3D_->SmoothingLength(support), h, 1e-14);
  }

  void test_NormalizationConstant()
  {
    const double h = 0.4;

    const double normalizationconstant_1D = 2.0 / (3.0 * h);
    const double normalizationconstant_2D = 10.0 * M_1_PI / (7.0 * UTILS::pow<2>(h));
    const double normalizationconstant_3D = M_1_PI / UTILS::pow<3>(h);

    TS_ASSERT_DELTA(kernel_1D_->NormalizationConstant(h), normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(kernel_2D_->NormalizationConstant(h), normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(kernel_3D_->NormalizationConstant(h), normalizationconstant_3D, 1e-14);
  }

  void test_W()
  {
    const double support = 0.8;
    const double h = 0.4;

    const double normalizationconstant_1D = 2.0 / (3.0 * h);
    const double normalizationconstant_2D = 10.0 * M_1_PI / (7.0 * UTILS::pow<2>(h));
    const double normalizationconstant_3D = M_1_PI / UTILS::pow<3>(h);

    double rij = 0.0;
    double q = rij / h;
    double w_unnormalized = 1.0;
    TS_ASSERT_DELTA(kernel_1D_->W(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(kernel_2D_->W(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(kernel_3D_->W(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.2;
    q = rij / h;
    w_unnormalized = 1.0 - 1.5 * UTILS::pow<2>(q) + 0.75 * UTILS::pow<3>(q);
    TS_ASSERT_DELTA(kernel_1D_->W(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(kernel_2D_->W(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(kernel_3D_->W(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.6;
    q = rij / h;
    w_unnormalized = UTILS::pow<3>(2.0 - q) / 4.0;
    TS_ASSERT_DELTA(kernel_1D_->W(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(kernel_2D_->W(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(kernel_3D_->W(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.8;
    q = rij / h;
    w_unnormalized = 0.0;
    TS_ASSERT_DELTA(kernel_1D_->W(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(kernel_2D_->W(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(kernel_3D_->W(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);
  }

  void test_dWdrij()
  {
    const double support = 0.8;
    const double h = 0.4;

    const double normalizationconstant_1D = 2.0 / (3.0 * h);
    const double normalizationconstant_2D = 10.0 * M_1_PI / (7.0 * UTILS::pow<2>(h));
    const double normalizationconstant_3D = M_1_PI / UTILS::pow<3>(h);

    double rij = 0.0;
    double q = rij / h;
    double w_unnormalized = 0.0;
    TS_ASSERT_DELTA(
        kernel_1D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_2D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_3D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.2;
    q = rij / h;
    w_unnormalized = (-3.0 * q + 2.25 * UTILS::pow<2>(q)) * (1.0 / h);
    TS_ASSERT_DELTA(
        kernel_1D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_2D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_3D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.6;
    q = rij / h;
    w_unnormalized = (-0.75 * UTILS::pow<2>(2.0 - q)) * (1.0 / h);
    TS_ASSERT_DELTA(
        kernel_1D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_2D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_3D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.8;
    q = rij / h;
    w_unnormalized = 0.0;
    TS_ASSERT_DELTA(
        kernel_1D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_2D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_3D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);
  }

  void test_d2Wdrij2()
  {
    const double support = 0.8;
    const double h = 0.4;

    const double normalizationconstant_1D = 2.0 / (3.0 * h);
    const double normalizationconstant_2D = 10.0 * M_1_PI / (7.0 * UTILS::pow<2>(h));
    const double normalizationconstant_3D = M_1_PI / UTILS::pow<3>(h);

    double rij = 0.0;
    double q = rij / h;
    double w_unnormalized = -3.0 * (1.0 / UTILS::pow<2>(h));
    TS_ASSERT_DELTA(
        kernel_1D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_2D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_3D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.2;
    q = rij / h;
    w_unnormalized = (-3.0 + 4.5 * q) * (1.0 / UTILS::pow<2>(h));
    TS_ASSERT_DELTA(
        kernel_1D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_2D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_3D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.6;
    q = rij / h;
    w_unnormalized = (1.5 * (2.0 - q)) * (1.0 / UTILS::pow<2>(h));
    TS_ASSERT_DELTA(
        kernel_1D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_2D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_3D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.8;
    q = rij / h;
    w_unnormalized = 0.0;
    TS_ASSERT_DELTA(
        kernel_1D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_2D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_3D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);
  }

  void test_GradWij()
  {
    const double rij = 0.2;
    const double support = 0.8;
    double eij[3];
    eij[0] = 0.5;
    eij[1] = std::sqrt(3.0) / 2.0;
    eij[2] = 0.0;

    const double h = 0.4;
    const double normalizationconstant_3D = M_1_PI / UTILS::pow<3>(h);
    const double q = rij / h;
    const double w_unnormalized = (-3.0 * q + 2.25 * UTILS::pow<2>(q)) * (1.0 / h);

    double gradWij_reference[3];
    for (int i = 0; i < 3; ++i)
      gradWij_reference[i] = w_unnormalized * normalizationconstant_3D * eij[i];


    double gradWij[3];
    kernel_3D_->GradWij(rij, support, eij, gradWij);

    for (int i = 0; i < 3; ++i) TS_ASSERT_DELTA(gradWij[i], gradWij_reference[i], 1e-14);
  }
};

class PARTICLEINTERACTION::SPHKernelQuinticSpline_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::SPHKernelQuinticSpline> kernel_1D_;
  std::unique_ptr<PARTICLEINTERACTION::SPHKernelQuinticSpline> kernel_2D_;
  std::unique_ptr<PARTICLEINTERACTION::SPHKernelQuinticSpline> kernel_3D_;

 public:
  void setUp()
  {
    // create a parameter list
    Teuchos::ParameterList params_1D;
    Teuchos::ParameterList& params_1D_sph = params_1D.sublist("SPH");

    Teuchos::setStringToIntegralParameter<int>("KERNEL_SPACE_DIM", "Kernel1D",
        "kernel space dimension number", Teuchos::tuple<std::string>("Kernel1D"),
        Teuchos::tuple<int>(INPAR::PARTICLE::Kernel1D), &params_1D_sph);

    Teuchos::ParameterList params_2D;
    Teuchos::ParameterList& params_2D_sph = params_2D.sublist("SPH");

    Teuchos::setStringToIntegralParameter<int>("KERNEL_SPACE_DIM", "Kernel2D",
        "kernel space dimension number", Teuchos::tuple<std::string>("Kernel2D"),
        Teuchos::tuple<int>(INPAR::PARTICLE::Kernel2D), &params_2D_sph);

    Teuchos::ParameterList params_3D;
    Teuchos::ParameterList& params_3D_sph = params_3D.sublist("SPH");

    Teuchos::setStringToIntegralParameter<int>("KERNEL_SPACE_DIM", "Kernel3D",
        "kernel space dimension number", Teuchos::tuple<std::string>("Kernel3D"),
        Teuchos::tuple<int>(INPAR::PARTICLE::Kernel3D), &params_3D_sph);

    // create contact law handler
    kernel_1D_ = std::unique_ptr<PARTICLEINTERACTION::SPHKernelQuinticSpline>(
        new PARTICLEINTERACTION::SPHKernelQuinticSpline(params_1D));
    kernel_2D_ = std::unique_ptr<PARTICLEINTERACTION::SPHKernelQuinticSpline>(
        new PARTICLEINTERACTION::SPHKernelQuinticSpline(params_2D));
    kernel_3D_ = std::unique_ptr<PARTICLEINTERACTION::SPHKernelQuinticSpline>(
        new PARTICLEINTERACTION::SPHKernelQuinticSpline(params_3D));

    // init contact law handler
    kernel_1D_->Init();
    kernel_2D_->Init();
    kernel_3D_->Init();

    // setup contact law handler
    kernel_1D_->Setup();
    kernel_2D_->Setup();
    kernel_3D_->Setup();
  }

  void tearDown()
  {
    kernel_1D_ = nullptr;
    kernel_2D_ = nullptr;
    kernel_3D_ = nullptr;
  }

  // note: the public functions Init() and Setup() of class SPHKernelQuinticSpline are called in
  // setUp()
  //       and thus implicitly tested by all following unittests

  void test_KernelSpaceDimension()
  {
    int dim = 0;

    kernel_1D_->KernelSpaceDimension(dim);
    TS_ASSERT_EQUALS(dim, 1);

    kernel_2D_->KernelSpaceDimension(dim);
    TS_ASSERT_EQUALS(dim, 2);

    kernel_3D_->KernelSpaceDimension(dim);
    TS_ASSERT_EQUALS(dim, 3);
  }

  void test_SmoothingLength()
  {
    const double support = 0.9;
    const double h = 0.3;

    TS_ASSERT_DELTA(kernel_1D_->SmoothingLength(support), h, 1e-14);
    TS_ASSERT_DELTA(kernel_2D_->SmoothingLength(support), h, 1e-14);
    TS_ASSERT_DELTA(kernel_3D_->SmoothingLength(support), h, 1e-14);
  }

  void test_NormalizationConstant()
  {
    const double h = 0.3;

    const double normalizationconstant_1D = 1.0 / (120.0 * h);
    const double normalizationconstant_2D = 7.0 * M_1_PI / (478.0 * UTILS::pow<2>(h));
    const double normalizationconstant_3D = 3.0 * M_1_PI / (359.0 * UTILS::pow<3>(h));

    TS_ASSERT_DELTA(kernel_1D_->NormalizationConstant(h), normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(kernel_2D_->NormalizationConstant(h), normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(kernel_3D_->NormalizationConstant(h), normalizationconstant_3D, 1e-14);
  }

  void test_W()
  {
    const double support = 0.9;
    const double h = 0.3;

    const double normalizationconstant_1D = 1.0 / (120.0 * h);
    const double normalizationconstant_2D = 7.0 * M_1_PI / (478.0 * UTILS::pow<2>(h));
    const double normalizationconstant_3D = 3.0 * M_1_PI / (359.0 * UTILS::pow<3>(h));

    double rij = 0.0;
    double q = rij / h;
    double w_unnormalized = 66.0;
    TS_ASSERT_DELTA(kernel_1D_->W(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(kernel_2D_->W(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(kernel_3D_->W(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.2;
    q = rij / h;
    w_unnormalized =
        UTILS::pow<5>(3.0 - q) - 6.0 * UTILS::pow<5>(2.0 - q) + 15.0 * UTILS::pow<5>(1.0 - q);
    TS_ASSERT_DELTA(kernel_1D_->W(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(kernel_2D_->W(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(kernel_3D_->W(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.5;
    q = rij / h;
    w_unnormalized = UTILS::pow<5>(3.0 - q) - 6.0 * UTILS::pow<5>(2.0 - q);
    TS_ASSERT_DELTA(kernel_1D_->W(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(kernel_2D_->W(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(kernel_3D_->W(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.8;
    q = rij / h;
    w_unnormalized = UTILS::pow<5>(3.0 - q);
    TS_ASSERT_DELTA(kernel_1D_->W(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(kernel_2D_->W(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(kernel_3D_->W(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.9;
    q = rij / h;
    w_unnormalized = 0.0;
    TS_ASSERT_DELTA(kernel_1D_->W(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(kernel_2D_->W(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(kernel_3D_->W(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);
  }

  void test_dWdrij()
  {
    const double support = 0.9;
    const double h = 0.3;

    const double normalizationconstant_1D = 1.0 / (120.0 * h);
    const double normalizationconstant_2D = 7.0 * M_1_PI / (478.0 * UTILS::pow<2>(h));
    const double normalizationconstant_3D = 3.0 * M_1_PI / (359.0 * UTILS::pow<3>(h));

    double rij = 0.0;
    double q = rij / h;
    double w_unnormalized = 0.0;
    TS_ASSERT_DELTA(
        kernel_1D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_2D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_3D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.2;
    q = rij / h;
    w_unnormalized = (-5.0 * UTILS::pow<4>(3.0 - q) + 30.0 * UTILS::pow<4>(2.0 - q) -
                         75.0 * UTILS::pow<4>(1.0 - q)) *
                     (1.0 / h);
    TS_ASSERT_DELTA(
        kernel_1D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_2D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_3D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.5;
    q = rij / h;
    w_unnormalized = (-5.0 * UTILS::pow<4>(3.0 - q) + 30.0 * UTILS::pow<4>(2.0 - q)) * (1.0 / h);
    TS_ASSERT_DELTA(
        kernel_1D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_2D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_3D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.8;
    q = rij / h;
    w_unnormalized = (-5.0 * UTILS::pow<4>(3.0 - q)) * (1.0 / h);
    TS_ASSERT_DELTA(
        kernel_1D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_2D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_3D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.9;
    q = rij / h;
    w_unnormalized = 0.0;
    TS_ASSERT_DELTA(
        kernel_1D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_2D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_3D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);
  }

  void test_d2Wdrij2()
  {
    const double support = 0.9;
    const double h = 0.3;

    const double normalizationconstant_1D = 1.0 / (120.0 * h);
    const double normalizationconstant_2D = 7.0 * M_1_PI / (478.0 * UTILS::pow<2>(h));
    const double normalizationconstant_3D = 3.0 * M_1_PI / (359.0 * UTILS::pow<3>(h));

    double rij = 0.0;
    double q = rij / h;
    double w_unnormalized = -4000.0 / 3.0;
    TS_ASSERT_DELTA(
        kernel_1D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_2D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_3D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.2;
    q = rij / h;
    w_unnormalized = (20.0 * UTILS::pow<3>(3.0 - q) - 120.0 * UTILS::pow<3>(2.0 - q) +
                         300.0 * UTILS::pow<3>(1.0 - q)) *
                     (1.0 / UTILS::pow<2>(h));
    TS_ASSERT_DELTA(
        kernel_1D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_2D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_3D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.5;
    q = rij / h;
    w_unnormalized =
        (20.0 * UTILS::pow<3>(3.0 - q) - 120.0 * UTILS::pow<3>(2.0 - q)) * (1.0 / UTILS::pow<2>(h));
    TS_ASSERT_DELTA(
        kernel_1D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_2D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_3D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.8;
    q = rij / h;
    w_unnormalized = (20.0 * UTILS::pow<3>(3.0 - q)) * (1.0 / UTILS::pow<2>(h));
    TS_ASSERT_DELTA(
        kernel_1D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_2D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_3D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);

    rij = 0.9;
    q = rij / h;
    w_unnormalized = 0.0;
    TS_ASSERT_DELTA(
        kernel_1D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_1D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_2D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_2D, 1e-14);
    TS_ASSERT_DELTA(
        kernel_3D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_3D, 1e-14);
  }

  void test_GradWij()
  {
    const double rij = 0.2;
    const double support = 0.9;
    double eij[3];
    eij[0] = 0.5;
    eij[1] = std::sqrt(3.0) / 2.0;
    eij[2] = 0.0;

    const double h = 0.3;
    const double normalizationconstant_3D = 3.0 * M_1_PI / (359.0 * UTILS::pow<3>(h));
    const double q = rij / h;
    const double w_unnormalized = (-5.0 * UTILS::pow<4>(3.0 - q) + 30.0 * UTILS::pow<4>(2.0 - q) -
                                      75.0 * UTILS::pow<4>(1.0 - q)) *
                                  (1.0 / h);

    double gradWij_reference[3];
    for (int i = 0; i < 3; ++i)
      gradWij_reference[i] = w_unnormalized * normalizationconstant_3D * eij[i];


    double gradWij[3];
    kernel_3D_->GradWij(rij, support, eij, gradWij);

    for (int i = 0; i < 3; ++i) TS_ASSERT_DELTA(gradWij[i], gradWij_reference[i], 1e-14);
  }
};

/*---------------------------------------------------------------------------*/
#endif
