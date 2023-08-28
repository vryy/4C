/*---------------------------------------------------------------------------*/
/*! \file
\brief unittests for kernel handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "baci_particle_interaction_sph_kernel.H"

#include "baci_inpar_validparameters.H"
#include "baci_particle_interaction_utils.H"
#include "baci_unittest_utils_assertions_test.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

namespace
{
  class SPHKernelCubicSplineTest : public ::testing::Test
  {
   protected:
    std::unique_ptr<PARTICLEINTERACTION::SPHKernelCubicSpline> kernel_1D_;
    std::unique_ptr<PARTICLEINTERACTION::SPHKernelCubicSpline> kernel_2D_;
    std::unique_ptr<PARTICLEINTERACTION::SPHKernelCubicSpline> kernel_3D_;

    SPHKernelCubicSplineTest()
    {
      // create a parameter list
      Teuchos::ParameterList params_sph_1D;
      Teuchos::setStringToIntegralParameter<int>("KERNEL_SPACE_DIM", "Kernel1D",
          "kernel space dimension number", Teuchos::tuple<std::string>("Kernel1D"),
          Teuchos::tuple<int>(INPAR::PARTICLE::Kernel1D), &params_sph_1D);

      Teuchos::ParameterList params_sph_2D;
      Teuchos::setStringToIntegralParameter<int>("KERNEL_SPACE_DIM", "Kernel2D",
          "kernel space dimension number", Teuchos::tuple<std::string>("Kernel2D"),
          Teuchos::tuple<int>(INPAR::PARTICLE::Kernel2D), &params_sph_2D);

      Teuchos::ParameterList params_sph_3D;
      Teuchos::setStringToIntegralParameter<int>("KERNEL_SPACE_DIM", "Kernel3D",
          "kernel space dimension number", Teuchos::tuple<std::string>("Kernel3D"),
          Teuchos::tuple<int>(INPAR::PARTICLE::Kernel3D), &params_sph_3D);

      // create kernel handler
      kernel_1D_ = std::make_unique<PARTICLEINTERACTION::SPHKernelCubicSpline>(params_sph_1D);
      kernel_2D_ = std::make_unique<PARTICLEINTERACTION::SPHKernelCubicSpline>(params_sph_2D);
      kernel_3D_ = std::make_unique<PARTICLEINTERACTION::SPHKernelCubicSpline>(params_sph_3D);

      // init kernel handler
      kernel_1D_->Init();
      kernel_2D_->Init();
      kernel_3D_->Init();

      // setup kernel handler
      kernel_1D_->Setup();
      kernel_2D_->Setup();
      kernel_3D_->Setup();
    }
    // note: the public functions Init() and Setup() of class SPHKernelCubicSpline are called in
    // Setup() and thus implicitly tested by all following unittests
  };

  TEST_F(SPHKernelCubicSplineTest, KernelSpaceDimension)
  {
    int dim = 0;

    kernel_1D_->KernelSpaceDimension(dim);
    EXPECT_EQ(dim, 1);

    kernel_2D_->KernelSpaceDimension(dim);
    EXPECT_EQ(dim, 2);

    kernel_3D_->KernelSpaceDimension(dim);
    EXPECT_EQ(dim, 3);
  }

  TEST_F(SPHKernelCubicSplineTest, SmoothingLength)
  {
    const double support = 0.8;
    const double h = 0.4;

    EXPECT_NEAR(kernel_1D_->SmoothingLength(support), h, 1.0e-10);
    EXPECT_NEAR(kernel_2D_->SmoothingLength(support), h, 1.0e-10);
    EXPECT_NEAR(kernel_3D_->SmoothingLength(support), h, 1.0e-10);
  }

  TEST_F(SPHKernelCubicSplineTest, NormalizationConstant)
  {
    const double h = 0.4;
    const double inv_h = 1.0 / h;

    const double normalizationconstant_1D = 2.0 / (3.0 * h);
    const double normalizationconstant_2D =
        10.0 * M_1_PI / (7.0 * PARTICLEINTERACTION::UTILS::Pow<2>(h));
    const double normalizationconstant_3D = M_1_PI / PARTICLEINTERACTION::UTILS::Pow<3>(h);

    EXPECT_NEAR(kernel_1D_->NormalizationConstant(inv_h), normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(kernel_2D_->NormalizationConstant(inv_h), normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(kernel_3D_->NormalizationConstant(inv_h), normalizationconstant_3D, 1.0e-10);
  }

  TEST_F(SPHKernelCubicSplineTest, W0)
  {
    const double support = 0.8;
    const double h = 0.4;

    const double normalizationconstant_1D = 2.0 / (3.0 * h);
    const double normalizationconstant_2D =
        10.0 * M_1_PI / (7.0 * PARTICLEINTERACTION::UTILS::Pow<2>(h));
    const double normalizationconstant_3D = M_1_PI / PARTICLEINTERACTION::UTILS::Pow<3>(h);

    double w_unnormalized = 1.0;
    EXPECT_NEAR(kernel_1D_->W0(support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(kernel_2D_->W0(support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(kernel_3D_->W0(support), w_unnormalized * normalizationconstant_3D, 1.0e-10);
  }

  TEST_F(SPHKernelCubicSplineTest, W)
  {
    const double support = 0.8;
    const double h = 0.4;

    const double normalizationconstant_1D = 2.0 / (3.0 * h);
    const double normalizationconstant_2D =
        10.0 * M_1_PI / (7.0 * PARTICLEINTERACTION::UTILS::Pow<2>(h));
    const double normalizationconstant_3D = M_1_PI / PARTICLEINTERACTION::UTILS::Pow<3>(h);

    double rij = 0.0;
    double q = rij / h;
    double w_unnormalized = 1.0;
    EXPECT_NEAR(kernel_1D_->W(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(kernel_2D_->W(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(kernel_3D_->W(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.2;
    q = rij / h;
    w_unnormalized = 1.0 - 1.5 * PARTICLEINTERACTION::UTILS::Pow<2>(q) +
                     0.75 * PARTICLEINTERACTION::UTILS::Pow<3>(q);
    EXPECT_NEAR(kernel_1D_->W(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(kernel_2D_->W(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(kernel_3D_->W(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.6;
    q = rij / h;
    w_unnormalized = PARTICLEINTERACTION::UTILS::Pow<3>(2.0 - q) / 4.0;
    EXPECT_NEAR(kernel_1D_->W(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(kernel_2D_->W(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(kernel_3D_->W(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.8;
    q = rij / h;
    w_unnormalized = 0.0;
    EXPECT_NEAR(kernel_1D_->W(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(kernel_2D_->W(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(kernel_3D_->W(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);
  }

  TEST_F(SPHKernelCubicSplineTest, dWdrij)
  {
    const double support = 0.8;
    const double h = 0.4;

    const double normalizationconstant_1D = 2.0 / (3.0 * h);
    const double normalizationconstant_2D =
        10.0 * M_1_PI / (7.0 * PARTICLEINTERACTION::UTILS::Pow<2>(h));
    const double normalizationconstant_3D = M_1_PI / PARTICLEINTERACTION::UTILS::Pow<3>(h);

    double rij = 0.0;
    double q = rij / h;
    double w_unnormalized = 0.0;
    EXPECT_NEAR(
        kernel_1D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(
        kernel_2D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(
        kernel_3D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.2;
    q = rij / h;
    w_unnormalized = (-3.0 * q + 2.25 * PARTICLEINTERACTION::UTILS::Pow<2>(q)) * (1.0 / h);
    EXPECT_NEAR(
        kernel_1D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(
        kernel_2D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(
        kernel_3D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.6;
    q = rij / h;
    w_unnormalized = (-0.75 * PARTICLEINTERACTION::UTILS::Pow<2>(2.0 - q)) * (1.0 / h);
    EXPECT_NEAR(
        kernel_1D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(
        kernel_2D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(
        kernel_3D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.8;
    q = rij / h;
    w_unnormalized = 0.0;
    EXPECT_NEAR(
        kernel_1D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(
        kernel_2D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(
        kernel_3D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);
  }

  TEST_F(SPHKernelCubicSplineTest, d2Wdrij2)
  {
    const double support = 0.8;
    const double h = 0.4;

    const double normalizationconstant_1D = 2.0 / (3.0 * h);
    const double normalizationconstant_2D =
        10.0 * M_1_PI / (7.0 * PARTICLEINTERACTION::UTILS::Pow<2>(h));
    const double normalizationconstant_3D = M_1_PI / PARTICLEINTERACTION::UTILS::Pow<3>(h);

    double rij = 0.0;
    double q = rij / h;
    double w_unnormalized = -3.0 * (1.0 / PARTICLEINTERACTION::UTILS::Pow<2>(h));
    EXPECT_NEAR(
        kernel_1D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(
        kernel_2D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(
        kernel_3D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.2;
    q = rij / h;
    w_unnormalized = (-3.0 + 4.5 * q) * (1.0 / PARTICLEINTERACTION::UTILS::Pow<2>(h));
    EXPECT_NEAR(
        kernel_1D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(
        kernel_2D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(
        kernel_3D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.6;
    q = rij / h;
    w_unnormalized = (1.5 * (2.0 - q)) * (1.0 / PARTICLEINTERACTION::UTILS::Pow<2>(h));
    EXPECT_NEAR(
        kernel_1D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(
        kernel_2D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(
        kernel_3D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.8;
    q = rij / h;
    w_unnormalized = 0.0;
    EXPECT_NEAR(
        kernel_1D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(
        kernel_2D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(
        kernel_3D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);
  }

  TEST_F(SPHKernelCubicSplineTest, GradWij)
  {
    const double rij = 0.2;
    const double support = 0.8;
    double eij[3];
    eij[0] = 0.5;
    eij[1] = std::sqrt(3.0) / 2.0;
    eij[2] = 0.0;

    const double h = 0.4;
    const double normalizationconstant_3D = M_1_PI / PARTICLEINTERACTION::UTILS::Pow<3>(h);
    const double q = rij / h;
    const double w_unnormalized =
        (-3.0 * q + 2.25 * PARTICLEINTERACTION::UTILS::Pow<2>(q)) * (1.0 / h);

    double gradWij_reference[3];
    for (int i = 0; i < 3; ++i)
      gradWij_reference[i] = w_unnormalized * normalizationconstant_3D * eij[i];


    double gradWij[3];
    kernel_3D_->GradWij(rij, support, eij, gradWij);

    BACI_EXPECT_ITERABLE_NEAR(gradWij, gradWij_reference, 3, 1.0e-10);
  }


  class SPHKernelQuinticSplineTest : public ::testing::Test
  {
   protected:
    std::unique_ptr<PARTICLEINTERACTION::SPHKernelQuinticSpline> kernel_1D_;
    std::unique_ptr<PARTICLEINTERACTION::SPHKernelQuinticSpline> kernel_2D_;
    std::unique_ptr<PARTICLEINTERACTION::SPHKernelQuinticSpline> kernel_3D_;

    SPHKernelQuinticSplineTest()
    {
      // create a parameter list
      Teuchos::ParameterList params_sph_1D;
      Teuchos::setStringToIntegralParameter<int>("KERNEL_SPACE_DIM", "Kernel1D",
          "kernel space dimension number", Teuchos::tuple<std::string>("Kernel1D"),
          Teuchos::tuple<int>(INPAR::PARTICLE::Kernel1D), &params_sph_1D);

      Teuchos::ParameterList params_sph_2D;
      Teuchos::setStringToIntegralParameter<int>("KERNEL_SPACE_DIM", "Kernel2D",
          "kernel space dimension number", Teuchos::tuple<std::string>("Kernel2D"),
          Teuchos::tuple<int>(INPAR::PARTICLE::Kernel2D), &params_sph_2D);

      Teuchos::ParameterList params_sph_3D;
      Teuchos::setStringToIntegralParameter<int>("KERNEL_SPACE_DIM", "Kernel3D",
          "kernel space dimension number", Teuchos::tuple<std::string>("Kernel3D"),
          Teuchos::tuple<int>(INPAR::PARTICLE::Kernel3D), &params_sph_3D);

      // create kernel handler
      kernel_1D_ = std::make_unique<PARTICLEINTERACTION::SPHKernelQuinticSpline>(params_sph_1D);
      kernel_2D_ = std::make_unique<PARTICLEINTERACTION::SPHKernelQuinticSpline>(params_sph_2D);
      kernel_3D_ = std::make_unique<PARTICLEINTERACTION::SPHKernelQuinticSpline>(params_sph_3D);

      // init kernel handler
      kernel_1D_->Init();
      kernel_2D_->Init();
      kernel_3D_->Init();

      // setup kernel handler
      kernel_1D_->Setup();
      kernel_2D_->Setup();
      kernel_3D_->Setup();
    }
    // note: the public functions Init() and Setup() of class SPHKernelQuinticSpline are called in
    // Setup() and thus implicitly tested by all following unittests
  };

  TEST_F(SPHKernelQuinticSplineTest, KernelSpaceDimension)
  {
    int dim = 0;

    kernel_1D_->KernelSpaceDimension(dim);
    EXPECT_EQ(dim, 1);

    kernel_2D_->KernelSpaceDimension(dim);
    EXPECT_EQ(dim, 2);

    kernel_3D_->KernelSpaceDimension(dim);
    EXPECT_EQ(dim, 3);
  }

  TEST_F(SPHKernelQuinticSplineTest, SmoothingLength)
  {
    const double support = 0.9;
    const double h = 0.3;

    EXPECT_NEAR(kernel_1D_->SmoothingLength(support), h, 1.0e-10);
    EXPECT_NEAR(kernel_2D_->SmoothingLength(support), h, 1.0e-10);
    EXPECT_NEAR(kernel_3D_->SmoothingLength(support), h, 1.0e-10);
  }

  TEST_F(SPHKernelQuinticSplineTest, NormalizationConstant)
  {
    const double h = 0.3;
    const double inv_h = 1.0 / h;

    const double normalizationconstant_1D = 1.0 / (120.0 * h);
    const double normalizationconstant_2D =
        7.0 * M_1_PI / (478.0 * PARTICLEINTERACTION::UTILS::Pow<2>(h));
    const double normalizationconstant_3D =
        3.0 * M_1_PI / (359.0 * PARTICLEINTERACTION::UTILS::Pow<3>(h));

    EXPECT_NEAR(kernel_1D_->NormalizationConstant(inv_h), normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(kernel_2D_->NormalizationConstant(inv_h), normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(kernel_3D_->NormalizationConstant(inv_h), normalizationconstant_3D, 1.0e-10);
  }

  TEST_F(SPHKernelQuinticSplineTest, W0)
  {
    const double support = 0.9;
    const double h = 0.3;

    const double normalizationconstant_1D = 1.0 / (120.0 * h);
    const double normalizationconstant_2D =
        7.0 * M_1_PI / (478.0 * PARTICLEINTERACTION::UTILS::Pow<2>(h));
    const double normalizationconstant_3D =
        3.0 * M_1_PI / (359.0 * PARTICLEINTERACTION::UTILS::Pow<3>(h));

    double w_unnormalized = 66.0;
    EXPECT_NEAR(kernel_1D_->W0(support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(kernel_2D_->W0(support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(kernel_3D_->W0(support), w_unnormalized * normalizationconstant_3D, 1.0e-10);
  }

  TEST_F(SPHKernelQuinticSplineTest, W)
  {
    const double support = 0.9;
    const double h = 0.3;

    const double normalizationconstant_1D = 1.0 / (120.0 * h);
    const double normalizationconstant_2D =
        7.0 * M_1_PI / (478.0 * PARTICLEINTERACTION::UTILS::Pow<2>(h));
    const double normalizationconstant_3D =
        3.0 * M_1_PI / (359.0 * PARTICLEINTERACTION::UTILS::Pow<3>(h));

    double rij = 0.0;
    double q = rij / h;
    double w_unnormalized = 66.0;
    EXPECT_NEAR(kernel_1D_->W(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(kernel_2D_->W(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(kernel_3D_->W(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.2;
    q = rij / h;
    w_unnormalized = PARTICLEINTERACTION::UTILS::Pow<5>(3.0 - q) -
                     6.0 * PARTICLEINTERACTION::UTILS::Pow<5>(2.0 - q) +
                     15.0 * PARTICLEINTERACTION::UTILS::Pow<5>(1.0 - q);
    EXPECT_NEAR(kernel_1D_->W(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(kernel_2D_->W(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(kernel_3D_->W(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.5;
    q = rij / h;
    w_unnormalized = PARTICLEINTERACTION::UTILS::Pow<5>(3.0 - q) -
                     6.0 * PARTICLEINTERACTION::UTILS::Pow<5>(2.0 - q);
    EXPECT_NEAR(kernel_1D_->W(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(kernel_2D_->W(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(kernel_3D_->W(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.8;
    q = rij / h;
    w_unnormalized = PARTICLEINTERACTION::UTILS::Pow<5>(3.0 - q);
    EXPECT_NEAR(kernel_1D_->W(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(kernel_2D_->W(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(kernel_3D_->W(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.9;
    q = rij / h;
    w_unnormalized = 0.0;
    EXPECT_NEAR(kernel_1D_->W(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(kernel_2D_->W(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(kernel_3D_->W(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);
  }

  TEST_F(SPHKernelQuinticSplineTest, dWdrij)
  {
    const double support = 0.9;
    const double h = 0.3;

    const double normalizationconstant_1D = 1.0 / (120.0 * h);
    const double normalizationconstant_2D =
        7.0 * M_1_PI / (478.0 * PARTICLEINTERACTION::UTILS::Pow<2>(h));
    const double normalizationconstant_3D =
        3.0 * M_1_PI / (359.0 * PARTICLEINTERACTION::UTILS::Pow<3>(h));

    double rij = 0.0;
    double q = rij / h;
    double w_unnormalized = 0.0;
    EXPECT_NEAR(
        kernel_1D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(
        kernel_2D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(
        kernel_3D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.2;
    q = rij / h;
    w_unnormalized = (-5.0 * PARTICLEINTERACTION::UTILS::Pow<4>(3.0 - q) +
                         30.0 * PARTICLEINTERACTION::UTILS::Pow<4>(2.0 - q) -
                         75.0 * PARTICLEINTERACTION::UTILS::Pow<4>(1.0 - q)) *
                     (1.0 / h);
    EXPECT_NEAR(
        kernel_1D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(
        kernel_2D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(
        kernel_3D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.5;
    q = rij / h;
    w_unnormalized = (-5.0 * PARTICLEINTERACTION::UTILS::Pow<4>(3.0 - q) +
                         30.0 * PARTICLEINTERACTION::UTILS::Pow<4>(2.0 - q)) *
                     (1.0 / h);
    EXPECT_NEAR(
        kernel_1D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(
        kernel_2D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(
        kernel_3D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.8;
    q = rij / h;
    w_unnormalized = (-5.0 * PARTICLEINTERACTION::UTILS::Pow<4>(3.0 - q)) * (1.0 / h);
    EXPECT_NEAR(
        kernel_1D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(
        kernel_2D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(
        kernel_3D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.9;
    q = rij / h;
    w_unnormalized = 0.0;
    EXPECT_NEAR(
        kernel_1D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(
        kernel_2D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(
        kernel_3D_->dWdrij(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);
  }

  TEST_F(SPHKernelQuinticSplineTest, d2Wdrij2)
  {
    const double support = 0.9;
    const double h = 0.3;

    const double normalizationconstant_1D = 1.0 / (120.0 * h);
    const double normalizationconstant_2D =
        7.0 * M_1_PI / (478.0 * PARTICLEINTERACTION::UTILS::Pow<2>(h));
    const double normalizationconstant_3D =
        3.0 * M_1_PI / (359.0 * PARTICLEINTERACTION::UTILS::Pow<3>(h));

    double rij = 0.0;
    double q = rij / h;
    double w_unnormalized = -4000.0 / 3.0;
    EXPECT_NEAR(
        kernel_1D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(
        kernel_2D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(
        kernel_3D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.2;
    q = rij / h;
    w_unnormalized = (20.0 * PARTICLEINTERACTION::UTILS::Pow<3>(3.0 - q) -
                         120.0 * PARTICLEINTERACTION::UTILS::Pow<3>(2.0 - q) +
                         300.0 * PARTICLEINTERACTION::UTILS::Pow<3>(1.0 - q)) *
                     (1.0 / PARTICLEINTERACTION::UTILS::Pow<2>(h));
    EXPECT_NEAR(
        kernel_1D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(
        kernel_2D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(
        kernel_3D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.5;
    q = rij / h;
    w_unnormalized = (20.0 * PARTICLEINTERACTION::UTILS::Pow<3>(3.0 - q) -
                         120.0 * PARTICLEINTERACTION::UTILS::Pow<3>(2.0 - q)) *
                     (1.0 / PARTICLEINTERACTION::UTILS::Pow<2>(h));
    EXPECT_NEAR(
        kernel_1D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(
        kernel_2D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(
        kernel_3D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.8;
    q = rij / h;
    w_unnormalized = (20.0 * PARTICLEINTERACTION::UTILS::Pow<3>(3.0 - q)) *
                     (1.0 / PARTICLEINTERACTION::UTILS::Pow<2>(h));
    EXPECT_NEAR(
        kernel_1D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(
        kernel_2D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(
        kernel_3D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);

    rij = 0.9;
    q = rij / h;
    w_unnormalized = 0.0;
    EXPECT_NEAR(
        kernel_1D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_1D, 1.0e-10);
    EXPECT_NEAR(
        kernel_2D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_2D, 1.0e-10);
    EXPECT_NEAR(
        kernel_3D_->d2Wdrij2(rij, support), w_unnormalized * normalizationconstant_3D, 1.0e-10);
  }

  TEST_F(SPHKernelQuinticSplineTest, GradWij)
  {
    const double rij = 0.2;
    const double support = 0.9;
    double eij[3];
    eij[0] = 0.5;
    eij[1] = std::sqrt(3.0) / 2.0;
    eij[2] = 0.0;

    const double h = 0.3;
    const double normalizationconstant_3D =
        3.0 * M_1_PI / (359.0 * PARTICLEINTERACTION::UTILS::Pow<3>(h));
    const double q = rij / h;
    const double w_unnormalized = (-5.0 * PARTICLEINTERACTION::UTILS::Pow<4>(3.0 - q) +
                                      30.0 * PARTICLEINTERACTION::UTILS::Pow<4>(2.0 - q) -
                                      75.0 * PARTICLEINTERACTION::UTILS::Pow<4>(1.0 - q)) *
                                  (1.0 / h);

    double gradWij_reference[3];
    for (int i = 0; i < 3; ++i)
      gradWij_reference[i] = w_unnormalized * normalizationconstant_3D * eij[i];


    double gradWij[3];
    kernel_3D_->GradWij(rij, support, eij, gradWij);

    BACI_EXPECT_ITERABLE_NEAR(gradWij, gradWij_reference, 3, 1.0e-10);
  }
}  // namespace
