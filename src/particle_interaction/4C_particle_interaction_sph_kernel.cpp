/*---------------------------------------------------------------------------*/
/*! \file
\brief kernel handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_interaction_sph_kernel.hpp"

#include "4C_particle_interaction_utils.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
ParticleInteraction::SPHKernelBase::SPHKernelBase(const Teuchos::ParameterList& params)
    : kernelspacedim_(Core::UTILS::IntegralValue<Inpar::PARTICLE::KernelSpaceDimension>(
          params, "KERNEL_SPACE_DIM"))
{
  // empty constructor
}

void ParticleInteraction::SPHKernelBase::init()
{
  // nothing to do
}

void ParticleInteraction::SPHKernelBase::setup()
{
  // nothing to do
}

void ParticleInteraction::SPHKernelBase::kernel_space_dimension(int& dim) const
{
  switch (kernelspacedim_)
  {
    case Inpar::PARTICLE::Kernel1D:
    {
      dim = 1;
      break;
    }
    case Inpar::PARTICLE::Kernel2D:
    {
      dim = 2;
      break;
    }
    case Inpar::PARTICLE::Kernel3D:
    {
      dim = 3;
      break;
    }
    default:
    {
      FOUR_C_THROW("unknown kernel space dimension!");
      break;
    }
  }
}

void ParticleInteraction::SPHKernelBase::grad_wij(
    const double& rij, const double& support, const double* eij, double* gradWij) const
{
  UTILS::VecSetScale(gradWij, this->d_wdrij(rij, support), eij);
}

ParticleInteraction::SPHKernelCubicSpline::SPHKernelCubicSpline(
    const Teuchos::ParameterList& params)
    : ParticleInteraction::SPHKernelBase(params)
{
  // empty constructor
}

double ParticleInteraction::SPHKernelCubicSpline::smoothing_length(const double& support) const
{
  return (0.5 * support);
}

double ParticleInteraction::SPHKernelCubicSpline::normalization_constant(const double& inv_h) const
{
  switch (kernelspacedim_)
  {
    case Inpar::PARTICLE::Kernel1D:
    {
      // (2.0 / 3.0) * inv_h
      return 0.6666666666666666 * inv_h;
    }
    case Inpar::PARTICLE::Kernel2D:
    {
      // (10.0 / 7.0) * M_1_PI * inv_h * inv_h
      return 0.4547284088339866 * UTILS::Pow<2>(inv_h);
    }
    case Inpar::PARTICLE::Kernel3D:
    {
      return M_1_PI * UTILS::Pow<3>(inv_h);
    }
    default:
    {
      FOUR_C_THROW("unknown kernel space dimension!");
      break;
    }
  }

  return 0.0;
}

double ParticleInteraction::SPHKernelCubicSpline::w0(const double& support) const
{
  return normalization_constant(2.0 / support);
}

double ParticleInteraction::SPHKernelCubicSpline::w(const double& rij, const double& support) const
{
  const double inv_h = 2.0 / support;
  const double q = rij * inv_h;

  if (q < 1.0)
    return (1.0 - 1.5 * UTILS::Pow<2>(q) + 0.75 * UTILS::Pow<3>(q)) * normalization_constant(inv_h);
  else if (q < 2.0)
    return (0.25 * UTILS::Pow<3>(2.0 - q)) * normalization_constant(inv_h);
  else
    return 0.0;
}

double ParticleInteraction::SPHKernelCubicSpline::d_wdrij(
    const double& rij, const double& support) const
{
  const double inv_h = 2.0 / support;
  const double q = rij * inv_h;

  if (q < 1.0)
    return (-3.0 * q + 2.25 * UTILS::Pow<2>(q)) * inv_h * normalization_constant(inv_h);
  else if (q < 2.0)
    return (-0.75 * UTILS::Pow<2>(2.0 - q)) * inv_h * normalization_constant(inv_h);
  else
    return 0.0;
}

double ParticleInteraction::SPHKernelCubicSpline::d2_wdrij2(
    const double& rij, const double& support) const
{
  const double inv_h = 2.0 / support;
  const double q = rij * inv_h;

  if (q < 1.0)
    return (-3.0 + 4.5 * q) * UTILS::Pow<2>(inv_h) * normalization_constant(inv_h);
  else if (q < 2.0)
    return (1.5 * (2.0 - q)) * UTILS::Pow<2>(inv_h) * normalization_constant(inv_h);
  else
    return 0.0;
}

ParticleInteraction::SPHKernelQuinticSpline::SPHKernelQuinticSpline(
    const Teuchos::ParameterList& params)
    : ParticleInteraction::SPHKernelBase(params)
{
  // empty constructor
}

double ParticleInteraction::SPHKernelQuinticSpline::smoothing_length(const double& support) const
{
  // (support / 3.0)
  return 0.3333333333333333 * support;
}

double ParticleInteraction::SPHKernelQuinticSpline::normalization_constant(
    const double& inv_h) const
{
  switch (kernelspacedim_)
  {
    case Inpar::PARTICLE::Kernel1D:
    {
      // (inv_h / 120.0)
      return 0.0083333333333333 * inv_h;
    }
    case Inpar::PARTICLE::Kernel2D:
    {
      // (7.0 / 478.0) * M_1_PI * inv_h * inv_h
      return 0.0046614418478797 * UTILS::Pow<2>(inv_h);
    }
    case Inpar::PARTICLE::Kernel3D:
    {
      // (3.0 / 359.0) * M_1_PI * inv_h * inv_h * inv_h
      return 0.0026599711937364 * UTILS::Pow<3>(inv_h);
    }
    default:
    {
      FOUR_C_THROW("unknown kernel space dimension!");
      break;
    }
  }

  return 0.0;
}

double ParticleInteraction::SPHKernelQuinticSpline::w0(const double& support) const
{
  return 66.0 * normalization_constant(3.0 / support);
}

double ParticleInteraction::SPHKernelQuinticSpline::w(
    const double& rij, const double& support) const
{
  const double inv_h = 3.0 / support;
  const double q = rij * inv_h;

  if (q < 1.0)
    return (UTILS::Pow<5>(3.0 - q) - 6.0 * UTILS::Pow<5>(2.0 - q) + 15.0 * UTILS::Pow<5>(1.0 - q)) *
           normalization_constant(inv_h);
  else if (q < 2.0)
    return (UTILS::Pow<5>(3.0 - q) - 6.0 * UTILS::Pow<5>(2.0 - q)) * normalization_constant(inv_h);
  else if (q < 3.0)
    return UTILS::Pow<5>(3.0 - q) * normalization_constant(inv_h);
  else
    return 0.0;
}

double ParticleInteraction::SPHKernelQuinticSpline::d_wdrij(
    const double& rij, const double& support) const
{
  const double inv_h = 3.0 / support;
  const double q = rij * inv_h;

  if (q < 1.0)
    return (-5.0 * UTILS::Pow<4>(3.0 - q) + 30.0 * UTILS::Pow<4>(2.0 - q) -
               75.0 * UTILS::Pow<4>(1.0 - q)) *
           inv_h * normalization_constant(inv_h);
  else if (q < 2.0)
    return (-5.0 * UTILS::Pow<4>(3.0 - q) + 30.0 * UTILS::Pow<4>(2.0 - q)) * inv_h *
           normalization_constant(inv_h);
  else if (q < 3.0)
    return (-5.0 * UTILS::Pow<4>(3.0 - q)) * inv_h * normalization_constant(inv_h);
  else
    return 0.0;
}

double ParticleInteraction::SPHKernelQuinticSpline::d2_wdrij2(
    const double& rij, const double& support) const
{
  const double inv_h = 3.0 / support;
  const double q = rij * inv_h;

  if (q < 1.0)
    return (20.0 * UTILS::Pow<3>(3.0 - q) - 120.0 * UTILS::Pow<3>(2.0 - q) +
               300.0 * UTILS::Pow<3>(1.0 - q)) *
           UTILS::Pow<2>(inv_h) * normalization_constant(inv_h);
  else if (q < 2.0)
    return (20.0 * UTILS::Pow<3>(3.0 - q) - 120.0 * UTILS::Pow<3>(2.0 - q)) * UTILS::Pow<2>(inv_h) *
           normalization_constant(inv_h);
  else if (q < 3.0)
    return (20.0 * UTILS::Pow<3>(3.0 - q)) * UTILS::Pow<2>(inv_h) * normalization_constant(inv_h);
  else
    return 0.0;
}

FOUR_C_NAMESPACE_CLOSE
