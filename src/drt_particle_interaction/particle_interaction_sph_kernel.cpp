/*---------------------------------------------------------------------------*/
/*! \file
\brief kernel handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_kernel.H"
#include "particle_interaction_utils.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHKernelBase::SPHKernelBase(const Teuchos::ParameterList& params)
    : kernelspacedim_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::KernelSpaceDimension>(
          params, "KERNEL_SPACE_DIM"))
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHKernelBase::Init()
{
  // nothing to do
}

void PARTICLEINTERACTION::SPHKernelBase::Setup()
{
  // nothing to do
}

void PARTICLEINTERACTION::SPHKernelBase::WriteRestart(const int step, const double time) const
{
  // nothing to do
}

void PARTICLEINTERACTION::SPHKernelBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

void PARTICLEINTERACTION::SPHKernelBase::KernelSpaceDimension(int& dim) const
{
  switch (kernelspacedim_)
  {
    case INPAR::PARTICLE::Kernel1D:
    {
      dim = 1;
      break;
    }
    case INPAR::PARTICLE::Kernel2D:
    {
      dim = 2;
      break;
    }
    case INPAR::PARTICLE::Kernel3D:
    {
      dim = 3;
      break;
    }
    default:
    {
      dserror("unknown kernel space dimension!");
      break;
    }
  }
}

void PARTICLEINTERACTION::SPHKernelBase::GradWij(
    const double& rij, const double& support, const double* eij, double* gradWij) const
{
  UTILS::vec_setscale(gradWij, this->dWdrij(rij, support), eij);
}

PARTICLEINTERACTION::SPHKernelCubicSpline::SPHKernelCubicSpline(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::SPHKernelBase(params)
{
  // empty constructor
}

double PARTICLEINTERACTION::SPHKernelCubicSpline::SmoothingLength(const double& support) const
{
  return (0.5 * support);
}

double PARTICLEINTERACTION::SPHKernelCubicSpline::NormalizationConstant(const double& inv_h) const
{
  switch (kernelspacedim_)
  {
    case INPAR::PARTICLE::Kernel1D:
    {
      // (2.0 / 3.0) * inv_h
      return 0.6666666666666666 * inv_h;
    }
    case INPAR::PARTICLE::Kernel2D:
    {
      // (10.0 / 7.0) * M_1_PI * inv_h * inv_h
      return 0.4547284088339866 * UTILS::pow<2>(inv_h);
    }
    case INPAR::PARTICLE::Kernel3D:
    {
      return M_1_PI * UTILS::pow<3>(inv_h);
    }
    default:
    {
      dserror("unknown kernel space dimension!");
      break;
    }
  }

  return 0.0;
}

double PARTICLEINTERACTION::SPHKernelCubicSpline::W0(const double& support) const
{
  return NormalizationConstant(2.0 / support);
}

double PARTICLEINTERACTION::SPHKernelCubicSpline::W(const double& rij, const double& support) const
{
  const double inv_h = 2.0 / support;
  const double q = rij * inv_h;

  if (q < 1.0)
    return (1.0 - 1.5 * UTILS::pow<2>(q) + 0.75 * UTILS::pow<3>(q)) * NormalizationConstant(inv_h);
  else if (q < 2.0)
    return (0.25 * UTILS::pow<3>(2.0 - q)) * NormalizationConstant(inv_h);
  else
    return 0.0;
}

double PARTICLEINTERACTION::SPHKernelCubicSpline::dWdrij(
    const double& rij, const double& support) const
{
  const double inv_h = 2.0 / support;
  const double q = rij * inv_h;

  if (q < 1.0)
    return (-3.0 * q + 2.25 * UTILS::pow<2>(q)) * inv_h * NormalizationConstant(inv_h);
  else if (q < 2.0)
    return (-0.75 * UTILS::pow<2>(2.0 - q)) * inv_h * NormalizationConstant(inv_h);
  else
    return 0.0;
}

double PARTICLEINTERACTION::SPHKernelCubicSpline::d2Wdrij2(
    const double& rij, const double& support) const
{
  const double inv_h = 2.0 / support;
  const double q = rij * inv_h;

  if (q < 1.0)
    return (-3.0 + 4.5 * q) * UTILS::pow<2>(inv_h) * NormalizationConstant(inv_h);
  else if (q < 2.0)
    return (1.5 * (2.0 - q)) * UTILS::pow<2>(inv_h) * NormalizationConstant(inv_h);
  else
    return 0.0;
}

PARTICLEINTERACTION::SPHKernelQuinticSpline::SPHKernelQuinticSpline(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::SPHKernelBase(params)
{
  // empty constructor
}

double PARTICLEINTERACTION::SPHKernelQuinticSpline::SmoothingLength(const double& support) const
{
  // (support / 3.0)
  return 0.3333333333333333 * support;
}

double PARTICLEINTERACTION::SPHKernelQuinticSpline::NormalizationConstant(const double& inv_h) const
{
  switch (kernelspacedim_)
  {
    case INPAR::PARTICLE::Kernel1D:
    {
      // (inv_h / 120.0)
      return 0.0083333333333333 * inv_h;
    }
    case INPAR::PARTICLE::Kernel2D:
    {
      // (7.0 / 478.0) * M_1_PI * inv_h * inv_h
      return 0.0046614418478797 * UTILS::pow<2>(inv_h);
    }
    case INPAR::PARTICLE::Kernel3D:
    {
      // (3.0 / 359.0) * M_1_PI * inv_h * inv_h * inv_h
      return 0.0026599711937364 * UTILS::pow<3>(inv_h);
    }
    default:
    {
      dserror("unknown kernel space dimension!");
      break;
    }
  }

  return 0.0;
}

double PARTICLEINTERACTION::SPHKernelQuinticSpline::W0(const double& support) const
{
  return 66.0 * NormalizationConstant(3.0 / support);
}

double PARTICLEINTERACTION::SPHKernelQuinticSpline::W(
    const double& rij, const double& support) const
{
  const double inv_h = 3.0 / support;
  const double q = rij * inv_h;

  if (q < 1.0)
    return (UTILS::pow<5>(3.0 - q) - 6.0 * UTILS::pow<5>(2.0 - q) + 15.0 * UTILS::pow<5>(1.0 - q)) *
           NormalizationConstant(inv_h);
  else if (q < 2.0)
    return (UTILS::pow<5>(3.0 - q) - 6.0 * UTILS::pow<5>(2.0 - q)) * NormalizationConstant(inv_h);
  else if (q < 3.0)
    return UTILS::pow<5>(3.0 - q) * NormalizationConstant(inv_h);
  else
    return 0.0;
}

double PARTICLEINTERACTION::SPHKernelQuinticSpline::dWdrij(
    const double& rij, const double& support) const
{
  const double inv_h = 3.0 / support;
  const double q = rij * inv_h;

  if (q < 1.0)
    return (-5.0 * UTILS::pow<4>(3.0 - q) + 30.0 * UTILS::pow<4>(2.0 - q) -
               75.0 * UTILS::pow<4>(1.0 - q)) *
           inv_h * NormalizationConstant(inv_h);
  else if (q < 2.0)
    return (-5.0 * UTILS::pow<4>(3.0 - q) + 30.0 * UTILS::pow<4>(2.0 - q)) * inv_h *
           NormalizationConstant(inv_h);
  else if (q < 3.0)
    return (-5.0 * UTILS::pow<4>(3.0 - q)) * inv_h * NormalizationConstant(inv_h);
  else
    return 0.0;
}

double PARTICLEINTERACTION::SPHKernelQuinticSpline::d2Wdrij2(
    const double& rij, const double& support) const
{
  const double inv_h = 3.0 / support;
  const double q = rij * inv_h;

  if (q < 1.0)
    return (20.0 * UTILS::pow<3>(3.0 - q) - 120.0 * UTILS::pow<3>(2.0 - q) +
               300.0 * UTILS::pow<3>(1.0 - q)) *
           UTILS::pow<2>(inv_h) * NormalizationConstant(inv_h);
  else if (q < 2.0)
    return (20.0 * UTILS::pow<3>(3.0 - q) - 120.0 * UTILS::pow<3>(2.0 - q)) * UTILS::pow<2>(inv_h) *
           NormalizationConstant(inv_h);
  else if (q < 3.0)
    return (20.0 * UTILS::pow<3>(3.0 - q)) * UTILS::pow<2>(inv_h) * NormalizationConstant(inv_h);
  else
    return 0.0;
}
