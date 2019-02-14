/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_kernel.cpp

\brief kernel handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_kernel.H"
#include "particle_interaction_utils.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHKernelBase::SPHKernelBase(const Teuchos::ParameterList& params)
    : kernelspacedim_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::KernelSpaceDimension>(
          params.sublist("SPH"), "KERNEL_SPACE_DIM"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init kernel handler                                        sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHKernelBase::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup kernel handler                                       sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHKernelBase::Setup()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | write restart of kernel handler                            sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHKernelBase::WriteRestart(const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of kernel handler                             sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHKernelBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | get kernel space dimension                                 sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
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

/*---------------------------------------------------------------------------*
 | evaluate gradient of kernel                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHKernelBase::GradWij(
    const double& rij, const double& support, const double* eij, double* gradWij) const
{
  UTILS::vec_setscale(gradWij, this->dWdrij(rij, support), eij);
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHKernelCubicSpline::SPHKernelCubicSpline(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::SPHKernelBase(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | get smoothing length from kernel support radius            sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::SPHKernelCubicSpline::SmoothingLength(const double& support) const
{
  return (0.5 * support);
}

/*---------------------------------------------------------------------------*
 | get normalization constant from smoothing length           sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::SPHKernelCubicSpline::NormalizationConstant(const double& h) const
{
  switch (kernelspacedim_)
  {
    case INPAR::PARTICLE::Kernel1D:
    {
      return (2.0 / (3.0 * h));
    }
    case INPAR::PARTICLE::Kernel2D:
    {
      return (10.0 * M_1_PI / (7.0 * UTILS::pow<2>(h)));
    }
    case INPAR::PARTICLE::Kernel3D:
    {
      return M_1_PI / UTILS::pow<3>(h);
    }
    default:
    {
      dserror("unknown kernel space dimension!");
      break;
    }
  }

  return 0.0;
}

/*---------------------------------------------------------------------------*
 | evaluate kernel (self-interaction)                         sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::SPHKernelCubicSpline::W0(const double& support) const
{
  return NormalizationConstant(0.5 * support);
}

/*---------------------------------------------------------------------------*
 | evaluate kernel                                            sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::SPHKernelCubicSpline::W(const double& rij, const double& support) const
{
  const double h = 0.5 * support;
  const double q = rij / h;

  if (q < 1.0)
    return (1.0 - 1.5 * UTILS::pow<2>(q) + 0.75 * UTILS::pow<3>(q)) * NormalizationConstant(h);
  else if (q < 2.0)
    return (0.25 * UTILS::pow<3>(2.0 - q)) * NormalizationConstant(h);
  else
    return 0.0;
}

/*---------------------------------------------------------------------------*
 | evaluate first derivative of kernel                        sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::SPHKernelCubicSpline::dWdrij(
    const double& rij, const double& support) const
{
  const double h = 0.5 * support;
  const double q = rij / h;

  if (q < 1.0)
    return (-3.0 * q + 2.25 * UTILS::pow<2>(q)) * (1.0 / h) * NormalizationConstant(h);
  else if (q < 2.0)
    return (-0.75 * UTILS::pow<2>(2.0 - q)) * (1.0 / h) * NormalizationConstant(h);
  else
    return 0.0;
}

/*---------------------------------------------------------------------------*
 | evaluate second derivative of kernel                       sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::SPHKernelCubicSpline::d2Wdrij2(
    const double& rij, const double& support) const
{
  const double h = 0.5 * support;
  const double q = rij / h;

  if (q < 1.0)
    return (-3.0 + 4.5 * q) * (1.0 / UTILS::pow<2>(h)) * NormalizationConstant(h);
  else if (q < 2.0)
    return (1.5 * (2.0 - q)) * (1.0 / UTILS::pow<2>(h)) * NormalizationConstant(h);
  else
    return 0.0;
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHKernelQuinticSpline::SPHKernelQuinticSpline(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::SPHKernelBase(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | get smoothing length from kernel support radius            sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::SPHKernelQuinticSpline::SmoothingLength(const double& support) const
{
  return (support / 3.0);
}

/*---------------------------------------------------------------------------*
 | get normalization constant from smoothing length           sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::SPHKernelQuinticSpline::NormalizationConstant(const double& h) const
{
  switch (kernelspacedim_)
  {
    case INPAR::PARTICLE::Kernel1D:
    {
      return 1.0 / (120.0 * h);
    }
    case INPAR::PARTICLE::Kernel2D:
    {
      return (7.0 * M_1_PI / (478.0 * UTILS::pow<2>(h)));
    }
    case INPAR::PARTICLE::Kernel3D:
    {
      return (3.0 * M_1_PI / (359.0 * UTILS::pow<3>(h)));
    }
    default:
    {
      dserror("unknown kernel space dimension!");
      break;
    }
  }

  return 0.0;
}

/*---------------------------------------------------------------------------*
 | evaluate kernel (self-interaction)                         sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::SPHKernelQuinticSpline::W0(const double& support) const
{
  return 66.0 * NormalizationConstant(support / 3.0);
}

/*---------------------------------------------------------------------------*
 | evaluate kernel                                            sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::SPHKernelQuinticSpline::W(
    const double& rij, const double& support) const
{
  const double h = support / 3.0;
  const double q = rij / h;

  if (q < 1.0)
    return (UTILS::pow<5>(3.0 - q) - 6.0 * UTILS::pow<5>(2.0 - q) + 15.0 * UTILS::pow<5>(1.0 - q)) *
           NormalizationConstant(h);
  else if (q < 2.0)
    return (UTILS::pow<5>(3.0 - q) - 6.0 * UTILS::pow<5>(2.0 - q)) * NormalizationConstant(h);
  else if (q < 3.0)
    return UTILS::pow<5>(3.0 - q) * NormalizationConstant(h);
  else
    return 0.0;
}

/*---------------------------------------------------------------------------*
 | evaluate first derivative of kernel                        sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::SPHKernelQuinticSpline::dWdrij(
    const double& rij, const double& support) const
{
  const double h = support / 3.0;
  const double q = rij / h;

  if (q < 1.0)
    return (-5.0 * UTILS::pow<4>(3.0 - q) + 30.0 * UTILS::pow<4>(2.0 - q) -
               75.0 * UTILS::pow<4>(1.0 - q)) *
           (1.0 / h) * NormalizationConstant(h);
  else if (q < 2.0)
    return (-5.0 * UTILS::pow<4>(3.0 - q) + 30.0 * UTILS::pow<4>(2.0 - q)) * (1.0 / h) *
           NormalizationConstant(h);
  else if (q < 3.0)
    return (-5.0 * UTILS::pow<4>(3.0 - q)) * (1.0 / h) * NormalizationConstant(h);
  else
    return 0.0;
}

/*---------------------------------------------------------------------------*
 | evaluate second derivative of kernel                       sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::SPHKernelQuinticSpline::d2Wdrij2(
    const double& rij, const double& support) const
{
  const double h = support / 3.0;
  const double q = rij / h;

  if (q < 1.0)
    return (20.0 * UTILS::pow<3>(3.0 - q) - 120.0 * UTILS::pow<3>(2.0 - q) +
               300.0 * UTILS::pow<3>(1.0 - q)) *
           (1.0 / UTILS::pow<2>(h)) * NormalizationConstant(h);
  else if (q < 2.0)
    return (20.0 * UTILS::pow<3>(3.0 - q) - 120.0 * UTILS::pow<3>(2.0 - q)) *
           (1.0 / UTILS::pow<2>(h)) * NormalizationConstant(h);
  else if (q < 3.0)
    return (20.0 * UTILS::pow<3>(3.0 - q)) * (1.0 / UTILS::pow<2>(h)) * NormalizationConstant(h);
  else
    return 0.0;
}
