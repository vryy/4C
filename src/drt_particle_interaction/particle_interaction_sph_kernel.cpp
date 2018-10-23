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

#include "../drt_inpar/inpar_particle.H"

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
  const double dWdrij = this->dWdrij(rij, support);

  for (int i = 0; i < 3; ++i) gradWij[i] = eij[i] * dWdrij;
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
      return (10.0 * M_1_PI / (7.0 * square(h)));
    }
    case INPAR::PARTICLE::Kernel3D:
    {
      return M_1_PI / cube(h);
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
 | evaluate kernel                                            sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::SPHKernelCubicSpline::W(const double& rij, const double& support) const
{
  const double h = 0.5 * support;
  const double q = rij / h;

  if (q < 1.0)
    return (1.0 - 1.5 * square(q) + 0.75 * cube(q)) * NormalizationConstant(h);
  else if (q < 2.0)
    return (0.25 * cube(2.0 - q)) * NormalizationConstant(h);
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
    return (-3.0 * q + 2.25 * square(q)) * (1.0 / h) * NormalizationConstant(h);
  else if (q < 2.0)
    return (-0.75 * square(2.0 - q)) * (1.0 / h) * NormalizationConstant(h);
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
    return (-3.0 + 4.5 * q) * (1.0 / square(h)) * NormalizationConstant(h);
  else if (q < 2.0)
    return (1.5 * (2.0 - q)) * (1.0 / square(h)) * NormalizationConstant(h);
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
      return (7.0 * M_1_PI / (478.0 * square(h)));
    }
    case INPAR::PARTICLE::Kernel3D:
    {
      return (3.0 * M_1_PI / (359.0 * cube(h)));
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
 | evaluate kernel                                            sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::SPHKernelQuinticSpline::W(
    const double& rij, const double& support) const
{
  const double h = support / 3.0;
  const double q = rij / h;

  if (q < 1.0)
    return (powint((3.0 - q), 5) - 6.0 * powint((2.0 - q), 5) + 15.0 * powint((1.0 - q), 5)) *
           NormalizationConstant(h);
  else if (q < 2.0)
    return (powint((3.0 - q), 5) - 6.0 * powint((2.0 - q), 5)) * NormalizationConstant(h);
  else if (q < 3.0)
    return powint((3.0 - q), 5) * NormalizationConstant(h);
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
    return (-5.0 * powint((3.0 - q), 4) + 30.0 * powint((2.0 - q), 4) -
               75.0 * powint((1.0 - q), 4)) *
           (1.0 / h) * NormalizationConstant(h);
  else if (q < 2.0)
    return (-5.0 * powint((3.0 - q), 4) + 30.0 * powint((2.0 - q), 4)) * (1.0 / h) *
           NormalizationConstant(h);
  else if (q < 3.0)
    return (-5.0 * powint((3.0 - q), 4)) * (1.0 / h) * NormalizationConstant(h);
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
    return (20.0 * powint((3.0 - q), 3) - 120.0 * powint((2.0 - q), 3) +
               300.0 * powint((1.0 - q), 3)) *
           (1.0 / square(h)) * NormalizationConstant(h);
  else if (q < 2.0)
    return (20.0 * powint((3.0 - q), 3) - 120.0 * powint((2.0 - q), 3)) * (1.0 / square(h)) *
           NormalizationConstant(h);
  else if (q < 3.0)
    return (20.0 * powint((3.0 - q), 3)) * (1.0 / square(h)) * NormalizationConstant(h);
  else
    return 0.0;
}
