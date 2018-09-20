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
  return (support / 2.0);
}

/*---------------------------------------------------------------------------*
 | get normalization constant from smoothing length           sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::SPHKernelCubicSpline::NormalizationConstant(const double& h) const
{
  double normalizationconstant = 0.0;

  switch (kernelspacedim_)
  {
    case INPAR::PARTICLE::Kernel1D:
    {
      normalizationconstant = (2.0 / (3.0 * h));
      break;
    }
    case INPAR::PARTICLE::Kernel2D:
    {
      normalizationconstant = (10.0 * M_1_PI / (7.0 * std::pow(h, 2)));
      break;
    }
    case INPAR::PARTICLE::Kernel3D:
    {
      normalizationconstant = M_1_PI / std::pow(h, 3);
      break;
    }
    default:
    {
      dserror("unknown kernel space dimension!");
      break;
    }
  }

  return normalizationconstant;
}

/*---------------------------------------------------------------------------*
 | evaluate kernel                                            sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::SPHKernelCubicSpline::W(const double& rij, const double& support) const
{
  const double h = support / 2.0;
  const double q = rij / h;

  if (q < 1.0)
    return (1.0 - 1.5 * std::pow(q, 2) + 0.75 * std::pow(q, 3)) * NormalizationConstant(h);
  else if (q < 2.0)
    return (std::pow((2.0 - q), 3) / 4.0) * NormalizationConstant(h);
  else
    return 0.0;
}

/*---------------------------------------------------------------------------*
 | evaluate first derivative of kernel                        sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::SPHKernelCubicSpline::dWdrij(
    const double& rij, const double& support) const
{
  const double h = support / 2.0;
  const double q = rij / h;

  if (q < 1.0)
    return (-3.0 * q + (9.0 / 4.0) * std::pow(q, 2)) * (1.0 / h) * NormalizationConstant(h);
  else if (q < 2.0)
    return (-(3.0 / 4.0) * std::pow((2.0 - q), 2)) * (1.0 / h) * NormalizationConstant(h);
  else
    return 0.0;
}

/*---------------------------------------------------------------------------*
 | evaluate second derivative of kernel                       sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::SPHKernelCubicSpline::d2Wdrij2(
    const double& rij, const double& support) const
{
  const double h = support / 2.0;
  const double q = rij / h;

  if (q < 1.0)
    return (-3.0 + (9.0 / 2.0) * q) * (1.0 / std::pow(h, 2)) * NormalizationConstant(h);
  else if (q < 2.0)
    return ((3.0 / 2.0) * (2.0 - q)) * (1.0 / std::pow(h, 2)) * NormalizationConstant(h);
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
  double normalizationconstant = 0.0;

  switch (kernelspacedim_)
  {
    case INPAR::PARTICLE::Kernel1D:
    {
      normalizationconstant = 1.0 / (120.0 * h);
      break;
    }
    case INPAR::PARTICLE::Kernel2D:
    {
      normalizationconstant = (7.0 * M_1_PI / (478.0 * std::pow(h, 2)));
      break;
    }
    case INPAR::PARTICLE::Kernel3D:
    {
      normalizationconstant = (3.0 * M_1_PI / (359.0 * std::pow(h, 3)));
      break;
    }
    default:
    {
      dserror("unknown kernel space dimension!");
      break;
    }
  }

  return normalizationconstant;
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
    return (std::pow((3.0 - q), 5) - 6.0 * std::pow((2.0 - q), 5) + 15.0 * std::pow((1.0 - q), 5)) *
           NormalizationConstant(h);
  else if (q < 2.0)
    return (std::pow((3.0 - q), 5) - 6.0 * std::pow((2.0 - q), 5)) * NormalizationConstant(h);
  else if (q < 3.0)
    return (std::pow((3.0 - q), 5)) * NormalizationConstant(h);
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
    return (-5.0 * std::pow((3.0 - q), 4) + 30.0 * std::pow((2.0 - q), 4) -
               75.0 * std::pow((1.0 - q), 4)) *
           (1.0 / h) * NormalizationConstant(h);
  else if (q < 2.0)
    return (-5.0 * std::pow((3.0 - q), 4) + 30.0 * std::pow((2.0 - q), 4)) * (1.0 / h) *
           NormalizationConstant(h);
  else if (q < 3.0)
    return (-5.0 * std::pow((3.0 - q), 4)) * (1.0 / h) * NormalizationConstant(h);
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
    return (20.0 * std::pow((3.0 - q), 3) - 120.0 * std::pow((2.0 - q), 3) +
               300.0 * std::pow((1.0 - q), 3)) *
           (1.0 / std::pow(h, 2)) * NormalizationConstant(h);
  else if (q < 2.0)
    return (20.0 * std::pow((3.0 - q), 3) - 120.0 * std::pow((2.0 - q), 3)) *
           (1.0 / std::pow(h, 2)) * NormalizationConstant(h);
  else if (q < 3.0)
    return (20.0 * std::pow((3.0 - q), 3)) * (1.0 / std::pow(h, 2)) * NormalizationConstant(h);
  else
    return 0.0;
}
