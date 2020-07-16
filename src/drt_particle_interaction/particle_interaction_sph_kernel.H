/*---------------------------------------------------------------------------*/
/*! \file
\brief kernel handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_KERNEL_H
#define PARTICLE_INTERACTION_SPH_KERNEL_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include <memory>

#include "../drt_inpar/inpar_particle.H"

#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHKernelBase
  {
   public:
    //! constructor
    explicit SPHKernelBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~SPHKernelBase() = default;

    //! init kernel handler
    virtual void Init();

    //! setup kernel handler
    virtual void Setup();

    //! write restart of kernel handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of kernel handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! get kernel space dimension
    virtual void KernelSpaceDimension(int& dim) const final;

    //! get smoothing length from kernel support radius
    virtual double SmoothingLength(const double& support) const = 0;

    //! get normalization constant from inverse smoothing length
    virtual double NormalizationConstant(const double& inv_h) const = 0;

    //! evaluate kernel (self-interaction)
    virtual double W0(const double& support) const = 0;

    //! evaluate kernel
    virtual double W(const double& rij, const double& support) const = 0;

    //! evaluate first derivative of kernel
    virtual double dWdrij(const double& rij, const double& support) const = 0;

    //! evaluate second derivative of kernel
    virtual double d2Wdrij2(const double& rij, const double& support) const = 0;

    //! evaluate gradient of kernel
    virtual void GradWij(
        const double& rij, const double& support, const double* eij, double* gradWij) const final;

   protected:
    // store problem dimension required by weight functions
    INPAR::PARTICLE::KernelSpaceDimension kernelspacedim_;
  };

  class SPHKernelCubicSpline final : public SPHKernelBase
  {
   public:
    //! constructor
    explicit SPHKernelCubicSpline(const Teuchos::ParameterList& params);

    //! get smoothing length from kernel support radius
    double SmoothingLength(const double& support) const override;

    //! get normalization constant from inverse smoothing length
    double NormalizationConstant(const double& inv_h) const override;

    //! evaluate kernel (self-interaction)
    double W0(const double& support) const override;

    //! evaluate kernel
    double W(const double& rij, const double& support) const override;

    //! evaluate first derivative of kernel
    double dWdrij(const double& rij, const double& support) const override;

    //! evaluate second derivative of kernel
    double d2Wdrij2(const double& rij, const double& support) const override;
  };

  class SPHKernelQuinticSpline final : public SPHKernelBase
  {
   public:
    //! constructor
    explicit SPHKernelQuinticSpline(const Teuchos::ParameterList& params);

    //! get smoothing length from kernel support radius
    double SmoothingLength(const double& support) const override;

    //! get normalization constant from inverse smoothing length
    double NormalizationConstant(const double& inv_h) const override;

    //! evaluate kernel (self-interaction)
    double W0(const double& support) const override;

    //! evaluate kernel
    double W(const double& rij, const double& support) const override;

    //! evaluate first derivative of kernel
    double dWdrij(const double& rij, const double& support) const override;

    //! evaluate second derivative of kernel
    double d2Wdrij2(const double& rij, const double& support) const override;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
