/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_kernel.H

\brief kernel handler for smoothed particle hydrodynamics (SPH) interactions

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
#ifndef PARTICLE_INTERACTION_SPH_KERNEL_H
#define PARTICLE_INTERACTION_SPH_KERNEL_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_inpar/inpar_particle.H"

#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  //! square of a double
  static inline double square(const double& x) { return x * x; };

  //! cube of a double
  static inline double cube(const double& x) { return x * x * x; };

  //! optimized version of std::pow with integer as exponent
  static inline double powint(const double& x, const int n)
  {
    double yy, ww;

    if (x == 0.0) return 0.0;
    int nn = (n > 0) ? n : -n;
    ww = x;

    for (yy = 1.0; nn != 0; nn >>= 1, ww *= ww)
      if (nn & 1) yy *= ww;

    return (n > 0) ? yy : 1.0 / yy;
  }

  class SPHKernelBase
  {
   public:
    //! constructor
    explicit SPHKernelBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~SPHKernelBase(){};

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

    //! get normalization constant from smoothing length
    virtual double NormalizationConstant(const double& h) const = 0;

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

  class SPHKernelCubicSpline : public SPHKernelBase
  {
   public:
    //! constructor
    explicit SPHKernelCubicSpline(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~SPHKernelCubicSpline(){};

    //! get smoothing length from kernel support radius
    virtual double SmoothingLength(const double& support) const override;

    //! get normalization constant from smoothing length
    virtual double NormalizationConstant(const double& h) const override;

    //! evaluate kernel
    virtual double W(const double& rij, const double& support) const override;

    //! evaluate first derivative of kernel
    virtual double dWdrij(const double& rij, const double& support) const override;

    //! evaluate second derivative of kernel
    virtual double d2Wdrij2(const double& rij, const double& support) const override;
  };

  class SPHKernelQuinticSpline : public SPHKernelBase
  {
   public:
    //! constructor
    explicit SPHKernelQuinticSpline(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~SPHKernelQuinticSpline(){};

    //! get smoothing length from kernel support radius
    virtual double SmoothingLength(const double& support) const override;

    //! get normalization constant from smoothing length
    virtual double NormalizationConstant(const double& h) const override;

    //! evaluate kernel
    virtual double W(const double& rij, const double& support) const override;

    //! evaluate first derivative of kernel
    virtual double dWdrij(const double& rij, const double& support) const override;

    //! evaluate second derivative of kernel
    virtual double d2Wdrij2(const double& rij, const double& support) const override;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
