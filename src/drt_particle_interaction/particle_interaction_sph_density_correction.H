/*---------------------------------------------------------------------------*/
/*! \file
\brief density correction handler in smoothed particle hydrodynamics (SPH)
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_DENSITY_CORRECTION_H
#define PARTICLE_INTERACTION_SPH_DENSITY_CORRECTION_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include <memory>

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
  class SPHDensityCorrectionBase
  {
   public:
    //! constructor
    explicit SPHDensityCorrectionBase();

    //! virtual destructor
    virtual ~SPHDensityCorrectionBase() = default;

    //! init density correction handler
    virtual void Init();

    //! setup density correction handler
    virtual void Setup();

    //! density boundary condition is needed
    virtual bool ComputeDensityBC() const = 0;

    //! set corrected density of interior particles
    virtual void CorrectedDensityInterior(const double* denssum, double* dens) const;

    //! set corrected density of free surface particles
    virtual void CorrectedDensityFreeSurface(const double* denssum, const double* colorfield,
        const double* dens_bc, double* dens) const = 0;
  };

  class SPHDensityCorrectionInterior : public SPHDensityCorrectionBase
  {
   public:
    //! constructor
    explicit SPHDensityCorrectionInterior();

    //! density boundary condition is needed
    bool ComputeDensityBC() const override;

    //! set corrected density of free surface particles
    void CorrectedDensityFreeSurface(const double* denssum, const double* colorfield,
        const double* dens_bc, double* dens) const override;
  };

  class SPHDensityCorrectionNormalized : public SPHDensityCorrectionBase
  {
   public:
    //! constructor
    explicit SPHDensityCorrectionNormalized();

    //! density boundary condition is needed
    bool ComputeDensityBC() const override;

    //! set corrected density of free surface particles
    void CorrectedDensityFreeSurface(const double* denssum, const double* colorfield,
        const double* dens_bc, double* dens) const override;
  };

  class SPHDensityCorrectionRandles : public SPHDensityCorrectionBase
  {
   public:
    //! constructor
    explicit SPHDensityCorrectionRandles();

    //! density boundary condition is needed
    bool ComputeDensityBC() const override;

    //! set corrected density of free surface particles
    void CorrectedDensityFreeSurface(const double* denssum, const double* colorfield,
        const double* dens_bc, double* dens) const override;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
