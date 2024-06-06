/*---------------------------------------------------------------------------*/
/*! \file
\brief density correction handler in smoothed particle hydrodynamics (SPH)
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_INTERACTION_SPH_DENSITY_CORRECTION_HPP
#define FOUR_C_PARTICLE_INTERACTION_SPH_DENSITY_CORRECTION_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace ParticleInteraction
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
    virtual void corrected_density_interior(const double* denssum, double* dens) const;

    //! set corrected density of free surface particles
    virtual void corrected_density_free_surface(const double* denssum, const double* colorfield,
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
    void corrected_density_free_surface(const double* denssum, const double* colorfield,
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
    void corrected_density_free_surface(const double* denssum, const double* colorfield,
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
    void corrected_density_free_surface(const double* denssum, const double* colorfield,
        const double* dens_bc, double* dens) const override;
  };

}  // namespace ParticleInteraction

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
