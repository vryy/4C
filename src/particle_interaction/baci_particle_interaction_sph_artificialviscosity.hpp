/*---------------------------------------------------------------------------*/
/*! \file
\brief artificial viscosity handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef BACI_PARTICLE_INTERACTION_SPH_ARTIFICIALVISCOSITY_HPP
#define BACI_PARTICLE_INTERACTION_SPH_ARTIFICIALVISCOSITY_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_config.hpp"

#include <memory>

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHArtificialViscosity final
  {
   public:
    //! constructor
    explicit SPHArtificialViscosity();

    //! init artificial viscosity handler
    void Init();

    //! setup artificial viscosity handler
    void Setup();

    //! evaluate artificial viscosity
    void ArtificialViscosity(const double* vel_i, const double* vel_j, const double* mass_i,
        const double* mass_j, const double& artvisc_i, const double& artvisc_j,
        const double& dWdrij, const double& dWdrji, const double& dens_ij, const double& h_ij,
        const double& c_ij, const double& abs_rij, const double* e_ij, double* acc_i,
        double* acc_j) const;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
