/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_artificialviscosity.H

\brief artificial viscosity handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_ARTIFICIALVISCOSITY_H
#define PARTICLE_INTERACTION_SPH_ARTIFICIALVISCOSITY_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include <memory>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHArtificialViscosity
  {
   public:
    //! constructor
    explicit SPHArtificialViscosity();

    //! virtual destructor
    virtual ~SPHArtificialViscosity(){};

    //! init artificial viscosity handler
    virtual void Init();

    //! setup artificial viscosity handler
    virtual void Setup();

    //! write restart of artificial viscosity handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of artificial viscosity handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! evaluate artificial viscosity
    virtual void ArtificialViscosity(const double* dens_i, const double* dens_j,
        const double* vel_i, const double* vel_j, const double* mass_j,
        const double& artificialviscosity, const double& dWdrij, const double& h_i,
        const double& h_j, const double& c_i, const double& c_j, const double& abs_rij,
        const double* e_ij, double* acc_i) const;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
