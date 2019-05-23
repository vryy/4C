/*---------------------------------------------------------------------------*/
/*!

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
    virtual void ArtificialViscosity(const double* vel_i, const double* vel_j, const double* mass_i,
        const double* mass_j, const double& artvisc_i, const double& artvisc_j,
        const double& dWdrij, const double& dWdrji, const double& dens_ij, const double& h_ij,
        const double& c_ij, const double& abs_rij, const double* e_ij, double* acc_i,
        double* acc_j) const;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
