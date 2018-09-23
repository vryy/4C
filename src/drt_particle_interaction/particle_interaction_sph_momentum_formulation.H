/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_momentum_formulation.H

\brief momentum formulation handler for smoothed particle hydrodynamics (SPH) interactions

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
#ifndef PARTICLE_INTERACTION_SPH_MOMENTUM_FORMULATION_H
#define PARTICLE_INTERACTION_SPH_MOMENTUM_FORMULATION_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2018 |
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
  class SPHMomentumFormulationBase
  {
   public:
    //! constructor
    explicit SPHMomentumFormulationBase();

    //! virtual destructor
    virtual ~SPHMomentumFormulationBase(){};

    //! init momentum formulation handler
    virtual void Init();

    //! setup momentum formulation handler
    virtual void Setup();

    //! write restart of momentum formulation handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of momentum formulation handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! evaluate specific coefficient
    virtual void SpecificCoefficient(const double* dens_i, const double* dens_j,
        const double* mass_i, const double* mass_j, const double& dWdrij,
        double& specificcoefficient_ij) const = 0;

    //! evaluate pressure gradient
    virtual void PressureGradient(const double* dens_i, const double* dens_j, const double* press_i,
        const double* press_j, const double& specificcoefficient_ij, const double* e_ij,
        double* acc_i) const = 0;

    //! evaluate shear forces
    virtual void ShearForces(const double* dens_i, const double* dens_j, const double* vel_i,
        const double* vel_j, const double& kernelfac, const double& visc_i, const double& visc_j,
        const double& bulk_visc_i, const double& bulk_visc_j, const double& abs_rij,
        const double& specificcoefficient_ij, const double* e_ij, double* acc_i) const = 0;

    //! evaluate background pressure (standard formulation)
    virtual void StandardBackgroundPressure(const double* dens_i, const double* dens_j,
        const double& backgroundpressure, const double& specificcoefficient_ij, const double* e_ij,
        double* acc_i) const = 0;

    //! evaluate background pressure  (generalized formulation)
    virtual void GeneralizedBackgroundPressure(const double* dens_i, const double* mass_i,
        const double* mass_j, const double& backgroundpressure_tilde, const double& dWdrij_tilde,
        const double* e_ij, double* acc_i) const = 0;

    //! evaluate modified velocity contribution
    virtual void ModifiedVelocityContribution(const double* dens_i, const double* dens_j,
        const double* vel_i, const double* vel_j, const double* mod_vel_i, const double* mod_vel_j,
        const double& specificcoefficient_ij, const double* e_ij, double* acc_i) const = 0;
  };

  class SPHMomentumFormulationMonaghan : public SPHMomentumFormulationBase
  {
   public:
    //! constructor
    explicit SPHMomentumFormulationMonaghan();

    //! virtual destructor
    virtual ~SPHMomentumFormulationMonaghan(){};

    //! evaluate specific coefficient
    virtual void SpecificCoefficient(const double* dens_i, const double* dens_j,
        const double* mass_i, const double* mass_j, const double& dWdrij,
        double& specificcoefficient_ij) const override;

    //! evaluate pressure gradient
    virtual void PressureGradient(const double* dens_i, const double* dens_j, const double* press_i,
        const double* press_j, const double& specificcoefficient_ij, const double* e_ij,
        double* acc_i) const override;

    //! evaluate shear forces
    virtual void ShearForces(const double* dens_i, const double* dens_j, const double* vel_i,
        const double* vel_j, const double& kernelfac, const double& visc_i, const double& visc_j,
        const double& bulk_visc_i, const double& bulk_visc_j, const double& abs_rij,
        const double& specificcoefficient_ij, const double* e_ij, double* acc_i) const override;

    //! evaluate background pressure (standard formulation)
    virtual void StandardBackgroundPressure(const double* dens_i, const double* dens_j,
        const double& backgroundpressure, const double& specificcoefficient_ij, const double* e_ij,
        double* acc_i) const override;

    //! evaluate background pressure (generalized formulation)
    virtual void GeneralizedBackgroundPressure(const double* dens_i, const double* mass_i,
        const double* mass_j, const double& backgroundpressure_tilde, const double& dWdrij_tilde,
        const double* e_ij, double* acc_i) const override;

    //! evaluate modified velocity contribution
    virtual void ModifiedVelocityContribution(const double* dens_i, const double* dens_j,
        const double* vel_i, const double* vel_j, const double* mod_vel_i, const double* mod_vel_j,
        const double& specificcoefficient_ij, const double* e_ij, double* acc_i) const override;
  };

  class SPHMomentumFormulationAdami : public SPHMomentumFormulationBase
  {
   public:
    //! constructor
    explicit SPHMomentumFormulationAdami();

    //! virtual destructor
    virtual ~SPHMomentumFormulationAdami(){};

    //! evaluate specific coefficient
    virtual void SpecificCoefficient(const double* dens_i, const double* dens_j,
        const double* mass_i, const double* mass_j, const double& dWdrij,
        double& specificcoefficient_ij) const override;

    //! evaluate pressure gradient
    virtual void PressureGradient(const double* dens_i, const double* dens_j, const double* press_i,
        const double* press_j, const double& specificcoefficient_ij, const double* e_ij,
        double* acc_i) const override;

    //! evaluate shear forces
    virtual void ShearForces(const double* dens_i, const double* dens_j, const double* vel_i,
        const double* vel_j, const double& kernelfac, const double& visc_i, const double& visc_j,
        const double& bulk_visc_i, const double& bulk_visc_j, const double& abs_rij,
        const double& specificcoefficient_ij, const double* e_ij, double* acc_i) const override;

    //! evaluate background pressure (standard formulation)
    virtual void StandardBackgroundPressure(const double* dens_i, const double* dens_j,
        const double& backgroundpressure, const double& specificcoefficient_ij, const double* e_ij,
        double* acc_i) const override;

    //! evaluate background pressure (generalized formulation)
    virtual void GeneralizedBackgroundPressure(const double* dens_i, const double* mass_i,
        const double* mass_j, const double& backgroundpressure_tilde, const double& dWdrij_tilde,
        const double* e_ij, double* acc_i) const override;

    //! evaluate modified velocity contribution
    virtual void ModifiedVelocityContribution(const double* dens_i, const double* dens_j,
        const double* vel_i, const double* vel_j, const double* mod_vel_i, const double* mod_vel_j,
        const double& specificcoefficient_ij, const double* e_ij, double* acc_i) const override;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
