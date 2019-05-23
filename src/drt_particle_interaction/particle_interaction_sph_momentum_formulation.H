/*---------------------------------------------------------------------------*/
/*!

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
        const double* mass_i, const double* mass_j, const double& dWdrij, const double& dWdrji,
        double& speccoeff_ij, double& speccoeff_ji) const = 0;

    //! evaluate pressure gradient
    virtual void PressureGradient(const double* dens_i, const double* dens_j, const double* press_i,
        const double* press_j, const double& speccoeff_ij, const double& speccoeff_ji,
        const double* e_ij, double* acc_i, double* acc_j) const = 0;

    //! evaluate shear forces
    virtual void ShearForces(const double* dens_i, const double* dens_j, const double* vel_i,
        const double* vel_j, const double& kernelfac, const double& visc_i, const double& visc_j,
        const double& bulk_visc_i, const double& bulk_visc_j, const double& abs_rij,
        const double& speccoeff_ij, const double& speccoeff_ji, const double* e_ij, double* acc_i,
        double* acc_j) const = 0;

    //! evaluate background pressure (standard formulation)
    virtual void StandardBackgroundPressure(const double* dens_i, const double* dens_j,
        const double& bg_press_i, const double& bg_press_j, const double& speccoeff_ij,
        const double& speccoeff_ji, const double* e_ij, double* mod_acc_i,
        double* mod_acc_j) const = 0;

    //! evaluate background pressure  (generalized formulation)
    virtual void GeneralizedBackgroundPressure(const double* dens_i, const double* dens_j,
        const double* mass_i, const double* mass_j, const double& mod_bg_press_i,
        const double& mod_bg_press_j, const double& mod_dWdrij, const double& mod_dWdrji,
        const double* e_ij, double* mod_acc_i, double* mod_acc_j) const = 0;

    //! evaluate modified velocity contribution
    virtual void ModifiedVelocityContribution(const double* dens_i, const double* dens_j,
        const double* vel_i, const double* vel_j, const double* mod_vel_i, const double* mod_vel_j,
        const double& speccoeff_ij, const double& speccoeff_ji, const double* e_ij, double* acc_i,
        double* acc_j) const = 0;
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
        const double* mass_i, const double* mass_j, const double& dWdrij, const double& dWdrji,
        double& speccoeff_ij, double& speccoeff_ji) const override;

    //! evaluate pressure gradient
    virtual void PressureGradient(const double* dens_i, const double* dens_j, const double* press_i,
        const double* press_j, const double& speccoeff_ij, const double& speccoeff_ji,
        const double* e_ij, double* acc_i, double* acc_j) const override;

    //! evaluate shear forces
    virtual void ShearForces(const double* dens_i, const double* dens_j, const double* vel_i,
        const double* vel_j, const double& kernelfac, const double& visc_i, const double& visc_j,
        const double& bulk_visc_i, const double& bulk_visc_j, const double& abs_rij,
        const double& speccoeff_ij, const double& speccoeff_ji, const double* e_ij, double* acc_i,
        double* acc_j) const override;

    //! evaluate background pressure (standard formulation)
    virtual void StandardBackgroundPressure(const double* dens_i, const double* dens_j,
        const double& bg_press_i, const double& bg_press_j, const double& speccoeff_ij,
        const double& speccoeff_ji, const double* e_ij, double* mod_acc_i,
        double* mod_acc_j) const override;

    //! evaluate background pressure (generalized formulation)
    virtual void GeneralizedBackgroundPressure(const double* dens_i, const double* dens_j,
        const double* mass_i, const double* mass_j, const double& mod_bg_press_i,
        const double& mod_bg_press_j, const double& mod_dWdrij, const double& mod_dWdrji,
        const double* e_ij, double* mod_acc_i, double* mod_acc_j) const override;

    //! evaluate modified velocity contribution
    virtual void ModifiedVelocityContribution(const double* dens_i, const double* dens_j,
        const double* vel_i, const double* vel_j, const double* mod_vel_i, const double* mod_vel_j,
        const double& speccoeff_ij, const double& speccoeff_ji, const double* e_ij, double* acc_i,
        double* acc_j) const override;
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
        const double* mass_i, const double* mass_j, const double& dWdrij, const double& dWdrji,
        double& speccoeff_ij, double& speccoeff_ji) const override;

    //! evaluate pressure gradient
    virtual void PressureGradient(const double* dens_i, const double* dens_j, const double* press_i,
        const double* press_j, const double& speccoeff_ij, const double& speccoeff_ji,
        const double* e_ij, double* acc_i, double* acc_j) const override;

    //! evaluate shear forces
    virtual void ShearForces(const double* dens_i, const double* dens_j, const double* vel_i,
        const double* vel_j, const double& kernelfac, const double& visc_i, const double& visc_j,
        const double& bulk_visc_i, const double& bulk_visc_j, const double& abs_rij,
        const double& speccoeff_ij, const double& speccoeff_ji, const double* e_ij, double* acc_i,
        double* acc_j) const override;

    //! evaluate background pressure (standard formulation)
    virtual void StandardBackgroundPressure(const double* dens_i, const double* dens_j,
        const double& bg_press_i, const double& bg_press_j, const double& speccoeff_ij,
        const double& speccoeff_ji, const double* e_ij, double* mod_acc_i,
        double* mod_acc_j) const override;

    //! evaluate background pressure (generalized formulation)
    virtual void GeneralizedBackgroundPressure(const double* dens_i, const double* dens_j,
        const double* mass_i, const double* mass_j, const double& mod_bg_press_i,
        const double& mod_bg_press_j, const double& mod_dWdrij, const double& mod_dWdrji,
        const double* e_ij, double* mod_acc_i, double* mod_acc_j) const override;

    //! evaluate modified velocity contribution
    virtual void ModifiedVelocityContribution(const double* dens_i, const double* dens_j,
        const double* vel_i, const double* vel_j, const double* mod_vel_i, const double* mod_vel_j,
        const double& speccoeff_ij, const double& speccoeff_ji, const double* e_ij, double* acc_i,
        double* acc_j) const override;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
