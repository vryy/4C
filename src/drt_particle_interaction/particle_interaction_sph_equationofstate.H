/*---------------------------------------------------------------------------*/
/*! \file
\brief equation of state handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_EQUATIONOFSTATE_H
#define PARTICLE_INTERACTION_SPH_EQUATIONOFSTATE_H

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
  class SPHEquationOfStateBase
  {
   public:
    //! constructor
    explicit SPHEquationOfStateBase();

    //! virtual destructor
    virtual ~SPHEquationOfStateBase() = default;

    //! init equation of state handler
    virtual void Init();

    //! setup equation of state handler
    virtual void Setup();

    //! determine the pressure
    virtual double DensityToPressure(const double& density, const double& density0) const = 0;

    //! determine the density
    virtual double PressureToDensity(const double& pressure, const double& density0) const = 0;

    //! determine the energy
    virtual double DensityToEnergy(
        const double& density, const double& mass, const double& density0) const = 0;
  };

  class SPHEquationOfStateGenTait final : public SPHEquationOfStateBase
  {
   public:
    //! constructor
    explicit SPHEquationOfStateGenTait(
        const double& speedofsound, const double& refdensfac, const double& exponent);

    //! determine the pressure
    double DensityToPressure(const double& density, const double& density0) const override;

    //! determine the density
    double PressureToDensity(const double& pressure, const double& density0) const override;

    //! determine the energy
    double DensityToEnergy(
        const double& density, const double& mass, const double& density0) const override;

   private:
    //! speed of sound
    const double speedofsound_;

    //! reference density factor
    const double refdensfac_;

    //! exponent
    const double exponent_;
  };

  class SPHEquationOfStateIdealGas final : public SPHEquationOfStateBase
  {
   public:
    //! constructor
    explicit SPHEquationOfStateIdealGas(const double& speedofsound);

    //! determine the pressure
    double DensityToPressure(const double& density, const double& density0) const override;

    //! determine the density
    double PressureToDensity(const double& pressure, const double& density0) const override;

    //! determine the energy
    double DensityToEnergy(
        const double& density, const double& mass, const double& density0) const override;

   private:
    //! speed of sound
    const double speedofsound_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
