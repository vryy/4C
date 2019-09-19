/*---------------------------------------------------------------------------*/
/*! \file
\brief equation of state handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer Sebastian Fuchs
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
 | class definitions                                                         |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHEquationOfStateBase
  {
   public:
    //! constructor
    explicit SPHEquationOfStateBase();

    //! virtual destructor
    virtual ~SPHEquationOfStateBase(){};

    //! init equation of state handler
    virtual void Init();

    //! setup equation of state handler
    virtual void Setup();

    //! write restart of equation of state handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of equation of state handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

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

    //! destructor
    ~SPHEquationOfStateGenTait() override{};

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

    //! destructor
    ~SPHEquationOfStateIdealGas() override{};

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
