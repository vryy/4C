/*----------------------------------------------------------------------*/
/*!
\file particle_sph_eos.H

\brief equation of state handler for SPH methods

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                             sfuchs 11/17 |
 *----------------------------------------------------------------------*/
#ifndef PARTICLE_SPH_EOS_H_
#define PARTICLE_SPH_EOS_H_

/*----------------------------------------------------------------------*
 | equation of state handler                               sfuchs 11/17 |
 *----------------------------------------------------------------------*/
namespace PARTICLE
{

//! equation of state base class
class EquationOfState_Base
{

public :

  //! constructor
  explicit EquationOfState_Base(){};

  //! virtual destructor
  virtual ~EquationOfState_Base(){};

  //! determine the pressure
  virtual double DensityToPressure(const double& density) = 0;

  //! determine the density
  virtual double PressureToDensity(const double& pressure) = 0;

  //! determine the energy
  virtual double DensityToEnergy(const double& density, const double& mass) = 0;

};

//! generalized tait's equation of state for water-like fluids (cf. Adami 2012, equation (3))
class EquationOfState_GenTait : public EquationOfState_Base
{

public :

  //! constructor
  explicit EquationOfState_GenTait(
      const double& speedOfSound,
      const double& initDensity,
      const double& refDensFac,
      const double& exponent
      );

  //! virtual destructor
  virtual ~EquationOfState_GenTait(){};

  //! determine the pressure
  virtual double DensityToPressure(const double& density);

  //! determine the density
  virtual double PressureToDensity(const double& pressure);

  //! determine the energy
  virtual double DensityToEnergy(const double& density, const double& mass);

private :

  //! speed of sound
  const double speedOfSound_;

  //! initial density
  const double initDensity_;

  //! reference density factor
  const double refDensFac_;

  //! exponent
  const double exponent_;

  //! rest pressure
  double initPressure_;

};

//! ideal gas law (cf. Hu and Adams 2006, equation (3) or similar Tartakovsky 2016, equation (50))
class EquationOfState_IdealGas : public EquationOfState_Base
{

public :

  //! constructor
  explicit EquationOfState_IdealGas(
      const double& speedOfSound,
      const double& initDensity
      );

  //! virtual destructor
  virtual ~EquationOfState_IdealGas(){};

  //! determine the pressure
  virtual double DensityToPressure(const double& density);

  //! determine the density
  virtual double PressureToDensity(const double& pressure);

  //! determine the energy
  virtual double DensityToEnergy(const double& density, const double& mass);

private :

  //! speed of sound
  const double speedOfSound_;

  //! initial density
  const double initDensity_;

};

}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif /* PARTICLE_SPH_EOS_H_ */
