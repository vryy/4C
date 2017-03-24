/*---------------------------------------------------------------------*/
/*!
\file particle_particle.H

\brief Small classes to locally store particle data. General purpose

\level 3

\maintainer Alessandro Cattabiani
*----------------------------------------------------------------------*/

#ifndef SRC_DRT_PARTICLE_PARTICLE_PARTICLE_H_
#define SRC_DRT_PARTICLE_PARTICLE_PARTICLE_H_

/*----------------------------------------------------------------------*/
/* headers */
#include "../linalg/linalg_fixedsizematrix.H"

#include <Sacado.hpp>
typedef Sacado::Fad::DFad<double> FAD;

namespace PARTICLE
{
// ---------------------------------------------------------- //
  /// base class for particles
  class ParticleBase
  {
  public:

  /// constructors

  ParticleBase(
        const int gid,
        const int owner,
        const std::vector<int> lm,
        const LINALG::Matrix<3,1> dis,
        const LINALG::Matrix<3,1> vel,
        const double radius,
        const double mass,
        const double density
        ) :
          gid_(gid),
          owner_(owner),
          lm_(lm),
          dis_(dis),
          vel_(vel),
          radius_(radius),
          mass_(mass),
          density_(density) {}

  ParticleBase() :
    gid_(-1),
    owner_(-1),
    radius_(-1),
    mass_(-1),
    density_(-1) {}

  ParticleBase(
        const int gid,
        const int owner,
        const std::vector<int> lm
        ) :
          gid_(gid),
          owner_(owner),
          lm_(lm),
          radius_(-1),
          mass_(-1),
          density_(-1) {}

  /// internal variables

  int gid_;
  int owner_;
  std::vector<int> lm_;
  LINALG::Matrix<3,1> dis_;
  LINALG::Matrix<3,1> vel_;
  double radius_;
  double mass_;
  double density_;

  }; // ParticleBase

// ---------------------------------------------------------- //
  /// base class for themodynamics
  class ParticleThermoBase
  {
  public:

  /// constructors

  ParticleThermoBase(const double specEnthalpy) : specEnthalpy_(specEnthalpy) {}

  ParticleThermoBase() : specEnthalpy_(-1) {}

  /// internal variables

  double specEnthalpy_;

  }; // ParticleThermoBase


// ---------------------------------------------------------- //
  /// base class for surface tension dynamics
  class ParticleSurfTensionBase
  {
  public:

  /// constructor

  ParticleSurfTensionBase(const LINALG::Matrix<3,1> colorFieldGradient) : colorFieldGradient_(colorFieldGradient) {}

  ParticleSurfTensionBase() {}

  /// internal variables

  LINALG::Matrix<3,1> colorFieldGradient_;

  }; // ParticleSurfTensionBase

// ---------------------------------------------------------- //
  /// base class for surface tension dynamics
  class ParticleMFmGradWbase
  {
  public:

  /// constructor

    ParticleMFmGradWbase(const LINALG::Matrix<3,1> mGradW) :
      mGradW_(mGradW) {}

    ParticleMFmGradWbase() {}

  /// internal variables

    LINALG::Matrix<3,1> mGradW_;

  }; // ParticleMFmGradHessWbase

// ---------------------------------------------------------- //
  /// class for material particles (no meshfree - GH stuff)
  class Particle : public ParticleBase
  {
  public:

  /// constructors

  Particle(
      const int gid,
      const int owner,
      const std::vector<int> lm,
      const LINALG::Matrix<3,1> dis,
      const LINALG::Matrix<3,1> vel,
      const double radius,
      const double mass,
      const double density,
      const LINALG::Matrix<3,1> angvel,
      const double ddt
      ) :
        ParticleBase(gid, owner, lm, dis, vel, radius, mass, density),
        angvel_(angvel),
        ddt_(ddt) {}

  Particle() : ParticleBase(), ddt_(-1) {}

  Particle(
      const int gid,
      const int owner,
      const std::vector<int> lm
      ) :
        ParticleBase(gid, owner, lm),
        ddt_(-1) {}

  /// internal variables

  LINALG::Matrix<3,1> angvel_;
  double ddt_;

  }; // Particle

// ---------------------------------------------------------- //
  /// particles with thermodynamics
  class ParticleThermo : public Particle, public ParticleThermoBase
  {
  public:

  /// constructor

  ParticleThermo(
      const int gid,
      const int owner,
      const std::vector<int> lm,
      const LINALG::Matrix<3,1> dis,
      const LINALG::Matrix<3,1> vel,
      const double radius,
      const double mass,
      const double density,
      const LINALG::Matrix<3,1> angvel,
      const double ddt,
      const double specEnthalpy
      ) :
        Particle(gid, owner, lm, dis, vel, radius, mass, density,  angvel, ddt),
        ParticleThermoBase(specEnthalpy) {}

  ParticleThermo() : Particle(), ParticleThermoBase() {}

  ParticleThermo(
      const int gid,
      const int owner,
      const std::vector<int> lm
      ) :
        Particle(gid, owner, lm),
        ParticleThermoBase() {}

  }; // ParticleThermo

// ---------------------------------------------------------- //
  /// Class containing boundary data required in ParticleMF
  class BoundaryData
  {
  public:

  /// constructor: empty
  BoundaryData(bool boundaryparticle) :
  boundaryparticle_(boundaryparticle),
  acc_(LINALG::Matrix<3,1>(true)),
  velMod_(LINALG::Matrix<3,1>(true)),
  pressureMod_(-1000),
  densityMod_(-1000)
  {
    //set distinct value as long as not set from outside
    for(int i=0;i<3;i++)
    {
      acc_(i)=-1000;
      velMod_(i)=-1000;
    }
  }

  /// internal variables
  bool boundaryparticle_;         //bool to check if boundary particle or not?
  LINALG::Matrix<3,1> acc_;       //acceleration vector of boundary particles
  LINALG::Matrix<3,1> velMod_;    //modified velocity vector of boundary particle (see Ref. Adami2012)
  double pressureMod_;            //modified pressure of boundary particle (see Ref. Adami2012)
  double densityMod_;            //modified density of boundary particle consistently derived from modified pressure via equation of state

  }; // BoundaryData

// ---------------------------------------------------------- //
  /// particles for meshfree dynamics (thermodynamics and surface tensions incuded)
  class ParticleMF : public ParticleBase, public ParticleThermoBase, public ParticleSurfTensionBase, public ParticleMFmGradWbase
  {
  public:

  /// constructor: full
  ParticleMF(
      const int gid,
      const int owner,
      const std::vector<int> lm,
      const LINALG::Matrix<3,1> dis,
      const LINALG::Matrix<3,1> vel,
      const double radius,
      const double mass,
      const double density,
      const double specEnthalpy,
      const LINALG::Matrix<3,1> colorFieldGradient,
      const double pressure,
      const double temperature,
      const double densityDot,
      const double alpha,
      const LINALG::Matrix<3,1> mGradW
      ) :
        ParticleBase(gid, owner, lm, dis, vel, radius, mass, density),
        ParticleThermoBase(specEnthalpy),
        ParticleSurfTensionBase(colorFieldGradient),
        ParticleMFmGradWbase(mGradW),
        pressure_(pressure),
        temperature_(temperature),
        densityDot_(densityDot),
        alpha_(alpha),
        boundarydata_(BoundaryData(false)) {}

  /// constructor: empty
  ParticleMF() : ParticleBase(), ParticleThermoBase(), ParticleSurfTensionBase(), ParticleMFmGradWbase(),
      pressure_(-1), temperature_(-1), densityDot_(-1), alpha_(-1), boundarydata_(BoundaryData(false)) {}

  // contructor: only identifiers
  ParticleMF(
      const int gid,
      const int owner,
      const std::vector<int> lm,
      bool boundaryparticle) :
        ParticleBase(gid, owner, lm),
        ParticleThermoBase(),
        ParticleSurfTensionBase(),
        ParticleMFmGradWbase(),
        pressure_(-1),
        temperature_(-1),
        densityDot_(-1),
        alpha_(-1),
        boundarydata_(BoundaryData(boundaryparticle)){}

  /// constructor: only the first round
  ParticleMF(
      const int gid,
      const int owner,
      const std::vector<int> lm,
      const LINALG::Matrix<3,1> dis,
      const LINALG::Matrix<3,1> vel,
      const double radius,
      const double mass,
      const double density,
      const double specEnthalpy,
      const double pressure,
      const double temperature
      ) :
        ParticleBase(gid, owner, lm, dis, vel, radius, mass, density),
        ParticleThermoBase(specEnthalpy),
        ParticleSurfTensionBase(),
        ParticleMFmGradWbase(),
        pressure_(pressure),
        temperature_(temperature),
        densityDot_(-1),
        alpha_(-1),
        boundarydata_(BoundaryData(false)){}

  /// internal variables
  double pressure_;
  double temperature_;
  double densityDot_;
  double alpha_;
  BoundaryData boundarydata_;
  }; // ParticleMF
}

/// class for particles based on FAD (for debugging)
class ParticleFAD
{
public:

/// constructors

  ParticleFAD(
      const int gid
      ) :
        gid_(gid),
        dis_(LINALG::TMatrix<FAD,3,1>(0.0)),
        density_(0.0) {}

/// internal variables

int gid_;
LINALG::TMatrix<FAD,3,1> dis_;
FAD density_;
}; // ParticleFAD

#endif /* SRC_DRT_PARTICLE_PARTICLE_PARTICLE_H_ */
