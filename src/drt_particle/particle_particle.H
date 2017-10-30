/*----------------------------------------------------------------------*/
/*!
\file particle_particle.H

\brief Small classes to locally store particle data. General purpose

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/

#ifndef SRC_DRT_PARTICLE_PARTICLE_PARTICLE_H_
#define SRC_DRT_PARTICLE_PARTICLE_PARTICLE_H_

/*----------------------------------------------------------------------*/
/* headers */
#include "../linalg/linalg_fixedsizematrix.H"
#include "particle_enums.H"

#include <Sacado.hpp>
typedef Sacado::Fad::DFad<double> FAD;

namespace MAT
{
  namespace PAR
  {
    class ExtParticleMat;
  }
}

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
  velModVisc_(LINALG::Matrix<3,1>(true)),
  pressureMod_(-1000)
  {
    //set distinct value as long as not set from outside
    for(int i=0;i<3;i++)
    {
      acc_(i)=-1000;
      velModVisc_(i)=-1000;
    }
  }

  /// internal variables
  bool boundaryparticle_;         //bool to check if boundary particle or not?
  LINALG::Matrix<3,1> acc_;       //acceleration vector of boundary particles
  LINALG::Matrix<3,1> velModVisc_;    //modified velocity vector of boundary particle required for fiscous force calculation(see Ref. Adami2012)
  double pressureMod_;            //modified pressure of boundary particle (see Ref. Adami2012)

  }; // BoundaryData

  /// Class containing free-surface data required in ParticleMF
  class FreeSurfaceData
  {
  public:

  /// constructor: empty
  FreeSurfaceData(FreeSurfaceParticleType freesurfaceparticletype, int phase_color) :
    freesurfaceparticletype_(freesurfaceparticletype),
    phase_color_(phase_color),
    color_field_(-1000),
    colorFieldGrad_(LINALG::Matrix<3,1>(true)),
    smoothedColorFieldGrad_(LINALG::Matrix<3,1>(true)),
    density_sum_(-1000),
    validNormal_(false)
  {
  }

  /// internal variables
  FreeSurfaceParticleType freesurfaceparticletype_;         //which type of free surface particle do we have (none, direct or indirect)?
  int phase_color_;                       //every phase is characterized by a discrete value of the phase color (e.g. 0 / 1  for water / air)
  double color_field_;                      //smooth interpolation of the phase_color_: small values of the color_function \in [0;\approx 1] indicate (direct) surface particles
  LINALG::Matrix<3,1> colorFieldGrad_;       //gradient of color field above (required as surface normal)
  LINALG::Matrix<3,1> smoothedColorFieldGrad_;    //gradient of smoothed color field (alternative for surface normal)
  double density_sum_;                      //particle density determined via naiv summation (also for free-surface particles)
  bool validNormal_;

  }; // FreeSurfaceData

// ---------------------------------------------------------- //
  /// particles for meshfree dynamics (thermodynamics included)
  class ParticleMF : public ParticleBase, public ParticleThermoBase
  {
  public:

  /// constructor: empty
  ParticleMF() : ParticleBase(), ParticleThermoBase(),
      pressure_(-1), temperature_(-1), densityDot_(-1), boundarydata_(BoundaryData(false)), freesurfacedata_(FreeSurfaceData(FS_NONE, 0)), velConv_(LINALG::Matrix<3,1>(true)), extParticleMat_(NULL) {}

  // contructor: only identifiers and material
  ParticleMF(
      const int gid,
      const int owner,
      const std::vector<int> lm,
      bool boundaryparticle,
      int phase_color,
      const MAT::PAR::ExtParticleMat* extParticleMat) :
        ParticleBase(gid, owner, lm),
        ParticleThermoBase(),
        pressure_(-1),
        temperature_(-1),
        densityDot_(-1),
        boundarydata_(BoundaryData(boundaryparticle)),
        freesurfacedata_(FreeSurfaceData(FS_NONE, phase_color)),
        velConv_(LINALG::Matrix<3,1>(true)),
        extParticleMat_(extParticleMat){}

  /// internal variables
  double pressure_;
  double temperature_;
  double densityDot_;
  BoundaryData boundarydata_;
  FreeSurfaceData freesurfacedata_;
  LINALG::Matrix<3,1> velConv_;    //particle convection velocity that may differ from the particle momentum velocity in general (see Ref. Adami2013)
  const MAT::PAR::ExtParticleMat* extParticleMat_;  //pointer to the extParticleMat material
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
