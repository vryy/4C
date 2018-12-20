/*---------------------------------------------------------------------------*/
/*!
\file unit_particle_container.H

\brief unittests for smart particle container class

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#ifndef UNIT_PARTICLE_CONTAINER_H
#define UNIT_PARTICLE_CONTAINER_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include <cxxtest/TestSuite.h>
#include "src/drt_particle_engine/particle_container.H"

/*---------------------------------------------------------------------------*
 | forward declaration                                        sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleContainer_TestSuite;
}

/*---------------------------------------------------------------------------*
 | particle container test suite                              sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
class PARTICLEENGINE::ParticleContainer_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEENGINE::ParticleContainer> container_;

 public:
  void setUp()
  {
    int size = 7;
    std::set<StateEnum> stateEnumSet = {
        PARTICLEENGINE::Position, PARTICLEENGINE::Velocity, PARTICLEENGINE::Mass};

    // create, init and setup container
    container_ =
        std::unique_ptr<PARTICLEENGINE::ParticleContainer>(new PARTICLEENGINE::ParticleContainer());
    container_->Init();
    container_->Setup(size, stateEnumSet);

    // init some particles
    int index(0);
    int globalid(0);
    ParticleStates particle;
    std::vector<double> pos(3);
    std::vector<double> vel(3);
    std::vector<double> mass(1);

    // first particle
    {
      particle.clear();
      globalid = 1;
      pos[0] = 1.20;
      pos[1] = 0.70;
      pos[2] = 2.10;
      vel[0] = 0.23;
      vel[1] = 1.76;
      vel[2] = 3.89;
      mass[0] = 0.12;
      particle.insert(std::make_pair(PARTICLEENGINE::Position, pos));
      particle.insert(std::make_pair(PARTICLEENGINE::Velocity, vel));
      particle.insert(std::make_pair(PARTICLEENGINE::Mass, mass));
      container_->AddParticle(index, globalid, particle);
    }

    // second particle
    {
      particle.clear();
      globalid = 2;
      pos[0] = -1.05;
      pos[1] = 12.6;
      pos[2] = -8.54;
      vel[0] = 0.25;
      vel[1] = -21.5;
      vel[2] = 1.0;
      mass[0] = 12.34;
      particle.insert(std::make_pair(PARTICLEENGINE::Position, pos));
      particle.insert(std::make_pair(PARTICLEENGINE::Velocity, vel));
      particle.insert(std::make_pair(PARTICLEENGINE::Mass, mass));
      container_->AddParticle(index, globalid, particle);
    }

    // third particle
    {
      particle.clear();
      globalid = 3;
      pos[0] = 61.0;
      pos[1] = -2.63;
      pos[2] = 0.11;
      vel[0] = -7.35;
      vel[1] = -5.98;
      vel[2] = 1.11;
      mass[0] = 0.5;
      particle.insert(std::make_pair(PARTICLEENGINE::Position, pos));
      particle.insert(std::make_pair(PARTICLEENGINE::Velocity, vel));
      particle.insert(std::make_pair(PARTICLEENGINE::Mass, mass));
      container_->AddParticle(index, globalid, particle);
    }
  }

  void tearDown() { container_ = nullptr; }

  // note: the public functions Init(), Setup() and AddParticle() of class ParticleContainer are
  // called in setUp() and thus implicitly tested by all following unittests

  void test_IncreaseContainerSize()
  {
    container_->IncreaseContainerSize();
    TS_ASSERT_EQUALS(container_->ParticlesStored(), 3);
    TS_ASSERT_EQUALS(container_->ContainerSize(), 14);
  }

  void test_DecreaseContainerSize()
  {
    container_->DecreaseContainerSize();
    TS_ASSERT_EQUALS(container_->ParticlesStored(), 3);
    TS_ASSERT_EQUALS(container_->ContainerSize(), 3);
  }

  void test_ClearContainer()
  {
    container_->ClearContainer();
    TS_ASSERT_EQUALS(container_->ParticlesStored(), 0);
  }

  void test_AddParticle()
  {
    // init a particle
    int globalid(0);
    ParticleStates particle;
    std::vector<double> pos(3);
    std::vector<double> vel(3);
    std::vector<double> mass(1);

    globalid = 4;
    pos[0] = -1.23;
    pos[1] = 1.70;
    pos[2] = 9.10;
    vel[0] = 6.23;
    vel[1] = 2.3;
    vel[2] = 6.9;
    mass[0] = 5.12;
    particle.insert(std::make_pair(PARTICLEENGINE::Position, pos));
    particle.insert(std::make_pair(PARTICLEENGINE::Velocity, vel));
    particle.insert(std::make_pair(PARTICLEENGINE::Mass, mass));

    int index(0);
    container_->AddParticle(index, globalid, particle);
    TS_ASSERT_EQUALS(index, 3);
    TS_ASSERT_EQUALS(container_->ParticlesStored(), 4);
  }

  void test_ReplaceParticle()
  {
    // init a particle
    int globalid(0);
    ParticleStates particle;
    ParticleStates particle_reference;
    std::vector<double> pos(3);
    std::vector<double> vel(3);
    std::vector<double> mass(1);

    pos[0] = -1.23;
    pos[1] = 1.70;
    pos[2] = 9.10;
    vel[0] = 6.23;
    vel[1] = 2.3;
    vel[2] = 6.9;
    mass[0] = 5.12;
    particle_reference.insert(std::make_pair(PARTICLEENGINE::Position, pos));
    particle_reference.insert(std::make_pair(PARTICLEENGINE::Velocity, vel));
    particle_reference.insert(std::make_pair(PARTICLEENGINE::Mass, mass));

    int index = 0;

    // replace only states and leave global id untouched
    container_->ReplaceParticle(index, -1, particle_reference);
    TS_ASSERT_EQUALS(container_->ParticlesStored(), 3);
    container_->GetParticle(index, globalid, particle);
    TS_ASSERT_EQUALS(globalid, 1);
    compareParticleStates(particle_reference, particle);

    // also replace global id
    container_->ReplaceParticle(index, 4, particle_reference);
    TS_ASSERT_EQUALS(container_->ParticlesStored(), 3);
    container_->GetParticle(index, globalid, particle);
    TS_ASSERT_EQUALS(globalid, 4);
    compareParticleStates(particle_reference, particle);
  }

  void test_GetParticle()
  {
    int globalid(0);
    int globalid_reference(0);
    ParticleStates particle;
    ParticleStates particle_reference;
    std::vector<double> pos(3);
    std::vector<double> vel(3);
    std::vector<double> mass(1);

    for (int index = 0; index < 3; ++index)
    {
      particle_reference.clear();

      if (index == 0)
      {
        globalid_reference = 1;
        pos[0] = 1.20;
        pos[1] = 0.70;
        pos[2] = 2.10;
        vel[0] = 0.23;
        vel[1] = 1.76;
        vel[2] = 3.89;
        mass[0] = 0.12;
      }
      else if (index == 1)
      {
        globalid_reference = 2;
        pos[0] = -1.05;
        pos[1] = 12.6;
        pos[2] = -8.54;
        vel[0] = 0.25;
        vel[1] = -21.5;
        vel[2] = 1.0;
        mass[0] = 12.34;
      }
      else if (index == 2)
      {
        globalid_reference = 3;
        pos[0] = 61.0;
        pos[1] = -2.63;
        pos[2] = 0.11;
        vel[0] = -7.35;
        vel[1] = -5.98;
        vel[2] = 1.11;
        mass[0] = 0.5;
      }

      particle_reference.insert(std::make_pair(PARTICLEENGINE::Position, pos));
      particle_reference.insert(std::make_pair(PARTICLEENGINE::Velocity, vel));
      particle_reference.insert(std::make_pair(PARTICLEENGINE::Mass, mass));

      container_->GetParticle(index, globalid, particle);
      TS_ASSERT_EQUALS(globalid_reference, globalid);
      compareParticleStates(particle_reference, particle);
    }
  }

  void test_RemoveParticle()
  {
    int globalid(0);
    int globalid_reference(0);
    ParticleStates particle;
    ParticleStates particle_reference;
    std::vector<double> pos(3);
    std::vector<double> vel(3);
    std::vector<double> mass(1);

    container_->RemoveParticle(0);
    TS_ASSERT_EQUALS(container_->ParticlesStored(), 2);

    for (int index = 0; index < 2; ++index)
    {
      particle_reference.clear();

      if (index == 0)
      {
        globalid_reference = 3;
        pos[0] = 61.0;
        pos[1] = -2.63;
        pos[2] = 0.11;
        vel[0] = -7.35;
        vel[1] = -5.98;
        vel[2] = 1.11;
        mass[0] = 0.5;
      }
      else if (index == 1)
      {
        globalid_reference = 2;
        pos[0] = -1.05;
        pos[1] = 12.6;
        pos[2] = -8.54;
        vel[0] = 0.25;
        vel[1] = -21.5;
        vel[2] = 1.0;
        mass[0] = 12.34;
      }

      particle_reference.insert(std::make_pair(PARTICLEENGINE::Position, pos));
      particle_reference.insert(std::make_pair(PARTICLEENGINE::Velocity, vel));
      particle_reference.insert(std::make_pair(PARTICLEENGINE::Mass, mass));

      container_->GetParticle(index, globalid, particle);
      TS_ASSERT_EQUALS(globalid_reference, globalid);
      compareParticleStates(particle_reference, particle);
    }
  }

  void test_GetPtrToParticleState()
  {
    std::vector<double> pos(3);
    std::vector<double> vel(3);
    std::vector<double> mass(1);

    for (int index = 0; index < 3; ++index)
    {
      if (index == 0)
      {
        pos[0] = 1.20;
        pos[1] = 0.70;
        pos[2] = 2.10;
        vel[0] = 0.23;
        vel[1] = 1.76;
        vel[2] = 3.89;
        mass[0] = 0.12;
      }
      else if (index == 1)
      {
        pos[0] = -1.05;
        pos[1] = 12.6;
        pos[2] = -8.54;
        vel[0] = 0.25;
        vel[1] = -21.5;
        vel[2] = 1.0;
        mass[0] = 12.34;
      }
      else if (index == 2)
      {
        pos[0] = 61.0;
        pos[1] = -2.63;
        pos[2] = 0.11;
        vel[0] = -7.35;
        vel[1] = -5.98;
        vel[2] = 1.11;
        mass[0] = 0.5;
      }

      double* currpos = container_->GetPtrToParticleState(PARTICLEENGINE::Position, index);
      for (int i = 0; i < (int)pos.size(); ++i) TS_ASSERT_DELTA(currpos[i], pos[i], 1e-14);

      double* currvel = container_->GetPtrToParticleState(PARTICLEENGINE::Velocity, index);
      for (int i = 0; i < (int)pos.size(); ++i) TS_ASSERT_DELTA(currvel[i], vel[i], 1e-14);

      double* currmass = container_->GetPtrToParticleState(PARTICLEENGINE::Mass, index);
      TS_ASSERT_DELTA(currmass[0], mass[0], 1e-14);
    }
  }

  void test_GetParticleState()
  {
    std::vector<double> pos(3);
    std::vector<double> vel(3);
    std::vector<double> mass(1);

    std::vector<double> currpos(3);
    std::vector<double> currvel(3);
    std::vector<double> currmass(1);

    for (int index = 0; index < 3; ++index)
    {
      if (index == 0)
      {
        pos[0] = 1.20;
        pos[1] = 0.70;
        pos[2] = 2.10;
        vel[0] = 0.23;
        vel[1] = 1.76;
        vel[2] = 3.89;
        mass[0] = 0.12;
      }
      else if (index == 1)
      {
        pos[0] = -1.05;
        pos[1] = 12.6;
        pos[2] = -8.54;
        vel[0] = 0.25;
        vel[1] = -21.5;
        vel[2] = 1.0;
        mass[0] = 12.34;
      }
      else if (index == 2)
      {
        pos[0] = 61.0;
        pos[1] = -2.63;
        pos[2] = 0.11;
        vel[0] = -7.35;
        vel[1] = -5.98;
        vel[2] = 1.11;
        mass[0] = 0.5;
      }

      currpos = container_->GetParticleState(PARTICLEENGINE::Position, index);
      for (int i = 0; i < (int)pos.size(); ++i) TS_ASSERT_DELTA(currpos[i], pos[i], 1e-14);

      currvel = container_->GetParticleState(PARTICLEENGINE::Velocity, index);
      for (int i = 0; i < (int)pos.size(); ++i) TS_ASSERT_DELTA(currvel[i], vel[i], 1e-14);

      currmass = container_->GetParticleState(PARTICLEENGINE::Mass, index);
      TS_ASSERT_DELTA(currmass[0], mass[0], 1e-14);
    }
  }

  void test_ScaleState()
  {
    int globalid(0);
    ParticleStates particle;
    ParticleStates particle_reference;
    std::vector<double> pos(3);
    std::vector<double> vel(3);
    std::vector<double> mass(1);

    container_->ScaleState(1.5, PARTICLEENGINE::Position);
    container_->ScaleState(3.25, PARTICLEENGINE::Velocity);
    container_->ScaleState(0.95, PARTICLEENGINE::Mass);

    for (int index = 0; index < 3; ++index)
    {
      particle_reference.clear();

      if (index == 0)
      {
        pos[0] = 1.8;
        pos[1] = 1.05;
        pos[2] = 3.15;
        vel[0] = 0.7475;
        vel[1] = 5.72;
        vel[2] = 12.6425;
        mass[0] = 0.114;
      }
      else if (index == 1)
      {
        pos[0] = -1.575;
        pos[1] = 18.9;
        pos[2] = -12.81;
        vel[0] = 0.8125;
        vel[1] = -69.875;
        vel[2] = 3.25;
        mass[0] = 11.723;
      }
      else if (index == 2)
      {
        pos[0] = 91.5;
        pos[1] = -3.945;
        pos[2] = 0.165;
        vel[0] = -23.8875;
        vel[1] = -19.435;
        vel[2] = 3.6075;
        mass[0] = 0.475;
      }

      particle_reference.insert(std::make_pair(PARTICLEENGINE::Position, pos));
      particle_reference.insert(std::make_pair(PARTICLEENGINE::Velocity, vel));
      particle_reference.insert(std::make_pair(PARTICLEENGINE::Mass, mass));

      container_->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }
  }

  void test_UpdateState()
  {
    int globalid(0);
    ParticleStates particle;
    ParticleStates particle_reference;
    std::vector<double> pos(3);
    std::vector<double> vel(3);
    std::vector<double> mass(1);

    container_->UpdateState(1.0, PARTICLEENGINE::Position, 0.5, PARTICLEENGINE::Velocity);

    for (int index = 0; index < 3; ++index)
    {
      particle_reference.clear();

      if (index == 0)
      {
        pos[0] = 1.315;
        pos[1] = 1.58;
        pos[2] = 4.045;
        vel[0] = 0.23;
        vel[1] = 1.76;
        vel[2] = 3.89;
        mass[0] = 0.12;
      }
      else if (index == 1)
      {
        pos[0] = -0.925;
        pos[1] = 1.85;
        pos[2] = -8.04;
        vel[0] = 0.25;
        vel[1] = -21.5;
        vel[2] = 1.0;
        mass[0] = 12.34;
      }
      else if (index == 2)
      {
        pos[0] = 57.325;
        pos[1] = -5.62;
        pos[2] = 0.665;
        vel[0] = -7.35;
        vel[1] = -5.98;
        vel[2] = 1.11;
        mass[0] = 0.5;
      }

      particle_reference.insert(std::make_pair(PARTICLEENGINE::Position, pos));
      particle_reference.insert(std::make_pair(PARTICLEENGINE::Velocity, vel));
      particle_reference.insert(std::make_pair(PARTICLEENGINE::Mass, mass));

      container_->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }
  }

  void test_SetState()
  {
    int globalid(0);
    ParticleStates particle;
    ParticleStates particle_reference;

    std::vector<double> pos(3);
    std::vector<double> vel(3);
    std::vector<double> mass(1);

    pos[0] = 3.15;
    pos[1] = -1.45;
    pos[2] = 9.5;
    vel[0] = -21.30;
    vel[1] = -4.33;
    vel[2] = 0.933;
    mass[0] = 1.234;

    particle_reference.insert(std::make_pair(PARTICLEENGINE::Position, pos));
    particle_reference.insert(std::make_pair(PARTICLEENGINE::Velocity, vel));
    particle_reference.insert(std::make_pair(PARTICLEENGINE::Mass, mass));

    container_->SetState(pos, PARTICLEENGINE::Position);
    container_->SetState(vel, PARTICLEENGINE::Velocity);
    container_->SetState(mass, PARTICLEENGINE::Mass);

    for (int index = 0; index < 3; ++index)
    {
      container_->GetParticle(index, globalid, particle);
      compareParticleStates(particle_reference, particle);
    }
  }

  void test_ClearState()
  {
    int globalid(0);
    ParticleStates particle;
    ParticleStates particle_reference;

    std::vector<double> pos(3, 0.0);
    std::vector<double> vel(3, 0.0);
    std::vector<double> mass(1, 0.0);

    particle_reference.insert(std::make_pair(PARTICLEENGINE::Position, pos));
    particle_reference.insert(std::make_pair(PARTICLEENGINE::Velocity, vel));
    particle_reference.insert(std::make_pair(PARTICLEENGINE::Mass, mass));

    container_->ClearState(PARTICLEENGINE::Position);
    container_->ClearState(PARTICLEENGINE::Velocity);
    container_->ClearState(PARTICLEENGINE::Mass);

    for (int index = 0; index < 3; ++index)
    {
      container_->GetParticle(index, globalid, particle);
      compareParticleStates(particle_reference, particle);
    }
  }

  void test_GetMinValueOfState()
  {
    TS_ASSERT_DELTA(container_->GetMinValueOfState(PARTICLEENGINE::Mass), 0.12, 1e-14);
    TS_ASSERT_DELTA(container_->GetMinValueOfState(PARTICLEENGINE::Position), -8.54, 1e-14);
  }

  void test_GetMaxValueOfState()
  {
    TS_ASSERT_DELTA(container_->GetMaxValueOfState(PARTICLEENGINE::Mass), 12.34, 1e-14);
    TS_ASSERT_DELTA(container_->GetMaxValueOfState(PARTICLEENGINE::Position), 61.0, 1e-14);
  }

  void test_GetStoredStates()
  {
    std::set<StateEnum> particleStates = container_->GetStoredStates();
    TS_ASSERT_EQUALS(particleStates.size(), 3);
    TS_ASSERT(particleStates.find(PARTICLEENGINE::Position) != particleStates.end());
    TS_ASSERT(particleStates.find(PARTICLEENGINE::Velocity) != particleStates.end());
    TS_ASSERT(particleStates.find(PARTICLEENGINE::Mass) != particleStates.end());
  }

  void test_ContainerSize() { TS_ASSERT_EQUALS(container_->ContainerSize(), 7); }

  void test_ParticlesStored() { TS_ASSERT_EQUALS(container_->ParticlesStored(), 3); }

  void compareParticleStates(ParticleStates& particle_reference, ParticleStates& particle)
  {
    TS_ASSERT_EQUALS(particle_reference.size(), particle.size());

    for (auto& refStateIt : particle_reference)
    {
      auto stateIt = particle.find(refStateIt.first);

      TS_ASSERT(stateIt != particle.end());

      std::vector<double>& state_reference = refStateIt.second;
      std::vector<double>& state = stateIt->second;

      TS_ASSERT_EQUALS(state_reference.size(), state.size());

      for (int i = 0; i < (int)state_reference.size(); ++i)
        TS_ASSERT_DELTA(state_reference[i], state[i], 1e-14);
    }
  }
};

/*---------------------------------------------------------------------------*/
#endif
