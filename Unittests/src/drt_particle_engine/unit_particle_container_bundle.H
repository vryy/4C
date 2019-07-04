/*---------------------------------------------------------------------------*/
/*!
\brief unittests for particle container bundle class

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#ifndef UNIT_PARTICLE_CONTAINER_BUNDLE_H
#define UNIT_PARTICLE_CONTAINER_BUNDLE_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#include <cxxtest/TestSuite.h>
#include "src/drt_particle_engine/particle_container_bundle.H"

/*---------------------------------------------------------------------------*
 | forward declaration                                        sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleContainerBundle_TestSuite;
}

/*---------------------------------------------------------------------------*
 | particle container bundle test suite                       sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
class PARTICLEENGINE::ParticleContainerBundle_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEENGINE::ParticleContainerBundle> particlecontainerbundle_;

  int statesvectorsize_;

 public:
  void setUp()
  {
    // set dummy rank
    int dummyrank = 0;

    // create and init particle container bundle
    particlecontainerbundle_ = std::unique_ptr<PARTICLEENGINE::ParticleContainerBundle>(
        new PARTICLEENGINE::ParticleContainerBundle(dummyrank));
    particlecontainerbundle_->Init();

    // init two phases with different particle states
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>> particlestatestotypes;
    std::set<PARTICLEENGINE::StateEnum> stateEnumSet = {
        PARTICLEENGINE::Position, PARTICLEENGINE::Mass, PARTICLEENGINE::Radius};
    particlestatestotypes.insert(std::make_pair(PARTICLEENGINE::Phase1, stateEnumSet));
    particlestatestotypes.insert(std::make_pair(PARTICLEENGINE::Phase2, stateEnumSet));

    // setup particle container bundle
    particlecontainerbundle_->Setup(particlestatestotypes);

    statesvectorsize_ = *(--stateEnumSet.end()) + 1;

    // init some particles
    int index(0);
    int globalid(0);

    ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>(0));

    std::vector<double> pos(3);
    std::vector<double> mass(1);
    std::vector<double> rad(1);

    // owned particles for phase 1
    {
      ParticleContainer* container = particlecontainerbundle_->GetSpecificContainer(
          PARTICLEENGINE::Phase1, PARTICLEENGINE::Owned);

      // first particle
      globalid = 1;
      pos[0] = 1.20;
      pos[1] = 0.70;
      pos[2] = 2.10;
      mass[0] = 0.1;
      rad[0] = 0.12;
      particle[PARTICLEENGINE::Position] = pos;
      particle[PARTICLEENGINE::Mass] = mass;
      particle[PARTICLEENGINE::Radius] = rad;
      container->AddParticle(index, globalid, particle);

      // second particle
      globalid = 2;
      pos[0] = -1.05;
      pos[1] = 12.6;
      pos[2] = -8.54;
      mass[0] = 0.5;
      rad[0] = 12.34;
      particle[PARTICLEENGINE::Position] = pos;
      particle[PARTICLEENGINE::Mass] = mass;
      particle[PARTICLEENGINE::Radius] = rad;
      container->AddParticle(index, globalid, particle);

      // third particle
      globalid = 3;
      pos[0] = -5.02;
      pos[1] = 2.26;
      pos[2] = -7.4;
      mass[0] = 0.2;
      rad[0] = 2.9;
      particle[PARTICLEENGINE::Position] = pos;
      particle[PARTICLEENGINE::Mass] = mass;
      particle[PARTICLEENGINE::Radius] = rad;
      container->AddParticle(index, globalid, particle);
    }

    // ghosted particles for phase 1
    {
      ParticleContainer* container = particlecontainerbundle_->GetSpecificContainer(
          PARTICLEENGINE::Phase1, PARTICLEENGINE::Ghosted);

      // first particle
      globalid = 4;
      pos[0] = 2.20;
      pos[1] = -0.52;
      pos[2] = 1.10;
      mass[0] = 0.8;
      rad[0] = 3.12;
      particle[PARTICLEENGINE::Position] = pos;
      particle[PARTICLEENGINE::Mass] = mass;
      particle[PARTICLEENGINE::Radius] = rad;
      container->AddParticle(index, globalid, particle);

      // second particle
      globalid = 5;
      pos[0] = -16.08;
      pos[1] = 1.46;
      pos[2] = -3.54;
      mass[0] = 1.4;
      rad[0] = 1.4;
      particle[PARTICLEENGINE::Position] = pos;
      particle[PARTICLEENGINE::Mass] = mass;
      particle[PARTICLEENGINE::Radius] = rad;
      container->AddParticle(index, globalid, particle);
    }

    // owned particles for phase 2
    {
      ParticleContainer* container = particlecontainerbundle_->GetSpecificContainer(
          PARTICLEENGINE::Phase2, PARTICLEENGINE::Owned);

      // first particle
      globalid = 6;
      pos[0] = 0.24;
      pos[1] = -1.71;
      pos[2] = -2.15;
      mass[0] = 1.91;
      rad[0] = 2.2;
      particle[PARTICLEENGINE::Position] = pos;
      particle[PARTICLEENGINE::Mass] = mass;
      particle[PARTICLEENGINE::Radius] = rad;
      container->AddParticle(index, globalid, particle);

      // second particle
      globalid = 7;
      pos[0] = -1.15;
      pos[1] = 2.6;
      pos[2] = 7.24;
      mass[0] = 0.4;
      rad[0] = 1.2;
      particle[PARTICLEENGINE::Position] = pos;
      particle[PARTICLEENGINE::Mass] = mass;
      particle[PARTICLEENGINE::Radius] = rad;
      container->AddParticle(index, globalid, particle);

      // third particle
      globalid = 8;
      pos[0] = 5.12;
      pos[1] = 4.26;
      pos[2] = -3.4;
      mass[0] = 1.1;
      rad[0] = 0.2;
      particle[PARTICLEENGINE::Position] = pos;
      particle[PARTICLEENGINE::Mass] = mass;
      particle[PARTICLEENGINE::Radius] = rad;
      container->AddParticle(index, globalid, particle);
    }
  }

  void tearDown() { particlecontainerbundle_ = nullptr; }

  // note: the public functions Init(), Setup() and GetSpecificContainer() of class
  // ParticleContainerBundle are called in setUp() and thus implicitly tested by all following
  // unittests

  void test_ScaleStateSpecificContainer()
  {
    particlecontainerbundle_->ScaleStateSpecificContainer(
        2.0, PARTICLEENGINE::Radius, PARTICLEENGINE::Phase1);

    ParticleContainer* container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase1, PARTICLEENGINE::Owned);

    TS_ASSERT_EQUALS(container->ParticlesStored(), 3);

    int globalid(0);

    ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>(0));
    ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>(0));

    std::vector<double> pos(3);
    std::vector<double> mass(1);
    std::vector<double> rad(1);

    for (int index = 0; index < 3; ++index)
    {
      if (index == 0)
      {
        pos[0] = 1.20;
        pos[1] = 0.70;
        pos[2] = 2.10;
        mass[0] = 0.1;
        rad[0] = 0.24;
      }
      else if (index == 1)
      {
        pos[0] = -1.05;
        pos[1] = 12.6;
        pos[2] = -8.54;
        mass[0] = 0.5;
        rad[0] = 24.68;
      }
      else if (index == 2)
      {
        pos[0] = -5.02;
        pos[1] = 2.26;
        pos[2] = -7.4;
        mass[0] = 0.2;
        rad[0] = 5.8;
      }

      particle_reference[PARTICLEENGINE::Position] = pos;
      particle_reference[PARTICLEENGINE::Mass] = mass;
      particle_reference[PARTICLEENGINE::Radius] = rad;

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }
  }

  void test_UpdateStateSpecificContainer()
  {
    particlecontainerbundle_->UpdateStateSpecificContainer(
        2.0, PARTICLEENGINE::Radius, 1.0, PARTICLEENGINE::Mass, PARTICLEENGINE::Phase1);

    ParticleContainer* container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase1, PARTICLEENGINE::Owned);

    TS_ASSERT_EQUALS(container->ParticlesStored(), 3);

    int globalid(0);

    ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>(0));
    ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>(0));

    std::vector<double> pos(3);
    std::vector<double> mass(1);
    std::vector<double> rad(1);

    for (int index = 0; index < 3; ++index)
    {
      if (index == 0)
      {
        pos[0] = 1.20;
        pos[1] = 0.70;
        pos[2] = 2.10;
        mass[0] = 0.1;
        rad[0] = 0.34;
      }
      else if (index == 1)
      {
        pos[0] = -1.05;
        pos[1] = 12.6;
        pos[2] = -8.54;
        mass[0] = 0.5;
        rad[0] = 25.18;
      }
      else if (index == 2)
      {
        pos[0] = -5.02;
        pos[1] = 2.26;
        pos[2] = -7.4;
        mass[0] = 0.2;
        rad[0] = 6.0;
      }

      particle_reference[PARTICLEENGINE::Position] = pos;
      particle_reference[PARTICLEENGINE::Mass] = mass;
      particle_reference[PARTICLEENGINE::Radius] = rad;

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }
  }

  void test_SetStateSpecificContainer()
  {
    std::vector<double> mass(1);
    mass[0] = 1.1;

    particlecontainerbundle_->SetStateSpecificContainer(
        mass, PARTICLEENGINE::Mass, PARTICLEENGINE::Phase2);

    ParticleContainer* container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase2, PARTICLEENGINE::Owned);

    TS_ASSERT_EQUALS(container->ParticlesStored(), 3);

    int globalid(0);

    ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>(0));
    ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>(0));

    std::vector<double> pos(3);
    std::vector<double> rad(1);

    for (int index = 0; index < 3; ++index)
    {
      if (index == 0)
      {
        pos[0] = 0.24;
        pos[1] = -1.71;
        pos[2] = -2.15;
        rad[0] = 2.2;
      }
      else if (index == 1)
      {
        pos[0] = -1.15;
        pos[1] = 2.6;
        pos[2] = 7.24;
        rad[0] = 1.2;
      }
      else if (index == 2)
      {
        pos[0] = 5.12;
        pos[1] = 4.26;
        pos[2] = -3.4;
        rad[0] = 0.2;
      }

      particle_reference[PARTICLEENGINE::Position] = pos;
      particle_reference[PARTICLEENGINE::Mass] = mass;
      particle_reference[PARTICLEENGINE::Radius] = rad;

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }
  }

  void test_ClearStateSpecificContainer()
  {
    std::vector<double> mass(1);
    mass[0] = 0.0;

    particlecontainerbundle_->ClearStateSpecificContainer(
        PARTICLEENGINE::Mass, PARTICLEENGINE::Phase2);

    ParticleContainer* container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase2, PARTICLEENGINE::Owned);

    TS_ASSERT_EQUALS(container->ParticlesStored(), 3);

    int globalid(0);

    ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>(0));
    ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>(0));

    std::vector<double> pos(3);
    std::vector<double> rad(1);

    for (int index = 0; index < 3; ++index)
    {
      if (index == 0)
      {
        pos[0] = 0.24;
        pos[1] = -1.71;
        pos[2] = -2.15;
        rad[0] = 2.2;
      }
      else if (index == 1)
      {
        pos[0] = -1.15;
        pos[1] = 2.6;
        pos[2] = 7.24;
        rad[0] = 1.2;
      }
      else if (index == 2)
      {
        pos[0] = 5.12;
        pos[1] = 4.26;
        pos[2] = -3.4;
        rad[0] = 0.2;
      }

      particle_reference[PARTICLEENGINE::Position] = pos;
      particle_reference[PARTICLEENGINE::Mass] = mass;
      particle_reference[PARTICLEENGINE::Radius] = rad;

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }
  }

  void test_ScaleStateAllContainers()
  {
    particlecontainerbundle_->ScaleStateAllContainers(2.0, PARTICLEENGINE::Mass);

    ParticleContainer* container = nullptr;
    int globalid(0);

    ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>(0));
    ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>(0));

    std::vector<double> pos(3);
    std::vector<double> mass(1);
    std::vector<double> rad(1);

    container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase1, PARTICLEENGINE::Owned);

    TS_ASSERT_EQUALS(container->ParticlesStored(), 3);

    for (int index = 0; index < 3; ++index)
    {
      if (index == 0)
      {
        pos[0] = 1.20;
        pos[1] = 0.70;
        pos[2] = 2.10;
        mass[0] = 0.2;
        rad[0] = 0.12;
      }
      else if (index == 1)
      {
        pos[0] = -1.05;
        pos[1] = 12.6;
        pos[2] = -8.54;
        mass[0] = 1.0;
        rad[0] = 12.34;
      }
      else if (index == 2)
      {
        pos[0] = -5.02;
        pos[1] = 2.26;
        pos[2] = -7.4;
        mass[0] = 0.4;
        rad[0] = 2.9;
      }

      particle_reference[PARTICLEENGINE::Position] = pos;
      particle_reference[PARTICLEENGINE::Mass] = mass;
      particle_reference[PARTICLEENGINE::Radius] = rad;

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }

    container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase2, PARTICLEENGINE::Owned);

    TS_ASSERT_EQUALS(container->ParticlesStored(), 3);

    for (int index = 0; index < 3; ++index)
    {
      if (index == 0)
      {
        pos[0] = 0.24;
        pos[1] = -1.71;
        pos[2] = -2.15;
        mass[0] = 3.82;
        rad[0] = 2.2;
      }
      else if (index == 1)
      {
        pos[0] = -1.15;
        pos[1] = 2.6;
        pos[2] = 7.24;
        mass[0] = 0.8;
        rad[0] = 1.2;
      }
      else if (index == 2)
      {
        pos[0] = 5.12;
        pos[1] = 4.26;
        pos[2] = -3.4;
        mass[0] = 2.2;
        rad[0] = 0.2;
      }

      particle_reference[PARTICLEENGINE::Position] = pos;
      particle_reference[PARTICLEENGINE::Mass] = mass;
      particle_reference[PARTICLEENGINE::Radius] = rad;

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }
  }

  void test_UpdateStateAllContainers()
  {
    particlecontainerbundle_->UpdateStateAllContainers(
        2.0, PARTICLEENGINE::Mass, 1.0, PARTICLEENGINE::Radius);

    ParticleContainer* container = nullptr;
    int globalid(0);

    ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>(0));
    ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>(0));

    std::vector<double> pos(3);
    std::vector<double> mass(1);
    std::vector<double> rad(1);

    container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase1, PARTICLEENGINE::Owned);

    TS_ASSERT_EQUALS(container->ParticlesStored(), 3);

    for (int index = 0; index < 3; ++index)
    {
      if (index == 0)
      {
        pos[0] = 1.20;
        pos[1] = 0.70;
        pos[2] = 2.10;
        mass[0] = 0.32;
        rad[0] = 0.12;
      }
      else if (index == 1)
      {
        pos[0] = -1.05;
        pos[1] = 12.6;
        pos[2] = -8.54;
        mass[0] = 13.34;
        rad[0] = 12.34;
      }
      else if (index == 2)
      {
        pos[0] = -5.02;
        pos[1] = 2.26;
        pos[2] = -7.4;
        mass[0] = 3.3;
        rad[0] = 2.9;
      }

      particle_reference[PARTICLEENGINE::Position] = pos;
      particle_reference[PARTICLEENGINE::Mass] = mass;
      particle_reference[PARTICLEENGINE::Radius] = rad;

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }

    container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase2, PARTICLEENGINE::Owned);

    TS_ASSERT_EQUALS(container->ParticlesStored(), 3);

    for (int index = 0; index < 3; ++index)
    {
      if (index == 0)
      {
        pos[0] = 0.24;
        pos[1] = -1.71;
        pos[2] = -2.15;
        mass[0] = 6.02;
        rad[0] = 2.2;
      }
      else if (index == 1)
      {
        pos[0] = -1.15;
        pos[1] = 2.6;
        pos[2] = 7.24;
        mass[0] = 2.0;
        rad[0] = 1.2;
      }
      else if (index == 2)
      {
        pos[0] = 5.12;
        pos[1] = 4.26;
        pos[2] = -3.4;
        mass[0] = 2.4;
        rad[0] = 0.2;
      }

      particle_reference[PARTICLEENGINE::Position] = pos;
      particle_reference[PARTICLEENGINE::Mass] = mass;
      particle_reference[PARTICLEENGINE::Radius] = rad;

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }
  }

  void test_SetStateAllContainers()
  {
    std::vector<double> mass(1);
    mass[0] = 1.1;

    particlecontainerbundle_->SetStateAllContainers(mass, PARTICLEENGINE::Mass);

    ParticleContainer* container = nullptr;
    int globalid(0);

    ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>(0));
    ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>(0));

    std::vector<double> pos(3);
    std::vector<double> rad(1);

    container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase1, PARTICLEENGINE::Owned);

    TS_ASSERT_EQUALS(container->ParticlesStored(), 3);

    for (int index = 0; index < 3; ++index)
    {
      if (index == 0)
      {
        pos[0] = 1.20;
        pos[1] = 0.70;
        pos[2] = 2.10;
        rad[0] = 0.12;
      }
      else if (index == 1)
      {
        pos[0] = -1.05;
        pos[1] = 12.6;
        pos[2] = -8.54;
        rad[0] = 12.34;
      }
      else if (index == 2)
      {
        pos[0] = -5.02;
        pos[1] = 2.26;
        pos[2] = -7.4;
        rad[0] = 2.9;
      }

      particle_reference[PARTICLEENGINE::Position] = pos;
      particle_reference[PARTICLEENGINE::Mass] = mass;
      particle_reference[PARTICLEENGINE::Radius] = rad;

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }

    container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase2, PARTICLEENGINE::Owned);

    TS_ASSERT_EQUALS(container->ParticlesStored(), 3);

    for (int index = 0; index < 3; ++index)
    {
      if (index == 0)
      {
        pos[0] = 0.24;
        pos[1] = -1.71;
        pos[2] = -2.15;
        rad[0] = 2.2;
      }
      else if (index == 1)
      {
        pos[0] = -1.15;
        pos[1] = 2.6;
        pos[2] = 7.24;
        rad[0] = 1.2;
      }
      else if (index == 2)
      {
        pos[0] = 5.12;
        pos[1] = 4.26;
        pos[2] = -3.4;
        rad[0] = 0.2;
      }

      particle_reference[PARTICLEENGINE::Position] = pos;
      particle_reference[PARTICLEENGINE::Mass] = mass;
      particle_reference[PARTICLEENGINE::Radius] = rad;

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }
  }

  void test_ClearStateAllContainers()
  {
    std::vector<double> mass(1);
    mass[0] = 0.0;

    particlecontainerbundle_->ClearStateAllContainers(PARTICLEENGINE::Mass);

    ParticleContainer* container = nullptr;
    int globalid(0);

    ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>(0));
    ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>(0));

    std::vector<double> pos(3);
    std::vector<double> rad(1);

    container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase1, PARTICLEENGINE::Owned);

    TS_ASSERT_EQUALS(container->ParticlesStored(), 3);

    for (int index = 0; index < 3; ++index)
    {
      if (index == 0)
      {
        pos[0] = 1.20;
        pos[1] = 0.70;
        pos[2] = 2.10;
        rad[0] = 0.12;
      }
      else if (index == 1)
      {
        pos[0] = -1.05;
        pos[1] = 12.6;
        pos[2] = -8.54;
        rad[0] = 12.34;
      }
      else if (index == 2)
      {
        pos[0] = -5.02;
        pos[1] = 2.26;
        pos[2] = -7.4;
        rad[0] = 2.9;
      }

      particle_reference[PARTICLEENGINE::Position] = pos;
      particle_reference[PARTICLEENGINE::Mass] = mass;
      particle_reference[PARTICLEENGINE::Radius] = rad;

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }

    container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase2, PARTICLEENGINE::Owned);

    TS_ASSERT_EQUALS(container->ParticlesStored(), 3);

    for (int index = 0; index < 3; ++index)
    {
      if (index == 0)
      {
        pos[0] = 0.24;
        pos[1] = -1.71;
        pos[2] = -2.15;
        rad[0] = 2.2;
      }
      else if (index == 1)
      {
        pos[0] = -1.15;
        pos[1] = 2.6;
        pos[2] = 7.24;
        rad[0] = 1.2;
      }
      else if (index == 2)
      {
        pos[0] = 5.12;
        pos[1] = 4.26;
        pos[2] = -3.4;
        rad[0] = 0.2;
      }

      particle_reference[PARTICLEENGINE::Position] = pos;
      particle_reference[PARTICLEENGINE::Mass] = mass;
      particle_reference[PARTICLEENGINE::Radius] = rad;

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }
  }

  void test_CheckAndDecreaseSizeAllContainersOfSpecificStatus()
  {
    ParticleContainer* container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase1, PARTICLEENGINE::Owned);

    TS_ASSERT_EQUALS(container->ParticlesStored(), 3);
    TS_ASSERT_EQUALS(container->ContainerSize(), 4);

    container->RemoveParticle(0);
    container->RemoveParticle(0);

    particlecontainerbundle_->CheckAndDecreaseSizeAllContainersOfSpecificStatus(
        PARTICLEENGINE::Owned);

    TS_ASSERT_EQUALS(container->ParticlesStored(), 1);
    TS_ASSERT_EQUALS(container->ContainerSize(), 2);
  }

  void test_ClearAllContainersOfSpecificStatus()
  {
    particlecontainerbundle_->ClearAllContainersOfSpecificStatus(PARTICLEENGINE::Ghosted);

    ParticleContainer* container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase1, PARTICLEENGINE::Ghosted);

    TS_ASSERT_EQUALS(container->ParticlesStored(), 0);
  }

  void test_GetVectorOfParticleObjectsOfAllContainers()
  {
    std::vector<ParticleObjShrdPtr> particlesstored;

    particlecontainerbundle_->GetVectorOfParticleObjectsOfAllContainers(particlesstored);

    TS_ASSERT_EQUALS(particlesstored.size(), 6);
  }

  void compareParticleStates(ParticleStates& particle_reference, ParticleStates& particle)
  {
    TS_ASSERT_EQUALS(particle_reference.size(), particle.size());

    for (int i = 0; i < (int)particle.size(); ++i)
    {
      std::vector<double>& state_reference = particle_reference[i];
      std::vector<double>& state = particle[i];

      TS_ASSERT_EQUALS(state_reference.size(), state.size());

      for (int i = 0; i < (int)state_reference.size(); ++i)
        TS_ASSERT_DELTA(state_reference[i], state[i], 1e-14);
    }
  }
};

/*---------------------------------------------------------------------------*/
#endif
