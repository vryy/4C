/*---------------------------------------------------------------------------*/
/*! \file
\brief unittests for particle container bundle class
\level 3
*/
/*---------------------------------------------------------------------------*/

#include <gtest/gtest.h>
#include "particle_engine_container_bundle.H"


namespace
{
  class ParticleContainerBundleTest : public ::testing::Test
  {
   protected:
    std::unique_ptr<PARTICLEENGINE::ParticleContainerBundle> particlecontainerbundle_;

    int statesvectorsize_;

    ParticleContainerBundleTest()
    {
      // create and init particle container bundle
      particlecontainerbundle_ = std::make_unique<PARTICLEENGINE::ParticleContainerBundle>();
      particlecontainerbundle_->Init();

      // init two phases with different particle states
      std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>> particlestatestotypes;
      std::set<PARTICLEENGINE::StateEnum> stateEnumSet = {
          PARTICLEENGINE::Position, PARTICLEENGINE::Mass, PARTICLEENGINE::Radius};
      particlestatestotypes.insert(std::make_pair(PARTICLEENGINE::Phase1, stateEnumSet));
      particlestatestotypes.insert(std::make_pair(PARTICLEENGINE::Phase2, stateEnumSet));

      // setup particle container bundle
      particlecontainerbundle_->Setup(particlestatestotypes);

      const auto GetMaximumStoredStateEnumSetValue = [&stateEnumSet]()
      { return *(--stateEnumSet.end()); };
      statesvectorsize_ = GetMaximumStoredStateEnumSetValue() + 1;

      // init some particles
      int index(0);
      int globalid(0);

      PARTICLEENGINE::ParticleStates particle;
      particle.assign(statesvectorsize_, std::vector<double>{});

      // owned particles for phase 1
      {
        PARTICLEENGINE::ParticleContainer* container =
            particlecontainerbundle_->GetSpecificContainer(
                PARTICLEENGINE::Phase1, PARTICLEENGINE::Owned);

        // first particle
        globalid = 1;
        particle = createTestParticle({1.20, 0.70, 2.10}, {0.1}, {0.12});
        container->AddParticle(index, globalid, particle);

        // second particle
        globalid = 2;
        particle = createTestParticle({-1.05, 12.6, -8.54}, {0.5}, {12.34});
        container->AddParticle(index, globalid, particle);

        // third particle
        globalid = 3;
        particle = createTestParticle({-5.02, 2.26, -7.4}, {0.2}, {2.9});
        container->AddParticle(index, globalid, particle);
      }

      // ghosted particles for phase 1
      {
        PARTICLEENGINE::ParticleContainer* container =
            particlecontainerbundle_->GetSpecificContainer(
                PARTICLEENGINE::Phase1, PARTICLEENGINE::Ghosted);

        // first particle
        globalid = 4;
        particle = createTestParticle({2.20, -0.52, 1.10}, {0.8}, {3.12});
        container->AddParticle(index, globalid, particle);

        // second particle
        globalid = 5;
        particle = createTestParticle({-16.08, 1.46, -3.54}, {1.4}, {1.4});
        container->AddParticle(index, globalid, particle);
      }

      // owned particles for phase 2
      {
        PARTICLEENGINE::ParticleContainer* container =
            particlecontainerbundle_->GetSpecificContainer(
                PARTICLEENGINE::Phase2, PARTICLEENGINE::Owned);

        // first particle
        globalid = 6;
        particle = createTestParticle({0.24, -1.71, -2.15}, {1.91}, {2.2});
        container->AddParticle(index, globalid, particle);

        // second particle
        globalid = 7;
        particle = createTestParticle({-1.15, 2.6, 7.24}, {0.4}, {1.2});
        container->AddParticle(index, globalid, particle);

        // third particle
        globalid = 8;
        particle = createTestParticle({5.12, 4.26, -3.4}, {1.1}, {0.2});
        container->AddParticle(index, globalid, particle);
      }
    }

    PARTICLEENGINE::ParticleStates createTestParticle(
        std::vector<double> pos, std::vector<double> mass, std::vector<double> rad)
    {
      PARTICLEENGINE::ParticleStates particle;
      particle.assign(statesvectorsize_, std::vector<double>{});

      particle[PARTICLEENGINE::Position] = pos;
      particle[PARTICLEENGINE::Mass] = mass;
      particle[PARTICLEENGINE::Radius] = rad;

      return particle;
    }

    // note: the public functions Init(), Setup() and GetSpecificContainer() of class
    // ParticleContainerBundle are called in the constructor and thus implicitly tested by all
    // following unittests
  };

  void compareParticleStates(
      PARTICLEENGINE::ParticleStates& particle_reference, PARTICLEENGINE::ParticleStates& particle)
  {
    ASSERT_EQ(particle_reference.size(), particle.size());

    for (std::size_t i = 0; i < particle.size(); ++i)
    {
      std::vector<double>& state_reference = particle_reference[i];
      std::vector<double>& state = particle[i];

      ASSERT_EQ(state_reference.size(), state.size());

      for (std::size_t j = 0; j < state_reference.size(); ++j)
        EXPECT_NEAR(state_reference[j], state[j], 1e-14)
            << "state '"
            << PARTICLEENGINE::EnumToStateName(static_cast<PARTICLEENGINE::ParticleState>(i))
            << "' j = " << j;
    }
  }

  TEST_F(ParticleContainerBundleTest, ScaleStateSpecificContainer)
  {
    particlecontainerbundle_->ScaleStateSpecificContainer(
        2.0, PARTICLEENGINE::Radius, PARTICLEENGINE::Phase1);

    PARTICLEENGINE::ParticleContainer* container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase1, PARTICLEENGINE::Owned);

    ASSERT_EQ(container->ParticlesStored(), 3);

    int globalid(0);

    PARTICLEENGINE::ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>{});
    PARTICLEENGINE::ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>{});

    for (int index = 0; index < 3; ++index)
    {
      SCOPED_TRACE("Particle " + std::to_string(index));
      if (index == 0)
      {
        particle_reference = createTestParticle({1.20, 0.70, 2.10}, {0.1}, {0.24});
      }
      else if (index == 1)
      {
        particle_reference = createTestParticle({-1.05, 12.6, -8.54}, {0.5}, {24.68});
      }
      else if (index == 2)
      {
        particle_reference = createTestParticle({-5.02, 2.26, -7.4}, {0.2}, {5.8});
      }

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }
  }

  TEST_F(ParticleContainerBundleTest, UpdateStateSpecificContainer)
  {
    particlecontainerbundle_->UpdateStateSpecificContainer(
        2.0, PARTICLEENGINE::Radius, 1.0, PARTICLEENGINE::Mass, PARTICLEENGINE::Phase1);

    PARTICLEENGINE::ParticleContainer* container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase1, PARTICLEENGINE::Owned);

    ASSERT_EQ(container->ParticlesStored(), 3);

    int globalid(0);

    PARTICLEENGINE::ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>{});
    PARTICLEENGINE::ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>{});

    for (int index = 0; index < 3; ++index)
    {
      SCOPED_TRACE("Particle " + std::to_string(index));
      if (index == 0)
      {
        particle_reference = createTestParticle({1.20, 0.70, 2.10}, {0.1}, {0.34});
      }
      else if (index == 1)
      {
        particle_reference = createTestParticle({-1.05, 12.6, -8.54}, {0.5}, {25.18});
      }
      else if (index == 2)
      {
        particle_reference = createTestParticle({-5.02, 2.26, -7.4}, {0.2}, {6.0});
      }

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }
  }

  TEST_F(ParticleContainerBundleTest, SetStateSpecificContainer)
  {
    std::vector<double> mass{1.1};

    particlecontainerbundle_->SetStateSpecificContainer(
        mass, PARTICLEENGINE::Mass, PARTICLEENGINE::Phase2);

    PARTICLEENGINE::ParticleContainer* container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase2, PARTICLEENGINE::Owned);

    ASSERT_EQ(container->ParticlesStored(), 3);

    int globalid(0);

    PARTICLEENGINE::ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>{});
    PARTICLEENGINE::ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>{});

    for (int index = 0; index < 3; ++index)
    {
      SCOPED_TRACE("Particle " + std::to_string(index));
      if (index == 0)
      {
        particle_reference = createTestParticle({0.24, -1.71, -2.15}, mass, {2.2});
      }
      else if (index == 1)
      {
        particle_reference = createTestParticle({-1.15, 2.6, 7.24}, mass, {1.2});
      }
      else if (index == 2)
      {
        particle_reference = createTestParticle({5.12, 4.26, -3.4}, mass, {0.2});
      }

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }
  }

  TEST_F(ParticleContainerBundleTest, ClearStateSpecificContainer)
  {
    std::vector<double> mass{0.0};

    particlecontainerbundle_->ClearStateSpecificContainer(
        PARTICLEENGINE::Mass, PARTICLEENGINE::Phase2);

    PARTICLEENGINE::ParticleContainer* container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase2, PARTICLEENGINE::Owned);

    ASSERT_EQ(container->ParticlesStored(), 3);

    int globalid(0);

    PARTICLEENGINE::ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>{});
    PARTICLEENGINE::ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>{});

    for (int index = 0; index < 3; ++index)
    {
      SCOPED_TRACE("Particle " + std::to_string(index));
      if (index == 0)
      {
        particle_reference = createTestParticle({0.24, -1.71, -2.15}, mass, {2.2});
      }
      else if (index == 1)
      {
        particle_reference = createTestParticle({-1.15, 2.6, 7.24}, mass, {1.2});
      }
      else if (index == 2)
      {
        particle_reference = createTestParticle({5.12, 4.26, -3.4}, mass, {0.2});
      }

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }
  }

  TEST_F(ParticleContainerBundleTest, ScaleStateAllContainers)
  {
    particlecontainerbundle_->ScaleStateAllContainers(2.0, PARTICLEENGINE::Mass);

    PARTICLEENGINE::ParticleContainer* container = nullptr;
    int globalid(0);

    PARTICLEENGINE::ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>{});
    PARTICLEENGINE::ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>{});

    container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase1, PARTICLEENGINE::Owned);

    ASSERT_EQ(container->ParticlesStored(), 3);

    for (int index = 0; index < 3; ++index)
    {
      SCOPED_TRACE("Phase1, Particle " + std::to_string(index));
      if (index == 0)
      {
        particle_reference = createTestParticle({1.20, 0.70, 2.10}, {0.2}, {0.12});
      }
      else if (index == 1)
      {
        particle_reference = createTestParticle({-1.05, 12.6, -8.54}, {1.0}, {12.34});
      }
      else if (index == 2)
      {
        particle_reference = createTestParticle({-5.02, 2.26, -7.4}, {0.4}, {2.9});
      }

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }

    container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase2, PARTICLEENGINE::Owned);

    ASSERT_EQ(container->ParticlesStored(), 3);

    for (int index = 0; index < 3; ++index)
    {
      SCOPED_TRACE("Phase2, Particle " + std::to_string(index));
      if (index == 0)
      {
        particle_reference = createTestParticle({0.24, -1.71, -2.15}, {3.82}, {2.2});
      }
      else if (index == 1)
      {
        particle_reference = createTestParticle({-1.15, 2.6, 7.24}, {0.8}, {1.2});
      }
      else if (index == 2)
      {
        particle_reference = createTestParticle({5.12, 4.26, -3.4}, {2.2}, {0.2});
      }

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }
  }

  TEST_F(ParticleContainerBundleTest, UpdateStateAllContainers)
  {
    particlecontainerbundle_->UpdateStateAllContainers(
        2.0, PARTICLEENGINE::Mass, 1.0, PARTICLEENGINE::Radius);

    PARTICLEENGINE::ParticleContainer* container = nullptr;
    int globalid(0);

    PARTICLEENGINE::ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>{});
    PARTICLEENGINE::ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>{});

    container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase1, PARTICLEENGINE::Owned);

    ASSERT_EQ(container->ParticlesStored(), 3);

    for (int index = 0; index < 3; ++index)
    {
      SCOPED_TRACE("Phase1, Particle " + std::to_string(index));
      if (index == 0)
      {
        particle_reference = createTestParticle({1.20, 0.70, 2.10}, {0.32}, {0.12});
      }
      else if (index == 1)
      {
        particle_reference = createTestParticle({-1.05, 12.6, -8.54}, {13.34}, {12.34});
      }
      else if (index == 2)
      {
        particle_reference = createTestParticle({-5.02, 2.26, -7.4}, {3.3}, {2.9});
      }

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }

    container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase2, PARTICLEENGINE::Owned);

    ASSERT_EQ(container->ParticlesStored(), 3);

    for (int index = 0; index < 3; ++index)
    {
      SCOPED_TRACE("Phase2, Particle " + std::to_string(index));
      if (index == 0)
      {
        particle_reference = createTestParticle({0.24, -1.71, -2.15}, {6.02}, {2.2});
      }
      else if (index == 1)
      {
        particle_reference = createTestParticle({-1.15, 2.6, 7.24}, {2.0}, {1.2});
      }
      else if (index == 2)
      {
        particle_reference = createTestParticle({5.12, 4.26, -3.4}, {2.4}, {0.2});
      }

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }
  }

  TEST_F(ParticleContainerBundleTest, SetStateAllContainers)
  {
    std::vector<double> mass{1.1};

    particlecontainerbundle_->SetStateAllContainers(mass, PARTICLEENGINE::Mass);

    PARTICLEENGINE::ParticleContainer* container = nullptr;
    int globalid(0);

    PARTICLEENGINE::ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>{});
    PARTICLEENGINE::ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>{});

    container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase1, PARTICLEENGINE::Owned);

    ASSERT_EQ(container->ParticlesStored(), 3);

    for (int index = 0; index < 3; ++index)
    {
      SCOPED_TRACE("Phase1, Particle " + std::to_string(index));
      if (index == 0)
      {
        particle_reference = createTestParticle({1.20, 0.70, 2.10}, mass, {0.12});
      }
      else if (index == 1)
      {
        particle_reference = createTestParticle({-1.05, 12.6, -8.54}, mass, {12.34});
      }
      else if (index == 2)
      {
        particle_reference = createTestParticle({-5.02, 2.26, -7.4}, mass, {2.9});
      }

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }

    container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase2, PARTICLEENGINE::Owned);

    ASSERT_EQ(container->ParticlesStored(), 3);

    for (int index = 0; index < 3; ++index)
    {
      SCOPED_TRACE("Phase2, Particle " + std::to_string(index));
      if (index == 0)
      {
        particle_reference = createTestParticle({0.24, -1.71, -2.15}, mass, {2.2});
      }
      else if (index == 1)
      {
        particle_reference = createTestParticle({-1.15, 2.6, 7.24}, mass, {1.2});
      }
      else if (index == 2)
      {
        particle_reference = createTestParticle({5.12, 4.26, -3.4}, mass, {0.2});
      }

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }
  }

  TEST_F(ParticleContainerBundleTest, ClearStateAllContainers)
  {
    std::vector<double> mass{0.0};

    particlecontainerbundle_->ClearStateAllContainers(PARTICLEENGINE::Mass);

    PARTICLEENGINE::ParticleContainer* container = nullptr;
    int globalid(0);

    PARTICLEENGINE::ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>{});
    PARTICLEENGINE::ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>{});

    container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase1, PARTICLEENGINE::Owned);

    ASSERT_EQ(container->ParticlesStored(), 3);

    for (int index = 0; index < 3; ++index)
    {
      SCOPED_TRACE("Phase1, Particle " + std::to_string(index));
      if (index == 0)
      {
        particle_reference = createTestParticle({1.20, 0.70, 2.10}, mass, {0.12});
      }
      else if (index == 1)
      {
        particle_reference = createTestParticle({-1.05, 12.6, -8.54}, mass, {12.34});
      }
      else if (index == 2)
      {
        particle_reference = createTestParticle({-5.02, 2.26, -7.4}, mass, {2.9});
      }

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }

    container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase2, PARTICLEENGINE::Owned);

    ASSERT_EQ(container->ParticlesStored(), 3);

    for (int index = 0; index < 3; ++index)
    {
      SCOPED_TRACE("Phase2, Particle " + std::to_string(index));
      if (index == 0)
      {
        particle_reference = createTestParticle({0.24, -1.71, -2.15}, mass, {2.2});
      }
      else if (index == 1)
      {
        particle_reference = createTestParticle({-1.15, 2.6, 7.24}, mass, {1.2});
      }
      else if (index == 2)
      {
        particle_reference = createTestParticle({5.12, 4.26, -3.4}, mass, {0.2});
      }

      container->GetParticle(index, globalid, particle);

      compareParticleStates(particle_reference, particle);
    }
  }

  TEST_F(ParticleContainerBundleTest, CheckAndDecreaseSizeAllContainersOfSpecificStatus)
  {
    PARTICLEENGINE::ParticleContainer* container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase1, PARTICLEENGINE::Owned);

    ASSERT_EQ(container->ParticlesStored(), 3);
    ASSERT_EQ(container->ContainerSize(), 4);

    container->RemoveParticle(0);
    container->RemoveParticle(0);

    particlecontainerbundle_->CheckAndDecreaseSizeAllContainersOfSpecificStatus(
        PARTICLEENGINE::Owned);

    EXPECT_EQ(container->ParticlesStored(), 1);
    EXPECT_EQ(container->ContainerSize(), 2);
  }

  TEST_F(ParticleContainerBundleTest, ClearAllContainersOfSpecificStatus)
  {
    particlecontainerbundle_->ClearAllContainersOfSpecificStatus(PARTICLEENGINE::Ghosted);

    PARTICLEENGINE::ParticleContainer* container = particlecontainerbundle_->GetSpecificContainer(
        PARTICLEENGINE::Phase1, PARTICLEENGINE::Ghosted);

    EXPECT_EQ(container->ParticlesStored(), 0);
  }

  TEST_F(ParticleContainerBundleTest, GetVectorOfParticleObjectsOfAllContainers)
  {
    std::vector<PARTICLEENGINE::ParticleObjShrdPtr> particlesstored;

    particlecontainerbundle_->GetVectorOfParticleObjectsOfAllContainers(particlesstored);

    EXPECT_EQ(particlesstored.size(), 6);
  }
}  // namespace
