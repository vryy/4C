/*---------------------------------------------------------------------------*/
/*! \file
\brief manage bundle of particle containers

\level 1

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_container_bundle.H"

#include "particle_object.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEENGINE::ParticleContainerBundle::ParticleContainerBundle()
{
  // empty constructor
}

void PARTICLEENGINE::ParticleContainerBundle::Init()
{
  // nothing to do
}

void PARTICLEENGINE::ParticleContainerBundle::Setup(
    const std::map<TypeEnum, std::set<StateEnum>>& particlestatestotypes)
{
  std::shared_ptr<PARTICLEENGINE::ParticleContainer> container;

  // determine necessary size of vector for particle types
  const int typevectorsize = ((--particlestatestotypes.end())->first) + 1;

  // allocate memory to hold particle types
  containers_.resize(typevectorsize);

  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // get particle type
    TypeEnum typeEnum = typeIt.first;

    // insert particle type into set of stored containers
    storedtypes_.insert(typeEnum);

    // allocate memory for container of owned and ghosted particles
    (containers_[typeEnum]).resize(2);

    // set of particle state enums of current particle type (equal for owned and ghosted particles)
    const std::set<StateEnum>& stateEnumSet = typeIt.second;

    // initial size of particle container
    int initialsize = 1;

    // create and init container of owned particles
    container = std::make_shared<PARTICLEENGINE::ParticleContainer>();
    container->Init();
    // setup container of owned particles
    container->Setup(initialsize, stateEnumSet);
    // set container of owned particles
    (containers_[typeEnum])[PARTICLEENGINE::Owned] = container;

    // create and init container of ghosted particles
    container = std::make_shared<PARTICLEENGINE::ParticleContainer>();
    container->Init();
    // setup container of ghosted particles
    container->Setup(initialsize, stateEnumSet);
    // set container of ghosted particles
    (containers_[typeEnum])[PARTICLEENGINE::Ghosted] = container;
  }
}

void PARTICLEENGINE::ParticleContainerBundle::GetPackedParticleObjectsOfAllContainers(
    std::shared_ptr<std::vector<char>>& particlebuffer) const
{
  // iterate over particle types
  for (auto& typeEnum : storedtypes_)
  {
    // get container of owned particles
    PARTICLEENGINE::ParticleContainer* container =
        (containers_[typeEnum])[PARTICLEENGINE::Owned].get();

    // loop over particles in container
    for (int index = 0; index < container->ParticlesStored(); ++index)
    {
      int globalid(0);
      ParticleStates particleStates;
      container->GetParticle(index, globalid, particleStates);

      ParticleObjShrdPtr particleobject =
          std::make_shared<PARTICLEENGINE::ParticleObject>(typeEnum, globalid, particleStates);

      // pack data for writing
      DRT::PackBuffer data;
      particleobject->Pack(data);
      data.StartPacking();
      particleobject->Pack(data);
      particlebuffer->insert(particlebuffer->end(), data().begin(), data().end());
    }
  }
}

void PARTICLEENGINE::ParticleContainerBundle::GetVectorOfParticleObjectsOfAllContainers(
    std::vector<ParticleObjShrdPtr>& particlesstored) const
{
  // iterate over particle types
  for (auto& typeEnum : storedtypes_)
  {
    // get container of owned particles
    PARTICLEENGINE::ParticleContainer* container =
        (containers_[typeEnum])[PARTICLEENGINE::Owned].get();

    // loop over particles in container
    for (int index = 0; index < container->ParticlesStored(); ++index)
    {
      int globalid(0);
      ParticleStates particleStates;
      container->GetParticle(index, globalid, particleStates);

      particlesstored.emplace_back(
          std::make_shared<PARTICLEENGINE::ParticleObject>(typeEnum, globalid, particleStates));
    }
  }
}
