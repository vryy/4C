/*---------------------------------------------------------------------------*/
/*!
\file particle_container_bundle.cpp

\brief class holding all particle containers

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_container_bundle.H"

#include "particle_container.H"
#include "particle_object.H"

#include "../drt_lib/drt_dserror.H"

#include <iostream>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEENGINE::ParticleContainerBundle::ParticleContainerBundle(const int myrank) : myrank_(myrank)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init particle container bundle                             sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup particle container bundle                            sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::Setup(
    const std::map<TypeEnum, std::set<StateEnum>>& particlestatestotypes)
{
  ParticleContainerShrdPtr container;

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

/*---------------------------------------------------------------------------*
 | scale state of owned particles in specific containers      sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::ScaleStateSpecificContainer(
    double fac, StateEnum stateEnum, TypeEnum typeEnum) const
{
  ((containers_[typeEnum])[PARTICLEENGINE::Owned])->ScaleState(fac, stateEnum);
}

/*---------------------------------------------------------------------------*
 | update state of owned particles in specific container      sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::UpdateStateSpecificContainer(
    double facA, StateEnum stateEnumA, double facB, StateEnum stateEnumB, TypeEnum typeEnum) const
{
  ((containers_[typeEnum])[PARTICLEENGINE::Owned])->UpdateState(facA, stateEnumA, facB, stateEnumB);
}

/*---------------------------------------------------------------------------*
 | set state of owned particles in specific container         sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::SetStateSpecificContainer(
    std::vector<double> val, StateEnum stateEnum, TypeEnum typeEnum) const
{
  ((containers_[typeEnum])[PARTICLEENGINE::Owned])->SetState(val, stateEnum);
}

/*---------------------------------------------------------------------------*
 | clear state of owned particles in specific container       sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::ClearStateSpecificContainer(
    StateEnum stateEnum, TypeEnum typeEnum) const
{
  ((containers_[typeEnum])[PARTICLEENGINE::Owned])->ClearState(stateEnum);
}

/*---------------------------------------------------------------------------*
 | scale state of owned particles in all containers           sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::ScaleStateAllContainers(
    double fac, StateEnum stateEnum) const
{
  // iterate over particle types
  for (auto& typeEnum : storedtypes_)
    ((containers_[typeEnum])[PARTICLEENGINE::Owned])->ScaleState(fac, stateEnum);
}

/*---------------------------------------------------------------------------*
 | update state of owned particles in all containers          sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::UpdateStateAllContainers(
    double facA, StateEnum stateEnumA, double facB, StateEnum stateEnumB) const
{
  // iterate over particle types
  for (auto& typeEnum : storedtypes_)
    ((containers_[typeEnum])[PARTICLEENGINE::Owned])
        ->UpdateState(facA, stateEnumA, facB, stateEnumB);
}

/*---------------------------------------------------------------------------*
 | set state of owned particles in all containers             sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::SetStateAllContainers(
    std::vector<double> val, StateEnum stateEnum) const
{
  // iterate over particle types
  for (auto& typeEnum : storedtypes_)
    ((containers_[typeEnum])[PARTICLEENGINE::Owned])->SetState(val, stateEnum);
}

/*---------------------------------------------------------------------------*
 | clear state of owned particles in all containers           sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::ClearStateAllContainers(StateEnum stateEnum) const
{
  // iterate over particle types
  for (auto& typeEnum : storedtypes_)
    ((containers_[typeEnum])[PARTICLEENGINE::Owned])->ClearState(stateEnum);
}

/*---------------------------------------------------------------------------*
 | get specific particle container                            sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEENGINE::ParticleContainerShrdPtr
PARTICLEENGINE::ParticleContainerBundle::GetSpecificContainer(
    TypeEnum typeEnum, StatusEnum statusEnum) const
{
  return ((containers_[typeEnum])[statusEnum]);
}

/*---------------------------------------------------------------------------*
 | clear all containers of specific status                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::ClearAllContainersOfSpecificStatus(
    StatusEnum statusEnum) const
{
  // iterate over particle types
  for (auto& typeEnum : storedtypes_) ((containers_[typeEnum])[statusEnum])->ClearContainer();
}

/*---------------------------------------------------------------------------*
 | pack all particle containers                               sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::PackParticleContainerBundle(
    Teuchos::RCP<std::vector<char>>& particlebuffer) const
{
  // iterate over particle types
  for (auto& typeEnum : storedtypes_)
  {
    // get container of owned particles
    PARTICLEENGINE::ParticleContainerShrdPtr container =
        ((containers_[typeEnum])[PARTICLEENGINE::Owned]);

    // loop over particles in container
    for (int index = 0; index < container->ParticlesStored(); ++index)
    {
      int globalid(0);
      ParticleStates particleStates;
      container->GetParticle(index, globalid, particleStates);

      ParticleObjShrdPtr particleobject = std::make_shared<PARTICLEENGINE::ParticleObject>();
      particleobject->Init(typeEnum, globalid, particleStates);

      // pack data for writing
      DRT::PackBuffer data;
      particleobject->Pack(data);
      data.StartPacking();
      particleobject->Pack(data);
      particlebuffer->insert(particlebuffer->end(), data().begin(), data().end());
    }
  }
}

/*---------------------------------------------------------------------------*
 | get vector of particle objects of all containers           sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::GetVectorOfParticleObjectsOfAllContainers(
    std::vector<ParticleObjShrdPtr>& particlesstored) const
{
  // iterate over particle types
  for (auto& typeEnum : storedtypes_)
  {
    // get container of owned particles
    PARTICLEENGINE::ParticleContainerShrdPtr container =
        ((containers_[typeEnum])[PARTICLEENGINE::Owned]);

    // loop over particles in container
    for (int index = 0; index < container->ParticlesStored(); ++index)
    {
      int globalid(0);
      ParticleStates particleStates;
      container->GetParticle(index, globalid, particleStates);

      ParticleObjShrdPtr particleobject = std::make_shared<PARTICLEENGINE::ParticleObject>();
      particleobject->Init(typeEnum, globalid, particleStates);

      particlesstored.push_back(particleobject);
    }
  }
}
