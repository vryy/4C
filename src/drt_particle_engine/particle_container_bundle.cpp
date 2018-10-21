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
  std::map<StatusEnum, ParticleContainerShrdPtr> statusMap;

  // iterate over particle types
  for (auto& typeEnum : particlestatestotypes)
  {
    // clear particle status map
    statusMap.clear();

    // get set of particle state enums of current particle type (equal for owned and ghosted
    // particles)
    const std::set<StateEnum>& stateEnumSet = typeEnum.second;

    // initial size of particle container
    int initialsize = 1;

    // create and init container of owned particles
    container = std::make_shared<PARTICLEENGINE::ParticleContainer>();
    container->Init();
    // setup container of owned particles
    container->Setup(initialsize, stateEnumSet);
    statusMap.insert(std::make_pair(PARTICLEENGINE::Owned, container));

    // create and init container of ghosted particles
    container = std::make_shared<PARTICLEENGINE::ParticleContainer>();
    container->Init();
    // setup container of ghosted particles
    container->Setup(initialsize, stateEnumSet);
    statusMap.insert(std::make_pair(PARTICLEENGINE::Ghosted, container));

    // insert into map that holds particle containers of all types and statuses
    containers_.insert(std::make_pair(typeEnum.first, statusMap));
  }
}

/*---------------------------------------------------------------------------*
 | scale state of owned particles in specific containers      sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::ScaleStateSpecificContainer(
    double fac, StateEnum stateEnum, TypeEnum particleType) const
{
  // get container of owned particles of current particle type
  ParticleContainerShrdPtr container = GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

  container->ScaleState(fac, stateEnum);
}

/*---------------------------------------------------------------------------*
 | update state of owned particles in specific container      sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::UpdateStateSpecificContainer(double facA,
    StateEnum stateEnumA, double facB, StateEnum stateEnumB, TypeEnum particleType) const
{
  // get container of owned particles of current particle type
  ParticleContainerShrdPtr container = GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

  container->UpdateState(facA, stateEnumA, facB, stateEnumB);
}

/*---------------------------------------------------------------------------*
 | set state of owned particles in specific container         sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::SetStateSpecificContainer(
    std::vector<double> val, StateEnum stateEnum, TypeEnum particleType) const
{
  // get container of owned particles of current particle type
  ParticleContainerShrdPtr container = GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

  container->SetState(val, stateEnum);
}

/*---------------------------------------------------------------------------*
 | clear state of owned particles in specific container       sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::ClearStateSpecificContainer(
    StateEnum stateEnum, TypeEnum particleType) const
{
  // get container of owned particles of current particle type
  ParticleContainerShrdPtr container = GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

  container->ClearState(stateEnum);
}

/*---------------------------------------------------------------------------*
 | scale state of owned particles in all containers           sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::ScaleStateAllContainers(
    double fac, StateEnum stateEnum) const
{
  // iterate over particle types
  for (auto& typeIt : containers_)
  {
    // iterate over particle statuses
    for (auto& statusIt : typeIt.second)
    {
      // check status of particles is owned
      if (statusIt.first == PARTICLEENGINE::Owned) (statusIt.second)->ScaleState(fac, stateEnum);
    }
  }
}

/*---------------------------------------------------------------------------*
 | update state of owned particles in all containers          sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::UpdateStateAllContainers(
    double facA, StateEnum stateEnumA, double facB, StateEnum stateEnumB) const
{
  // iterate over particle types
  for (auto& typeIt : containers_)
  {
    // iterate over particle statuses
    for (auto& statusIt : typeIt.second)
    {
      // check status of particles is owned
      if (statusIt.first == PARTICLEENGINE::Owned)
        (statusIt.second)->UpdateState(facA, stateEnumA, facB, stateEnumB);
    }
  }
}

/*---------------------------------------------------------------------------*
 | set state of owned particles in all containers             sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::SetStateAllContainers(
    std::vector<double> val, StateEnum stateEnum) const
{
  // iterate over particle types
  for (auto& typeIt : containers_)
  {
    // iterate over particle statuses
    for (auto& statusIt : typeIt.second)
    {
      // check status of particles is owned
      if (statusIt.first == PARTICLEENGINE::Owned) (statusIt.second)->SetState(val, stateEnum);
    }
  }
}

/*---------------------------------------------------------------------------*
 | clear state of owned particles in all containers           sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::ClearStateAllContainers(StateEnum stateEnum) const
{
  // iterate over particle types
  for (auto& typeIt : containers_)
  {
    // iterate over particle statuses
    for (auto& statusIt : typeIt.second)
    {
      // check status of particles is owned
      if (statusIt.first == PARTICLEENGINE::Owned) (statusIt.second)->ClearState(stateEnum);
    }
  }
}

/*---------------------------------------------------------------------------*
 | get specific particle container                            sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEENGINE::ParticleContainerShrdPtr
PARTICLEENGINE::ParticleContainerBundle::GetSpecificContainer(
    TypeEnum particleType, StatusEnum particleStatus) const
{
  auto typeIt = containers_.find(particleType);
  if (typeIt == containers_.end())
    dserror("particle type '%s' not found!", PARTICLEENGINE::EnumToTypeName(particleType).c_str());

  auto statusIt = (typeIt->second).find(particleStatus);
  if (statusIt == (typeIt->second).end())
    dserror("particle status '%s' not found!",
        PARTICLEENGINE::EnumToStatusName(particleStatus).c_str());

  return statusIt->second;
}

/*---------------------------------------------------------------------------*
 | clear all containers of specific status                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::ClearAllContainersOfSpecificStatus(
    StatusEnum particleStatus) const
{
  // iterate over particle types
  for (auto& typeIt : containers_)
  {
    // iterate over particle statuses
    for (auto& statusIt : typeIt.second)
    {
      // check status of current particle container
      if (statusIt.first == particleStatus) (statusIt.second)->ClearContainer();
    }
  }
}

/*---------------------------------------------------------------------------*
 | pack all particle containers                               sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainerBundle::PackParticleContainerBundle(
    Teuchos::RCP<std::vector<char>>& particlebuffer) const
{
  // iterate over particle types
  for (auto& typeIt : containers_)
  {
    // get type of particles
    TypeEnum particleType = typeIt.first;

    // get container of owned particles
    auto statusIt = (typeIt.second).find(PARTICLEENGINE::Owned);
    if (statusIt == (typeIt.second).end())
      dserror("particle status '%s' not found!",
          PARTICLEENGINE::EnumToStatusName(PARTICLEENGINE::Owned).c_str());
    ParticleContainerShrdPtr container = statusIt->second;

    // loop over particles in container
    for (int index = 0; index < container->ParticlesStored(); ++index)
    {
      int globalid(0);
      ParticleStates particleStates;
      container->GetParticle(index, globalid, particleStates);

      ParticleObjShrdPtr particleobject = std::make_shared<PARTICLEENGINE::ParticleObject>();
      particleobject->Init(particleType, globalid, particleStates);

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
  for (auto& typeIt : containers_)
  {
    // get type of particles
    TypeEnum particleType = typeIt.first;

    // get container of owned particles
    auto statusIt = (typeIt.second).find(PARTICLEENGINE::Owned);
    if (statusIt == (typeIt.second).end())
      dserror("particle status '%s' not found!",
          PARTICLEENGINE::EnumToStatusName(PARTICLEENGINE::Owned).c_str());
    ParticleContainerShrdPtr container = statusIt->second;

    // loop over particles in container
    for (int index = 0; index < container->ParticlesStored(); ++index)
    {
      int globalid(0);
      ParticleStates particleStates;
      container->GetParticle(index, globalid, particleStates);

      ParticleObjShrdPtr particleobject = std::make_shared<PARTICLEENGINE::ParticleObject>();
      particleobject->Init(particleType, globalid, particleStates);

      particlesstored.push_back(particleobject);
    }
  }
}
