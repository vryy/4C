/*---------------------------------------------------------------------------*/
/*!
\file particle_container.cpp

\brief smart particle container class

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_container.H"

#include "../drt_lib/drt_dserror.H"

#include <iostream>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEENGINE::ParticleContainer::ParticleContainer()
    : containersize_(0), particlestored_(0), globalids_(0, -1)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init particle container                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainer::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup particle container                                   sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainer::Setup(
    int numContainerSize, const std::set<StateEnum>& stateEnumSet)
{
  // set size of particle container (at least one)
  containersize_ = (numContainerSize > 0) ? numContainerSize : 1;

  // set of stored particle states
  storedstates_ = stateEnumSet;

  // allocate memory for global ids
  globalids_.resize(containersize_, -1);

  // determine necessary size of vector for states
  const int statesvectorsize = *(--storedstates_.end()) + 1;

  // allocate memory to hold particle states and dimension
  states_.resize(statesvectorsize);
  statedim_.resize(statesvectorsize);

  // iterate over states to be stored in container
  for (auto& stateEnum : storedstates_)
  {
    // set particle state dimension for current state
    statedim_[stateEnum] = EnumToStateDim(stateEnum);

    // allocate memory for current state in particle container
    states_[stateEnum].resize(containersize_ * statedim_[stateEnum]);
  }
}

/*---------------------------------------------------------------------------*
 | increase the container size                                sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainer::IncreaseContainerSize()
{
  // size of container is doubled
  containersize_ *= 2;

  // resize vector of global ids
  globalids_.resize(containersize_);

  // iterate over states stored in container
  for (auto& stateEnum : storedstates_)
  {
    // resize vector of current state
    states_[stateEnum].resize(containersize_ * statedim_[stateEnum]);
  }
}

/*---------------------------------------------------------------------------*
 | decrease the container size                                sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainer::DecreaseContainerSize()
{
  // size of container is halved
  int newsize = static_cast<int>(0.5 * containersize_);

  // set size of particle container (at least one)
  containersize_ = (newsize > 0) ? newsize : 1;

  // safety check
  if (particlestored_ > containersize_)
    dserror(
        "decreasing size of container not possible: particles stored %d > new container size %d!",
        particlestored_, containersize_);

  // resize vector of global ids
  globalids_.resize(containersize_);

  // iterate over states stored in container
  for (auto& stateEnum : storedstates_)
  {
    // resize vector of current state
    states_[stateEnum].resize(containersize_ * statedim_[stateEnum]);
  }
}

/*---------------------------------------------------------------------------*
 | add particle to particle container and get index           sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainer::AddParticle(
    int& index, int globalid, const ParticleStates& particle)
{
  // increase size of container
  if (particlestored_ == containersize_) IncreaseContainerSize();

  // store global id
  globalids_[particlestored_] = globalid;

  // iterate over states stored in container
  for (auto& stateEnum : storedstates_)
  {
    auto particleStateIt = particle.find(stateEnum);
    // state not handed over
    if (particleStateIt == particle.end())
    {
      // initialize to zero
      for (int dim = 0; dim < statedim_[stateEnum]; ++dim)
        (states_[stateEnum])[particlestored_ * statedim_[stateEnum] + dim] = 0.0;
    }
    // state handed over
    else
    {
      const std::vector<double>& particleState = particleStateIt->second;

      // check dimensions
      if (static_cast<int>(particleState.size()) != statedim_[stateEnum])
        dserror("Cannot add particle: dimensions of state '%s' do not match!",
            PARTICLEENGINE::EnumToStateName(stateEnum).c_str());

      // store state in container
      for (int dim = 0; dim < statedim_[stateEnum]; ++dim)
        (states_[stateEnum])[particlestored_ * statedim_[stateEnum] + dim] = particleState[dim];
    }
  }

  // set index of added particle
  index = particlestored_;

  // increase counter of stored particles
  particlestored_++;
}

/*---------------------------------------------------------------------------*
 | replace particle in particle container at given index      sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainer::ReplaceParticle(
    int index, int globalid, const ParticleStates& particle)
{
  if (index < 0 or index > (particlestored_ - 1))
    dserror("can not replace particle as index %d out of bounds!", index);

  // replace global id in container
  if (globalid >= 0) globalids_[index] = globalid;

  // iterate over states stored in container
  for (auto& stateEnum : storedstates_)
  {
    auto particleStateIt = particle.find(stateEnum);
    // state not handed over
    if (particleStateIt == particle.end())
    {
      // leave state untouched
    }
    // state handed over
    else
    {
      const std::vector<double>& particleState = particleStateIt->second;

      // check dimensions
      if (static_cast<int>(particleState.size()) != statedim_[stateEnum])
        dserror("Cannot add particle: dimensions of state '%s' do not match!",
            PARTICLEENGINE::EnumToStateName(stateEnum).c_str());

      // replace state in container
      for (int dim = 0; dim < statedim_[stateEnum]; ++dim)
        (states_[stateEnum])[index * statedim_[stateEnum] + dim] = particleState[dim];
    }
  }
}

/*---------------------------------------------------------------------------*
 | get particle at index from particle container              sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainer::GetParticle(
    int index, int& globalid, ParticleStates& particle) const
{
  if (index < 0 or index > (particlestored_ - 1))
    dserror("can not return particle as index %d out of bounds!", index);

  // get global id from container
  globalid = globalids_[index];

  // clear particle
  particle.clear();

  // iterate over states stored in container
  for (auto& stateEnum : storedstates_)
  {
    std::vector<double> particleState(statedim_[stateEnum]);

    // store current state in vector
    for (int dim = 0; dim < statedim_[stateEnum]; ++dim)
      particleState[dim] = (states_[stateEnum])[index * statedim_[stateEnum] + dim];

    // get state from container
    particle.insert(std::make_pair(stateEnum, particleState));
  }
}

/*---------------------------------------------------------------------------*
 | remove particle from particle container                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainer::RemoveParticle(int index)
{
  if (index < 0 or index > (particlestored_ - 1))
    dserror("can not remove particle as index %d out of bounds!", index);

  // decrease counter of stored particles
  --particlestored_;

  // overwrite global id in container
  globalids_[index] = globalids_[particlestored_];

  // iterate over states stored in container
  for (auto& stateEnum : storedstates_)
  {
    // overwrite state in container
    for (int dim = 0; dim < statedim_[stateEnum]; ++dim)
      (states_[stateEnum])[index * statedim_[stateEnum] + dim] =
          (states_[stateEnum])[particlestored_ * statedim_[stateEnum] + dim];
  }

  // decrease size of container
  if (particlestored_ < 0.45 * containersize_) DecreaseContainerSize();
}

/*---------------------------------------------------------------------------*
 | return pointer to state of a particle at index             sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
double* PARTICLEENGINE::ParticleContainer::GetPtrToParticleState(StateEnum stateEnum, int index)
{
  if (index < 0 or index > (particlestored_ - 1))
    dserror("can not return pointer to state of particle as index %d out of bounds!", index);

  return &((states_[stateEnum])[index * statedim_[stateEnum]]);
}

/*---------------------------------------------------------------------------*
 | return pointer to global id of a particle at index         sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
int* PARTICLEENGINE::ParticleContainer::GetPtrToParticleGlobalID(int index)
{
  if (index < 0 or index > (particlestored_ - 1))
    dserror("can not return pointer to global id of particle as index %d out of bounds!", index);

  return &(globalids_[index]);
}

/*---------------------------------------------------------------------------*
 | get state of a particle at index                           sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
std::vector<double> PARTICLEENGINE::ParticleContainer::GetParticleState(
    StateEnum stateEnum, int index) const
{
  if (index < 0 or index > (particlestored_ - 1))
    dserror("can not return state of particle as index %d out of bounds!", index);

  std::vector<double> particleState(statedim_[stateEnum]);
  for (int dim = 0; dim < statedim_[stateEnum]; ++dim)
    particleState[dim] = (states_[stateEnum])[index * statedim_[stateEnum] + dim];

  return particleState;
}

/*---------------------------------------------------------------------------*
 | scale state of particles                                   sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainer::ScaleState(double fac, StateEnum stateEnum)
{
  for (int i = 0; i < (particlestored_ * statedim_[stateEnum]); ++i) (states_[stateEnum])[i] *= fac;
}

/*---------------------------------------------------------------------------*
 | scale and add states to update state of particles          sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainer::UpdateState(
    double facA, StateEnum stateEnumA, double facB, StateEnum stateEnumB)
{
  if (statedim_[stateEnumA] != statedim_[stateEnumB]) dserror("dimensions of states do not match!");

  for (int i = 0; i < (particlestored_ * statedim_[stateEnumA]); ++i)
    (states_[stateEnumA])[i] = facA * (states_[stateEnumA])[i] + facB * (states_[stateEnumB])[i];
}

/*---------------------------------------------------------------------------*
 | set state to particles                                     sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainer::SetState(std::vector<double> val, StateEnum stateEnum)
{
  if (statedim_[stateEnum] != static_cast<int>(val.size()))
    dserror("dimensions of states do not match!");

  for (int i = 0; i < particlestored_; ++i)
    for (int dim = 0; dim < statedim_[stateEnum]; ++dim)
      (states_[stateEnum])[i * statedim_[stateEnum] + dim] = val[dim];
}

/*---------------------------------------------------------------------------*
 | clear state of particles                                   sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainer::ClearState(StateEnum stateEnum)
{
  for (int i = 0; i < (particlestored_ * statedim_[stateEnum]); ++i) (states_[stateEnum])[i] = 0.0;
}

/*---------------------------------------------------------------------------*
 | get minimum stored value of state                          sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEENGINE::ParticleContainer::GetMinValueOfState(StateEnum stateEnum) const
{
  if (particlestored_ <= 0) return 0.0;

  double min = (states_[stateEnum])[0];

  for (int i = 0; i < (particlestored_ * statedim_[stateEnum]); ++i)
    if ((states_[stateEnum])[i] < min) min = (states_[stateEnum])[i];

  return min;
}

/*---------------------------------------------------------------------------*
 | get maximum stored value of state                          sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEENGINE::ParticleContainer::GetMaxValueOfState(StateEnum stateEnum) const
{
  if (particlestored_ <= 0) return 0.0;

  double max = (states_[stateEnum])[0];

  for (int i = 0; i < (particlestored_ * statedim_[stateEnum]); ++i)
    if ((states_[stateEnum])[i] > max) max = (states_[stateEnum])[i];

  return max;
}
