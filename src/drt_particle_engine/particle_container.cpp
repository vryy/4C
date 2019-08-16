/*---------------------------------------------------------------------------*/
/*!
\brief smart particle container class

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_container.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEENGINE::ParticleContainer::ParticleContainer()
    : containersize_(0), particlestored_(0), statesvectorsize_(0), globalids_(0, -1)
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

  // determine necessary size of vector for states
  statesvectorsize_ = *(--storedstates_.end()) + 1;

  // allocate memory for global ids
  globalids_.resize(containersize_, -1);

  // allocate memory to hold particle states and dimension
  states_.resize(statesvectorsize_);
  statedim_.resize(statesvectorsize_);

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

#ifdef DEBUG
  if (particlestored_ > containersize_)
    dserror(
        "decreasing size of container not possible: particles stored %d > new container size %d!",
        particlestored_, containersize_);
#endif

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
    // state not handed over
    if (particle.size() <= stateEnum or particle[stateEnum].size() == 0)
    {
      // initialize to zero
      for (int dim = 0; dim < statedim_[stateEnum]; ++dim)
        (states_[stateEnum])[particlestored_ * statedim_[stateEnum] + dim] = 0.0;
    }
    // state handed over
    else
    {
#ifdef DEBUG
      if (static_cast<int>(particle[stateEnum].size()) != statedim_[stateEnum])
        dserror("can not add particle: dimensions of state '%s' do not match!",
            PARTICLEENGINE::EnumToStateName(stateEnum).c_str());
#endif

      // store state in container
      for (int dim = 0; dim < statedim_[stateEnum]; ++dim)
        (states_[stateEnum])[particlestored_ * statedim_[stateEnum] + dim] =
            particle[stateEnum][dim];
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
#ifdef DEBUG
  if (index < 0 or index > (particlestored_ - 1))
    dserror("can not replace particle as index %d out of bounds!", index);
#endif

  // replace global id in container
  if (globalid >= 0) globalids_[index] = globalid;

  // iterate over states stored in container
  for (auto& stateEnum : storedstates_)
  {
    // state not handed over
    if (particle.size() <= stateEnum or particle[stateEnum].size() == 0)
    {
      // leave state untouched
    }
    // state handed over
    else
    {
#ifdef DEBUG
      if (static_cast<int>(particle[stateEnum].size()) != statedim_[stateEnum])
        dserror("can not replace particle: dimensions of state '%s' do not match!",
            PARTICLEENGINE::EnumToStateName(stateEnum).c_str());
#endif

      // replace state in container
      for (int dim = 0; dim < statedim_[stateEnum]; ++dim)
        (states_[stateEnum])[index * statedim_[stateEnum] + dim] = particle[stateEnum][dim];
    }
  }
}

/*---------------------------------------------------------------------------*
 | get particle at index from particle container              sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainer::GetParticle(
    int index, int& globalid, ParticleStates& particle) const
{
#ifdef DEBUG
  if (index < 0 or index > (particlestored_ - 1))
    dserror("can not return particle as index %d out of bounds!", index);
#endif

  // get global id from container
  globalid = globalids_[index];

  // allocate memory to hold particle states
  particle.assign(statesvectorsize_, std::vector<double>(0));

  // iterate over states stored in container
  for (auto& stateEnum : storedstates_)
  {
    // get pointer to particle state
    const double* state_ptr = &((states_[stateEnum])[index * statedim_[stateEnum]]);

    // fill particle state
    particle[stateEnum].assign(state_ptr, state_ptr + statedim_[stateEnum]);
  }
}

/*---------------------------------------------------------------------------*
 | remove particle from particle container                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleContainer::RemoveParticle(int index)
{
#ifdef DEBUG
  if (index < 0 or index > (particlestored_ - 1))
    dserror("can not remove particle as index %d out of bounds!", index);
#endif

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
}

/*---------------------------------------------------------------------------*
 | get minimum stored value of state                          sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEENGINE::ParticleContainer::GetMinValueOfState(StateEnum stateEnum) const
{
#ifdef DEBUG
  if (not storedstates_.count(stateEnum))
    dserror("particle state '%s' not stored in container!",
        PARTICLEENGINE::EnumToStateName(stateEnum).c_str());
#endif

  if (particlestored_ <= 0) return 0.0;

  double min = (states_[stateEnum])[0];

  for (int i = 0; i < (particlestored_ * statedim_[stateEnum]); ++i)
    min = std::min(min, states_[stateEnum][i]);

  return min;
}

/*---------------------------------------------------------------------------*
 | get maximum stored value of state                          sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEENGINE::ParticleContainer::GetMaxValueOfState(StateEnum stateEnum) const
{
#ifdef DEBUG
  if (not storedstates_.count(stateEnum))
    dserror("particle state '%s' not stored in container!",
        PARTICLEENGINE::EnumToStateName(stateEnum).c_str());
#endif

  if (particlestored_ <= 0) return 0.0;

  double max = (states_[stateEnum])[0];

  for (int i = 0; i < (particlestored_ * statedim_[stateEnum]); ++i)
    max = std::max(max, states_[stateEnum][i]);

  return max;
}
