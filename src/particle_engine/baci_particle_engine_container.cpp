/*---------------------------------------------------------------------------*/
/*! \file
\brief smart particle container
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_engine_container.H"

#include "baci_utils_exceptions.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEENGINE::ParticleContainer::ParticleContainer()
    : containersize_(0), particlestored_(0), statesvectorsize_(0), globalids_(0, -1)
{
  // empty constructor
}

void PARTICLEENGINE::ParticleContainer::Init()
{
  // nothing to do
}

void PARTICLEENGINE::ParticleContainer::Setup(
    int containersize, const std::set<ParticleState>& stateset)
{
  // set size of particle container (at least one)
  containersize_ = (containersize > 0) ? containersize : 1;

  // set of stored particle states
  storedstates_ = stateset;

  // determine necessary size of vector for states
  statesvectorsize_ = *(--storedstates_.end()) + 1;

  // allocate memory for global ids
  globalids_.resize(containersize_, -1);

  // allocate memory to hold particle states and dimension
  states_.resize(statesvectorsize_);
  statedim_.resize(statesvectorsize_);

  // iterate over states to be stored in container
  for (const auto& state : storedstates_)
  {
    // set particle state dimension for current state
    statedim_[state] = EnumToStateDim(state);

    // allocate memory for current state in particle container
    states_[state].resize(containersize_ * statedim_[state]);
  }
}

void PARTICLEENGINE::ParticleContainer::IncreaseContainerSize()
{
  // size of container is doubled
  containersize_ *= 2;

  // resize vector of global ids
  globalids_.resize(containersize_);

  // iterate over states stored in container
  for (const auto& state : storedstates_)
  {
    // resize vector of current state
    states_[state].resize(containersize_ * statedim_[state]);
  }
}

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
  for (const auto& state : storedstates_)
  {
    // resize vector of current state
    states_[state].resize(containersize_ * statedim_[state]);
  }
}

void PARTICLEENGINE::ParticleContainer::AddParticle(
    int& index, int globalid, const ParticleStates& states)
{
  // increase size of container
  if (particlestored_ == containersize_) IncreaseContainerSize();

  // store global id
  globalids_[particlestored_] = globalid;

  // iterate over states stored in container
  for (const auto& state : storedstates_)
  {
    // state not handed over
    if (states.size() <= state or states[state].empty())
    {
      // initialize to zero
      for (int dim = 0; dim < statedim_[state]; ++dim)
        (states_[state])[particlestored_ * statedim_[state] + dim] = 0.0;
    }
    // state handed over
    else
    {
#ifdef DEBUG
      if (static_cast<int>(states[state].size()) != statedim_[state])
        dserror("can not add particle: dimensions of state '%s' do not match!",
            EnumToStateName(state).c_str());
#endif

      // store state in container
      for (int dim = 0; dim < statedim_[state]; ++dim)
        (states_[state])[particlestored_ * statedim_[state] + dim] = states[state][dim];
    }
  }

  // set index of added particle
  index = particlestored_;

  // increase counter of stored particles
  particlestored_++;
}

void PARTICLEENGINE::ParticleContainer::ReplaceParticle(
    int index, int globalid, const ParticleStates& states)
{
#ifdef DEBUG
  if (index < 0 or index > (particlestored_ - 1))
    dserror("can not replace particle as index %d out of bounds!", index);
#endif

  // replace global id in container
  if (globalid >= 0) globalids_[index] = globalid;

  // iterate over states stored in container
  for (const auto& state : storedstates_)
  {
    // state not handed over
    if (states.size() <= state or states[state].empty())
    {
      // leave state untouched
    }
    // state handed over
    else
    {
#ifdef DEBUG
      if (static_cast<int>(states[state].size()) != statedim_[state])
        dserror("can not replace particle: dimensions of state '%s' do not match!",
            EnumToStateName(state).c_str());
#endif

      // replace state in container
      for (int dim = 0; dim < statedim_[state]; ++dim)
        (states_[state])[index * statedim_[state] + dim] = states[state][dim];
    }
  }
}

void PARTICLEENGINE::ParticleContainer::GetParticle(
    int index, int& globalid, ParticleStates& states) const
{
#ifdef DEBUG
  if (index < 0 or index > (particlestored_ - 1))
    dserror("can not return particle as index %d out of bounds!", index);
#endif

  // get global id from container
  globalid = globalids_[index];

  // allocate memory to hold particle states
  states.assign(statesvectorsize_, std::vector<double>(0));

  // iterate over states stored in container
  for (const auto& state : storedstates_)
  {
    // get pointer to particle state
    const double* state_ptr = &((states_[state])[index * statedim_[state]]);

    // fill particle state
    states[state].assign(state_ptr, state_ptr + statedim_[state]);
  }
}

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
  for (const auto& state : storedstates_)
  {
    // overwrite state in container
    for (int dim = 0; dim < statedim_[state]; ++dim)
      (states_[state])[index * statedim_[state] + dim] =
          (states_[state])[particlestored_ * statedim_[state] + dim];
  }
}

double PARTICLEENGINE::ParticleContainer::GetMinValueOfState(ParticleState state) const
{
#ifdef DEBUG
  if (not storedstates_.count(state))
    dserror("particle state '%s' not stored in container!", EnumToStateName(state).c_str());
#endif

  if (particlestored_ <= 0) return 0.0;

  double min = (states_[state])[0];

  for (int i = 0; i < (particlestored_ * statedim_[state]); ++i)
    min = std::min(min, states_[state][i]);

  return min;
}

double PARTICLEENGINE::ParticleContainer::GetMaxValueOfState(ParticleState state) const
{
#ifdef DEBUG
  if (not storedstates_.count(state))
    dserror("particle state '%s' not stored in container!", EnumToStateName(state).c_str());
#endif

  if (particlestored_ <= 0) return 0.0;

  double max = (states_[state])[0];

  for (int i = 0; i < (particlestored_ * statedim_[state]); ++i)
    max = std::max(max, states_[state][i]);

  return max;
}
