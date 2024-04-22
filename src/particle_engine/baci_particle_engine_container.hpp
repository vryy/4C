/*---------------------------------------------------------------------------*/
/*! \file
\brief smart particle container
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_ENGINE_CONTAINER_HPP
#define FOUR_C_PARTICLE_ENGINE_CONTAINER_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_particle_engine_enums.hpp"
#include "baci_particle_engine_typedefs.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  /*!
   * \brief smart particle container class
   *
   * A particle container class that allows for dynamic adding and removing of particles while
   * avoiding both expensive memory reallocations and memory fragmentation. Fast access to particle
   * states is provided.
   *
   * \author Sebastian Fuchs \date 03/2018
   */
  class ParticleContainer final
  {
   public:
    //! constructor
    explicit ParticleContainer();

    /*!
     * \brief init particle container
     *
     * \author Sebastian Fuchs \date 03/2018
     */
    void Init();

    /*!
     * \brief setup particle container
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] containersize size of particle container
     * \param[in] stateset      set of particle states to be stored
     */
    void Setup(int containersize, const std::set<ParticleState>& stateset);

    //! \name manipulate container size
    //! @{

    /*!
     * \brief increase the container size
     *
     * The size of the particle container is doubled to reduce memory (re-)allocation costs.
     *
     * \author Sebastian Fuchs \date 04/2018
     */
    void IncreaseContainerSize();

    /*!
     * \brief decrease the container size
     *
     * The size of the particle container is halved.
     *
     * \author Sebastian Fuchs \date 04/2018
     */
    void DecreaseContainerSize();

    /*!
     * \brief check and decrease the container size
     *
     * \author Sebastian Fuchs \date 04/2018
     */
    inline void CheckAndDecreaseContainerSize()
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (particlestored_ > containersize_)
        FOUR_C_THROW(
            "checking size of container not possible: particles stored %d > new container size %d!",
            particlestored_, containersize_);
#endif

      if (particlestored_ < 0.45 * containersize_) DecreaseContainerSize();
    };

    //! @}

    //! \name manage particles stored in container
    //! @{

    /*!
     * \brief clear particle container
     *
     * Clear the particle container by resetting the number of particles stored in the container.
     * The container size remains untouched. The global id and the particle states are not cleared.
     *
     * \author Sebastian Fuchs \date 04/2018
     */
    inline void ClearContainer() { particlestored_ = 0; };

    /*!
     * \brief add particle to particle container and get index
     *
     * Add a particle to the particle container and return the index of the particle in the
     * container. The size of the particle container is increased if necessary. A state that is not
     * handed over is initialized to zero.
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[out] index    index of particle in container
     * \param[in]  globalid global id of particle
     * \param[in]  states   states of particle
     */
    void AddParticle(int& index, int globalid, const ParticleStates& states);

    /*!
     * \brief replace particle in particle container at given index
     *
     * Replace a particle at the given index in the particle container, meaning the global id and
     * the particle states are overwritten. A negative global id is ignored. A state that is not
     * handed over is untouched.
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] index    index of particle in container
     * \param[in] globalid global id of particle
     * \param[in] states   states of particle
     */
    void ReplaceParticle(int index, int globalid, const ParticleStates& states);

    /*!
     * \brief get particle at index from particle container
     *
     * Get a particle at the given index from the particle container.
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in]  index    index of particle in container
     * \param[out] globalid global id of particle
     * \param[out] states   states of particle
     */
    void GetParticle(int index, int& globalid, ParticleStates& states) const;

    /*!
     * \brief remove particle from particle container
     *
     * Remove a particle at the given index from the particle container. The particle is swapped
     * with the particle at the end of the container and the number of particles stored in the
     * container is decreased.
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] index index of particle in container
     */
    void RemoveParticle(int index);

    //! @}

    /*!
     * \brief get particle state dimension
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] state particle state
     *
     * \return dimension of particle state
     */
    inline int GetStateDim(ParticleState state)
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (not storedstates_.count(state))
        FOUR_C_THROW(
            "particle state '%s' not stored in container!", EnumToStateName(state).c_str());
#endif

      return statedim_[state];
    };

    //! \name access global id and particle states
    //! @{

    /*!
     * \brief get pointer to state of a particle at index
     *
     * This is the default method to be used to get a pointer to the state of a particle at a
     * certain index.
     *
     * \note Throws an error in the debug version in case the requested state is not stored in the
     *       particle container.
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] state particle state
     * \param[in] index index of particle in container
     *
     * \return pointer to particle state
     */
    inline double* GetPtrToState(ParticleState state, int index)
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (not storedstates_.count(state))
        FOUR_C_THROW(
            "particle state '%s' not stored in container!", EnumToStateName(state).c_str());

      if (index < 0 or index > (particlestored_ - 1))
        FOUR_C_THROW(
            "can not return pointer to state of particle as index %d out of bounds!", index);
#endif

      return &((states_[state])[index * statedim_[state]]);
    };

    /*!
     * \brief conditionally get pointer to state of a particle at index
     *
     * This method to get a pointer to the state of a particle at a certain index is used in cases
     * when a state may not be stored in the particle container. Conditionally, a pointer is
     * returned in case the state is stored in the particle container, otherwise, a nullptr is
     * returned.
     *
     * \note The returned pointer may not be used to access memory without checking for a nullptr.
     *
     * \author Sebastian Fuchs \date 12/2020
     *
     * \param[in] state particle state
     * \param[in] index index of particle in container
     *
     * \return pointer to particle state or nullptr
     */
    inline double* CondGetPtrToState(ParticleState state, int index)
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (index < 0 or index > (particlestored_ - 1))
        FOUR_C_THROW(
            "can not return pointer to state of particle as index %d out of bounds!", index);
#endif

      if (storedstates_.count(state)) return &((states_[state])[index * statedim_[state]]);

      return nullptr;
    };

    /*!
     * \brief get pointer to global id of a particle at index
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] index index of particle in container
     *
     * \return pointer to global id
     */
    inline int* GetPtrToGlobalID(int index)
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (index < 0 or index > (particlestored_ - 1))
        FOUR_C_THROW(
            "can not return pointer to global id of particle as index %d out of bounds!", index);
#endif

      return &(globalids_[index]);
    };

    //! @}

    //! \name manipulate particle states
    //! @{

    /*!
     * \brief scale state of particles
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] fac   scale factor
     * \param[in] state particle state
     */
    inline void ScaleState(double fac, ParticleState state)
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (not storedstates_.count(state))
        FOUR_C_THROW(
            "particle state '%s' not stored in container!", EnumToStateName(state).c_str());
#endif

      for (int i = 0; i < (particlestored_ * statedim_[state]); ++i) (states_[state])[i] *= fac;
    };

    /*!
     * \brief add scaled states to first state of particles
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] facA   first scale factor
     * \param[in] stateA first particle state
     * \param[in] facB   second scale factor
     * \param[in] stateB second particle state
     */
    inline void UpdateState(double facA, ParticleState stateA, double facB, ParticleState stateB)
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (not storedstates_.count(stateA))
        FOUR_C_THROW(
            "particle state '%s' not stored in container!", EnumToStateName(stateA).c_str());

      if (not storedstates_.count(stateB))
        FOUR_C_THROW(
            "particle state '%s' not stored in container!", EnumToStateName(stateB).c_str());

      if (statedim_[stateA] != statedim_[stateB])
        FOUR_C_THROW("dimensions of states do not match!");
#endif

      for (int i = 0; i < (particlestored_ * statedim_[stateA]); ++i)
        (states_[stateA])[i] = facA * (states_[stateA])[i] + facB * (states_[stateB])[i];
    };

    /*!
     * \brief set given state to all particles
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] val   particle state
     * \param[in] state particle state
     */
    inline void SetState(std::vector<double> val, ParticleState state)
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (not storedstates_.count(state))
        FOUR_C_THROW(
            "particle state '%s' not stored in container!", EnumToStateName(state).c_str());

      if (statedim_[state] != static_cast<int>(val.size()))
        FOUR_C_THROW("dimensions of states do not match!");
#endif

      for (int i = 0; i < particlestored_; ++i)
        for (int dim = 0; dim < statedim_[state]; ++dim)
          (states_[state])[i * statedim_[state] + dim] = val[dim];
    };

    /*!
     * \brief clear state of all particles
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] state particle state
     */
    inline void ClearState(ParticleState state)
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (not storedstates_.count(state))
        FOUR_C_THROW(
            "particle state '%s' not stored in container!", EnumToStateName(state).c_str());
#endif

      for (int i = 0; i < (particlestored_ * statedim_[state]); ++i) (states_[state])[i] = 0.0;
    };

    //! @}

    /*!
     * \brief get stored particle states
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \return stored particle states
     */
    inline const std::set<ParticleState>& GetStoredStates() const { return storedstates_; };

    /*!
     * \brief get flag indicating stored state
     *
     * Get a flag that is indicating if a state is stored in the particle container.
     *
     * \author Sebastian Fuchs \date 11/2019
     *
     * \param[in] state particle state
     *
     * \return flag indicating stored state
     */
    inline bool HaveStoredState(ParticleState state) const { return storedstates_.count(state); };

    /*!
     * \brief get size of particle container
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \return size of particle container
     */
    inline int ContainerSize() const { return containersize_; };

    /*!
     * \brief get number of particles stored in container
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \return number of particles stored in container
     */
    inline int ParticlesStored() const { return particlestored_; };

    /*!
     * \brief get minimum stored value of state in container
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \param[in] state particle state
     *
     * \return minimum stored value of state in container
     */
    double GetMinValueOfState(ParticleState state) const;

    /*!
     * \brief get maximum stored value of state in container
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \param[in] state particle state
     *
     * \return maximum stored value of state in container
     */
    double GetMaxValueOfState(ParticleState state) const;

   private:
    //! size of particles container
    int containersize_;

    //! number of particles stored in container
    int particlestored_;

    //! set of stored particle states
    std::set<ParticleState> storedstates_;

    //! size of vector for states
    int statesvectorsize_;

    //! global ids of stored particles
    std::vector<int> globalids_;

    //! particle states in container indexed by particle state enum
    std::vector<std::vector<double>> states_;

    //! particle state dimension indexed by particle state enum
    std::vector<int> statedim_;
  };

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
