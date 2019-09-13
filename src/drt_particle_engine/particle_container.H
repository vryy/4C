/*---------------------------------------------------------------------------*/
/*! \file
\brief smart particle container

\level 1

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_CONTAINER_H
#define PARTICLE_CONTAINER_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_enums.H"
#include "particle_typedefs.H"

/*---------------------------------------------------------------------------*
 | class definitions                                                         |
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

    //! destructor
    ~ParticleContainer(){};

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
     * \param[in] numContainerSize size of particle container
     * \param[in] stateEnumSet     set of particle states to be stored
     */
    void Setup(int numContainerSize, const std::set<StateEnum>& stateEnumSet);

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
#ifdef DEBUG
      if (particlestored_ > containersize_)
        dserror(
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
     * \param[in]  particle states of particle
     */
    void AddParticle(int& index, int globalid, const ParticleStates& particle);

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
     * \param[in] particle states of particle
     */
    void ReplaceParticle(int index, int globalid, const ParticleStates& particle);

    /*!
     * \brief get particle at index from particle container
     *
     * Get a particle at the given index from the particle container.
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in]  index    index of particle in container
     * \param[out] globalid global id of particle
     * \param[out] particle states of particle
     */
    void GetParticle(int index, int& globalid, ParticleStates& particle) const;

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
     * \param[in] stateEnum particle state enum
     *
     * \return dimension of particle state
     */
    inline int GetParticleStateDim(StateEnum stateEnum)
    {
#ifdef DEBUG
      if (not storedstates_.count(stateEnum))
        dserror("particle state '%s' not stored in container!",
            PARTICLEENGINE::EnumToStateName(stateEnum).c_str());
#endif

      return statedim_[stateEnum];
    };

    //! \name access global id and particle states
    //! @{

    /*!
     * \brief get pointer to state of a particle at index
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] stateEnum particle state enum
     * \param[in] index     index of particle in container
     *
     * \return pointer to particle state
     */
    inline double* GetPtrToParticleState(StateEnum stateEnum, int index)
    {
#ifdef DEBUG
      if (not storedstates_.count(stateEnum))
        dserror("particle state '%s' not stored in container!",
            PARTICLEENGINE::EnumToStateName(stateEnum).c_str());

      if (index < 0 or index > (particlestored_ - 1))
        dserror("can not return pointer to state of particle as index %d out of bounds!", index);
#endif

      return &((states_[stateEnum])[index * statedim_[stateEnum]]);
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
    inline int* GetPtrToParticleGlobalID(int index)
    {
#ifdef DEBUG
      if (index < 0 or index > (particlestored_ - 1))
        dserror(
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
     * \param[in] fac       scale factor
     * \param[in] stateEnum particle state enum
     */
    inline void ScaleState(double fac, StateEnum stateEnum)
    {
#ifdef DEBUG
      if (not storedstates_.count(stateEnum))
        dserror("particle state '%s' not stored in container!",
            PARTICLEENGINE::EnumToStateName(stateEnum).c_str());
#endif

      for (int i = 0; i < (particlestored_ * statedim_[stateEnum]); ++i)
        (states_[stateEnum])[i] *= fac;
    };

    /*!
     * \brief add scaled states to first state of particles
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] facA       first scale factor
     * \param[in] stateEnumA first particle state enum
     * \param[in] facB       second scale factor
     * \param[in] stateEnumB second particle state enum
     */
    inline void UpdateState(double facA, StateEnum stateEnumA, double facB, StateEnum stateEnumB)
    {
#ifdef DEBUG
      if (not storedstates_.count(stateEnumA))
        dserror("particle state '%s' not stored in container!",
            PARTICLEENGINE::EnumToStateName(stateEnumA).c_str());

      if (not storedstates_.count(stateEnumB))
        dserror("particle state '%s' not stored in container!",
            PARTICLEENGINE::EnumToStateName(stateEnumB).c_str());

      if (statedim_[stateEnumA] != statedim_[stateEnumB])
        dserror("dimensions of states do not match!");
#endif

      for (int i = 0; i < (particlestored_ * statedim_[stateEnumA]); ++i)
        (states_[stateEnumA])[i] =
            facA * (states_[stateEnumA])[i] + facB * (states_[stateEnumB])[i];
    };

    /*!
     * \brief set given state to all particles
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] val       particle state
     * \param[in] stateEnum particle state enum
     */
    inline void SetState(std::vector<double> val, StateEnum stateEnum)
    {
#ifdef DEBUG
      if (not storedstates_.count(stateEnum))
        dserror("particle state '%s' not stored in container!",
            PARTICLEENGINE::EnumToStateName(stateEnum).c_str());

      if (statedim_[stateEnum] != static_cast<int>(val.size()))
        dserror("dimensions of states do not match!");
#endif

      for (int i = 0; i < particlestored_; ++i)
        for (int dim = 0; dim < statedim_[stateEnum]; ++dim)
          (states_[stateEnum])[i * statedim_[stateEnum] + dim] = val[dim];
    };

    /*!
     * \brief clear state of all particles
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] stateEnum particle state enum
     */
    inline void ClearState(StateEnum stateEnum)
    {
#ifdef DEBUG
      if (not storedstates_.count(stateEnum))
        dserror("particle state '%s' not stored in container!",
            PARTICLEENGINE::EnumToStateName(stateEnum).c_str());
#endif

      for (int i = 0; i < (particlestored_ * statedim_[stateEnum]); ++i)
        (states_[stateEnum])[i] = 0.0;
    };

    //! @}

    /*!
     * \brief get stored particle states
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \return stored particle state enums
     */
    inline const std::set<StateEnum>& GetStoredStates() const { return storedstates_; };

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
     * \param[in] stateEnum particle state enum
     *
     * \return minimum stored value of state in container
     */
    double GetMinValueOfState(StateEnum stateEnum) const;

    /*!
     * \brief get maximum stored value of state in container
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \param[in] stateEnum particle state enum
     *
     * \return maximum stored value of state in container
     */
    double GetMaxValueOfState(StateEnum stateEnum) const;

   private:
    //! size of particles container
    int containersize_;

    //! number of particles stored in container
    int particlestored_;

    //! set of stored particle states
    std::set<StateEnum> storedstates_;

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
#endif
