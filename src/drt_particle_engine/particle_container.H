/*---------------------------------------------------------------------------*/
/*!
\file particle_container.H

\brief smart particle container class

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_CONTAINER_H
#define PARTICLE_CONTAINER_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_enums.H"
#include "particle_typedefs.H"

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleContainer
  {
   public:
    //! constructor
    explicit ParticleContainer();

    //! destructor
    ~ParticleContainer(){};

    //! init particle container
    void Init();

    //! setup particle container
    void Setup(int numContainerSize, const std::set<StateEnum>& stateEnumSet);

    //! increase the container size
    void IncreaseContainerSize();

    //! decrease the container size
    void DecreaseContainerSize();

    //! clear particle container
    inline void ClearContainer() { particlestored_ = 0; };

    //! add particle to particle container and get index
    void AddParticle(int& index, int globalid, const ParticleStates& particle);

    //! replace particle in particle container at given index
    void ReplaceParticle(int index, int globalid, const ParticleStates& particle);

    //! get particle at index from particle container
    void GetParticle(int index, int& globalid, ParticleStates& particle) const;

    //! remove particle from particle container
    void RemoveParticle(int index);

    //! get state of a particle at index
    std::vector<double> GetParticleState(StateEnum stateEnum, int index) const;

    //! return pointer to state of a particle at index
    inline double* GetPtrToParticleState(StateEnum stateEnum, int index)
    {
      if (index < 0 or index > (particlestored_ - 1))
        dserror("can not return pointer to state of particle as index %d out of bounds!", index);

      return &((states_[stateEnum])[index * statedim_[stateEnum]]);
    };

    //! return pointer to global id of a particle at index
    inline int* GetPtrToParticleGlobalID(int index)
    {
      if (index < 0 or index > (particlestored_ - 1))
        dserror(
            "can not return pointer to global id of particle as index %d out of bounds!", index);

      return &(globalids_[index]);
    };

    //! scale state of particles
    inline void ScaleState(double fac, StateEnum stateEnum)
    {
      for (int i = 0; i < (particlestored_ * statedim_[stateEnum]); ++i)
        (states_[stateEnum])[i] *= fac;
    };

    //! scale and add states to update state of particles
    inline void UpdateState(double facA, StateEnum stateEnumA, double facB, StateEnum stateEnumB)
    {
      if (statedim_[stateEnumA] != statedim_[stateEnumB])
        dserror("dimensions of states do not match!");

      for (int i = 0; i < (particlestored_ * statedim_[stateEnumA]); ++i)
        (states_[stateEnumA])[i] =
            facA * (states_[stateEnumA])[i] + facB * (states_[stateEnumB])[i];
    };

    //! set state to particles
    inline void SetState(std::vector<double> val, StateEnum stateEnum)
    {
      if (statedim_[stateEnum] != static_cast<int>(val.size()))
        dserror("dimensions of states do not match!");

      for (int i = 0; i < particlestored_; ++i)
        for (int dim = 0; dim < statedim_[stateEnum]; ++dim)
          (states_[stateEnum])[i * statedim_[stateEnum] + dim] = val[dim];
    };

    //! clear state of particles
    inline void ClearState(StateEnum stateEnum)
    {
      for (int i = 0; i < (particlestored_ * statedim_[stateEnum]); ++i)
        (states_[stateEnum])[i] = 0.0;
    };

    //! get stored particle states
    inline const std::set<StateEnum> GetStoredStates() const { return storedstates_; };

    //! get size of particle container
    inline int ContainerSize() const { return containersize_; };

    //! get number of particles stored in container
    inline int ParticlesStored() const { return particlestored_; };

    //! get minimum stored value of state
    double GetMinValueOfState(StateEnum stateEnum) const;

    //! get maximum stored value of state
    double GetMaxValueOfState(StateEnum stateEnum) const;

   private:
    //! size of particles container
    int containersize_;

    //! number of particles stored in container
    int particlestored_;

    //! set of stored particle states
    std::set<StateEnum> storedstates_;

    //! global ids of stored particles
    std::vector<int> globalids_;

    //! particle states in container indexed by state enum
    std::vector<std::vector<double>> states_;

    //! particle state dimension indexed by state enum
    std::vector<int> statedim_;
  };

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
