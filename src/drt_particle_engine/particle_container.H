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

    //! return pointer to state of a particle at index
    double* GetPtrToParticleState(StateEnum stateEnum, int index);

    //! return pointer to global id of a particle at index
    int* GetPtrToParticleGlobalID(int index);

    //! get state of a particle at index
    std::vector<double> GetParticleState(StateEnum stateEnum, int index) const;

    //! scale state of particles
    void ScaleState(double fac, StateEnum stateEnum);

    //! scale and add states to update state of particles
    void UpdateState(double facA, StateEnum stateEnumA, double facB, StateEnum stateEnumB);

    //! set state to particles
    void SetState(std::vector<double> val, StateEnum stateEnum);

    //! clear state of particles
    void ClearState(StateEnum stateEnum);

    //! get minimum stored value of state
    double GetMinValueOfState(StateEnum stateEnum) const;

    //! get maximum stored value of state
    double GetMaxValueOfState(StateEnum stateEnum) const;

    //! get stored particle states
    inline const std::set<StateEnum> GetParticleStates() const { return storedstates_; };

    //! get size of particle container
    inline int ContainerSize() const { return containersize_; };

    //! get number of particles stored in container
    inline int ParticlesStored() const { return particlestored_; };

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
