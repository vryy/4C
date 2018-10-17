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
#include <map>
#include <set>
#include <vector>
#include <memory>

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
    void ClearContainer();

    //! add particle to particle container and get index
    void AddParticle(int& index, int globalid, const ParticleStates& particle);

    //! replace particle in particle container at given index
    void ReplaceParticle(int index, int globalid, const ParticleStates& particle);

    //! get particle at index from particle container
    void GetParticle(int index, int& globalid, ParticleStates& particle) const;

    //! remove particle from particle container
    void RemoveParticle(int index);

    //! swap particles at index a and b in particle container
    void SwapParticle(int index_a, int index_b);

    //! return pointer to state of a particle at index
    double* GetPtrToParticleState(StateEnum stateEnum, int index) const;

    //! get state of a particle at index
    std::vector<double> GetParticleState(StateEnum stateEnum, int index) const;

    //! scale state of particles
    void ScaleState(double fac, StateEnum stateEnum) const;

    //! scale and add states to update state of particles
    void UpdateState(double facA, StateEnum stateEnumA, double facB, StateEnum stateEnumB) const;

    //! set state to particles
    void SetState(std::vector<double> val, StateEnum stateEnum) const;

    //! clear state of particles
    void ClearState(StateEnum stateEnum) const;

    //! get stored particle states
    const std::set<StateEnum> GetParticleStates() const;

    //! get size of particle container
    int ContainerSize() const { return containersize_; };

    //! get number of particles stored in container
    int ParticlesStored() const { return particlestored_; };

   private:
    //! size of particles container
    int containersize_;

    //! number of particles stored in container
    int particlestored_;

    //! global ids of stored particles
    std::vector<int> globalids_;

    //! maps type of particle state to corresponding data type
    std::map<StateEnum, std::shared_ptr<std::vector<double>>> states_;
  };

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
