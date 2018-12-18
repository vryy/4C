/*---------------------------------------------------------------------------*/
/*!
\file particle_container_bundle.H

\brief class holding all particle containers

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_CONTAINER_BUNDLE_H
#define PARTICLE_CONTAINER_BUNDLE_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_enums.H"
#include "particle_typedefs.H"

#include <map>
#include <set>
#include <vector>
#include <memory>

#include <Teuchos_RCP.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleContainer;
  class ParticleObject;
}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleContainerBundle
  {
   public:
    //! constructor
    explicit ParticleContainerBundle(const int myrank);

    //! destructor
    ~ParticleContainerBundle(){};

    //! init particle container bundle
    void Init();

    //! setup particle container bundle
    void Setup(const std::map<TypeEnum, std::set<StateEnum>>& particlestatestotypes);

    //! scale state of owned particles in specific container
    void ScaleStateSpecificContainer(double fac, StateEnum stateEnum, TypeEnum typeEnum) const;

    //! update state of owned particles in specific container
    void UpdateStateSpecificContainer(double facA, StateEnum stateEnumA, double facB,
        StateEnum stateEnumB, TypeEnum typeEnum) const;

    //! set state of owned particles in specific container
    void SetStateSpecificContainer(
        std::vector<double> val, StateEnum stateEnum, TypeEnum typeEnum) const;

    //! clear state of owned particles in specific container
    void ClearStateSpecificContainer(StateEnum stateEnum, TypeEnum typeEnum) const;

    //! scale state of owned particles in all containers
    void ScaleStateAllContainers(double fac, StateEnum stateEnum) const;

    //! update state of owned particles in all containers
    void UpdateStateAllContainers(
        double facA, StateEnum stateEnumA, double facB, StateEnum stateEnumB) const;

    //! set state of owned particles in all containers
    void SetStateAllContainers(std::vector<double> val, StateEnum stateEnum) const;

    //! clear state of owned particles in all containers
    void ClearStateAllContainers(StateEnum stateEnum) const;

    //! get particle types of stored containers
    inline const std::set<TypeEnum> GetParticleTypes() const { return storedtypes_; };

    //! get specific particle container
    ParticleContainerShrdPtr GetSpecificContainer(TypeEnum typeEnum, StatusEnum statusEnum) const;

    //! clear all containers of specific status
    void ClearAllContainersOfSpecificStatus(StatusEnum statusEnum) const;

    //! pack all particle containers
    void PackParticleContainerBundle(Teuchos::RCP<std::vector<char>>& particlebuffer) const;

    //! get vector of particle objects of all containers
    void GetVectorOfParticleObjectsOfAllContainers(
        std::vector<ParticleObjShrdPtr>& particlesstored) const;

   private:
    //! processor id
    const int myrank_;

    //! set of particle types of stored containers
    std::set<TypeEnum> storedtypes_;

    //! data structure holding all particle containers indexed by type enum and status enum
    TypeStatusContainers containers_;
  };

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
