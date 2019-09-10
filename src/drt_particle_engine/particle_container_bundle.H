/*---------------------------------------------------------------------------*/
/*! \file
\brief class holding all particle containers

\level 3

\maintainer  Sebastian Fuchs
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

#include "particle_container.H"

#include <Teuchos_RCP.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
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
    explicit ParticleContainerBundle();

    //! destructor
    ~ParticleContainerBundle(){};

    //! init particle container bundle
    void Init();

    //! setup particle container bundle
    void Setup(const std::map<TypeEnum, std::set<StateEnum>>& particlestatestotypes);

    //! get particle types of stored containers
    inline const std::set<TypeEnum> GetParticleTypes() const { return storedtypes_; };

    //! get specific particle container
    inline ParticleContainer* GetSpecificContainer(TypeEnum typeEnum, StatusEnum statusEnum) const
    {
#ifdef DEBUG
      if (not storedtypes_.count(typeEnum))
        dserror("container for particle type '%s' is not stored!",
            PARTICLEENGINE::EnumToTypeName(typeEnum).c_str());
#endif

      return (containers_[typeEnum])[statusEnum].get();
    };

    //! scale state of owned particles in specific container
    inline void ScaleStateSpecificContainer(
        double fac, StateEnum stateEnum, TypeEnum typeEnum) const
    {
#ifdef DEBUG
      if (not storedtypes_.count(typeEnum))
        dserror("container for particle type '%s' is not stored!",
            PARTICLEENGINE::EnumToTypeName(typeEnum).c_str());
#endif

      ((containers_[typeEnum])[PARTICLEENGINE::Owned])->ScaleState(fac, stateEnum);
    };

    //! update state of owned particles in specific container
    inline void UpdateStateSpecificContainer(double facA, StateEnum stateEnumA, double facB,
        StateEnum stateEnumB, TypeEnum typeEnum) const
    {
#ifdef DEBUG
      if (not storedtypes_.count(typeEnum))
        dserror("container for particle type '%s' is not stored!",
            PARTICLEENGINE::EnumToTypeName(typeEnum).c_str());
#endif

      ((containers_[typeEnum])[PARTICLEENGINE::Owned])
          ->UpdateState(facA, stateEnumA, facB, stateEnumB);
    };

    //! set state of owned particles in specific container
    inline void SetStateSpecificContainer(
        std::vector<double> val, StateEnum stateEnum, TypeEnum typeEnum) const
    {
#ifdef DEBUG
      if (not storedtypes_.count(typeEnum))
        dserror("container for particle type '%s' is not stored!",
            PARTICLEENGINE::EnumToTypeName(typeEnum).c_str());
#endif

      ((containers_[typeEnum])[PARTICLEENGINE::Owned])->SetState(val, stateEnum);
    };

    //! clear state of owned particles in specific container
    inline void ClearStateSpecificContainer(StateEnum stateEnum, TypeEnum typeEnum) const
    {
#ifdef DEBUG
      if (not storedtypes_.count(typeEnum))
        dserror("container for particle type '%s' is not stored!",
            PARTICLEENGINE::EnumToTypeName(typeEnum).c_str());
#endif

      ((containers_[typeEnum])[PARTICLEENGINE::Owned])->ClearState(stateEnum);
    };

    //! scale state of owned particles in all containers
    inline void ScaleStateAllContainers(double fac, StateEnum stateEnum) const
    {
      for (auto& typeEnum : storedtypes_)
        ((containers_[typeEnum])[PARTICLEENGINE::Owned])->ScaleState(fac, stateEnum);
    };

    //! update state of owned particles in all containers
    inline void UpdateStateAllContainers(
        double facA, StateEnum stateEnumA, double facB, StateEnum stateEnumB) const
    {
      for (auto& typeEnum : storedtypes_)
        ((containers_[typeEnum])[PARTICLEENGINE::Owned])
            ->UpdateState(facA, stateEnumA, facB, stateEnumB);
    };

    //! set state of owned particles in all containers
    inline void SetStateAllContainers(std::vector<double> val, StateEnum stateEnum) const
    {
      for (auto& typeEnum : storedtypes_)
        ((containers_[typeEnum])[PARTICLEENGINE::Owned])->SetState(val, stateEnum);
    };

    //! clear state of owned particles in all containers
    inline void ClearStateAllContainers(StateEnum stateEnum) const
    {
      for (auto& typeEnum : storedtypes_)
        ((containers_[typeEnum])[PARTICLEENGINE::Owned])->ClearState(stateEnum);
    };

    //! check and decrease the size of all containers of specific status
    inline void CheckAndDecreaseSizeAllContainersOfSpecificStatus(StatusEnum statusEnum) const
    {
      for (auto& typeEnum : storedtypes_)
        ((containers_[typeEnum])[statusEnum])->CheckAndDecreaseContainerSize();
    }

    //! clear all containers of specific status
    inline void ClearAllContainersOfSpecificStatus(StatusEnum statusEnum) const
    {
      for (auto& typeEnum : storedtypes_) ((containers_[typeEnum])[statusEnum])->ClearContainer();
    };

    //! pack all particle containers
    void PackParticleContainerBundle(Teuchos::RCP<std::vector<char>>& particlebuffer) const;

    //! get vector of particle objects of all containers
    void GetVectorOfParticleObjectsOfAllContainers(
        std::vector<ParticleObjShrdPtr>& particlesstored) const;

   private:
    //! set of particle types of stored containers
    std::set<TypeEnum> storedtypes_;

    //! data structure holding all particle containers indexed by type enum and status enum
    TypeStatusContainers containers_;
  };

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
