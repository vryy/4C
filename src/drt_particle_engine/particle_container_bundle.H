/*---------------------------------------------------------------------------*/
/*! \file
\brief manage bundle of particle containers

\level 1

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_CONTAINER_BUNDLE_H
#define PARTICLE_CONTAINER_BUNDLE_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_enums.H"
#include "particle_typedefs.H"

#include "particle_container.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleObject;
}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*
 | class definitions                                                         |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  /*!
   * \brief handler managing bundle of particle containers
   *
   * A handler managing the access to the bundle of particle containers. For each particle type a
   * container for owned particles and a container for ghosted particles is initialized.
   *
   * \author Sebastian Fuchs \date 05/2018
   */
  class ParticleContainerBundle final
  {
   public:
    //! constructor
    explicit ParticleContainerBundle();

    //! destructor
    ~ParticleContainerBundle(){};

    /*!
     * \brief init particle container bundle
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void Init();

    /*!
     * \brief setup particle container bundle
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] particlestatestotypes particle types and corresponding states
     */
    void Setup(const std::map<TypeEnum, std::set<StateEnum>>& particlestatestotypes);

    /*!
     * \brief get particle types of stored containers
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \return reference to particle types of stored containers
     */
    inline const std::set<TypeEnum>& GetParticleTypes() const { return storedtypes_; };

    /*!
     * \brief get specific particle container
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] typeEnum   particle type enum
     * \param[in] statusEnum particle status enum
     *
     * @return pointer to particle container
     */
    inline ParticleContainer* GetSpecificContainer(TypeEnum typeEnum, StatusEnum statusEnum) const
    {
#ifdef DEBUG
      if (not storedtypes_.count(typeEnum))
        dserror("container for particle type '%s' is not stored!",
            PARTICLEENGINE::EnumToTypeName(typeEnum).c_str());
#endif

      return (containers_[typeEnum])[statusEnum].get();
    };

    //! \name manipulate particle states of owned particles of specific type
    //! @{

    /*!
     * \brief scale state of particles in container of owned particles of specific type
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] fac       scale factor
     * \param[in] stateEnum particle state enum
     * \param[in] typeEnum  particle type enum
     */
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

    /*!
     * \brief add scaled states to first state of particles in container of owned particles of
     *        specific type
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] facA       first scale factor
     * \param[in] stateEnumA first particle state enum
     * \param[in] facB       second scale factor
     * \param[in] stateEnumB second particle state enum
     * \param[in] typeEnum   particle type enum
     */
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

    /*!
     * \brief set given state to all particles in container of owned particles of specific type
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] val       particle state
     * \param[in] stateEnum particle state enum
     * \param[in] typeEnum  particle type enum
     */
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

    /*!
     * \brief clear state of all particles in container of owned particles of specific type
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] stateEnum particle state enum
     * \param[in] typeEnum  particle type enum
     */
    inline void ClearStateSpecificContainer(StateEnum stateEnum, TypeEnum typeEnum) const
    {
#ifdef DEBUG
      if (not storedtypes_.count(typeEnum))
        dserror("container for particle type '%s' is not stored!",
            PARTICLEENGINE::EnumToTypeName(typeEnum).c_str());
#endif

      ((containers_[typeEnum])[PARTICLEENGINE::Owned])->ClearState(stateEnum);
    };

    //! @}

    //! \name manipulate particle states of owned particles of all types
    //! @{

    /*!
     * \brief scale state of particles in container of owned particles of all types
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] fac       scale factor
     * \param[in] stateEnum particle state enum
     */
    inline void ScaleStateAllContainers(double fac, StateEnum stateEnum) const
    {
      for (auto& typeEnum : storedtypes_)
        ((containers_[typeEnum])[PARTICLEENGINE::Owned])->ScaleState(fac, stateEnum);
    };

    /*!
     * \brief add scaled states to first state of particles in container of owned particles of all
     *        types
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] facA       first scale factor
     * \param[in] stateEnumA first particle state enum
     * \param[in] facB       second scale factor
     * \param[in] stateEnumB second particle state enum
     */
    inline void UpdateStateAllContainers(
        double facA, StateEnum stateEnumA, double facB, StateEnum stateEnumB) const
    {
      for (auto& typeEnum : storedtypes_)
        ((containers_[typeEnum])[PARTICLEENGINE::Owned])
            ->UpdateState(facA, stateEnumA, facB, stateEnumB);
    };

    /*!
     * \brief set given state to all particles in container of owned particles of all types
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] val       particle state
     * \param[in] stateEnum particle state enum
     */
    inline void SetStateAllContainers(std::vector<double> val, StateEnum stateEnum) const
    {
      for (auto& typeEnum : storedtypes_)
        ((containers_[typeEnum])[PARTICLEENGINE::Owned])->SetState(val, stateEnum);
    };

    /*!
     * \brief clear state of all particles in container of owned particles of all types
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] stateEnum particle state enum
     */
    inline void ClearStateAllContainers(StateEnum stateEnum) const
    {
      for (auto& typeEnum : storedtypes_)
        ((containers_[typeEnum])[PARTICLEENGINE::Owned])->ClearState(stateEnum);
    };

    //! @}

    //! \name manipulate particle container of specific status
    //! @{

    /*!
     * \brief check and decrease the size of all containers of specific status
     *
     * \author Sebastian Fuchs \date 07/2019
     *
     * \param[in] statusEnum particle status enum
     */
    inline void CheckAndDecreaseSizeAllContainersOfSpecificStatus(StatusEnum statusEnum) const
    {
      for (auto& typeEnum : storedtypes_)
        ((containers_[typeEnum])[statusEnum])->CheckAndDecreaseContainerSize();
    }

    /*!
     * \brief clear all containers of specific status
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] statusEnum particle status enum
     */
    inline void ClearAllContainersOfSpecificStatus(StatusEnum statusEnum) const
    {
      for (auto& typeEnum : storedtypes_) ((containers_[typeEnum])[statusEnum])->ClearContainer();
    };

    //! @}

    //! \name get particle objects of all particle containers
    //! @{

    /*!
     * \brief get packed particle objects of all containers
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[out] particlebuffer buffer of packed particle objects of all containers
     */
    void GetPackedParticleObjectsOfAllContainers(
        std::shared_ptr<std::vector<char>>& particlebuffer) const;

    /*!
     * \brief get particle objects of all containers
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[out] particlesstored particle objects of all containers
     */
    void GetVectorOfParticleObjectsOfAllContainers(
        std::vector<ParticleObjShrdPtr>& particlesstored) const;

    //! @}

   private:
    //! set of particle types of stored containers
    std::set<TypeEnum> storedtypes_;

    //! collection of particle containers indexed by particle type enum and particle status enum
    TypeStatusContainers containers_;
  };

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
