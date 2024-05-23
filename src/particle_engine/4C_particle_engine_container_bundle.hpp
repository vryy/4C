/*---------------------------------------------------------------------------*/
/*! \file
\brief manage bundle of particle containers
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_ENGINE_CONTAINER_BUNDLE_HPP
#define FOUR_C_PARTICLE_ENGINE_CONTAINER_BUNDLE_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_particle_engine_container.hpp"
#include "4C_particle_engine_enums.hpp"
#include "4C_particle_engine_typedefs.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleObject;
}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
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
    void Setup(const std::map<ParticleType, std::set<ParticleState>>& particlestatestotypes);

    /*!
     * \brief get particle types of stored containers
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \return reference to particle types of stored containers
     */
    inline const std::set<ParticleType>& GetParticleTypes() const { return storedtypes_; };

    /*!
     * \brief get specific particle container
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] type   particle type
     * \param[in] status particle status
     *
     * @return pointer to particle container
     */
    inline ParticleContainer* get_specific_container(ParticleType type, ParticleStatus status) const
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (not storedtypes_.count(type))
        FOUR_C_THROW("container for particle type '%s' not stored!", EnumToTypeName(type).c_str());
#endif

      return (containers_[type])[status].get();
    };

    //! \name manipulate particle states of owned particles of specific type
    //! @{

    /*!
     * \brief scale state of particles in container of owned particles of specific type
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] fac       scale factor
     * \param[in] state particle state
     * \param[in] type  particle type
     */
    inline void scale_state_specific_container(
        double fac, ParticleState state, ParticleType type) const
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (not storedtypes_.count(type))
        FOUR_C_THROW("container for particle type '%s' not stored!", EnumToTypeName(type).c_str());
#endif

      ((containers_[type])[Owned])->ScaleState(fac, state);
    };

    /*!
     * \brief add scaled states to first state of particles in container of owned particles of
     *        specific type
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] facA   first scale factor
     * \param[in] stateA first particle state
     * \param[in] facB   second scale factor
     * \param[in] stateB second particle state
     * \param[in] type   particle type
     */
    inline void update_state_specific_container(double facA, ParticleState stateA, double facB,
        ParticleState stateB, ParticleType type) const
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (not storedtypes_.count(type))
        FOUR_C_THROW("container for particle type '%s' not stored!", EnumToTypeName(type).c_str());
#endif

      ((containers_[type])[Owned])->UpdateState(facA, stateA, facB, stateB);
    };

    /*!
     * \brief set given state to all particles in container of owned particles of specific type
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] val   particle state value
     * \param[in] state particle state
     * \param[in] type  particle type
     */
    inline void set_state_specific_container(
        std::vector<double> val, ParticleState state, ParticleType type) const
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (not storedtypes_.count(type))
        FOUR_C_THROW("container for particle type '%s' not stored!", EnumToTypeName(type).c_str());
#endif

      ((containers_[type])[Owned])->SetState(val, state);
    };

    /*!
     * \brief clear state of all particles in container of owned particles of specific type
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] state particle state
     * \param[in] type  particle type
     */
    inline void clear_state_specific_container(ParticleState state, ParticleType type) const
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (not storedtypes_.count(type))
        FOUR_C_THROW("container for particle type '%s' not stored!", EnumToTypeName(type).c_str());
#endif

      ((containers_[type])[Owned])->ClearState(state);
    };

    //! @}

    //! \name manipulate particle states of owned particles of all types
    //! @{

    /*!
     * \brief scale state of particles in container of owned particles of all types
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] fac   scale factor
     * \param[in] state particle state
     */
    inline void scale_state_all_containers(double fac, ParticleState state) const
    {
      for (const auto& type : storedtypes_) ((containers_[type])[Owned])->ScaleState(fac, state);
    };

    /*!
     * \brief add scaled states to first state of particles in container of owned particles of all
     *        types
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] facA   first scale factor
     * \param[in] stateA first particle state
     * \param[in] facB   second scale factor
     * \param[in] stateB second particle state
     */
    inline void update_state_all_containers(
        double facA, ParticleState stateA, double facB, ParticleState stateB) const
    {
      for (const auto& type : storedtypes_)
        ((containers_[type])[Owned])->UpdateState(facA, stateA, facB, stateB);
    };

    /*!
     * \brief set given state to all particles in container of owned particles of all types
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] val   particle state value
     * \param[in] state particle state
     */
    inline void set_state_all_containers(std::vector<double> val, ParticleState state) const
    {
      for (const auto& type : storedtypes_) ((containers_[type])[Owned])->SetState(val, state);
    };

    /*!
     * \brief clear state of all particles in container of owned particles of all types
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] state particle state
     */
    inline void clear_state_all_containers(ParticleState state) const
    {
      for (const auto& type : storedtypes_) ((containers_[type])[Owned])->ClearState(state);
    };

    //! @}

    //! \name manipulate particle container of specific status
    //! @{

    /*!
     * \brief check and decrease the size of all containers of specific status
     *
     * \author Sebastian Fuchs \date 07/2019
     *
     * \param[in] status particle status
     */
    inline void check_and_decrease_size_all_containers_of_specific_status(
        ParticleStatus status) const
    {
      for (const auto& type : storedtypes_)
        ((containers_[type])[status])->check_and_decrease_container_size();
    }

    /*!
     * \brief clear all containers of specific status
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] status particle status
     */
    inline void clear_all_containers_of_specific_status(ParticleStatus status) const
    {
      for (const auto& type : storedtypes_) ((containers_[type])[status])->ClearContainer();
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
    void get_packed_particle_objects_of_all_containers(
        std::shared_ptr<std::vector<char>>& particlebuffer) const;

    /*!
     * \brief get particle objects of all containers
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[out] particlesstored particle objects of all containers
     */
    void get_vector_of_particle_objects_of_all_containers(
        std::vector<ParticleObjShrdPtr>& particlesstored) const;

    //! @}

   private:
    //! set of particle types of stored containers
    std::set<ParticleType> storedtypes_;

    //! collection of particle containers indexed by particle type enum and particle status enum
    TypeStatusContainers containers_;
  };

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
