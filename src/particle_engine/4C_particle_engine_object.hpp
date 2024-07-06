/*---------------------------------------------------------------------------*/
/*! \file
\brief particle object for parallel communication
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_ENGINE_OBJECT_HPP
#define FOUR_C_PARTICLE_ENGINE_OBJECT_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_comm_parobjectfactory.hpp"
#include "4C_particle_engine_enums.hpp"
#include "4C_particle_engine_typedefs.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace Core::Communication
{
  class PackBuffer;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  /*!
   * \brief particle object type singleton
   *
   * \author Sebastian Fuchs \date 03/2018
   */
  class ParticleObjectType final : public Core::Communication::ParObjectType
  {
   public:
    std::string name() const override { return "ParticleObjectType"; };

    /*!
     * \brief get instance of particle object type
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \return reference to instance of particle object type
     */
    static ParticleObjectType& instance() { return instance_; };

    Core::Communication::ParObject* create(const std::vector<char>& data) override;

   private:
    //! particle object type instance
    static ParticleObjectType instance_;
  };

  /*!
   * \brief particle object for parallel communication of particle data
   *
   * The class members of the particle object are either initialized via an constructor with
   * initializer list or via unpacking of the data (after communication and receiving of packed
   * data).
   *
   * \note A separate method for initializing or modifying of class members is not provided due to
   *       performance reasons (avoid copy operations) and due to the fact that a particle object
   *       only needs to be packed once (without being modified afterwards) and is directly sent to
   *       another processor.
   *
   * \author Sebastian Fuchs \date 03/2018
   */
  class ParticleObject : public Core::Communication::ParObject
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 03/2018
     */
    ParticleObject();

    /*!
     * \brief constructor with initializer list
     *
     * Construct particle object and set class members via initializer list.
     *
     * \author Sebastian Fuchs \date 02/2019
     *
     * \param[in] type     particle type
     * \param[in] globalid global id of particle
     * \param[in] states   states of particle
     * \param[in] bingid   global id of bin the particle is located in
     *                     optional: set to -1 if omnitted
     * \param[in] index    index of particle in container
     *                     optional: set to -1 if omnitted
     */
    ParticleObject(ParticleType type, int globalid, const ParticleStates& states, int bingid = -1,
        int index = -1);

    int unique_par_object_id() const override
    {
      return ParticleObjectType::instance().unique_par_object_id();
    };

    void pack(Core::Communication::PackBuffer& data) const override;

    void unpack(const std::vector<char>& data) override;

    //! \name set particle object members
    //! @{

    /*!
     * \brief set global id of particle
     *
     * \author Sebastian Fuchs \date 11/2019
     *
     * \param[in] particleglobalid global id of particle
     */
    inline void set_particle_global_id(int globalid) { globalid_ = globalid; };

    //! @}

    //! \name get particle object members
    //! @{

    /*!
     * \brief get particle type
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \return particle type
     */
    inline ParticleType return_particle_type() const { return type_; };

    /*!
     * \brief get global id of particle
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \return global id of particle
     */
    inline int return_particle_global_id() const { return globalid_; };

    /*!
     * \brief get states of particle
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \return reference to states of particle
     */
    inline const ParticleStates& return_particle_states() const { return states_; };

    /*!
     * \brief get global id of bin the particle is located in
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \return global id of bin the particle is located in
     */
    inline int return_bin_gid() const { return bingid_; };

    /*!
     * \brief get index of particle in container
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \return index of particle in container
     */
    inline int return_container_index() const { return index_; };

    //! @}

   private:
    //! particle type
    ParticleType type_;

    //! global id of particle
    int globalid_;

    //! states of particle
    ParticleStates states_;

    //! global id of bin the particle is located in
    int bingid_;

    //! index of particle in container (owned or ghosted depending on case)
    int index_;
  };

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
