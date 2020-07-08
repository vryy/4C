/*---------------------------------------------------------------------------*/
/*! \file
\brief particle object for parallel communication
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_OBJECT_H
#define PARTICLE_OBJECT_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_enums.H"
#include "particle_typedefs.H"

#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_parobjectfactory.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace DRT
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
  class ParticleObjectType final : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "ParticleObjectType"; };

    /*!
     * \brief get instance of particle object type
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \return reference to instance of particle object type
     */
    static ParticleObjectType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

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
  class ParticleObject : public DRT::ParObject
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
     * \param[in] particletype     particle type enum
     * \param[in] particleglobalid global id of particle
     * \param[in] particlestates   states of particle
     * \param[in] bingid           global id of bin the particle is located in
     *                             optional: set to -1 if omnitted
     * \param[in] containerindex   index of particle in container
     *                             optional: set to -1 if omnitted
     */
    ParticleObject(TypeEnum particletype, int particleglobalid,
        const ParticleStates& particlestates, int bingid = -1, int containerindex = -1);

    //! destructor
    ~ParticleObject() override{};

    int UniqueParObjectId() const override
    {
      return ParticleObjectType::Instance().UniqueParObjectId();
    };

    void Pack(DRT::PackBuffer& data) const override;

    void Unpack(const std::vector<char>& data) override;

    //! \name set particle object members
    //! @{

    /*!
     * \brief set global id of particle
     *
     * \author Sebastian Fuchs \date 11/2019
     *
     * \param[in] particleglobalid global id of particle
     */
    inline void SetParticleGlobalID(int particleglobalid) { particleglobalid_ = particleglobalid; };

    //! @}

    //! \name get particle object members
    //! @{

    /*!
     * \brief get particle type enum
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \return particle type enum
     */
    inline TypeEnum ReturnParticleType() const { return particletype_; };

    /*!
     * \brief get global id of particle
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \return global id of particle
     */
    inline int ReturnParticleGlobalID() const { return particleglobalid_; };

    /*!
     * \brief get states of particle
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \return reference to states of particle
     */
    inline const ParticleStates& ReturnParticleStates() const { return particlestates_; };

    /*!
     * \brief get global id of bin the particle is located in
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \return global id of bin the particle is located in
     */
    inline int ReturnBinGid() const { return bingid_; };

    /*!
     * \brief get index of particle in container
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \return index of particle in container
     */
    inline int ReturnContainerIndex() const { return containerindex_; };

    //! @}

   private:
    //! particle type enum
    TypeEnum particletype_;

    //! global id of particle
    int particleglobalid_;

    //! states of particle
    ParticleStates particlestates_;

    //! global id of bin the particle is located in
    int bingid_;

    //! index of particle in container (owned or ghosted depending on case)
    int containerindex_;
  };

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
