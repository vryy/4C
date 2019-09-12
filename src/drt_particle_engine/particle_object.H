/*---------------------------------------------------------------------------*/
/*! \file
\brief particle object class for communication

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_OBJECT_H
#define PARTICLE_OBJECT_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_enums.H"
#include "particle_typedefs.H"

#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_parobjectfactory.H"

#include <map>
#include <vector>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
namespace DRT
{
  class PackBuffer;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleObjectType final : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "ParticleObjectType"; };

    static ParticleObjectType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ParticleObjectType instance_;
  };

  /**
   * \brief A particle object for parallel communication of particle data.
   *
   * The class members of the particle object are either initialized via an constructor with
   * initializer list or via unpacking of the data (after communication and receiving of packed
   * data).
   *
   * \note A separate method for initializing or modifying of class members is not provided due to
   * performance reasons (avoid copy operations) and due to the fact that a particle object only
   * needs to be packed once (without being modified afterwards) and is directly sent to another
   * processor.
   */
  class ParticleObject : public DRT::ParObject
  {
   public:
    //! constructor
    ParticleObject();

    //! standard constructor
    ParticleObject(TypeEnum particletype, int particleglobalid,
        const ParticleStates& particlestates, int bingid = -1, int containerindex = -1);

    //! destructor
    ~ParticleObject() override{};

    //! return unique ParObject id
    int UniqueParObjectId() const override
    {
      return ParticleObjectType::Instance().UniqueParObjectId();
    };

    //! pack this class so it can be communicated
    void Pack(DRT::PackBuffer& data) const override;

    //! unpack data from a char vector into this class
    void Unpack(const std::vector<char>& data) override;

    //! return particle type
    inline TypeEnum ReturnParticleType() const { return particletype_; };

    //! return particle global id
    inline int ReturnParticleGlobalID() const { return particleglobalid_; };

    //! return reference to particle states
    inline const ParticleStates& ReturnParticleStates() const { return particlestates_; };

    //! return bin gid
    inline int ReturnBinGid() const { return bingid_; };

    //! return index
    inline int ReturnContainerIndex() const { return containerindex_; };

   private:
    //! particle type enum
    TypeEnum particletype_;

    //! particle global id
    int particleglobalid_;

    //! maps particle state enum to state
    ParticleStates particlestates_;

    //! bin gid particle is located in
    int bingid_;

    //! index of particle in container (owned or ghosted depending on case)
    int containerindex_;
  };

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
