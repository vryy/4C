/*---------------------------------------------------------------------------*/
/*!
\file particle_object.H

\brief particle object class for communication

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
  class ParticleObjectType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "ParticleObjectType"; };

    static ParticleObjectType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ParticleObjectType instance_;
  };

  //! particle object
  class ParticleObject : public DRT::ParObject
  {
   public:
    //! constructor
    ParticleObject();

    //! virtual destructor
    virtual ~ParticleObject(){};

    //! init particle object
    void Init(TypeEnum particletype,
        const std::map<StateEnum, std::vector<double>>& particlestates_, int bingid = -1,
        int containerindex = -1);

    //! return unique ParObject id
    virtual int UniqueParObjectId() const
    {
      return ParticleObjectType::Instance().UniqueParObjectId();
    };

    //! pack this class so it can be communicated
    virtual void Pack(DRT::PackBuffer& data) const;

    //! unpack data from a char vector into this class
    virtual void Unpack(const std::vector<char>& data);

    //! return particle type
    TypeEnum ReturnParticleType() const { return particletype_; };

    //! return particle states
    std::map<StateEnum, std::vector<double>> ReturnParticleStates() const
    {
      return particlestates_;
    };

    //! return bin gid
    int ReturnBinGid() const { return bingid_; };

    //! return index
    int ReturnContainerIndex() const { return containerindex_; };

   private:
    //! particle type enum
    TypeEnum particletype_;

    //! maps particle state enum to state
    std::map<StateEnum, std::vector<double>> particlestates_;

    //! bin gid particle is located in
    int bingid_;

    //! index of particle in container (owned or ghosted depending on case)
    int containerindex_;
  };

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
