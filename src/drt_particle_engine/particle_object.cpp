/*---------------------------------------------------------------------------*/
/*!
\file particle_object.cpp

\brief particle object class for communication

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_object.H"

/*---------------------------------------------------------------------------*
 | define static class member                                 sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEENGINE::ParticleObjectType PARTICLEENGINE::ParticleObjectType::instance_;

/*---------------------------------------------------------------------------*
 | create particle object                                     sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
DRT::ParObject* PARTICLEENGINE::ParticleObjectType::Create(const std::vector<char>& data)
{
  PARTICLEENGINE::ParticleObject* my_particleobject = new PARTICLEENGINE::ParticleObject();
  my_particleobject->Unpack(data);
  return my_particleobject;
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEENGINE::ParticleObject::ParticleObject()
    : particletype_(PARTICLEENGINE::Phase1), bingid_(-1), containerindex_(-1)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init particle object                                       sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleObject::Init(TypeEnum particletype,
    const std::map<StateEnum, std::vector<double>>& particlestates, int bingid, int containerindex)
{
  particletype_ = particletype;
  particlestates_ = particlestates;
  bingid_ = bingid;
  containerindex_ = containerindex;
}

/*---------------------------------------------------------------------------*
 | pack this class so it can be communicated                  sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleObject::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // particletype_
  AddtoPack(data, particletype_);

  // particle_
  AddtoPack(data, particlestates_);

  // bingid_
  AddtoPack(data, bingid_);

  // containerindex_
  AddtoPack(data, containerindex_);
}

/*---------------------------------------------------------------------------*
 | unpack data from a char vector into this class             sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::ParticleObject::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // particletype_
  ExtractfromPack(position, data, particletype_);

  // particle_
  ExtractfromPack(position, data, particlestates_);

  // bingid_
  ExtractfromPack(position, data, bingid_);

  // containerindex_
  ExtractfromPack(position, data, containerindex_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}
