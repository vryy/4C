/*---------------------------------------------------------------------------*/
/*! \file
\brief particle object for parallel communication
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_object.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEENGINE::ParticleObjectType PARTICLEENGINE::ParticleObjectType::instance_;

DRT::ParObject* PARTICLEENGINE::ParticleObjectType::Create(const std::vector<char>& data)
{
  PARTICLEENGINE::ParticleObject* my_particleobject = new PARTICLEENGINE::ParticleObject();
  my_particleobject->Unpack(data);
  return my_particleobject;
}

PARTICLEENGINE::ParticleObject::ParticleObject()
    : particletype_(PARTICLEENGINE::Phase1), particleglobalid_(0), bingid_(-1), containerindex_(-1)
{
  // empty constructor
}

PARTICLEENGINE::ParticleObject::ParticleObject(TypeEnum particletype, int particleglobalid,
    const ParticleStates& particlestates, int bingid, int containerindex)
    : particletype_(particletype),
      particleglobalid_(particleglobalid),
      particlestates_(particlestates),
      bingid_(bingid),
      containerindex_(containerindex)
{
  // empty constructor
}

void PARTICLEENGINE::ParticleObject::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // particletype_
  AddtoPack(data, particletype_);

  // particleglobalid_
  AddtoPack(data, particleglobalid_);

  // particle states
  int numstates = particlestates_.size();
  AddtoPack(data, numstates);
  for (int i = 0; i < numstates; ++i) AddtoPack(data, particlestates_[i]);

  // bingid_
  AddtoPack(data, bingid_);

  // containerindex_
  AddtoPack(data, containerindex_);
}

void PARTICLEENGINE::ParticleObject::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // particletype_
  ExtractfromPack(position, data, particletype_);

  // particleglobalid_
  ExtractfromPack(position, data, particleglobalid_);

  // particle states
  int numstates = 0;
  ExtractfromPack(position, data, numstates);
  particlestates_.resize(numstates);
  for (int i = 0; i < numstates; ++i) ExtractfromPack(position, data, particlestates_[i]);

  // bingid_
  ExtractfromPack(position, data, bingid_);

  // containerindex_
  ExtractfromPack(position, data, containerindex_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", static_cast<int>(data.size()), position);
}
