/*---------------------------------------------------------------------------*/
/*! \file
\brief particle object for parallel communication
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_engine_object.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEENGINE::ParticleObjectType PARTICLEENGINE::ParticleObjectType::instance_;

CORE::COMM::ParObject* PARTICLEENGINE::ParticleObjectType::Create(const std::vector<char>& data)
{
  ParticleObject* my_particleobject = new ParticleObject();
  my_particleobject->Unpack(data);
  return my_particleobject;
}

PARTICLEENGINE::ParticleObject::ParticleObject()
    : type_(Phase1), globalid_(0), bingid_(-1), index_(-1)
{
  // empty constructor
}

PARTICLEENGINE::ParticleObject::ParticleObject(
    ParticleType type, int globalid, const ParticleStates& states, int bingid, int index)
    : type_(type), globalid_(globalid), states_(states), bingid_(bingid), index_(index)
{
  // empty constructor
}

void PARTICLEENGINE::ParticleObject::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // particletype_
  AddtoPack(data, type_);

  // particleglobalid_
  AddtoPack(data, globalid_);

  // particle states
  int numstates = states_.size();
  AddtoPack(data, numstates);
  for (int i = 0; i < numstates; ++i) AddtoPack(data, states_[i]);

  // bingid_
  AddtoPack(data, bingid_);

  // containerindex_
  AddtoPack(data, index_);
}

void PARTICLEENGINE::ParticleObject::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // particletype_
  ExtractfromPack(position, data, type_);

  // particleglobalid_
  ExtractfromPack(position, data, globalid_);

  // particle states
  int numstates = 0;
  ExtractfromPack(position, data, numstates);
  states_.resize(numstates);
  for (int i = 0; i < numstates; ++i) ExtractfromPack(position, data, states_[i]);

  // bingid_
  ExtractfromPack(position, data, bingid_);

  // containerindex_
  ExtractfromPack(position, data, index_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", static_cast<int>(data.size()), position);
}
