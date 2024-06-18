/*---------------------------------------------------------------------------*/
/*! \file
\brief particle object for parallel communication
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_engine_object.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEENGINE::ParticleObjectType PARTICLEENGINE::ParticleObjectType::instance_;

Core::Communication::ParObject* PARTICLEENGINE::ParticleObjectType::Create(
    const std::vector<char>& data)
{
  ParticleObject* my_particleobject = new ParticleObject();
  my_particleobject->unpack(data);
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

void PARTICLEENGINE::ParticleObject::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // particletype_
  add_to_pack(data, type_);

  // particleglobalid_
  add_to_pack(data, globalid_);

  // particle states
  int numstates = states_.size();
  add_to_pack(data, numstates);
  for (int i = 0; i < numstates; ++i) add_to_pack(data, states_[i]);

  // bingid_
  add_to_pack(data, bingid_);

  // containerindex_
  add_to_pack(data, index_);
}

void PARTICLEENGINE::ParticleObject::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // particletype_
  extract_from_pack(position, data, type_);

  // particleglobalid_
  extract_from_pack(position, data, globalid_);

  // particle states
  int numstates = 0;
  extract_from_pack(position, data, numstates);
  states_.resize(numstates);
  for (int i = 0; i < numstates; ++i) extract_from_pack(position, data, states_[i]);

  // bingid_
  extract_from_pack(position, data, bingid_);

  // containerindex_
  extract_from_pack(position, data, index_);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", static_cast<int>(data.size()), position);
}

FOUR_C_NAMESPACE_CLOSE
