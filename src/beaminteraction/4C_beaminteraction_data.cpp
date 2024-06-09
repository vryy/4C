/*-----------------------------------------------------------*/
/*! \file

\brief small data containers for beam interaction framework


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_beaminteraction_data.hpp"

#include "4C_comm_pack_buffer.hpp"
#include "4C_comm_parobject.hpp"
#include "4C_global_data.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamInteractionParams::BeamInteractionParams()
    : isinit_(false),
      issetup_(false),
      rep_strategy_(Inpar::BEAMINTERACTION::repstr_adaptive),
      search_strategy_(Inpar::BEAMINTERACTION::SearchStrategy::bruteforce_with_binning)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamInteractionParams::Init()
{
  issetup_ = false;

  Teuchos::ParameterList const& params_list =
      Global::Problem::Instance()->beam_interaction_params();

  rep_strategy_ = Core::UTILS::IntegralValue<Inpar::BEAMINTERACTION::RepartitionStrategy>(
      params_list, "REPARTITIONSTRATEGY");

  search_strategy_ = Teuchos::getIntegralValue<Inpar::BEAMINTERACTION::SearchStrategy>(
      params_list, "SEARCH_STRATEGY");

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamInteractionParams::Setup()
{
  check_init();

  // empty for now

  issetup_ = true;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::Data::CrosslinkerData::CrosslinkerData() : id_(-1), pos_(true), numbond_(0)
{
  bspots_ = std::vector<std::pair<int, int>>(2, std::pair<int, int>(-1, -1));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::Data::CrosslinkerData::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack id
  Core::Communication::ParObject::add_to_pack(data, id_);
  // pack position
  Core::Communication::ParObject::add_to_pack(data, pos_);
  // pack bspot status
  Core::Communication::ParObject::add_to_pack(data, bspots_);
  // pack number of bonds
  Core::Communication::ParObject::add_to_pack(data, numbond_);


  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::Data::CrosslinkerData::Unpack(std::vector<char> const& data)
{
  std::vector<char>::size_type position = 0;

  // id
  Core::Communication::ParObject::extract_from_pack(position, data, id_);
  // position
  Core::Communication::ParObject::extract_from_pack(position, data, pos_);
  // bspot status
  Core::Communication::ParObject::extract_from_pack(position, data, bspots_);
  // number of bonds
  Core::Communication::ParObject::extract_from_pack(position, data, numbond_);

  if (position != data.size())
    FOUR_C_THROW(
        "Mismatch in size of data %d and position %d", static_cast<int>(data.size()), position);

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::Data::BeamData::BeamData() : id_(-1)
{
  bspotpos_.clear();
  bspottriad_.clear();
  bspotstatus_.clear();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::Data::BeamData::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack id
  Core::Communication::ParObject::add_to_pack(data, id_);
  // pack bspotpos
  Core::Communication::ParObject::add_to_pack(data, bspotpos_);
  // pack bspottriad
  Core::Communication::ParObject::add_to_pack(data, bspottriad_);
  // pack bspotstatus_
  Core::Communication::ParObject::add_to_pack(data, bspotstatus_);
  //  // pack filamenttype
  //  Core::Communication::ParObject::add_to_pack( data, filamenttype_ );

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::Data::BeamData::Unpack(std::vector<char> const& data)
{
  std::vector<char>::size_type position = 0;

  // id
  Core::Communication::ParObject::extract_from_pack(position, data, id_);
  // bspotpos
  Core::Communication::ParObject::extract_from_pack(position, data, bspotpos_);
  // bspottriad
  Core::Communication::ParObject::extract_from_pack(position, data, bspottriad_);
  // bspotstatus
  Core::Communication::ParObject::extract_from_pack(position, data, bspotstatus_);
  //  // filamenttype
  //  filamenttype_ = static_cast<Inpar::BEAMINTERACTION::FilamentType>(
  //      Core::Communication::ParObject::ExtractInt(position,data) );


  if (position != data.size())
    FOUR_C_THROW(
        "Mismatch in size of data %d and position %d", static_cast<int>(data.size()), position);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::Data::BindEventData::BindEventData()
    : clgid_(-1), elegid_(-1), bspotlocn_(-1), requestproc_(-1), permission_(-1)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::Data::BindEventData::Init(
    int clgid, int elegid, int bspotlocn, int requestproc, int permission)
{
  clgid_ = clgid;
  elegid_ = elegid;
  bspotlocn_ = bspotlocn;
  requestproc_ = requestproc;
  permission_ = permission;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::Data::BindEventData::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  Core::Communication::ParObject::add_to_pack(data, clgid_);

  Core::Communication::ParObject::add_to_pack(data, elegid_);

  Core::Communication::ParObject::add_to_pack(data, bspotlocn_);

  Core::Communication::ParObject::add_to_pack(data, requestproc_);

  Core::Communication::ParObject::add_to_pack(data, permission_);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::Data::BindEventData::Unpack(std::vector<char> const& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ParObject::extract_from_pack(position, data, clgid_);

  Core::Communication::ParObject::extract_from_pack(position, data, elegid_);

  Core::Communication::ParObject::extract_from_pack(position, data, bspotlocn_);

  Core::Communication::ParObject::extract_from_pack(position, data, requestproc_);

  Core::Communication::ParObject::extract_from_pack(position, data, permission_);

  if (position != data.size())
    FOUR_C_THROW(
        "Mismatch in size of data %d and position %d", static_cast<int>(data.size()), position);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::Data::UnBindEventData::UnBindEventData()
    : clgid_(-1),
      eletoupdate_(std::make_pair(-1, -1)),
      linkertype_(Inpar::BEAMINTERACTION::linkertype_arbitrary)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::Data::UnBindEventData::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  Core::Communication::ParObject::add_to_pack(data, clgid_);

  Core::Communication::ParObject::add_to_pack(data, eletoupdate_);

  Core::Communication::ParObject::add_to_pack(data, linkertype_);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::Data::UnBindEventData::Unpack(std::vector<char> const& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ParObject::extract_from_pack(position, data, clgid_);

  Core::Communication::ParObject::extract_from_pack(position, data, eletoupdate_);

  linkertype_ = static_cast<Inpar::BEAMINTERACTION::CrosslinkerType>(
      Core::Communication::ParObject::ExtractInt(position, data));

  if (position != data.size())
    FOUR_C_THROW(
        "Mismatch in size of data %d and position %d", static_cast<int>(data.size()), position);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
BEAMINTERACTION::Data::BspotLinkerData::BspotLinkerData()
    : elegid_1_(-1),
      elegid_2_(-1),
      locbspot_1_(-1),
      locbspot_2_(-1),
      type_(-1),
      mat_id_(-1),
      number_of_bonds_1_(-1),
      number_of_bonds_2_(-1)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::Data::BspotLinkerData::SameAs(BspotLinkerData bspotlinker)
{
  if (bspotlinker.GetEleGid1() == elegid_1_ and bspotlinker.GetEleGid2() == elegid_2_ and
      bspotlinker.GetLocBspotId1() == locbspot_1_ and
      bspotlinker.GetLocBspotId2() == locbspot_2_ and bspotlinker.GetMatId() == mat_id_ and
      bspotlinker.GetType() == type_)
    return true;

  return false;
}

FOUR_C_NAMESPACE_CLOSE
