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
      rep_strategy_(INPAR::BEAMINTERACTION::repstr_adaptive),
      search_strategy_(INPAR::BEAMINTERACTION::SearchStrategy::bruteforce_with_binning)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamInteractionParams::Init()
{
  issetup_ = false;

  Teuchos::ParameterList const& params_list = GLOBAL::Problem::Instance()->BeamInteractionParams();

  rep_strategy_ = CORE::UTILS::IntegralValue<INPAR::BEAMINTERACTION::RepartitionStrategy>(
      params_list, "REPARTITIONSTRATEGY");

  search_strategy_ = Teuchos::getIntegralValue<INPAR::BEAMINTERACTION::SearchStrategy>(
      params_list, "SEARCH_STRATEGY");

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamInteractionParams::Setup()
{
  CheckInit();

  // empty for now

  issetup_ = true;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::DATA::CrosslinkerData::CrosslinkerData() : id_(-1), pos_(true), numbond_(0)
{
  bspots_ = std::vector<std::pair<int, int>>(2, std::pair<int, int>(-1, -1));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::DATA::CrosslinkerData::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack id
  CORE::COMM::ParObject::AddtoPack(data, id_);
  // pack position
  CORE::COMM::ParObject::AddtoPack(data, pos_);
  // pack bspot status
  CORE::COMM::ParObject::AddtoPack(data, bspots_);
  // pack number of bonds
  CORE::COMM::ParObject::AddtoPack(data, numbond_);


  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::DATA::CrosslinkerData::Unpack(std::vector<char> const& data)
{
  std::vector<char>::size_type position = 0;

  // id
  CORE::COMM::ParObject::ExtractfromPack(position, data, id_);
  // position
  CORE::COMM::ParObject::ExtractfromPack(position, data, pos_);
  // bspot status
  CORE::COMM::ParObject::ExtractfromPack(position, data, bspots_);
  // number of bonds
  CORE::COMM::ParObject::ExtractfromPack(position, data, numbond_);

  if (position != data.size())
    FOUR_C_THROW(
        "Mismatch in size of data %d and position %d", static_cast<int>(data.size()), position);

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::DATA::BeamData::BeamData() : id_(-1)
{
  bspotpos_.clear();
  bspottriad_.clear();
  bspotstatus_.clear();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::DATA::BeamData::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack id
  CORE::COMM::ParObject::AddtoPack(data, id_);
  // pack bspotpos
  CORE::COMM::ParObject::AddtoPack(data, bspotpos_);
  // pack bspottriad
  CORE::COMM::ParObject::AddtoPack(data, bspottriad_);
  // pack bspotstatus_
  CORE::COMM::ParObject::AddtoPack(data, bspotstatus_);
  //  // pack filamenttype
  //  CORE::COMM::ParObject::AddtoPack( data, filamenttype_ );

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::DATA::BeamData::Unpack(std::vector<char> const& data)
{
  std::vector<char>::size_type position = 0;

  // id
  CORE::COMM::ParObject::ExtractfromPack(position, data, id_);
  // bspotpos
  CORE::COMM::ParObject::ExtractfromPack(position, data, bspotpos_);
  // bspottriad
  CORE::COMM::ParObject::ExtractfromPack(position, data, bspottriad_);
  // bspotstatus
  CORE::COMM::ParObject::ExtractfromPack(position, data, bspotstatus_);
  //  // filamenttype
  //  filamenttype_ = static_cast<INPAR::BEAMINTERACTION::FilamentType>(
  //      CORE::COMM::ParObject::ExtractInt(position,data) );


  if (position != data.size())
    FOUR_C_THROW(
        "Mismatch in size of data %d and position %d", static_cast<int>(data.size()), position);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::DATA::BindEventData::BindEventData()
    : clgid_(-1), elegid_(-1), bspotlocn_(-1), requestproc_(-1), permission_(-1)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::DATA::BindEventData::Init(
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
void BEAMINTERACTION::DATA::BindEventData::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  CORE::COMM::ParObject::AddtoPack(data, clgid_);

  CORE::COMM::ParObject::AddtoPack(data, elegid_);

  CORE::COMM::ParObject::AddtoPack(data, bspotlocn_);

  CORE::COMM::ParObject::AddtoPack(data, requestproc_);

  CORE::COMM::ParObject::AddtoPack(data, permission_);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::DATA::BindEventData::Unpack(std::vector<char> const& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ParObject::ExtractfromPack(position, data, clgid_);

  CORE::COMM::ParObject::ExtractfromPack(position, data, elegid_);

  CORE::COMM::ParObject::ExtractfromPack(position, data, bspotlocn_);

  CORE::COMM::ParObject::ExtractfromPack(position, data, requestproc_);

  CORE::COMM::ParObject::ExtractfromPack(position, data, permission_);

  if (position != data.size())
    FOUR_C_THROW(
        "Mismatch in size of data %d and position %d", static_cast<int>(data.size()), position);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::DATA::UnBindEventData::UnBindEventData()
    : clgid_(-1),
      eletoupdate_(std::make_pair(-1, -1)),
      linkertype_(INPAR::BEAMINTERACTION::linkertype_arbitrary)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::DATA::UnBindEventData::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  CORE::COMM::ParObject::AddtoPack(data, clgid_);

  CORE::COMM::ParObject::AddtoPack(data, eletoupdate_);

  CORE::COMM::ParObject::AddtoPack(data, linkertype_);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::DATA::UnBindEventData::Unpack(std::vector<char> const& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ParObject::ExtractfromPack(position, data, clgid_);

  CORE::COMM::ParObject::ExtractfromPack(position, data, eletoupdate_);

  linkertype_ = static_cast<INPAR::BEAMINTERACTION::CrosslinkerType>(
      CORE::COMM::ParObject::ExtractInt(position, data));

  if (position != data.size())
    FOUR_C_THROW(
        "Mismatch in size of data %d and position %d", static_cast<int>(data.size()), position);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
BEAMINTERACTION::DATA::BspotLinkerData::BspotLinkerData()
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
bool BEAMINTERACTION::DATA::BspotLinkerData::SameAs(BspotLinkerData bspotlinker)
{
  if (bspotlinker.GetEleGid1() == elegid_1_ and bspotlinker.GetEleGid2() == elegid_2_ and
      bspotlinker.GetLocBspotId1() == locbspot_1_ and
      bspotlinker.GetLocBspotId2() == locbspot_2_ and bspotlinker.GetMatId() == mat_id_ and
      bspotlinker.GetType() == type_)
    return true;

  return false;
}

FOUR_C_NAMESPACE_CLOSE
