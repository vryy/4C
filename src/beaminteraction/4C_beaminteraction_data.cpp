/*-----------------------------------------------------------*/
/*! \file

\brief small data containers for beam interaction framework


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_beaminteraction_data.hpp"

#include "4C_comm_pack_buffer.hpp"
#include "4C_comm_pack_helpers.hpp"
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
void BEAMINTERACTION::BeamInteractionParams::init()
{
  issetup_ = false;

  Teuchos::ParameterList const& params_list =
      Global::Problem::instance()->beam_interaction_params();

  rep_strategy_ = Teuchos::getIntegralValue<Inpar::BEAMINTERACTION::RepartitionStrategy>(
      params_list, "REPARTITIONSTRATEGY");

  search_strategy_ = Teuchos::getIntegralValue<Inpar::BEAMINTERACTION::SearchStrategy>(
      params_list, "SEARCH_STRATEGY");

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamInteractionParams::setup()
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
void BEAMINTERACTION::Data::CrosslinkerData::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack id
  add_to_pack(data, id_);
  // pack position
  add_to_pack(data, pos_);
  // pack bspot status
  add_to_pack(data, bspots_);
  // pack number of bonds
  add_to_pack(data, numbond_);


  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::Data::CrosslinkerData::unpack(Core::Communication::UnpackBuffer& buffer)
{
  // id
  extract_from_pack(buffer, id_);
  // position
  extract_from_pack(buffer, pos_);
  // bspot status
  extract_from_pack(buffer, bspots_);
  // number of bonds
  extract_from_pack(buffer, numbond_);

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");

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
void BEAMINTERACTION::Data::BeamData::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack id
  add_to_pack(data, id_);
  // pack bspotpos
  add_to_pack(data, bspotpos_);
  // pack bspottriad
  add_to_pack(data, bspottriad_);
  // pack bspotstatus_
  add_to_pack(data, bspotstatus_);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::Data::BeamData::unpack(Core::Communication::UnpackBuffer& buffer)
{
  // id
  extract_from_pack(buffer, id_);
  // bspotpos
  extract_from_pack(buffer, bspotpos_);
  // bspottriad
  extract_from_pack(buffer, bspottriad_);
  // bspotstatus
  extract_from_pack(buffer, bspotstatus_);

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");

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
void BEAMINTERACTION::Data::BindEventData::init(
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
void BEAMINTERACTION::Data::BindEventData::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  add_to_pack(data, clgid_);

  add_to_pack(data, elegid_);

  add_to_pack(data, bspotlocn_);

  add_to_pack(data, requestproc_);

  add_to_pack(data, permission_);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::Data::BindEventData::unpack(Core::Communication::UnpackBuffer& buffer)
{
  extract_from_pack(buffer, clgid_);

  extract_from_pack(buffer, elegid_);

  extract_from_pack(buffer, bspotlocn_);

  extract_from_pack(buffer, requestproc_);

  extract_from_pack(buffer, permission_);

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");

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
void BEAMINTERACTION::Data::UnBindEventData::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  add_to_pack(data, clgid_);

  add_to_pack(data, eletoupdate_);

  add_to_pack(data, linkertype_);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::Data::UnBindEventData::unpack(Core::Communication::UnpackBuffer& buffer)
{
  extract_from_pack(buffer, clgid_);

  extract_from_pack(buffer, eletoupdate_);

  linkertype_ = static_cast<Inpar::BEAMINTERACTION::CrosslinkerType>(extract_int(buffer));

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");

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
bool BEAMINTERACTION::Data::BspotLinkerData::same_as(BspotLinkerData bspotlinker)
{
  if (bspotlinker.get_ele_gid1() == elegid_1_ and bspotlinker.get_ele_gid2() == elegid_2_ and
      bspotlinker.get_loc_bspot_id1() == locbspot_1_ and
      bspotlinker.get_loc_bspot_id2() == locbspot_2_ and bspotlinker.get_mat_id() == mat_id_ and
      bspotlinker.get_type() == type_)
    return true;

  return false;
}

FOUR_C_NAMESPACE_CLOSE
