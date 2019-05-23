/*-----------------------------------------------------------*/
/*!

\brief small data containers for beam interaction framework

\maintainer Jonas Eichinger

\level 3

*/
/*-----------------------------------------------------------*/

#include "beaminteraction_data.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_pack_buffer.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamInteractionParams::BeamInteractionParams()
    : isinit_(false), issetup_(false), rep_strategy_(INPAR::BEAMINTERACTION::repstr_adaptive)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamInteractionParams::Init()
{
  issetup_ = false;

  Teuchos::ParameterList const& params_list = DRT::Problem::Instance()->BeamInteractionParams();

  rep_strategy_ = DRT::INPUT::IntegralValue<INPAR::BEAMINTERACTION::RepartitionStrategy>(
      params_list, "REPARTITIONSTRATEGY");

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
void BEAMINTERACTION::DATA::CrosslinkerData::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack id
  DRT::ParObject::AddtoPack(data, id_);
  // pack position
  DRT::ParObject::AddtoPack(data, pos_);
  // pack bspot status
  DRT::ParObject::AddtoPack(data, bspots_);
  // pack number of bonds
  DRT::ParObject::AddtoPack(data, numbond_);


  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::DATA::CrosslinkerData::Unpack(std::vector<char> const& data)
{
  std::vector<char>::size_type position = 0;

  // id
  DRT::ParObject::ExtractfromPack(position, data, id_);
  // position
  DRT::ParObject::ExtractfromPack(position, data, pos_);
  // bspot status
  DRT::ParObject::ExtractfromPack(position, data, bspots_);
  // number of bonds
  DRT::ParObject::ExtractfromPack(position, data, numbond_);

  if (position != data.size())
    dserror("Mismatch in size of data %d and position %d", static_cast<int>(data.size()), position);

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
void BEAMINTERACTION::DATA::BeamData::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack id
  DRT::ParObject::AddtoPack(data, id_);
  // pack bspotpos
  DRT::ParObject::AddtoPack(data, bspotpos_);
  // pack bspottriad
  DRT::ParObject::AddtoPack(data, bspottriad_);
  // pack bspotstatus_
  DRT::ParObject::AddtoPack(data, bspotstatus_);
  //  // pack filamenttype
  //  DRT::ParObject::AddtoPack( data, filamenttype_ );

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::DATA::BeamData::Unpack(std::vector<char> const& data)
{
  std::vector<char>::size_type position = 0;

  // id
  DRT::ParObject::ExtractfromPack(position, data, id_);
  // bspotpos
  DRT::ParObject::ExtractfromPack(position, data, bspotpos_);
  // bspottriad
  DRT::ParObject::ExtractfromPack(position, data, bspottriad_);
  // bspotstatus
  DRT::ParObject::ExtractfromPack(position, data, bspotstatus_);
  //  // filamenttype
  //  filamenttype_ = static_cast<INPAR::BEAMINTERACTION::FilamentType>(
  //      DRT::ParObject::ExtractInt(position,data) );


  if (position != data.size())
    dserror("Mismatch in size of data %d and position %d", static_cast<int>(data.size()), position);

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
void BEAMINTERACTION::DATA::BindEventData::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  DRT::ParObject::AddtoPack(data, clgid_);

  DRT::ParObject::AddtoPack(data, elegid_);

  DRT::ParObject::AddtoPack(data, bspotlocn_);

  DRT::ParObject::AddtoPack(data, requestproc_);

  DRT::ParObject::AddtoPack(data, permission_);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::DATA::BindEventData::Unpack(std::vector<char> const& data)
{
  std::vector<char>::size_type position = 0;

  DRT::ParObject::ExtractfromPack(position, data, clgid_);

  DRT::ParObject::ExtractfromPack(position, data, elegid_);

  DRT::ParObject::ExtractfromPack(position, data, bspotlocn_);

  DRT::ParObject::ExtractfromPack(position, data, requestproc_);

  DRT::ParObject::ExtractfromPack(position, data, permission_);

  if (position != data.size())
    dserror("Mismatch in size of data %d and position %d", static_cast<int>(data.size()), position);

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
void BEAMINTERACTION::DATA::UnBindEventData::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  DRT::ParObject::AddtoPack(data, clgid_);

  DRT::ParObject::AddtoPack(data, eletoupdate_);

  DRT::ParObject::AddtoPack(data, linkertype_);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::DATA::UnBindEventData::Unpack(std::vector<char> const& data)
{
  std::vector<char>::size_type position = 0;

  DRT::ParObject::ExtractfromPack(position, data, clgid_);

  DRT::ParObject::ExtractfromPack(position, data, eletoupdate_);

  linkertype_ = static_cast<INPAR::BEAMINTERACTION::CrosslinkerType>(
      DRT::ParObject::ExtractInt(position, data));

  if (position != data.size())
    dserror("Mismatch in size of data %d and position %d", static_cast<int>(data.size()), position);

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
