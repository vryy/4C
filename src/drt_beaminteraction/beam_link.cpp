/*----------------------------------------------------------------------*/
/*!
\file beam_link.cpp

\brief One beam-to-beam pair (two beam elements) connected by a mechanical link

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/

#include "../drt_fem_general/largerotations.H"

#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_RCP.hpp>

#include "beam_link.H"

BEAMINTERACTION::BeamLinkType BEAMINTERACTION::BeamLinkType::instance_;


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamLink::BeamLink() :
    ParObject(),
    isinit_(false),
    issetup_(false),
    id_(-1),
    bspotpos1_(true),
    bspotpos2_(true),
    linkertype_(INPAR::BEAMINTERACTION::linkertype_arbitrary),
    timelinkwasset_(-1.0),
    reflength_(-1.0)
{
  bspotIds_.clear();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
BEAMINTERACTION::BeamLink::BeamLink( const BEAMINTERACTION::BeamLink & old) :
    ParObject(old),
    isinit_(old.isinit_),
    issetup_(old.issetup_),
    id_(old.id_),
    bspotIds_(old.bspotIds_),
    bspotpos1_(old.bspotpos1_),
    bspotpos2_(old.bspotpos2_),
    linkertype_(old.linkertype_),
    timelinkwasset_(old.timelinkwasset_),
    reflength_(old.reflength_)
{
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLink::Init(
    const int id,
    const std::vector<std::pair<int, int> >& eleids,
    const std::vector<LINALG::Matrix< 3, 1 > >& initpos,
    const std::vector<LINALG::Matrix< 3, 3 > >& inittriad,
    INPAR::BEAMINTERACTION::CrosslinkerType linkertype,
    double timelinkwasset )
{
  issetup_ = false;

  id_ = id;
  bspotIds_ = eleids;

  bspotpos1_ = initpos[0];
  bspotpos2_ = initpos[1];

  linkertype_ = linkertype;

  timelinkwasset_ = timelinkwasset;

  reflength_ = 0.0;
  for ( unsigned int i = 0; i < 3; ++i )
    reflength_ += ( initpos[1](i) - initpos[0](i) ) * ( initpos[1](i) - initpos[0](i) );
  reflength_ = sqrt(reflength_);

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLink::Setup( const int matnum )
{
  CheckInit();

  // the flag issetup_ will be set in the derived method!
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLink::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // isinit_
  AddtoPack(data,isinit_);
  // issetup
  AddtoPack(data,issetup_);
  // add id
  AddtoPack(data,id_);

  // add eleids_
  AddtoPack(data,bspotIds_);
  // bspotpos1_
  AddtoPack(data,bspotpos1_);
  // bspotpos2_
  AddtoPack(data,bspotpos2_);
  // linkertype
  AddtoPack(data,linkertype_);
  // timelinkwasset
  AddtoPack(data,timelinkwasset_);
  // reflength
  AddtoPack(data,reflength_);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLink::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // isinit_
  isinit_ = DRT::ParObject::ExtractInt(position,data);
  // issetup
  issetup_ = DRT::ParObject::ExtractInt(position,data);
  // id_
  ExtractfromPack(position,data,id_);

  // eleids_
  ExtractfromPack(position,data,bspotIds_);
  // bspotpos1
  ExtractfromPack(position,data,bspotpos1_);
  // bspotpos2
  ExtractfromPack(position,data,bspotpos2_);
  // linkertype
  linkertype_ = static_cast<INPAR::BEAMINTERACTION::CrosslinkerType>( ExtractInt(position,data) );
  // timelinkwasset
  ExtractfromPack(position,data,timelinkwasset_);
  // reflength
  ExtractfromPack(position,data,reflength_);

  if ( position != data.size() )
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(), position );

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLink::ResetState(
    std::vector<LINALG::Matrix<3,1> >& bspotpos,
    std::vector<LINALG::Matrix<3,3> >& bspottriad)
{
  CheckInitSetup();

  /* the two positions of the linkage element coincide with the positions of the
   * binding spots on the parent elements */
  bspotpos1_ = bspotpos[0];
  bspotpos2_ = bspotpos[1];
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLink::Print( std::ostream & out ) const
{
  CheckInit();

  out << "\nBeamLinkRigidJointed (ID " << Id() << "):";
  out << "\nbspotIds_[0] = ";
  out << "EleGID " << GetEleGid(0) << " locbspotnum " << GetLocBSpotNum(0);
  out << "\nbspotIds_[1] = ";
  out << "EleGID " << GetEleGid(1)  << " locbspotnum " << GetLocBSpotNum(1);
  out << "\n";
  out << "\nbspotpos1_ = ";
  GetBindSpotPos1().Print(out);
  out << "\nbspotpos2_ = ";
  GetBindSpotPos2().Print(out);

  out << "\n";
}
