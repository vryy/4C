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
    bspotpos2_(true)
{
  bspotIds_.clear();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLink::Init(
    const int id,
    const std::vector<std::pair<int, int> >& eleids,
    const std::vector<LINALG::Matrix<3,1> >& initpos)
{
  issetup_ = false;

  id_ = id;
  bspotIds_ = eleids;

  bspotpos1_ = initpos[0];
  bspotpos2_ = initpos[1];

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
  // bspotpos1_
  ExtractfromPack(position,data,bspotpos1_);
  // bspotpos2_
  ExtractfromPack(position,data,bspotpos2_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLink::ResetState( std::vector<LINALG::Matrix<3,1> >& bspotpos )
{
  CheckInitSetup();

  /* the two positions of the linkage element coincide with the positions of the
   * binding spots on the parent elements */
  bspotpos1_ = bspotpos[0];
  bspotpos2_ = bspotpos[1];
}
