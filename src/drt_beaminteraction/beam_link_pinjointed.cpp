/*----------------------------------------------------------------------*/
/*!
\file beam_link_pinjointed.cpp

\brief connecting beam linked by pin joint

\level 3

\maintainer Jonas Eichinger
*/
/*----------------------------------------------------------------------*/

#include "../drt_fem_general/largerotations.H"

#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_RCP.hpp>

#include "beam_link.H"

#include "beam_link_pinjointed.H"

#include "beam_link_beam3r_lin2_pinjointed.H"


BEAMINTERACTION::BeamLinkPinJointedType BEAMINTERACTION::BeamLinkPinJointedType::instance_;


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamLinkPinJointed::BeamLinkPinJointed() :
    BeamLink()
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkPinJointed::Init(
    int id,
    const std::vector<std::pair<int, int> >& eleids,
    const std::vector<LINALG::Matrix<3,1> >& initpos)
{
  issetup_ = false;

  BeamLink::Init( id, eleids, initpos );

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkPinJointed::Setup( const int matnum )
{
  CheckInit();

  // the flag issetup_ will be set in the derived method!
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkPinJointed::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  BeamLink::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkPinJointed::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  BeamLink::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkPinJointed::ResetState(
    std::vector<LINALG::Matrix<3,1> >& bspotpos)
{
  CheckInitSetup();

  BeamLink::ResetState(bspotpos);

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed>
BEAMINTERACTION::BeamLinkPinJointed::Create()
{
  // for now, we always use a 2-noded linear Reissner element
  return Teuchos::rcp( new BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed() );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkPinJointed::Print(std::ostream& out) const
{
  CheckInit();

  out << "\nBeamLinkPinJointed (ID " << Id() << "):";
  out << "\nbspotIds_[0] = ";
  out << "EleGID " << GetEleGid(0) << " locbspotnum " << GetLocBSpotNum(0);
  out << "\nbspotIds_[1] = ";
  out << "EleGID " << GetEleGid(1)  << " locbspotnum " << GetLocBSpotNum(1);
  out << "\n";
  out << "\nbspotpos1_ = ";
  GetBindSpotPos1().Print(out);
  out << "\nbspotpos2_ = ";
  GetBindSpotPos1().Print(out);
  out << "\n";
}
