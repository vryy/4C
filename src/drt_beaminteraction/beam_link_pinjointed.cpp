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

#include "../drt_inpar/inpar_beaminteraction.H"

#include <Teuchos_RCP.hpp>

#include "beam_link.H"

#include "beam_link_pinjointed.H"
#include "beam_link_truss.H"

#include "beam_link_beam3r_lin2_pinjointed.H"


BEAMINTERACTION::BeamLinkPinJointedType BEAMINTERACTION::BeamLinkPinJointedType::instance_;


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamLinkPinJointed::BeamLinkPinJointed() :
    BeamLink()
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
BEAMINTERACTION::BeamLinkPinJointed::BeamLinkPinJointed( const BEAMINTERACTION::BeamLinkPinJointed & old) :
    BEAMINTERACTION::BeamLink(old)
{
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkPinJointed::Init(
    int id,
    const std::vector<std::pair<int, int> >& eleids,
    const std::vector<LINALG::Matrix<3,1> >& initpos,
    const std::vector<LINALG::Matrix<3,3> >& inittriad,
    double timelinkwasset)
{
  issetup_ = false;

  BeamLink::Init( id, eleids, initpos, inittriad, timelinkwasset );

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
    std::vector<LINALG::Matrix<3,1> >& bspotpos,
    std::vector<LINALG::Matrix<3,3> >& bspottriad)
{
  CheckInitSetup();

  BeamLink::ResetState( bspotpos, bspottriad );

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed>
BEAMINTERACTION::BeamLinkPinJointed::Create(
    INPAR::BEAMINTERACTION::JointType type)
{
  if ( type == INPAR::BEAMINTERACTION::beam3r_lin2_pin )
    return Teuchos::rcp( new BEAMINTERACTION::BeamLinkBeam3rLin2PinJointed() );
  else if ( type == INPAR::BEAMINTERACTION::truss )
    return Teuchos::rcp( new BEAMINTERACTION::BeamLinkTruss() );
  else
    dserror("instantiation of new BeamLinkPinJointed object failed due to "
        "unknown type of linker");

  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkPinJointed::Print(std::ostream& out) const
{
  CheckInit();

  BeamLink::Print( out );

  out << "\n";
}
